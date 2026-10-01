// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <cmath>
#include <Eigen/LU>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <mutex>

#include "BaseLib/Logging.h"
#include "ComputeMicroPorosity.h"
#include "ConstitutiveRelations/ConstitutiveModels.h"
#include "ConstitutiveRelations/PotentialExchange.h"
#include "IntegrationPointData.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/PorosityFromMassBalance.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Fem/Interpolation.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ProcessLib/Utils/TransposeInPlace.h"
#include "RichardsMechanicsFEM.h"

namespace ProcessLib
{
namespace RichardsMechanics
{
inline bool isPotentialExchangeEnabled(
    PotentialExchangeParameters const* const potential_exchange_parameters)
{
    return potential_exchange_parameters &&
           potential_exchange_parameters->enabled;
}

inline bool isPotentialExchangeEnabled(
    std::optional<PotentialExchangeParameters> const&
        potential_exchange_parameters)
{
    return isPotentialExchangeEnabled(
        potential_exchange_parameters ? &*potential_exchange_parameters
                                         : nullptr);
}

// Film-pressure coupling (maxwell sec.5) requires the exchange to be enabled
// AND the film_pressure_coupling master flag set. Default OFF -> false, so the
// old vdW/eigenstress path runs unchanged bit-for-bit.
inline bool isFilmPressureCouplingEnabled(
    PotentialExchangeParameters const* const potential_exchange_parameters)
{
    return isPotentialExchangeEnabled(potential_exchange_parameters) &&
           potential_exchange_parameters->film_pressure_coupling;
}

// Drained bulk modulus K [Pa] from a Kelvin stiffness tensor: for any isotropic
// C, m^T C m = 9K with m = identity2, so K = identity2 . (C identity2) / 9. This
// is exactly d sigma'_m / d eps_v = -dp_conf/deps_v, the MECHANICAL stiffness the
// film-pressure coupling needs for the p-u tangent (mu_lR film delta, exchange
// equation) AND its one-Psi transpose, the swelling-stress u-eps tangent
// (+(1-phi_M)*b*K). Using ONE formula keeps the two sides of the Maxwell pair on
// the SAME K. NOTE (2026-06-06): this is the MECHANICAL drained K only; it is no
// longer used as a swelling-stress modulus (the K_sw eigenstress form is retired
// -- the swelling stress is now the transmitted pressure -(1-phi_M)*p_film).
template <int DisplacementDim>
inline double drainedBulkModulusFromStiffness(
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C)
{
    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;
    return identity2.dot(C * identity2) / 9.0;
}

inline double getPotentialPressureTolerance(
    PotentialExchangeParameters const* const potential_exchange_parameters)
{
    if (!isPotentialExchangeEnabled(potential_exchange_parameters))
    {
        return 0.0;
    }

    return potential_exchange_parameters->pressure_tolerance;
}

inline double getPotentialPressureTolerance(
    std::optional<PotentialExchangeParameters> const&
        potential_exchange_parameters)
{
    return getPotentialPressureTolerance(
        potential_exchange_parameters ? &*potential_exchange_parameters
                                         : nullptr);
}

inline void requirePositiveViscosity(char const* caller, double const mu)
{
    if (!(std::isfinite(mu) && mu > 0.0))
    {
        OGS_FATAL("{} requires finite mu > 0, got {:g}.", caller, mu);
    }
}

struct PotentialExchangeUpdateData
{
    YoungLaplaceMacroPotentialData macro_potential;
    PotentialDrivenMassExchangeData exchange;

    double alpha_M_effective = 0.0;
    double mu_LR_active = 0.0;
    double mu_lR_exchange_input = 0.0;
    bool use_macro_potential_for_active_exchange = false;
    bool use_vdw_micro_potential_for_active_exchange = false;
    bool use_fd_jacobian_for_direct_macro_derivative = false;
    double fd_jacobian_perturbation = 0.0;

    // Direct macro derivative (with density dependence through rho_LR), while
    // keeping the microscale state lagged.
    double drho_L_hat_dpL_direct = 0.0;
};

inline PotentialExchangeUpdateData computePotentialExchangeUpdate(
    double const alpha_bar, double const mu, double const p_L_ip,
    double const p_L_m, double const rho_LR, double const beta_LR,
    double const rho_lR_exchange_input = std::numeric_limits<double>::quiet_NaN(),
    double const drho_lR_exchange_input_dpL =
        std::numeric_limits<double>::quiet_NaN(),
    double const pressure_tolerance = 0.0,
    bool const use_macro_potential_for_active_exchange = false,
    bool const use_vdw_micro_potential_for_active_exchange = false,
    double const mu_lR_vdw = 0.0,
    double const dmu_lR_vdw_drho_lR = 0.0,
    bool const use_custom_dmu_lR_vdw_dpL = false,
    double const dmu_lR_vdw_dpL = 0.0,
    bool const use_fd_jacobian_for_direct_macro_derivative = false,
    double const fd_jacobian_perturbation = 1e-8)
{
    requirePositiveViscosity("computePotentialExchangeUpdate", mu);

    PotentialExchangeUpdateData out;

    // Keep the exchange coefficient scaling in mass-density units.
    out.alpha_M_effective = alpha_bar * rho_LR / mu;

    out.macro_potential =
        computeYoungLaplaceMacroPotential(p_L_ip, rho_LR, pressure_tolerance);
    out.use_macro_potential_for_active_exchange =
        use_macro_potential_for_active_exchange;
    out.use_vdw_micro_potential_for_active_exchange =
        use_vdw_micro_potential_for_active_exchange;
    out.use_fd_jacobian_for_direct_macro_derivative =
        use_fd_jacobian_for_direct_macro_derivative;
    out.fd_jacobian_perturbation = fd_jacobian_perturbation;

    // rho_LR depends on liquid pressure in RM through beta_LR = (1/rho) drho/dp.
    double const drho_LR_dpL = rho_LR * beta_LR;
    bool const use_custom_micro_density_for_exchange =
        std::isfinite(rho_lR_exchange_input) && rho_lR_exchange_input > 0.0;
    double const rho_lR_exchange =
        use_custom_micro_density_for_exchange ? rho_lR_exchange_input : rho_LR;
    double const drho_lR_exchange_dpL =
        use_custom_micro_density_for_exchange
            ? (std::isfinite(drho_lR_exchange_input_dpL)
                   ? drho_lR_exchange_input_dpL
                   : 0.0)
            : drho_LR_dpL;

    out.mu_lR_exchange_input = use_vdw_micro_potential_for_active_exchange
                                   ? mu_lR_vdw
                                   : p_L_m / rho_lR_exchange;
    out.mu_LR_active = use_macro_potential_for_active_exchange
                           ? out.macro_potential.mu_LR
                           : p_L_ip / rho_LR;

    out.exchange = computePotentialDrivenMassExchange(
        out.alpha_M_effective, out.mu_LR_active, out.mu_lR_exchange_input);

    if (use_fd_jacobian_for_direct_macro_derivative)
    {
        auto const compute_rho_L_hat = [&](double const p_L_ip_eval,
                                           double const rho_LR_eval)
        {
            auto const macro_potential_eval = computeYoungLaplaceMacroPotential(
                p_L_ip_eval, rho_LR_eval, pressure_tolerance);
            double const alpha_M_effective_eval =
                alpha_bar * rho_LR_eval / mu;
            double const mu_LR_active_eval = use_macro_potential_for_active_exchange
                                                 ? macro_potential_eval.mu_LR
                                                 : p_L_ip_eval / rho_LR_eval;
            double const rho_lR_eval =
                use_custom_micro_density_for_exchange ? rho_lR_exchange
                                                      : rho_LR_eval;
            double const mu_lR_active_eval =
                use_vdw_micro_potential_for_active_exchange
                    ? mu_lR_vdw
                    : p_L_m / rho_lR_eval;
            auto const exchange_eval = computePotentialDrivenMassExchange(
                alpha_M_effective_eval, mu_LR_active_eval, mu_lR_active_eval);
            return -exchange_eval.rho_l_hat;
        };

        double const h =
            fd_jacobian_perturbation * std::max(1.0, std::abs(p_L_ip));
        if (!(h > 0.0) || !std::isfinite(h))
        {
            OGS_FATAL(
                "computePotentialExchangeUpdate requires finite h > 0 for FD Jacobian, got {:g} (from fd_jacobian_perturbation={:g}, p_L_ip={:g}).",
                h, fd_jacobian_perturbation, p_L_ip);
        }

        constexpr double rho_floor = 1e-16;
        double const rho_plus = std::max(rho_floor, rho_LR + drho_LR_dpL * h);
        double const rho_minus = rho_LR - drho_LR_dpL * h;
        double const rho_L_hat_plus = compute_rho_L_hat(p_L_ip + h, rho_plus);
        if (rho_minus > rho_floor)
        {
            double const rho_L_hat_minus =
                compute_rho_L_hat(p_L_ip - h, rho_minus);
            out.drho_L_hat_dpL_direct =
                (rho_L_hat_plus - rho_L_hat_minus) / (2.0 * h);
        }
        else
        {
            double const rho_L_hat = -out.exchange.rho_l_hat;
            out.drho_L_hat_dpL_direct = (rho_L_hat_plus - rho_L_hat) / h;
        }

        return out;
    }

    // alpha_M_effective = alpha_bar * rho_LR / mu (mu dependence is lagged).
    double const dalpha_M_effective_dpL = alpha_bar / mu * drho_LR_dpL;

    double const dmu_LR_dpL = use_macro_potential_for_active_exchange
                                  ? out.macro_potential.dmu_LR_dpLR +
                                        out.macro_potential.dmu_LR_drho_LR *
                                            drho_LR_dpL
                                  : 1.0 / rho_LR -
                                        p_L_ip / (rho_LR * rho_LR) * drho_LR_dpL;

    double const dmu_lR_exchange_input_dpL =
        use_vdw_micro_potential_for_active_exchange
            ? (use_custom_dmu_lR_vdw_dpL
                   ? dmu_lR_vdw_dpL
                   : dmu_lR_vdw_drho_lR * drho_LR_dpL)
            : -p_L_m / (rho_lR_exchange * rho_lR_exchange) *
                  drho_lR_exchange_dpL;

    double const drho_l_hat_dpL_direct =
        out.exchange.drho_l_hat_dalpha_M * dalpha_M_effective_dpL +
        out.exchange.drho_l_hat_dmu_LR * dmu_LR_dpL +
        out.exchange.drho_l_hat_dmu_lR * dmu_lR_exchange_input_dpL;

    out.drho_L_hat_dpL_direct = -drho_l_hat_dpL_direct;
    return out;
}

struct ImplicitMicroWaterContentUpdateData
{
    double n_l = 0.0;
    // Converged micro liquid density (ScalarReferenceMassStorage mode only;
    // NaN otherwise). Threaded into computeImplicitNlDpL so the F1 REV-mass
    // tangent can evaluate the EOS partials at the converged state.
    double rho_lR = std::numeric_limits<double>::quiet_NaN();
    VanDerWaalsMicroPotentialData micro_potential;
    PotentialDrivenMassExchangeData exchange;
    bool converged = true;
};

struct CompatibilityMicroHydraulicOutputData
{
    double p_L_m = 0.0;
    double S_L_m = 0.0;
    double n_l_ref = 0.0;
    VanDerWaalsMicroPotentialData micro_potential;
};

struct PotentialExchangeLocalSolveContext;
inline CompatibilityMicroHydraulicOutputData
computeCompatibilityMicroHydraulicOutput(
    double const n_l, double const rho_LR,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params);

inline double microPotentialSignFactorFromParameters(
    PotentialExchangeParameters const& potential_exchange_params)
{
    return microPotentialSignFactor(potential_exchange_params.micro_potential_convention);
}

// NOTE: the only active overload of computeCompatibilityMicroHydraulicOutput is
// the 4-argument version below (with local_context, defined after
// computeActiveMicroPotential). The local_context overload derives the micro
// liquid density (rho_lR ~ 1100 kg/m³) from the EOS and uses it in both the
// vdW potential formula and in p_L_m = -rho_lR * mu_lR when
// use_micro_liquid_density_for_micro_pressure = true (set in all MS33 PRJs).
// A 3-argument overload without local_context existed here previously but was
// dead code and used bulk rho_LR (~1000 kg/m³) in the vdW denominator —
// a ~10% error — so it was removed (2026-05-22).

struct TransportPorosityUpdateData
{
    double phi_M = 0.0;
    double phi_M_prev = 0.0;
    double phi_m = 0.0;
    double phi_m_prev = 0.0;
};

struct PotentialExchangeLocalSolveContext
{
    double phi = std::numeric_limits<double>::infinity();
    double phi_M_prev = 0.0;
    double phi_m_prev = 0.0;
    double volumetric_strain = 0.0;
    double volumetric_strain_prev = 0.0;
    // Confining pressure p_conf = -tr(sigma_eff)/3 [Pa] (>0 in compression, OGS
    // tension-positive). Film-pressure coupling only (maxwell sec.5). NaN
    // sentinel = "stress not supplied" -> the film term self-disables, so a
    // default-constructed context (e.g. the eigenstress-difference driver, which
    // has no stress in scope) gets NO film term. Set explicitly at the call
    // sites that have EffectiveStressData available.
    double confining_pressure_p_conf = std::numeric_limits<double>::quiet_NaN();
    // Poroelastic Biot coefficient (= poroelastic Biot alpha, the SAME solid-
    // fluid volume partitioning) threaded from the medium MPL biot_coefficient
    // property. Used as the eigenstrain Biot b in the film-pressure coupling
    // (one-Psi consistency; previously a separate film_pressure_biot_b scalar).
    // Default 1.0 -> sensible incompressible-grain fallback for default-
    // constructed contexts (the GP eigenstress-difference driver and tests that
    // do not set it explicitly). Set at the assembly / micro-solve context sites
    // from the evaluated MPL biot_coefficient.
    double biot_coefficient = 1.0;
    // Drained bulk modulus K_drained = identity2:(C_el:identity2)/9 [Pa], the
    // MECHANICAL stiffness dp_conf/deps_v = -dsigma'_m/deps_v. Needed by the
    // INTEGRABLE Maxwell mechanical partner mu_lR_mech (film_pressure_coupling
    // ON) for the 0.5*b*K_drained*eps_v^2 term. NaN sentinel = "not supplied" ->
    // the partner's quadratic (p_conf-chain) piece is dropped, leaving the linear
    // -(Pi + n_l*Pi')*eps_v/rho_lR conjugate (the GP eigenstress-difference driver
    // and tests that have no stiffness in scope get this reduced, finite form).
    // Set at the assembly / micro-solve context sites from the elastic stiffness.
    double drained_bulk_modulus = std::numeric_limits<double>::quiet_NaN();
};

inline double boundedMicroWaterContentCeiling(
    PotentialExchangeLocalSolveContext const& local_context,
    double const n_l_floor)
{
    constexpr double porosity_upper = 1.0 - 1e-12;
    auto const compute_total_porosity_bound = [&]()
    {
        double const phi_prev_sum = std::clamp(
            std::max(0.0, local_context.phi_M_prev) +
                std::max(0.0, local_context.phi_m_prev),
            0.0, porosity_upper);
        if (std::isfinite(local_context.phi))
        {
            return std::clamp(std::max(0.0, local_context.phi), 0.0,
                              porosity_upper);
        }

        double const delta_eps_v =
            local_context.volumetric_strain - local_context.volumetric_strain_prev;
        double const denominator = 1.0 + delta_eps_v;
        if (std::isfinite(denominator) && std::abs(denominator) > 1e-12)
        {
            double const phi_from_kinematics =
                (phi_prev_sum + delta_eps_v) / denominator;
            if (std::isfinite(phi_from_kinematics))
            {
                return std::clamp(phi_from_kinematics, 0.0, porosity_upper);
            }
        }

        return phi_prev_sum;
    };

    return std::max(n_l_floor, compute_total_porosity_bound());
}

inline TransportPorosityUpdateData computeTransportPorosityUpdate(
    double const phi, double const phi_M_prev, double const phi_m_prev,
    double const n_l, double const volumetric_strain,
    double const volumetric_strain_prev,
    MacroPorosityUpdateMode const macro_porosity_update_mode)
{
    constexpr double porosity_upper = 1.0 - 1e-12;
    double const phi_prev_sum = std::clamp(
        std::max(0.0, phi_M_prev) + std::max(0.0, phi_m_prev), 0.0,
        porosity_upper);
    double const delta_eps_v = volumetric_strain - volumetric_strain_prev;
    double const denominator = 1.0 + delta_eps_v;
    double phi_safe = phi_prev_sum;
    if (std::isfinite(phi))
    {
        phi_safe = std::clamp(std::max(0.0, phi), 0.0, porosity_upper);
    }
    else if (std::isfinite(denominator) && std::abs(denominator) > 1e-12)
    {
        double const phi_from_kinematics =
            (phi_prev_sum + delta_eps_v) / denominator;
        if (std::isfinite(phi_from_kinematics))
        {
            phi_safe = std::clamp(phi_from_kinematics, 0.0, porosity_upper);
        }
    }
    double const phi_M_prev_safe = std::min(std::max(0.0, phi_M_prev), phi_safe);
    double const phi_m_prev_safe =
        std::min(std::max(0.0, phi_m_prev), std::max(0.0, phi_safe - phi_M_prev_safe));

    // Hierarchical split:
    //   phi = phi_M + (1 - phi_M) * n_l
    //   phi_M = (phi - n_l) / (1 - n_l)
    //   phi_m = (1 - phi_M) * n_l
    //
    // Keep the legacy "additive_macro_porosity_rate_mode" keyword as a config
    // alias, but use the hierarchical split law here as requested.
    if (macro_porosity_update_mode ==
        MacroPorosityUpdateMode::ReferenceAdditiveRate)
    {
        static std::once_flag once;
        std::call_once(once, []
        {
            INFO(
                "DSM: macro_porosity_update_mode='additive_macro_porosity_rate_mode' now evaluates the hierarchical porosity split.");
        });
    }

    double const n_l_safe = std::clamp(std::max(0.0, n_l), 0.0, phi_safe);
    double const one_minus_n_l = std::max(1e-12, 1.0 - n_l_safe);
    double const phi_M_candidate = (phi_safe - n_l_safe) / one_minus_n_l;
    double const phi_M = std::clamp(phi_M_candidate, 0.0, phi_safe);
    double const phi_m = std::clamp(
        (1.0 - phi_M) * n_l_safe, 0.0, std::max(0.0, phi_safe - phi_M));

    return {
        .phi_M = phi_M,
        .phi_M_prev = phi_M_prev_safe,
        .phi_m = phi_m,
        .phi_m_prev = phi_m_prev_safe,
    };
}

// DIAGNOSTIC (mass-strip A/B test, 2026-09-30; Vinay R-03; NOT adopted).
// True when the stored micro water content n_l sits on its ceiling n_l = phi,
// i.e. on the clamp of boundedMicroWaterContentCeiling. The micro solve stores
// n_l = clamp(., floor, phi_clamped) exactly, and the assembly-side phi is
// phi_M + phi_m recomputed from the same n_l, so the two agree to a few
// round-off units: the detection tolerance is 1e3 machine epsilons
// (2.2e-13, dimensionless), derived from double precision, not a physical value.
inline bool microWaterContentIsAtCeiling(double const n_l, double const phi)
{
    constexpr double porosity_upper = 1.0 - 1e-12;
    if (!std::isfinite(phi) || !std::isfinite(n_l))
    {
        return false;
    }
    double const ceiling =
        std::clamp(std::max(0.0, phi), 0.0, porosity_upper);  // [-]
    constexpr double tol =
        1e3 * std::numeric_limits<double>::epsilon();  // [-]
    return n_l >= ceiling - tol;
}

// [2026-05-26 PHYSICS FIX] The returned quantity is the aggregate SOLID
// fraction (V_solid/V_aggregate = 1 - n_l) used as the denominator of
// the gravimetric water content omega_l = n_l * rho_lR / (nS * rho_SR).
// Earlier this returned the aggregate VOLUME fraction in REV
// (1 - phi_M = (1-phi0)/(1-n_l)), which produces a non-standard
// omega_l that deviates from the dry-solid-mass-referenced gravimetric
// content by factor (1-n_l)^2/(1-phi0) (state-dependent: up to +80%
// at low n_l, down to -50% near saturation). See the OPEN section in
// agents_dsm_mfront_hierarchical.md (commit fc21a3dd1d) for the full
// algebra and numerical verification. The same fix is applied to the
// mfront bridge in RichardsMechanicsDSMMicroMacroBridge.mfront.
inline double computeActiveMicroSolidVolumeFraction(
    double const n_l, PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    if (potential_exchange_params.micro_solid_volume_fraction_mode ==
        MicroSolidVolumeFractionMode::Reference)
    {
        return std::max(1e-16, potential_exchange_params.micro_solid_volume_fraction_reference);
    }

    // Aggregate solid fraction = 1 - n_l (with clamps).
    double const n_l_safe = std::clamp(std::max(0.0, n_l), 0.0, 1.0 - 1e-12);
    return std::max(1e-16, 1.0 - n_l_safe);
}

inline double computePreviousMicroSolidVolumeFraction(
    double const n_l_prev,
    PotentialExchangeLocalSolveContext const& /*local_context*/,
    PotentialExchangeParameters const& potential_exchange_params)
{
    if (potential_exchange_params.micro_solid_volume_fraction_mode ==
        MicroSolidVolumeFractionMode::Reference)
    {
        return std::max(1e-16, potential_exchange_params.micro_solid_volume_fraction_reference);
    }

    // [2026-05-26 PHYSICS FIX] Previously returned 1 - phi_M_prev. Now
    // returns the previous aggregate solid fraction 1 - n_l_prev to match
    // the corrected active definition.
    double const n_l_prev_safe = std::clamp(std::max(0.0, n_l_prev), 0.0, 1.0 - 1e-12);
    return std::max(1e-16, 1.0 - n_l_prev_safe);
}

struct ReducedMicroLiquidDensityData
{
    double rho_lR = 0.0;
    double omega_l = 0.0;
    double drho_lR_dnl = 0.0;
    double drho_l_dn_l = 0.0;
    // KKT micro-water ceiling (DESIGN.md 3.1(a), C2): d rho_lR_EOS/d rho_LR at
    // fixed n_l = 1/dg_drho (implicit function of g = rho_lR - rho_LR -
    // rho_l0*exp(-a*omega(rho_lR)), dg/drho_LR = -1). = 1 to round-off when
    // micro_liquid_density_a = 1e-16. 0 where dg_drho is degenerate (same
    // fallback convention as drho_lR_dnl). Appended at the end: no existing
    // value changes.
    double drho_lR_drho_LR = 0.0;
};

inline ReducedMicroLiquidDensityData computeReducedMicroLiquidDensity(
    double const n_l, double const rho_LR, double const nS,
    PotentialExchangeParameters const& potential_exchange_params)
{
    double const n_l_safe = std::max(1e-16, n_l);
    double const nS_safe = std::max(1e-16, nS);
    double const rho_SR = std::max(1e-16, potential_exchange_params.micro_solid_density_reference);
    double const rho_l0 = std::max(1e-16, potential_exchange_params.micro_liquid_density_reference);
    double const a_rho = std::max(1e-16, potential_exchange_params.micro_liquid_density_a);
    double const b_rho = std::max(1e-16, potential_exchange_params.micro_liquid_density_b);
    double const denominator = nS_safe * rho_SR;

    auto const eval_rhs = [&](double const rho_lR)
    {
        double const omega_l =
            std::max(1e-16, n_l_safe * rho_lR / denominator);
        double const exp_term =
            std::exp(-a_rho * std::pow(omega_l, b_rho));
        return std::pair{omega_l, rho_l0 * exp_term + rho_LR};
    };

    double rho_lR = rho_LR +
                    rho_l0 *
                        std::exp(-a_rho *
                                 std::pow(std::max(1e-16, n_l_safe * rho_LR /
                                                              denominator),
                                          b_rho));
    constexpr int max_iterations = 30;
    constexpr double tolerance = 1e-14;
    bool converged = false;

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        auto const [omega_l, rhs] = eval_rhs(rho_lR);
        double const residual = rho_lR - rhs;
        if (std::abs(residual) <=
            tolerance * std::max(1.0, std::abs(rho_lR)))
        {
            converged = true;
            break;
        }

        double const common =
            (rhs - rho_LR) * a_rho * b_rho *
            std::pow(omega_l, b_rho - 1.0);
        double const jacobian =
            1.0 + common * (n_l_safe / denominator);
        if (!(std::isfinite(jacobian) && std::abs(jacobian) > 1e-20))
        {
            break;
        }

        double const rho_candidate =
            std::max(1e-16, rho_lR - residual / jacobian);
        if (std::abs(rho_candidate - rho_lR) <=
            tolerance * std::max(1.0, std::abs(rho_lR)))
        {
            rho_lR = rho_candidate;
            converged = true;
            break;
        }
        rho_lR = rho_candidate;
    }

    if (!converged)
    {
        static std::once_flag once;
        std::call_once(once, []
        {
            WARN(
                "DSM: reduced microscale liquid-density EOS did not converge at least once; using the last Newton iterate.");
        });
    }

    auto const [omega_l, rhs] = eval_rhs(rho_lR);
    (void)rhs;
    double const common =
        (rho_lR - rho_LR) * a_rho * b_rho *
        std::pow(omega_l, b_rho - 1.0);
    double const dg_drho =
        1.0 + common * (n_l_safe / denominator);
    // dg/dn_l with nS FROZEN: common * domega/dn_l|_{nS} = common*(rho_lR/denom).
    double dg_dn = common * (rho_lR / denominator);
    // F2 (2026-06-06, tangent-only): under current_porosity_split nS = 1 - n_l is
    // LIVE, so the EOS slaved drho_lR/dn_l picks up the domega/dnS*(dnS/dnl)
    // chain. domega/dnS = -omega/nS, dnS/dnl = -1 -> chain = +omega/nS, added to
    // domega/dn_l, i.e. dg_dn += common*(omega_l/nS). Reference mode: nS constant
    // -> NO change (exact). This makes the returned drho_lR_dnl the SAME slaved
    // derivative the forward 2x2 FD solve sees (active_nS recomputed per n_l),
    // which F1 relies on. Tangent-only: the forward EOS Newton uses an analytic
    // 1D self-Jacobian dg_drho (unchanged) and converges to the same rho_lR.
    if (potential_exchange_params.micro_solid_volume_fraction_mode ==
        MicroSolidVolumeFractionMode::CurrentPorositySplit)
    {
        dg_dn += common * (omega_l / nS_safe);
    }
    double const drho_lR_dnl =
        (std::isfinite(dg_drho) && std::abs(dg_drho) > 1e-20)
            ? -dg_dn / dg_drho
            : 0.0;

    double const drho_lR_drho_LR =
        (std::isfinite(dg_drho) && std::abs(dg_drho) > 1e-20)
            ? 1.0 / dg_drho
            : 0.0;  // [-], KKT active p-p tangent (DESIGN.md 3.2)

    return {
        .rho_lR = rho_lR,
        .omega_l = omega_l,
        .drho_lR_dnl = drho_lR_dnl,
        .drho_l_dn_l = rho_lR + n_l_safe * drho_lR_dnl,
        .drho_lR_drho_LR = drho_lR_drho_LR,
    };
}

inline ReducedMicroLiquidDensityData computeActiveMicroLiquidDensity(
    double const n_l, double const rho_LR,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    double const active_nS =
        computeActiveMicroSolidVolumeFraction(n_l, local_context, potential_exchange_params);
    return computeReducedMicroLiquidDensity(n_l, rho_LR, active_nS, potential_exchange_params);
}

// ── Film-pressure folding (maxwell sec.5), shared by every local micro solve ──
// Given a BARE van-der-Waals micro potential `out` (already evaluated at this
// n_l with the SAME rho_lR_used the vdW formula consumed), ADD the smoothly-
// gated film delta mu_lR -> -(Pi - b*p_conf)/rho_lR in place. Strictly gated on
// the film_pressure_coupling master flag AND a finite confining_pressure_p_conf
// sentinel, so flag OFF or stress-not-supplied (default-constructed context,
// e.g. the eigenstress-difference driver) leaves `out` BIT-FOR-BIT unchanged.
//
// This is the ONE evaluator factored out of computeActiveMicroPotential so the
// scalar/microstate local solve, the macro exchange assembly, AND the
// mass-storage 2x2 local solve (which builds rho_lR itself and so cannot route
// through computeActiveMicroPotential's internal density) all fold the IDENTICAL
// film term — equipresence across local-solve modes (increment E, 2026-06-06).
// ── Macro-porosity floor as a SMOOTH film-to-bulk cutoff (Vinay 2026-06-06) ──
// Fades the disjoining micro potential to bulk (mu_lR -> 0) as interlayer water
// n_l approaches n_l_cap = (phi - floor)/(1 - floor), over width w, so the
// exchange equilibrates at n_l ~ n_l_cap (phi_M ~ floor) WITHOUT a hard clamp.
// Gate g = t*(2 - t), t = (n_l_cap - n_l)/w in [0,1]: C1 at onset (t=1) but a
// NONZERO slope at the cutoff edge (t->0, dg/dt->2), so d mu_lR/d n_l stays
// nonzero where the saturated equilibrium sits -> the pressure-block diagonal
// stays conditioned (a smoothstep, dg/dt=0 at t=0, would re-singularise it).
// floor == 0 -> unchanged (bit-for-bit).
inline void applyMacroFloorCutoff(
    VanDerWaalsMicroPotentialData& out, double const n_l,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& params)
{
    double const floor = params.macro_porosity_floor;
    if (!(floor > 0.0 && std::isfinite(local_context.phi)))
    {
        return;
    }
    double const phi = std::clamp(local_context.phi, 0.0, 1.0 - 1e-12);
    if (!(phi > floor))
    {
        return;
    }
    double const n_l_cap = (phi - floor) / std::max(1e-12, 1.0 - floor);
    double const w = (params.macro_floor_cutoff_width > 0.0)
                         ? params.macro_floor_cutoff_width
                         : std::max(1e-8, 0.05 * n_l_cap);
    double const t = (n_l_cap - n_l) / w;
    if (t >= 1.0)
    {
        return;  // n_l <= n_l_cap - w: full disjoining, unchanged
    }
    double g_cut, dg_dt;
    if (t <= 0.0)
    {
        g_cut = 0.0;  // n_l >= n_l_cap: bulk
        dg_dt = 0.0;
    }
    else
    {
        g_cut = t * (2.0 - t);   // C1 at t=1; dg/dt=2 at t=0 (nonzero edge slope)
        dg_dt = 2.0 * (1.0 - t);
    }
    double const dg_dnl = dg_dt * (-1.0 / w);
    double const mu0 = out.mu_lR;
    out.mu_lR = g_cut * mu0;                                 // J/kg
    out.dmu_lR_dnl = g_cut * out.dmu_lR_dnl + mu0 * dg_dnl;  // J/kg
    out.dmu_lR_drho_lR *= g_cut;
    out.dmu_lR_dnS *= g_cut;
    out.dmu_lR_drho_SR *= g_cut;
}

inline void applyFilmPressureMicroPotential(
    VanDerWaalsMicroPotentialData& out, double const n_l,
    double const rho_lR_used, double const active_nS,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    applyMacroFloorCutoff(out, n_l, local_context, potential_exchange_params);
    if (!(potential_exchange_params.film_pressure_coupling &&
          std::isfinite(local_context.confining_pressure_p_conf)))
    {
        return;  // flag OFF or NaN p_conf -> unchanged (pure vdW), bit-for-bit.
    }

    // ── EXACT energy route (film_energy_route = exact; PI_OF_NL_EV §4.4) ───
    // REPLACES the operational mu assembly below: the bare adsorption part
    // stays evaluated at the TRUE n_l (already in `out`, cutoff-folded above —
    // mirror of the shipped Off-path structure), and the strain coupling is
    // the one-Psi partner mu_mech from computeStrainedFilmEnergyPair (closed-
    // form strain integrals; Maxwell-exact pair with the unchanged eigenstress
    // half). No +b*p_conf/rho_lR bolt-on. The macro-floor cutoff factor g is
    // recovered from post/pre bare values and product-ruled into the n_l
    // chain. nS chains FROZEN (B1), as in the operational route.
    if (potential_exchange_params.film_strain_coupling !=
            FilmStrainCouplingMode::Off &&
        potential_exchange_params.film_energy_route == FilmEnergyRoute::Exact)
    {
        if (potential_exchange_params.film_strain_coupling !=
            FilmStrainCouplingMode::Kinematic)
        {
            OGS_FATAL(
                "film_energy_route = 'exact' requires film_strain_coupling = "
                "'kinematic' (create-time validated; "
                "DSM/PI_OF_NL_EV_IMPLEMENTATION.md §3).");
        }
        double const eps_v_ex = local_context.volumetric_strain;
        double const sign_ex =
            microPotentialSignFactorFromParameters(potential_exchange_params);
        double const kappa_ex =
            potential_exchange_params.film_strain_kappa ==
                    FilmStrainKappaMode::Aggregate
                ? active_nS
                : 1.0;
        auto const pair = computeStrainedFilmEnergyPair(
            n_l, eps_v_ex, kappa_ex, local_context.biot_coefficient,
            local_context.drained_bulk_modulus, true /*include_S, route R3*/,
            rho_lR_used, active_nS,
            potential_exchange_params.micro_solid_density_reference,
            potential_exchange_params.hamaker_constant,
            potential_exchange_params.specific_surface, sign_ex,
            // Live K(rho_d): rho_d = rho_SR*(1-phi); off / phi sentinel ->
            // parse scalar (K_OF_RHO_D_LIVE.md). H2 fix (2026-06-14): use the
            // SAME effective K the bare out.mu_lR was built with (:1364), so
            // g_cut = out.mu_lR/pair.mu_bare_pre stays the macro-floor cutoff
            // factor under live K (scalar-mode bit-for-bit; matches :768/:2084).
            effectiveAugmentationPrefactor(potential_exchange_params,
                                           local_context.phi),  // K [J/kg]
            potential_exchange_params.potential_augmentation_exponent,
            0.0 /*dnS_dnl: frozen nS (B1)*/,
            potential_exchange_params.micro_water_content_floor);

        // Macro-floor cutoff factor g = mu_post/mu_pre (g == 1 when the
        // cutoff is inactive; g == 0 at full bulk -> film physics off).
        // |mu_bare_pre| > 0 is guaranteed by the bare law's OGS_FATALs.
        double const g_cut = out.mu_lR / pair.mu_bare_pre;            // [-]
        // L1 fix (2026-06-14): pair.dmu_bare_dnl_pre was computed with frozen
        // nS (dnS_dnl = 0, L738), but out.dmu_lR_dnl carries the caller's F2
        // chain (dnS_dnl = -1 under CurrentPorositySplit; see computeActive-
        // MicroPotential L1353). Mixing the two legs corrupts dg_dnl under
        // CurrentPorositySplit. Recompute the bare-pre derivative with the
        // caller's dnS_dnl so both legs of dg_dnl use the SAME nS chain.
        // mu_bare_pre (the value) is dnS-independent, so g_cut is unaffected;
        // only the derivative changes. Reference mode -> dnS_dnl = 0 ->
        // bit-for-bit identical.
        double const dnS_dnl_caller =
            potential_exchange_params.micro_solid_volume_fraction_mode ==
                    MicroSolidVolumeFractionMode::CurrentPorositySplit
                ? -1.0
                : 0.0;
        double const dmu_bare_dnl_pre_caller =
            computeVanDerWaalsMicroPotential(
                n_l, rho_lR_used, active_nS,
                potential_exchange_params.micro_solid_density_reference,
                potential_exchange_params.hamaker_constant,
                potential_exchange_params.specific_surface, sign_ex,
                effectiveAugmentationPrefactor(potential_exchange_params,
                                               local_context.phi),  // K [J/kg]
                potential_exchange_params.potential_augmentation_exponent,
                dnS_dnl_caller,
                potential_exchange_params.micro_water_content_floor)
                .dmu_lR_dnl;  // J/kg per n_l (caller's nS chain)
        double const dg_dnl =
            (out.dmu_lR_dnl - g_cut * dmu_bare_dnl_pre_caller) /
            pair.mu_bare_pre;  // [1 per n_l]

        out.mu_lR += g_cut * pair.mu_mech;  // J/kg, additive (never =)
        out.dmu_lR_dnl +=
            g_cut * pair.dmu_mech_dnl + pair.mu_mech * dg_dnl;  // J/kg per n_l
        out.dmu_lR_drho_lR += g_cut * pair.dmu_mech_drho_lR;  // (J/kg)/(kg/m^3)
        return;
    }

    // ── Strained-film modes (DSM/STRAINED_FILM_IMPLEMENTATION.md) ───────────
    // film_strain_coupling != Off REPLACES the frozen-geometry path below: the
    // bare law is evaluated at the strained film state w_eff and mu_lR carries
    // the Derjaguin load term +b*p_conf/rho_lR (squeezing confined liquid
    // raises its chemical potential). The shipped integrable partner is NOT
    // added on top — it is the frozen-h, O(eps_v) truncation of the same
    // coupling (no double counting; D3 provisional, see the shipped-limit unit
    // test). nS chains are FROZEN here (B1): the strained re-evaluation uses
    // dnS_dnl = 0 regardless of the caller's F2 mode.
    if (potential_exchange_params.film_strain_coupling !=
        FilmStrainCouplingMode::Off)
    {
        double const p_conf_sf = local_context.confining_pressure_p_conf;
        double const eps_v_sf = local_context.volumetric_strain;
        double const sign_sf =
            microPotentialSignFactorFromParameters(potential_exchange_params);
        // Live K(rho_d) (K_OF_RHO_D_LIVE.md): rho_d = rho_SR*(1-phi) from the
        // context's total porosity; off-mode / phi sentinel -> parse scalar.
        double const K_aug_sf = effectiveAugmentationPrefactor(
            potential_exchange_params, local_context.phi);  // K [J/kg]

        auto const film_state = computeStrainedFilmState(
            potential_exchange_params.film_strain_coupling,
            potential_exchange_params.film_strain_kappa, n_l, active_nS,
            eps_v_sf, p_conf_sf, rho_lR_used,
            potential_exchange_params.micro_solid_density_reference,
            potential_exchange_params.hamaker_constant,
            potential_exchange_params.specific_surface, sign_sf,
            K_aug_sf,
            potential_exchange_params.potential_augmentation_exponent,
            potential_exchange_params.micro_water_content_floor,
            rho_lR_used /*rho_pi: mirrors Pi = -rho_lR_used*mu_lR below*/);

        // Bare law at the strained state (frozen nS, B1).
        auto strained = computeVanDerWaalsMicroPotential(
            film_state.w_eff, rho_lR_used, active_nS,
            potential_exchange_params.micro_solid_density_reference,
            potential_exchange_params.hamaker_constant,
            potential_exchange_params.specific_surface, sign_sf,
            K_aug_sf,
            potential_exchange_params.potential_augmentation_exponent,
            0.0 /*dnS_dnl: frozen nS in strained modes*/,
            potential_exchange_params.micro_water_content_floor);

        // Chain the law's n_l-derivative through w_eff BEFORE the cutoff
        // product rule (the cutoff factor g is a function of the TRUE n_l).
        strained.dmu_lR_dnl *= film_state.dw_eff_dnl;
        strained.d2mu_lR_dnl2 *=
            film_state.dw_eff_dnl * film_state.dw_eff_dnl;
        applyMacroFloorCutoff(strained, n_l, local_context,
                              potential_exchange_params);

        // Derjaguin load term, UNCUT (mirrors the Off path, where the
        // mechanical partner is added after the cutoff):
        //   mu_load = +b*p_conf/rho_lR  [J/kg]  — compression raises mu_lR
        //   (expulsion channel); reversible because the disjoining half
        //   Pi(w_eff) stiffens through the SAME strained state (one Psi).
        double const b_sf = local_context.biot_coefficient;
        double const mu_load = b_sf * p_conf_sf / rho_lR_used;  // J/kg

        out.mu_lR = strained.mu_lR + mu_load;               // J/kg
        out.dmu_lR_dnl = strained.dmu_lR_dnl;               // J/kg per n_l
        out.d2mu_lR_dnl2 = strained.d2mu_lR_dnl2;           // J/kg per n_l^2
        out.dmu_lR_dnS = strained.dmu_lR_dnS;
        out.dmu_lR_drho_SR = strained.dmu_lR_drho_SR;
        out.dmu_lR_drho_lR =
            strained.dmu_lR_drho_lR - mu_load / rho_lR_used;  // (J/kg)/(kg/m^3)
        return;
    }
    (void)active_nS;
    // ── INTEGRABLE Maxwell mechanical partner (spec item 2; REPLACES the old
    // non-integrable +g*b*p_conf/rho_lR film delta) ─────────────────────────
    // mu_lR_mech = -[ (Pi + n_l*Pi')*eps_v + 0.5*b*K_drained*eps_v^2 ]/rho_lR,
    // additive to mu_lR_vdw (NEVER overwrites it). Pi, Pi', Pi'' are the BARE
    // van-der-Waals disjoining pressure and its n_l-derivatives at THIS state,
    // density-mirrored to rho_lR_used exactly as the eigenstress half does
    // (Pi = -rho_lR_used*mu_lR_vdw). This is the ONE equipresent fold point:
    // every micro local solve and the macro exchange route mu_lR through here,
    // so they all see the SAME mu_lR(n_l, eps_v). PARAMETER-FREE: eps_v and b
    // from the context, K_drained from the elastic stiffness (context), Pi/Pi'/
    // Pi'' from the vdW + augmentation potential.
    double const Pi = -rho_lR_used * out.mu_lR;            // Pa (disjoining)
    double const dPi_dnl = -rho_lR_used * out.dmu_lR_dnl;  // Pa per n_l = Pi'
    double const d2Pi_dnl2 =
        -rho_lR_used * out.d2mu_lR_dnl2;  // Pa per n_l^2 = Pi''
    double const eps_v = local_context.volumetric_strain;
    // K_drained from the context (mechanical drained bulk modulus). NaN sentinel
    // (default-constructed context / GP eigenstress-difference driver) -> drop
    // the 0.5*b*K_drained*eps_v^2 quadratic so the partner reduces to the finite
    // linear conjugate -(Pi + n_l*Pi')*eps_v/rho_lR.
    double const K_drained =
        std::isfinite(local_context.drained_bulk_modulus)
            ? local_context.drained_bulk_modulus
            : 0.0;
    auto const mech = computeIntegrableMechanicalMicroPotential(
        Pi, dPi_dnl, d2Pi_dnl2, n_l, eps_v, local_context.biot_coefficient,
        K_drained, rho_lR_used);
    out.mu_lR += mech.mu_lR_mech;              // J/kg, additive
    out.dmu_lR_dnl += mech.dmu_lR_mech_dnl;    // J/kg per n_l
    out.dmu_lR_drho_lR += mech.dmu_lR_mech_drho_lR;  // (J/kg)/(kg/m^3)
}

inline ReducedMicroLiquidDensityData computePreviousMicroLiquidDensity(
    double const n_l_prev, double const rho_LR,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    double const previous_nS =
        computePreviousMicroSolidVolumeFraction(n_l_prev, local_context,
                                                potential_exchange_params);
    return computeReducedMicroLiquidDensity(n_l_prev, rho_LR, previous_nS,
                                              potential_exchange_params);
}

struct MicroMacroMassStorageCoupledSolveData
{
    double n_l = 0.0;
    double rho_lR = 0.0;
    double phi_m = 0.0;
    double phi_M = 0.0;
    double p_L_m = 0.0;
    double S_L_m = 0.0;
    VanDerWaalsMicroPotentialData micro_potential;
    PotentialDrivenMassExchangeData exchange;
    bool converged = true;
};

inline MicroMacroMassStorageCoupledSolveData
solveReferenceMassStoragePredictorState(
    double const n_l_prev, double const rho_l_prev, double const rho_lR_prev,
    double const dt, double const rho_LR, double const alpha_bar,
    double const mu, YoungLaplaceMacroPotentialData const& macro_potential,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    requirePositiveViscosity("solveReferenceMassStoragePredictorState", mu);
    constexpr double n_l_floor = 1e-16;
    constexpr double rho_floor = 1e-16;
    // V2 (F3): +1 shipped sign of the volume-change term, -1 Eulerian
    // (micro_mass_strain_term_eulerian); +1 -> bit-identical.
    double const strain_sign =
        microMassStrainTermSign(potential_exchange_params);  // [-]
    double const dt_safe = std::isfinite(dt) && dt > 0.0 ? dt : 0.0;
    double const alpha_M_effective = alpha_bar * rho_LR / mu;
    double const volumetric_strain_rate =
        dt_safe > 0.0
            ? (local_context.volumetric_strain -
               local_context.volumetric_strain_prev) /
                  dt_safe
            : 0.0;
    double const n_l_ceiling =
        boundedMicroWaterContentCeiling(local_context, n_l_floor);

    auto evaluate = [&](double const n_l)
    {
        double const active_nS = computeActiveMicroSolidVolumeFraction(
            n_l, local_context, potential_exchange_params);
        auto const micro_liquid_density = computeReducedMicroLiquidDensity(
            n_l, rho_LR, active_nS, potential_exchange_params);
        auto micro_potential = computeVanDerWaalsMicroPotential(
            n_l, micro_liquid_density.rho_lR, active_nS,
            potential_exchange_params.micro_solid_density_reference, potential_exchange_params.hamaker_constant,
            potential_exchange_params.specific_surface,
            microPotentialSignFactorFromParameters(potential_exchange_params),
            // Live K(rho_d): rho_d = rho_SR*(1-phi); off / phi sentinel ->
            // parse scalar (K_OF_RHO_D_LIVE.md).
            effectiveAugmentationPrefactor(potential_exchange_params,
                                           local_context.phi),  // K [J/kg]
            potential_exchange_params.potential_augmentation_exponent,
            0.0 /*dnS_dnl*/,
            potential_exchange_params.micro_water_content_floor);
        // Increment E: fold the film delta into mu_lR (and dmu_lR_dnl, which the
        // analytic predictor Jacobian below reads through micro_potential) using
        // the SAME micro density the vdW formula consumed. No-op when the flag is
        // OFF / p_conf NaN -> predictor is bit-for-bit the bare-vdW path.
        applyFilmPressureMicroPotential(micro_potential, n_l,
                                        micro_liquid_density.rho_lR, active_nS,
                                        local_context, potential_exchange_params);
        double const mu_LR_active = macro_potential.mu_LR;
        double const mu_lR_active = micro_potential.mu_lR;
        auto const exchange = computePotentialDrivenMassExchange(
            alpha_M_effective, mu_LR_active, mu_lR_active);
        // REV-scale liquid apparent density: rho_l = phi_m * rho_lR
        // Hierarchical split: phi_m = (1-phi)/(1-n_l)*n_l.
        // Previously this was n_l*rho_lR (aggregate scale — missing (1-phi_M)).
        double const phi_h = std::isfinite(local_context.phi)
            ? std::clamp(local_context.phi, 0.0, 1.0 - 1e-12)
            : std::clamp(local_context.phi_M_prev + local_context.phi_m_prev,
                         0.0, 1.0 - 1e-12);
        double const one_minus_n_l_h = std::max(1e-12, 1.0 - n_l);
        double const rho_l =
            (1.0 - phi_h) / one_minus_n_l_h * n_l * micro_liquid_density.rho_lR;
        double const residual = rho_l - rho_l_prev -
                                dt_safe * exchange.rho_l_hat -
                                strain_sign * dt_safe * rho_l *
                                    volumetric_strain_rate;  // kg/m^3
        return std::tuple{residual, micro_potential, exchange,
                          micro_liquid_density};
    };

    MicroMacroMassStorageCoupledSolveData out;
    if (dt_safe <= 0.0)
    {
        out.n_l = std::clamp(n_l_prev, n_l_floor, n_l_ceiling);
        out.rho_lR = std::max(rho_floor, rho_lR_prev);
        auto const [residual, micro_potential, exchange, micro_density] =
            evaluate(out.n_l);
        (void)residual;
        (void)micro_density;
        out.micro_potential = micro_potential;
        out.exchange = exchange;
        return out;
    }

    double n_l = std::clamp(n_l_prev, n_l_floor, n_l_ceiling);
    constexpr int max_iterations = 40;
    constexpr double residual_tolerance = 1e-14;
    constexpr double increment_tolerance = 1e-14;

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        auto const [residual, micro_potential, exchange, micro_density] =
            evaluate(n_l);
        if (std::abs(residual) <=
            residual_tolerance * std::max(1.0, std::abs(rho_l_prev)))
        {
            out.n_l = n_l;
            out.rho_lR = micro_density.rho_lR;
            out.micro_potential = micro_potential;
            out.exchange = exchange;
            return out;
        }

        double const drho_l_hat_dn_l =
            exchange.drho_l_hat_dmu_lR * micro_potential.dmu_lR_dnl;
        // d(rho_l_REV)/dn_l where rho_l_REV = (1-phi)/(1-n_l)*n_l*rho_lR:
        //   = (1-phi)/(1-n_l)^2 * rho_lR  +  (1-phi)/(1-n_l) * drho_lR_dnl
        //   = (1-phi_M)/one_minus_n_l * rho_lR  +  (1-phi_M) * drho_lR_dnl
        double const phi_jac = std::isfinite(local_context.phi)
            ? std::clamp(local_context.phi, 0.0, 1.0 - 1e-12)
            : std::clamp(local_context.phi_M_prev + local_context.phi_m_prev,
                         0.0, 1.0 - 1e-12);
        double const one_minus_n_l_jac = std::max(1e-12, 1.0 - n_l);
        double const one_minus_phi_M_jac = (1.0 - phi_jac) / one_minus_n_l_jac;
        double const drho_l_REV_dn_l =
            one_minus_phi_M_jac / one_minus_n_l_jac * micro_density.rho_lR +
            one_minus_phi_M_jac * micro_density.drho_lR_dnl;
        double const jacobian = drho_l_REV_dn_l -
                                dt_safe * drho_l_hat_dn_l -
                                strain_sign * dt_safe * drho_l_REV_dn_l *
                                    volumetric_strain_rate;  // kg/m^3
        if (!(std::isfinite(jacobian) && std::abs(jacobian) > 1e-20))
        {
            break;
        }

        double delta_n_l = -residual / jacobian;
        double n_l_candidate =
            std::clamp(n_l + delta_n_l, n_l_floor, n_l_ceiling);
        auto const [candidate_residual_initial, candidate_micro_potential,
                    candidate_exchange, candidate_micro_density] =
            evaluate(n_l_candidate);
        double candidate_residual = candidate_residual_initial;
        int backtracking_steps = 0;
        while (std::abs(candidate_residual) > std::abs(residual) &&
               backtracking_steps < 12)
        {
            delta_n_l *= 0.5;
            n_l_candidate =
                std::clamp(n_l + delta_n_l, n_l_floor, n_l_ceiling);
            auto const [retry_residual, retry_micro_potential,
                        retry_exchange, retry_micro_density] =
                evaluate(n_l_candidate);
            (void)retry_micro_potential;
            (void)retry_exchange;
            (void)retry_micro_density;
            candidate_residual = retry_residual;
            ++backtracking_steps;
        }

        if (std::abs(n_l_candidate - n_l) <=
            increment_tolerance * std::max(1.0, std::abs(n_l)))
        {
            out.n_l = n_l_candidate;
            auto const [final_residual, final_micro_potential, final_exchange,
                        final_micro_density] =
                evaluate(n_l_candidate);
            (void)final_residual;
            out.rho_lR = final_micro_density.rho_lR;
            out.micro_potential = final_micro_potential;
            out.exchange = final_exchange;
            out.converged = true;
            return out;
        }

        n_l = n_l_candidate;
        out.n_l = n_l;
        out.rho_lR = micro_density.rho_lR;
        out.micro_potential = micro_potential;
        out.exchange = exchange;
    }

    auto const [residual, micro_potential, exchange, micro_density] =
        evaluate(n_l);
    (void)residual;
    out.n_l = n_l;
    out.rho_lR = micro_density.rho_lR;
    out.micro_potential = micro_potential;
    out.exchange = exchange;
    out.converged = false;
    return out;
}

// Localized floating-point-contraction guard (Phase A, 2026-06-09).
// The micro 2x2 local solve forms FD central differences (r(x+h) - r(x-h)) that
// are cancellation-prone; under the build default -ffp-contract=fast clang is
// free to fuse the subtractions and divisions into FMAs and reassociate them.
// On the dd1800 conditioning cliff (near-singular assembled tangent) that
// reassociation perturbs the last bits enough to tip the global Newton path —
// the if/else refactor boundary alone changed which fusions clang chose and
// broke dd1800. Pinning FP_CONTRACT OFF for this translation unit's evaluation
// of these differences removes that fragility and restores parent-identical
// numerics. Portable STDC pragma (file scope; gcc honours this) + clang-specific
// statement-scope `#pragma clang fp contract(off)` inside the body (clang honours
// that form most reliably).
#pragma STDC FP_CONTRACT OFF
inline MicroMacroMassStorageCoupledSolveData
solveReferenceMassStorageCoupledState(
    double const n_l_prev, double const rho_l_prev, double const rho_lR_prev,
    double const dt, double const rho_LR, double const alpha_bar,
    double const mu, YoungLaplaceMacroPotentialData const& macro_potential,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
    requirePositiveViscosity("solveReferenceMassStorageCoupledState", mu);
    constexpr double n_l_floor = 1e-16;
    constexpr double rho_floor = 1e-16;
    // V2 (F3): see solveReferenceMassStoragePredictorState.
    double const strain_sign =
        microMassStrainTermSign(potential_exchange_params);  // [-]
    double const dt_safe = std::isfinite(dt) && dt > 0.0 ? dt : 0.0;
    double const alpha_M_effective = alpha_bar * rho_LR / mu;
    double const volumetric_strain_rate =
        dt_safe > 0.0
            ? (local_context.volumetric_strain -
               local_context.volumetric_strain_prev) /
                  dt_safe
            : 0.0;
    double const n_l_ceiling =
        boundedMicroWaterContentCeiling(local_context, n_l_floor);

    auto evaluate = [&](double const n_l, double const rho_lR)
    {
        double const active_nS = computeActiveMicroSolidVolumeFraction(
            n_l, local_context, potential_exchange_params);
        auto micro_potential = computeVanDerWaalsMicroPotential(
            n_l, rho_lR, active_nS, potential_exchange_params.micro_solid_density_reference,
            potential_exchange_params.hamaker_constant, potential_exchange_params.specific_surface,
            microPotentialSignFactorFromParameters(potential_exchange_params),
            // Live K(rho_d): rho_d = rho_SR*(1-phi); off / phi sentinel ->
            // parse scalar (K_OF_RHO_D_LIVE.md).
            effectiveAugmentationPrefactor(potential_exchange_params,
                                           local_context.phi),  // K [J/kg]
            potential_exchange_params.potential_augmentation_exponent,
            0.0 /*dnS_dnl*/,
            potential_exchange_params.micro_water_content_floor);
        // Increment E: fold the film delta into mu_lR with the in-iteration micro
        // density rho_lR (the 2x2 unknown the vdW formula just consumed). The 2x2
        // local Jacobian is finite-difference over THIS lambda, so it picks up the
        // film term in both residuals automatically. No-op when the flag is OFF /
        // p_conf NaN -> bit-for-bit the bare-vdW path.
        applyFilmPressureMicroPotential(micro_potential, n_l, rho_lR, active_nS,
                                        local_context, potential_exchange_params);
        double const mu_LR_active = macro_potential.mu_LR;
        double const mu_lR_active = micro_potential.mu_lR;
        auto const exchange = computePotentialDrivenMassExchange(
            alpha_M_effective, mu_LR_active, mu_lR_active);
        // REV-scale liquid apparent density: phi_m * rho_lR (hierarchical split).
        // phi_m = (1-phi)/(1-n_l)*n_l. Previously n_l*rho_lR (aggregate scale).
        double const phi_cs = std::isfinite(local_context.phi)
            ? std::clamp(local_context.phi, 0.0, 1.0 - 1e-12)
            : std::clamp(local_context.phi_M_prev + local_context.phi_m_prev,
                         0.0, 1.0 - 1e-12);
        double const one_minus_n_l_cs = std::max(1e-12, 1.0 - n_l);
        double const rho_l = (1.0 - phi_cs) / one_minus_n_l_cs * n_l * rho_lR;
        double const mass_residual = rho_l - rho_l_prev -
                                     dt_safe * exchange.rho_l_hat -
                                     strain_sign * dt_safe * rho_l *
                                         volumetric_strain_rate;  // kg/m^3
        auto const density = computeReducedMicroLiquidDensity(
            n_l, rho_LR, active_nS, potential_exchange_params);
        double const density_residual = rho_lR - density.rho_lR;
        return std::tuple{mass_residual, density_residual, micro_potential,
                          exchange};
    };

    // ── Analytic micro 2x2 Jacobian (replaces the 4 FD evaluate() calls) ──────
    // J = d(mass_residual, density_residual)/d(n_l, rho_lR), assembled from the
    // closed-form residual structure. Mirrors `evaluate` term-for-term so the
    // converged (n_l, rho_lR) — and hence the global solution — is byte-for-byte
    // the FD-path result (tangent-only: residuals are untouched). The micro
    // potential's mu_lR-derivatives (which already fold the macro-floor cutoff,
    // the film-pressure coupling and the augmentation through the shared helper
    // chain) are reused; the ONLY piece `evaluate` does not carry is the live-nS
    // chain in d mu_lR/d n_l (it passes dnS_dnl=0), so the analytic evaluator
    // re-runs the vdW helper with the correct dnS_dnl and re-applies the film
    // coupling, reproducing the TOTAL n_l-derivative the FD of `evaluate` sees.
    //
    // Closed forms (phi fixed in the local solve; one_minus_n_l = max(1e-12,
    // 1-n_l) matching `evaluate`):
    //   rho_l            = (1-phi)/(1-n_l) * n_l * rho_lR                [kg/m^3]
    //   d rho_l/d rho_lR = (1-phi)/(1-n_l) * n_l        (= rho_l/rho_lR) [-]
    //   d rho_l/d n_l    = (1-phi) * rho_lR / (1-n_l)^2                  [kg/m^3]
    //   d(mass_res)/d n_l    = d rho_l/d n_l*(1 - dt*eps_dot)
    //                          + dt*alpha_M*d mu_lR/d n_l   (rho_l_hat =
    //                          alpha_M*(mu_LR-mu_lR), d rho_l_hat/d n_l =
    //                          -alpha_M*d mu_lR/d n_l; sign flips with the -dt*).
    //   d(mass_res)/d rho_lR = d rho_l/d rho_lR*(1 - dt*eps_dot)
    //                          + dt*alpha_M*d mu_lR/d rho_lR
    //   d(dens_res)/d n_l    = -d rho_lR_EOS/d n_l   (EOS drho_lR_dnl; carries the
    //                          live-nS chain in CurrentPorositySplit, 0 otherwise)
    //   d(dens_res)/d rho_lR = 1   (rho_lR_EOS depends on n_l & rho_LR only, not
    //                          on the outer 2x2 unknown rho_lR -> exact).
    auto evaluate_analytic_jacobian = [&](double const n_l, double const rho_lR)
    {
        double const active_nS = computeActiveMicroSolidVolumeFraction(
            n_l, local_context, potential_exchange_params);
        // Live-nS chain: nS = 1 - n_l under CurrentPorositySplit (dnS_dnl = -1);
        // constant in Reference mode (dnS_dnl = 0, exact — no change).
        double const dnS_dnl =
            potential_exchange_params.micro_solid_volume_fraction_mode ==
                    MicroSolidVolumeFractionMode::CurrentPorositySplit
                ? -1.0
                : 0.0;
        auto micro_potential = computeVanDerWaalsMicroPotential(
            n_l, rho_lR, active_nS,
            potential_exchange_params.micro_solid_density_reference,
            potential_exchange_params.hamaker_constant,
            potential_exchange_params.specific_surface,
            microPotentialSignFactorFromParameters(potential_exchange_params),
            // MERGE FIX 2026-08-11 (jac_parallel -> maxwell_conjugate): use the
            // SAME live K(rho_d) the residual `evaluate` uses. This lambda came
            // from dsm_maxwell_jac_parallel, which branched off d98f5f8324 —
            // i.e. BEFORE maxwell_conjugate converted every
            // computeVanDerWaalsMicroPotential call site to
            // effectiveAugmentationPrefactor. The three-way merge was textually
            // clean, so this one call site was left on the parse-time scalar:
            // with potential_augmentation_prefactor_live_dry_density=true the
            // analytic 2x2 would then be the derivative of a DIFFERENT potential
            // than the residual. That is tangent-only, but do NOT read it as
            // harmless: if the local 2x2 fails to converge in max_iterations
            // this routine silently returns the decoupled PREDICTOR state (see
            // the `return out.converged ? out : predictor;` below) and no caller
            // inspects `converged` — so a bad micro tangent CAN change the
            // computed state without any warning or time-step rejection.
            // No-op when live mode is off (the helper returns the parse scalar),
            // which is the default and the state of every MS33 suite PRJ.
            effectiveAugmentationPrefactor(potential_exchange_params,
                                           local_context.phi),  // K [J/kg]
            potential_exchange_params.potential_augmentation_exponent, dnS_dnl,
            potential_exchange_params.micro_water_content_floor);
        // Re-apply the film coupling so d mu_lR/d n_l and d mu_lR/d rho_lR pick up
        // the integrable Maxwell partner exactly as `evaluate` did.
        applyFilmPressureMicroPotential(micro_potential, n_l, rho_lR, active_nS,
                                        local_context, potential_exchange_params);
        double const dmu_lR_dnl = micro_potential.dmu_lR_dnl;       // J/kg per n_l
        double const dmu_lR_drho_lR = micro_potential.dmu_lR_drho_lR;  // (J/kg)/(kg/m^3)

        double const phi_cs = std::isfinite(local_context.phi)
            ? std::clamp(local_context.phi, 0.0, 1.0 - 1e-12)
            : std::clamp(local_context.phi_M_prev + local_context.phi_m_prev,
                         0.0, 1.0 - 1e-12);
        double const one_minus_n_l_cs = std::max(1e-12, 1.0 - n_l);
        // d rho_l/d n_l = (1-phi)*rho_lR/(1-n_l)^2   [kg/m^3];
        // d rho_l/d rho_lR = (1-phi)/(1-n_l)*n_l     [-].
        double const drho_l_dnl = (1.0 - phi_cs) * rho_lR /
                                  (one_minus_n_l_cs * one_minus_n_l_cs);  // kg/m^3
        double const drho_l_drho_lR =
            (1.0 - phi_cs) / one_minus_n_l_cs * n_l;  // [-]

        // d rho_l_hat/d n_l = -alpha_M*d mu_lR/d n_l;  d/d rho_lR likewise.
        // mass_residual = rho_l - rho_l_prev - dt*rho_l_hat - s*dt*rho_l*eps_dot,
        // s = +1 shipped, s = -1 Eulerian (V2 F3 switch).
        double const J11 =
            drho_l_dnl * (1.0 - strain_sign * dt_safe * volumetric_strain_rate) +
            dt_safe * alpha_M_effective * dmu_lR_dnl;  // kg/m^3 per n_l
        double const J12 =
            drho_l_drho_lR * (1.0 - strain_sign * dt_safe * volumetric_strain_rate) +
            dt_safe * alpha_M_effective * dmu_lR_drho_lR;  // [-] (kg/m^3 per kg/m^3)

        // density_residual = rho_lR - rho_lR_EOS(n_l): d/d n_l = -drho_lR_dnl,
        // d/d rho_lR = 1.
        auto const density = computeReducedMicroLiquidDensity(
            n_l, rho_LR, active_nS, potential_exchange_params);
        double const J21 = -density.drho_lR_dnl;  // kg/m^3 per n_l
        double const J22 = 1.0;                    // [-]
        return std::array{J11, J12, J21, J22};
    };

    auto const predictor = solveReferenceMassStoragePredictorState(
        n_l_prev, rho_l_prev, rho_lR_prev, dt, rho_LR, alpha_bar, mu,
        macro_potential, local_context, potential_exchange_params);
    if (!predictor.converged)
    {
        return predictor;
    }

    MicroMacroMassStorageCoupledSolveData out = predictor;
    if (dt_safe <= 0.0)
    {
        return out;
    }

    double n_l = std::clamp(n_l_prev, n_l_floor, n_l_ceiling);
    double rho_lR = std::max(rho_floor, rho_lR_prev);
    constexpr int max_iterations = 60;
    constexpr double residual_tolerance = 1e-10;
    constexpr double increment_tolerance = 1e-10;

    // Phase D (2026-06-10): the analytic micro 2x2 local Jacobian is now the
    // DEFAULT. It was verified correct in Phase B — analytic == FD central
    // difference to round-off on every MS33 case, INCLUDING the dense /
    // EOS-active dd1800 (the earlier suspected J11/J12 error did not exist; the
    // dd1800 break was a global linear-solver conditioning fragility, not a
    // micro-tangent error — see DSM/AUDIT_maxwell_local_jacobian_2026-06-09.md).
    // REQUIREMENT: this path needs <scaling>false</scaling> in the Eigen
    // <linear_solver> block of the DSM PRJs. With Eigen's IterScaling (Ruiz)
    // enabled, the row/column equilibration overflows on the intrinsically
    // near-singular pressure block (cond ~5.8e22 from micro_liquid_density_a=
    // 1e-16 + Bishop saturation cutoff + tiny intrinsic-k); with scaling off and
    // SparseLU the system factorizes cleanly. The fp-contract(off) pragma and
    // the u-side analytic blocks (still OFF) are left as-is.
    // Deliberately decoupled from use_fd_jacobian_for_exchange (the block #3
    // macro p-p tangent flag), which keeps its own default so block #3 stays
    // analytic exactly as parent; setting use_fd_jacobian_for_exchange still
    // forces the FD micro path here.
    //
    // 2026-08-11, merge of dsm_maxwell_jac_parallel into
    // dsm_native_maxwell_conjugate — DEFAULT REVERTED true -> false by Vinay.
    // Reason (MEASURED on this tree/binary, not predicted): the analytic path
    // changes the ADAPTIVE TIME-STEP PATH on two of the six gating MS33 models,
    // so they no longer reproduce the reference VTUs Vinay approved 2026-06-23:
    //   Model III  405 -> 376 steps (2 rejected steps -> 0);  3/11 fields pass
    //   Model VII  682 -> 675 steps;                          5/11 fields pass
    //   dd1400 / dd1600 / dd1800 / Model IV: step path unchanged, 11/11 pass
    // The differences are time-discretisation scale, NOT a physics change
    // (Model III max: displacement 1.73e-6 m, sigma 1.86e4 Pa / 0.54% rel,
    // swelling_stress 5.2% rel); the TIER-A tolerances (1e-9 abs on
    // displacement) were calibrated for a bit-identical step path and cannot
    // survive a changed one. Isolated by a 2x2 experiment: the micro-Jacobian
    // choice alone drives the step path (III = 376 steps under BOTH <scaling>
    // settings), so this is not the linear-solver scaling flag.
    // With this flag false the merged tree reproduces all 6 references
    // step-for-step identically to the pre-merge binary (VERIFIED 6/6, 66/66
    // field comparisons). The Phase-D 2026-06-10 "land it as default" decision
    // is therefore parked, NOT withdrawn: re-enabling is this one line, and
    // requires re-baselining the Model III + VII reference VTUs first
    // (CLAUDE.md §3 / §12.5 — Vinay's call).
    // RE-ENABLED 2026-08-12 (Vinay; references re-baselined in the same commit).
    // The parking rationale above is the dated record of 2026-08-11. Measured for
    // this tree/toolchain before re-enabling: unit suite 1418/1418; dd1400/dd1600/
    // dd1800/ModelIV reproduce their committed references 11/11 fields each under
    // the analytic path (step paths unchanged: 308/311/308/637); Model III (405 ->
    // 376 steps, 2 rejected -> 0) and Model VII (682 -> 675) change step path only,
    // and their references are re-baselined (ratified) in this commit. The new
    // III/VII baselines are run-to-run AND cross-build bit-identical (exact 0.0
    // between two independent builds; the file-scope FP_CONTRACT OFF pragma is
    // what pins this). u-side blocks remain OFF (unchanged, known-unsafe).
    // NOTE the Phase-B "analytic == FD to round-off" claim above was measured
    // against the jac-branch residual; mc has since changed that residual
    // (live K(rho_d), strained film), so it is NOT re-verified for this tree.
    constexpr bool use_analytic_micro_jacobian = true;

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        auto const [mass_residual, density_residual, micro_potential, exchange] =
            evaluate(n_l, rho_lR);

        double const residual_norm =
            std::abs(mass_residual) / std::max(1.0, std::abs(rho_l_prev)) +
            std::abs(density_residual) / std::max(1.0, std::abs(rho_lR));
        if (residual_norm <= residual_tolerance)
        {
            out.n_l = n_l;
            out.rho_lR = rho_lR;
            out.micro_potential = micro_potential;
            out.exchange = exchange;
            out.converged = true;
            return out;
        }

        double J11 = 0.0, J12 = 0.0, J21 = 0.0, J22 = 0.0;
        // Default (Phase A): FD micro 2x2 — parent-identical. The analytic micro
        // 2x2 is taken ONLY when explicitly opted in via the parked
        // use_analytic_micro_jacobian constexpr above; the existing
        // use_fd_jacobian_for_exchange flag still forces FD when set.
        if (potential_exchange_params.use_fd_jacobian_for_exchange ||
            !use_analytic_micro_jacobian)
        {
            // FD path (default; also opt-in via use_fd_jacobian_for_exchange):
            // central difference of the full `evaluate` over (n_l, rho_lR).
            double const h_n = 1e-8 * std::max(1.0, std::abs(n_l));
            double const h_rho = 1e-8 * std::max(1.0, std::abs(rho_lR));
            auto const [r1_n_plus, r2_n_plus] =
                [&]() {
                    auto const [r1, r2, _, __] = evaluate(n_l + h_n, rho_lR);
                    (void)_;
                    (void)__;
                    return std::pair{r1, r2};
                }();
            auto const [r1_n_minus, r2_n_minus] =
                [&]() {
                    auto const [r1, r2, _, __] =
                        evaluate(std::max(n_l_floor, n_l - h_n), rho_lR);
                    (void)_;
                    (void)__;
                    return std::pair{r1, r2};
                }();
            auto const [r1_rho_plus, r2_rho_plus] =
                [&]() {
                    auto const [r1, r2, _, __] = evaluate(n_l, rho_lR + h_rho);
                    (void)_;
                    (void)__;
                    return std::pair{r1, r2};
                }();
            auto const [r1_rho_minus, r2_rho_minus] =
                [&]() {
                    auto const [r1, r2, _, __] =
                        evaluate(n_l, std::max(rho_floor, rho_lR - h_rho));
                    (void)_;
                    (void)__;
                    return std::pair{r1, r2};
                }();

            double const denom_n = (n_l + h_n) - std::max(n_l_floor, n_l - h_n);
            double const denom_rho =
                (rho_lR + h_rho) - std::max(rho_floor, rho_lR - h_rho);
            if (!(denom_n > 0.0 && denom_rho > 0.0))
            {
                break;
            }

            J11 = (r1_n_plus - r1_n_minus) / denom_n;
            J21 = (r2_n_plus - r2_n_minus) / denom_n;
            J12 = (r1_rho_plus - r1_rho_minus) / denom_rho;
            J22 = (r2_rho_plus - r2_rho_minus) / denom_rho;
        }
        else
        {
            // Analytic micro 2x2 Jacobian (default): no extra evaluate() calls.
            auto const [a11, a12, a21, a22] =
                evaluate_analytic_jacobian(n_l, rho_lR);
            J11 = a11;
            J12 = a12;
            J21 = a21;
            J22 = a22;
        }

        double const det = J11 * J22 - J12 * J21;
        if (!(std::isfinite(det) && std::abs(det) > 1e-24))
        {
            break;
        }

        double const delta_n = (-mass_residual * J22 +
                                density_residual * J12) /
                               det;
        double const delta_rho = (J21 * mass_residual -
                                  J11 * density_residual) /
                                 det;

        double step_scale = 1.0;
        bool accepted = false;
        for (int backtrack = 0; backtrack < 12; ++backtrack)
        {
            double const n_candidate =
                std::clamp(n_l + step_scale * delta_n, n_l_floor, n_l_ceiling);
            double const rho_candidate =
                std::max(rho_floor, rho_lR + step_scale * delta_rho);
            auto const [cand_mass_residual, cand_density_residual,
                        cand_micro_potential, cand_exchange] =
                evaluate(n_candidate, rho_candidate);
            double const current_norm =
                std::abs(mass_residual) /
                    std::max(1.0, std::abs(rho_l_prev)) +
                std::abs(density_residual) / std::max(1.0, std::abs(rho_lR));
            double const candidate_norm =
                std::abs(cand_mass_residual) /
                    std::max(1.0, std::abs(rho_l_prev)) +
                std::abs(cand_density_residual) /
                    std::max(1.0, std::abs(rho_candidate));
            if (candidate_norm <= current_norm || step_scale < 1e-3)
            {
                n_l = n_candidate;
                rho_lR = rho_candidate;
                out.n_l = n_l;
                out.rho_lR = rho_lR;
                out.micro_potential = cand_micro_potential;
                out.exchange = cand_exchange;
                accepted = true;
                break;
            }
            step_scale *= 0.5;
        }

        if (!accepted)
        {
            break;
        }

        if (std::abs(step_scale * delta_n) <=
                increment_tolerance * std::max(1.0, std::abs(n_l)) &&
            std::abs(step_scale * delta_rho) <=
                increment_tolerance * std::max(1.0, std::abs(rho_lR)))
        {
            out.converged = true;
            return out;
        }
    }

    auto const [mass_residual, density_residual, micro_potential, exchange] =
        evaluate(n_l, rho_lR);
    (void)mass_residual;
    (void)density_residual;
    out.n_l = n_l;
    out.rho_lR = rho_lR;
    out.micro_potential = micro_potential;
    out.exchange = exchange;
    out.converged = false;
    return out.converged ? out : predictor;
}

// ════════════════════════════════════════════════════════════════════════════
// KKT micro-water ceiling (branch dsm_mass_conservation_v3_kkt_ceiling_
// 2026-09-30; NOT adopted; reached only through potential_exchange
// micro_ceiling_treatment = kkt, default clamp = the code above, bitwise).
//
// Design, derivation and weak forms: DESIGN.md (3.2-3.5), DERIVATION.md (D-n),
// WEAK_FORMS.md (W-n) in ~/ogs-models/scratch/2026-09-30_kkt_ceiling_impl/.
// Symbols as in D-0: lambda (code: multiplier) is the multiplier of the ceiling
// n_l <= n_max(eps_v) = phi, r = alpha_M*lambda/rho_lR the rejected exchange,
// s the sign of the strain term (strain_sign), c_s = 1 - s*Deps.
//
// Nothing above is edited or shared: the solver below carries its own
// evaluation, written in the same expression order as the residual of
// solveReferenceMassStorageCoupledState (this file, `mass_residual`), because
// the file-scope FP_CONTRACT OFF above records that a refactor boundary alone
// once changed clang's fusions and broke dd1800 (DESIGN.md 3).
// ════════════════════════════════════════════════════════════════════════════

// Status of the local KKT solve at one evaluation (DESIGN.md 3.2). The value is
// stored in the state field MicroCeilingStatus and read by the assembler as the
// active flag (status == Active); no detector, no tolerance.
enum class MicroCeilingKktStatus : int
{
    Interior = 0,        // f(n_max) >= 0 (ties and NaN included): base solve
    Active = 1,          // f(n_max) < 0 and the scan found no interior root
    NonMonotone = 2,     // f(n_max) < 0 but an interior root exists (scan)
    PremiseViolated = 3  // f(n_floor) >= 0, or a non-finite value in the scan
};

struct MicroCeilingKktSolveData
{
    // n_l, rho_lR, micro_potential, exchange (= rhohat_pot at the RETURNED
    // state), converged. Interior / NonMonotone / PremiseViolated: the base
    // solve's result; Active: the wall state.
    MicroMacroMassStorageCoupledSolveData local;
    MicroCeilingKktStatus status = MicroCeilingKktStatus::Interior;
    double n_max = 0.0;              // [-], boundedMicroWaterContentCeiling
    double multiplier = 0.0;         // lambda [Pa], 0 unless Active
    double rejected_exchange = 0.0;  // r = rhohat_pot - S_s [kg/(m3 s)], 0 unless Active
    double exchange_received = 0.0;  // rhohat [kg/(m3 s)]: S_s Active, else rhohat_pot
    double booked_storage_rate = 0.0;  // S_s(n_max) [kg/(m3 s)], at every active candidate
    double f_at_n_max = 0.0;           // f(n_max) [kg/m3], diagnostic
};

// Bracketed unique-root test at an active candidate (DERIVATION.md 2.5, D-2.5).
// f and f_n are injected callables (n -> f(n), n -> df/dn) so that unit tests
// can supply synthetic functions. Precondition of the caller: f(n_max) < 0.
// The wall is the unique solution of the NCP iff f < 0 on all of
// [n_floor, n_max] (continuity, f(n_floor) < 0).
//   (i)   f(n_floor), f_n(n_floor) finite and f(n_floor) < 0, else PremiseViolated
//   (ii)  f(n_j) < 0 at every node of a log-uniform grid from n_floor to n_max
//         with nodes_per_decade nodes per decade, else NonMonotone
//   (iii) in every cell with f_n > 0 at the left node and f_n < 0 at the right
//         node, bisection on f_n (fixed number of halvings, no tolerance)
//         locates the local maximum; f < 0 is required there, else NonMonotone
// A non-finite value anywhere gives PremiseViolated (nothing can be certified).
// Limit (DERIVATION.md 2.5): two sign changes of f_n inside ONE cell escape.
template <typename F, typename Fn>
MicroCeilingKktStatus scanForInteriorCeilingRoot(F const& f, Fn const& f_n,
                                                 double const n_floor,
                                                 double const n_max,
                                                 int const nodes_per_decade)
{
    using Status = MicroCeilingKktStatus;
    double const f_floor = f(n_floor);
    double n_left = n_floor;
    double fn_left = f_n(n_floor);
    if (!(std::isfinite(f_floor) && std::isfinite(fn_left)))
    {
        return Status::PremiseViolated;
    }
    if (!(f_floor < 0.0))
    {
        return Status::PremiseViolated;
    }
    if (!(n_max > n_floor))
    {
        return Status::Active;  // degenerate bracket: the wall is the floor
    }

    double const decades = std::log10(n_max / n_floor);
    int const cells = std::max(
        1, static_cast<int>(std::ceil(static_cast<double>(nodes_per_decade) *
                                      decades)));
    for (int j = 1; j <= cells; ++j)
    {
        double const n_right =
            (j == cells) ? n_max
                         : n_floor * std::pow(n_max / n_floor,
                                              static_cast<double>(j) /
                                                  static_cast<double>(cells));
        double const f_right = f(n_right);
        double const fn_right = f_n(n_right);
        if (!(std::isfinite(f_right) && std::isfinite(fn_right)))
        {
            return Status::PremiseViolated;
        }
        if (!(f_right < 0.0))
        {
            return Status::NonMonotone;  // a node with f >= 0 proves a root
        }
        if (fn_left > 0.0 && fn_right < 0.0)
        {
            // Local maximum of f inside the cell: bisection on f_n in the
            // geometric mean (log n); 64 halvings exhaust the double
            // resolution of the bracket (count, not a tolerance).
            double a = n_left;
            double b = n_right;
            for (int it = 0; it < 64; ++it)
            {
                double const m = std::sqrt(a * b);
                double const fn_m = f_n(m);
                if (!std::isfinite(fn_m))
                {
                    return Status::PremiseViolated;
                }
                if (fn_m > 0.0)
                {
                    a = m;
                }
                else
                {
                    b = m;
                }
            }
            double const f_star = f(std::sqrt(a * b));
            if (!std::isfinite(f_star))
            {
                return Status::PremiseViolated;
            }
            if (!(f_star < 0.0))
            {
                return Status::NonMonotone;
            }
        }
        n_left = n_right;
        fn_left = fn_right;
    }
    return Status::Active;
}

// Drop-in for solveReferenceMassStorageCoupledState (same argument list) with
// the ceiling as a complementarity condition, DESIGN.md 3.2:
//   R_rho explicit: rho_lR = rho_lR_EOS(n) slaves the second unknown, so the
//   local problem is the scalar NCP  f(n) + dt*alpha_M*lambda/rho_lR = 0,
//   lambda >= 0, n <= n_max, lambda*(n_max - n) = 0,
//   f(n) = R_m(n, rho_lR_EOS(n), lambda = 0) (the residual of the base solve).
// Active candidate f(n_max) < 0 (no tolerance, no stored state read): the scan
// decides whether the wall is the unique solution. Active: n = n_max,
// rhohat = S_s(n_max) (the booked storage rate, received by micro AND macro),
// lambda = rho_lR*(rhohat_pot - S_s)/alpha_M >= 0, r = rhohat_pot - S_s.
// Interior (f(n_max) >= 0): the base solve, unchanged arithmetic, bitwise.
// NonMonotone / PremiseViolated: the KKT rule chooses nothing; the base result
// is returned with the status, the assembler treats the point like a clamp
// point. Memoryless: no argument and no field of the previous evaluation is read.
inline MicroCeilingKktSolveData solveReferenceMassStorageKktState(
    double const n_l_prev, double const rho_l_prev, double const rho_lR_prev,
    double const dt, double const rho_LR, double const alpha_bar,
    double const mu, YoungLaplaceMacroPotentialData const& macro_potential,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
    requirePositiveViscosity("solveReferenceMassStorageKktState", mu);
    constexpr double n_l_floor = 1e-16;  // as F:939, F:1141 (existing literal)
    double const strain_sign =
        microMassStrainTermSign(potential_exchange_params);  // s [-]
    double const dt_safe = std::isfinite(dt) && dt > 0.0 ? dt : 0.0;
    double const alpha_M_effective = alpha_bar * rho_LR / mu;
    double const volumetric_strain_rate =
        dt_safe > 0.0
            ? (local_context.volumetric_strain -
               local_context.volumetric_strain_prev) /
                  dt_safe
            : 0.0;
    double const n_max =
        boundedMicroWaterContentCeiling(local_context, n_l_floor);

    auto const base_result = [&](MicroCeilingKktStatus const status,
                                 double const storage_rate,
                                 double const f_at_n_max)
    {
        MicroCeilingKktSolveData out;
        out.local = solveReferenceMassStorageCoupledState(
            n_l_prev, rho_l_prev, rho_lR_prev, dt, rho_LR, alpha_bar, mu,
            macro_potential, local_context, potential_exchange_params);
        out.status = status;
        out.n_max = n_max;
        out.multiplier = 0.0;
        out.rejected_exchange = 0.0;
        out.exchange_received = out.local.exchange.rho_l_hat;
        out.booked_storage_rate = storage_rate;
        out.f_at_n_max = f_at_n_max;
        return out;
    };

    if (dt_safe <= 0.0)
    {
        return base_result(MicroCeilingKktStatus::Interior, 0.0, 0.0);
    }

    // One evaluation of the EOS-slaved local problem at micro water content n:
    // f, S_s, f_n and the intermediates. No side effect.
    struct Evaluation
    {
        double f = 0.0;    // [kg/m3] R_m(n, rho_lR_EOS(n), lambda = 0)
        double S_s = 0.0;  // [kg/(m3 s)] (c_s rho_l - rho_l_prev)/dt
        double f_n = 0.0;  // [kg/m3] df/dn, total derivative along the EOS
        double rho_lR = 0.0;
        VanDerWaalsMicroPotentialData micro_potential;
        PotentialDrivenMassExchangeData exchange;
    };
    auto const evaluate = [&](double const n_l)
    {
        double const active_nS = computeActiveMicroSolidVolumeFraction(
            n_l, local_context, potential_exchange_params);
        auto const density = computeReducedMicroLiquidDensity(
            n_l, rho_LR, active_nS, potential_exchange_params);
        double const rho_lR = density.rho_lR;
        // C5 (unit-test finding, 2026-10-01): the VALUE of mu_lR must be
        // evaluated exactly as the base residual `evaluate` of
        // solveReferenceMassStorageCoupledState does, with the nS chain FROZEN
        // (dnS_dnl = 0). The live-nS chain (F2, dnS_dnl = -1 under
        // CurrentPorositySplit) changes dmu_lR_dnl and, through
        // Pi' = -rho_lR*dmu_lR_dnl, ALSO the value of the integrable Maxwell
        // partner (Pi + n_l*Pi')*eps_v: the comment of C2 ("only the derivative,
        // not values") was wrong whenever the film term is on and eps_v != 0.
        // The tangent potential below carries the live chain, as the base
        // analytic Jacobian `evaluate_analytic_jacobian` does.
        auto micro_potential = computeVanDerWaalsMicroPotential(
            n_l, rho_lR, active_nS,
            potential_exchange_params.micro_solid_density_reference,
            potential_exchange_params.hamaker_constant,
            potential_exchange_params.specific_surface,
            microPotentialSignFactorFromParameters(potential_exchange_params),
            effectiveAugmentationPrefactor(potential_exchange_params,
                                           local_context.phi),  // K [J/kg]
            potential_exchange_params.potential_augmentation_exponent,
            0.0 /*dnS_dnl: frozen nS, as the base residual*/,
            potential_exchange_params.micro_water_content_floor);
        // Film coupling is always on (D-1.8): this folds the Maxwell partner
        // mu_m into mu_lR, exactly as the base residual does.
        applyFilmPressureMicroPotential(micro_potential, n_l, rho_lR,
                                        active_nS, local_context,
                                        potential_exchange_params);
        // Tangent potential (live-nS chain), for J11 and J12 only.
        double const dnS_dnl =
            potential_exchange_params.micro_solid_volume_fraction_mode ==
                    MicroSolidVolumeFractionMode::CurrentPorositySplit
                ? -1.0
                : 0.0;
        auto micro_potential_tan = computeVanDerWaalsMicroPotential(
            n_l, rho_lR, active_nS,
            potential_exchange_params.micro_solid_density_reference,
            potential_exchange_params.hamaker_constant,
            potential_exchange_params.specific_surface,
            microPotentialSignFactorFromParameters(potential_exchange_params),
            effectiveAugmentationPrefactor(potential_exchange_params,
                                           local_context.phi),  // K [J/kg]
            potential_exchange_params.potential_augmentation_exponent,
            dnS_dnl, potential_exchange_params.micro_water_content_floor);
        applyFilmPressureMicroPotential(micro_potential_tan, n_l, rho_lR,
                                        active_nS, local_context,
                                        potential_exchange_params);
        double const mu_LR_active = macro_potential.mu_LR;
        double const mu_lR_active = micro_potential.mu_lR;
        auto const exchange = computePotentialDrivenMassExchange(
            alpha_M_effective, mu_LR_active, mu_lR_active);
        double const phi_cs =
            std::isfinite(local_context.phi)
                ? std::clamp(local_context.phi, 0.0, 1.0 - 1e-12)
                : std::clamp(
                      local_context.phi_M_prev + local_context.phi_m_prev, 0.0,
                      1.0 - 1e-12);
        double const one_minus_n_l_cs = std::max(1e-12, 1.0 - n_l);
        double const rho_l = (1.0 - phi_cs) / one_minus_n_l_cs * n_l * rho_lR;
        // Same expression order as `mass_residual` of the base solve.
        double const f = rho_l - rho_l_prev - dt_safe * exchange.rho_l_hat -
                         strain_sign * dt_safe * rho_l * volumetric_strain_rate;
        double const storage_term =
            rho_l - strain_sign * dt_safe * rho_l * volumetric_strain_rate -
            rho_l_prev;
        // df/dn = J11 + J12 * d rho_lR_EOS/dn (D-2.4), J11/J12 as the analytic
        // micro Jacobian of the base solve.
        double const time_factor =
            1.0 - strain_sign * dt_safe * volumetric_strain_rate;  // c_s
        double const drho_l_dnl = (1.0 - phi_cs) * rho_lR /
                                  (one_minus_n_l_cs * one_minus_n_l_cs);
        double const drho_l_drho_lR = (1.0 - phi_cs) / one_minus_n_l_cs * n_l;
        double const J11 = drho_l_dnl * time_factor +
                           dt_safe * alpha_M_effective *
                               micro_potential_tan.dmu_lR_dnl;
        double const J12 = drho_l_drho_lR * time_factor +
                           dt_safe * alpha_M_effective *
                               micro_potential_tan.dmu_lR_drho_lR;
        Evaluation out;
        out.f = f;
        out.S_s = storage_term / dt_safe;
        out.f_n = J11 + J12 * density.drho_lR_dnl;
        out.rho_lR = rho_lR;
        out.micro_potential = micro_potential;
        out.exchange = exchange;
        return out;
    };

    // Active candidate: f(n_max) < 0, no tolerance, no stored state. Ties and
    // non-finite values are inactive (the base arithmetic, bitwise).
    Evaluation const at_wall = evaluate(n_max);
    if (!(at_wall.f < 0.0))
    {
        return base_result(MicroCeilingKktStatus::Interior, at_wall.S_s,
                           at_wall.f);
    }
    if (!(std::isfinite(alpha_M_effective) && alpha_M_effective > 0.0 &&
          std::isfinite(at_wall.S_s) && std::isfinite(at_wall.rho_lR)))
    {
        return base_result(MicroCeilingKktStatus::PremiseViolated,
                           at_wall.S_s, at_wall.f);
    }

    // One-entry memo: the scan asks for f and f_n at the same n back to back.
    double memo_n = std::numeric_limits<double>::quiet_NaN();
    Evaluation memo;
    auto const evaluate_memo = [&](double const n_l) -> Evaluation const&
    {
        if (!(n_l == memo_n))
        {
            memo = evaluate(n_l);
            memo_n = n_l;
        }
        return memo;
    };
    auto const scan_status = scanForInteriorCeilingRoot(
        [&](double const n_l) { return evaluate_memo(n_l).f; },
        [&](double const n_l) { return evaluate_memo(n_l).f_n; }, n_l_floor,
        n_max, potential_exchange_params.micro_ceiling_scan_nodes_per_decade);
    if (scan_status != MicroCeilingKktStatus::Active)
    {
        return base_result(scan_status, at_wall.S_s, at_wall.f);
    }

    // Active, unique: the wall. rhohat = S_s (computed directly, so the micro
    // residual holds exactly); r = rhohat_pot - S_s >= 0 by the branch
    // condition (clamped at 0 against round-off only), D-2.4, D-3.1.
    MicroCeilingKktSolveData out;
    out.local.n_l = n_max;
    out.local.rho_lR = at_wall.rho_lR;
    out.local.micro_potential = at_wall.micro_potential;
    out.local.exchange = at_wall.exchange;  // rhohat_pot at the wall (D-7, G6)
    out.local.converged = true;
    out.status = MicroCeilingKktStatus::Active;
    out.n_max = n_max;
    double const r =
        std::max(0.0, at_wall.exchange.rho_l_hat - at_wall.S_s);  // [kg/(m3 s)]
    out.rejected_exchange = r;
    out.multiplier = at_wall.rho_lR * r / alpha_M_effective;  // [Pa]
    out.exchange_received = at_wall.S_s;
    out.booked_storage_rate = at_wall.S_s;
    out.f_at_n_max = at_wall.f;
    return out;
}

// Tangents of the booked exchange rhohat = S_s on the ACTIVE branch (n = n_max =
// phi, rho_lR = rho_lR_EOS(phi; p_L)), D-5.2 / DESIGN.md 3.2. Free function so
// that a unit test can compare it with central differences of the solver.
//   c_s = 1 - s*d_eps, rho_l = phi*rho_lR
//   d rhohat/d p_L   = c_s*phi*(d rho_lR/d rho_LR)*(d rho_LR/d p_L)/dt
//   d rhohat/d eps_v = [ -s*rho_l + c_s*phi'*(rho_lR + phi*d rho_lR/dn) ]/dt
// (s = -1, alpha = 1, constant rho_lR: rho_lR/dt). Both independent of the
// Maxwell partner (the active exchange is S_s, not mu_lR-driven).
struct CeilingKktActiveExchangeTangents
{
    double drhohat_dpL = 0.0;      // [kg/(m3 s)/Pa]
    double drhohat_deps_v = 0.0;   // [kg/(m3 s)] per unit eps_v
};

inline CeilingKktActiveExchangeTangents computeCeilingKktActiveExchangeTangents(
    double const strain_sign, double const d_eps, double const phi,
    double const rho_lR, double const drho_lR_dnl,
    double const drho_lR_drho_LR, double const drho_LR_dpL,
    double const dphi_deps_v, double const dt)
{
    double const c_s = 1.0 - strain_sign * d_eps;  // [-]
    double const rho_l = phi * rho_lR;             // [kg/m3]
    return {
        .drhohat_dpL = c_s * phi * drho_lR_drho_LR * drho_LR_dpL / dt,
        .drhohat_deps_v =
            (-strain_sign * rho_l +
             c_s * dphi_deps_v * (rho_lR + phi * drho_lR_dnl)) /
            dt};
}

// d phi/d eps_v of the porosity law of PorosityFromMassBalance,
//   phi = (phi_prev + alpha*w)/(1 + w),  w = Deps + Dp_eff*beta_SR,
// = (alpha - phi)/(1 + w); 0 where the porosity clamp acts (same test as
// clamp_active_B of the variant-B block, existing literal 1e-12). The dValue of
// the property is not implemented, which is why the derivative is hand-coded.
inline double porosityDerivativeWrtVolumetricStrain(double const alpha,
                                                    double const phi,
                                                    double const phi_prev,
                                                    double const w_eps)
{
    double const phi_unclamped = (phi_prev + alpha * w_eps) / (1.0 + w_eps);
    bool const clamp_active = std::abs(phi_unclamped - phi) >
                              1e-12 * std::max(1.0, std::abs(phi));
    return clamp_active ? 0.0 : (alpha - phi) / (1.0 + w_eps);  // [-]
}

// Group-3 iteration diagnostics (DESIGN.md 3.3, M3). Memory fields and counters
// are updated at every evaluation of an active-or-not IP:
//   new attempt (t != attempt_t): attempt_t = t, inc_last = 0, eps_seen = e,
//     no counting;
//   same iterate (e == eps_seen, e.g. the output re-evaluation): nothing;
//   else inc = e - eps_seen; if inc_last != 0 and status is Active now and was
//     Active at the previous evaluation: inc_alt += sign(inc) != sign(inc_last),
//     inc_same += sign(inc) == sign(inc_last); then inc_last = inc, eps_seen = e.
inline void updateCeilingIterationDiagnostics(
    double& eps_seen, double& inc_last, double& attempt_t, double& inc_alt,
    double& inc_same, double const eps_v, double const t,
    bool const active_now, bool const active_before)
{
    if (t != attempt_t)
    {
        attempt_t = t;
        inc_last = 0.0;
        eps_seen = eps_v;
        return;
    }
    if (eps_v == eps_seen)
    {
        return;
    }
    double const inc = eps_v - eps_seen;
    if (inc_last != 0.0 && active_now && active_before)
    {
        if ((inc > 0.0) != (inc_last > 0.0))
        {
            inc_alt += 1.0;
        }
        else
        {
            inc_same += 1.0;
        }
    }
    inc_last = inc;
    eps_seen = eps_v;
}

// Iteration trace of the KKT micro-water ceiling (DESIGN.md 3.8, M3; the
// discriminator of H-A branch chatter against H-B missing p-u entry, 4.7).
// With micro_ceiling_trace_elements = id1 id2 ... (default empty = off, no
// cost) every evaluation at the integration points of the listed elements
// (every assembleWithJacobian evaluation and the output re-evaluation, flagged)
// appends one line to kkt_trace.csv in the working directory of the process
// (append mode, one mutex, flushed per line so that a run killed by a wall-clock
// budget keeps its trace). DEVIATION from DESIGN.md 3.8 ("<output_prefix>_kkt_
// trace.csv"): the output prefix is not reachable from the local assembler; the
// suite runs each deck in its own directory.
// Columns: t, dt, eval_index_in_attempt, element, ip, eps_v, eps_v_prev_step,
// status, S_s, rhohat_received, rhohat_pot, f_at_n_max, lambda, n_l, phi,
// is_output_reeval. eval_index_in_attempt counts the evaluations at this IP with
// the same t (an attempt key as in 3.3); a violation shows as a non-monotone
// index.
inline void writeMicroCeilingTraceLine(
    MicroCeilingTraceTag const& tag, MicroCeilingKktSolveData const& kkt,
    PotentialExchangeLocalSolveContext const& local_context)
{
    static std::mutex mutex;
    static std::ofstream file;
    static std::map<std::pair<std::size_t, std::size_t>,
                    std::pair<double, long>>
        attempts;
    std::lock_guard<std::mutex> const lock(mutex);
    if (!file.is_open())
    {
        char const* const name = "kkt_trace.csv";
        std::error_code ec;
        bool const has_content =
            std::filesystem::exists(name, ec) &&
            std::filesystem::file_size(name, ec) > 0;
        file.open(name, std::ios::app);
        if (!file)
        {
            OGS_FATAL("micro_ceiling_trace_elements: cannot open {}.", name);
        }
        if (!has_content)
        {
            file << "t,dt,eval_index_in_attempt,element,ip,eps_v,"
                    "eps_v_prev_step,status,S_s,rhohat_received,rhohat_pot,"
                    "f_at_n_max,lambda,n_l,phi,is_output_reeval\n";
        }
    }
    auto& attempt =
        attempts[{tag.element_id, tag.integration_point}];
    if (attempt.second == 0 || attempt.first != tag.t)
    {
        attempt.first = tag.t;
        attempt.second = 0;
    }
    long const index = attempt.second++;
    file << std::setprecision(17) << tag.t << ',' << tag.dt << ',' << index
         << ',' << tag.element_id << ',' << tag.integration_point << ','
         << local_context.volumetric_strain << ','
         << local_context.volumetric_strain_prev << ','
         << static_cast<int>(kkt.status) << ',' << kkt.booked_storage_rate
         << ',' << kkt.exchange_received << ','
         << kkt.local.exchange.rho_l_hat << ',' << kkt.f_at_n_max << ','
         << kkt.multiplier << ',' << kkt.local.n_l << ',' << local_context.phi
         << ',' << (tag.output_reevaluation ? 1 : 0) << std::endl;
}

template <int DisplacementDim>
inline void applyReferenceMassStorageLocalState(
    StatefulData<DisplacementDim>& state_current,
    StatefulDataPrev<DisplacementDim> const& state_previous,
    MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
    double const rho_LR, PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params,
    MicroMacroMassStorageCoupledSolveData const& coupled_update)
{
    auto const transport_porosity_update = computeTransportPorosityUpdate(
        local_context.phi, local_context.phi_M_prev, local_context.phi_m_prev,
        coupled_update.n_l, local_context.volumetric_strain,
        local_context.volumetric_strain_prev,
        potential_exchange_params.macro_porosity_update_mode);

    auto const compatibility_output =
        computeCompatibilityMicroHydraulicOutput(
            coupled_update.n_l, rho_LR, local_context,
            potential_exchange_params);

    auto& n_l = std::get<MicroWaterContent>(state_current);
    *n_l = coupled_update.n_l;

    auto& rho_lR = std::get<MicroLiquidDensity>(state_current);
    *rho_lR = coupled_update.rho_lR;

    auto& micro_porosity = std::get<MicroPorosity>(state_current);
    *micro_porosity = transport_porosity_update.phi_m;

    auto& transport_porosity =
        std::get<ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(
            state_current)
            .phi;
    transport_porosity = transport_porosity_update.phi_M;
    variables.transport_porosity = transport_porosity_update.phi_M;
    variables_prev.transport_porosity = transport_porosity_update.phi_M_prev;

    auto& porosity =
        std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(state_current).phi;
    porosity = transport_porosity_update.phi_M + transport_porosity_update.phi_m;
    variables.porosity = porosity;
    variables_prev.porosity =
        transport_porosity_update.phi_M_prev + transport_porosity_update.phi_m_prev;

    auto& p_L_m = std::get<MicroPressure>(state_current);
    auto& S_L_m = std::get<MicroSaturation>(state_current);
    auto& rho_l_hat = std::get<MicroExchangeSource>(state_current);
    *p_L_m = compatibility_output.p_L_m;
    *S_L_m = compatibility_output.S_L_m;
    rho_l_hat = MicroExchangeSource{coupled_update.exchange.rho_l_hat};

    (void)state_previous;
}

inline VanDerWaalsMicroPotentialData computeActiveMicroPotential(
    double const n_l, double const rho_lR,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    double const active_nS =
        computeActiveMicroSolidVolumeFraction(n_l, local_context, potential_exchange_params);
    double const rho_lR_effective =
        potential_exchange_params.local_nonlinear_solve_mode ==
                LocalNonlinearSolveMode::ScalarReferenceMassStorage
            ? computeReducedMicroLiquidDensity(n_l, rho_lR, active_nS, potential_exchange_params)
                  .rho_lR
            : rho_lR;
    // F2: under current_porosity_split active_nS = 1 - n_l is LIVE, so the vdW
    // dmu_lR/dnl total derivative picks up the dmu_lR_dnS*(dnS/dnl) chain with
    // dnS/dnl = -1. In reference mode nS is constant -> dnS_dnl = 0 (exact, no
    // change). Tangent-only: forward solves use FD Jacobians.
    double const dnS_dnl =
        potential_exchange_params.micro_solid_volume_fraction_mode ==
                MicroSolidVolumeFractionMode::CurrentPorositySplit
            ? -1.0
            : 0.0;
    auto out = computeVanDerWaalsMicroPotential(
        n_l, rho_lR_effective, active_nS, potential_exchange_params.micro_solid_density_reference,
        potential_exchange_params.hamaker_constant, potential_exchange_params.specific_surface,
        microPotentialSignFactorFromParameters(potential_exchange_params),
            // Live K(rho_d): rho_d = rho_SR*(1-phi); off / phi sentinel ->
            // parse scalar (K_OF_RHO_D_LIVE.md).
            effectiveAugmentationPrefactor(potential_exchange_params,
                                           local_context.phi),  // K [J/kg]
            potential_exchange_params.potential_augmentation_exponent, dnS_dnl,
            potential_exchange_params.micro_water_content_floor);

    // ── Film-pressure coupling (maxwell sec.5): ONE evaluator ────────────────
    // Fold the smoothly-gated film delta into mu_lR via the shared helper, so the
    // SAME mu_lR(p_film) propagates to EVERY consumer of computeActiveMicroPotential
    // (the scalar local n_l solve eval_at, the p_L_m writer
    // computeCompatibilityMicroHydraulicOutput, the residual/Jacobian assembly)
    // AND stays identical to the mass-storage 2x2 solve (increment E). Flag OFF or
    // NaN p_conf -> the helper is a no-op (pure vdW), bit-for-bit. Uses the SAME
    // density (rho_lR_effective) the vdW formula above consumed.
    applyFilmPressureMicroPotential(out, n_l, rho_lR_effective, active_nS,
                                    local_context, potential_exchange_params);
    return out;
}

inline CompatibilityMicroHydraulicOutputData
computeCompatibilityMicroHydraulicOutput(
    double const n_l, double const rho_LR,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    double const n_l_safe = std::max(1e-16, n_l);
    double const n_l_ref = std::max(
        1e-16, potential_exchange_params.initial_micro_water_content.value_or(
                   potential_exchange_params.micro_solid_volume_fraction_reference));

    // micro_liquid_density.rho_lR: EOS-derived confined water density (~1100 kg/m³).
    // This is the physically correct density for the vdW specific free energy
    // denominator (energy/area × area/REV / mass/REV) and for p_L_m = -rho*mu_lR.
    // computeActiveMicroPotential internally uses micro density for the vdW formula
    // when local_nonlinear_solve_mode == ScalarReferenceMassStorage (all MS33 cases).
    // use_micro_liquid_density_for_micro_pressure should be true in all PRJ files.
    auto const micro_liquid_density = computeActiveMicroLiquidDensity(
        n_l_safe, rho_LR, local_context, potential_exchange_params);
    auto const micro_potential = computeActiveMicroPotential(
        n_l_safe, rho_LR, local_context, potential_exchange_params);
    double const p_L_m_density =
        potential_exchange_params.use_micro_liquid_density_for_micro_pressure
            ? micro_liquid_density.rho_lR   // correct: confined water density
            : rho_LR;                        // fallback: bulk density (~10% error)

    return {
        .p_L_m = -p_L_m_density * micro_potential.mu_lR,
        .S_L_m = n_l_safe / n_l_ref,
        .n_l_ref = n_l_ref,
        .micro_potential = micro_potential,
    };
}

inline ImplicitMicroWaterContentUpdateData solveImplicitMicroWaterContent(
    double const n_l_prev, double const dt, double const rho_LR,
    double const alpha_bar, double const mu,
    YoungLaplaceMacroPotentialData const& macro_potential,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    requirePositiveViscosity("solveImplicitMicroWaterContent", mu);
    constexpr double n_l_floor = 1e-16;
    double const dt_safe = std::isfinite(dt) && dt > 0.0 ? dt : 0.0;
    double const alpha_M_effective = alpha_bar * rho_LR / mu;
    bool const use_microstate_storage_mode =
        potential_exchange_params.local_nonlinear_solve_mode !=
        LocalNonlinearSolveMode::ScalarExchange;
    bool const use_mass_storage =
        potential_exchange_params.local_nonlinear_solve_mode ==
        LocalNonlinearSolveMode::ScalarReferenceMassStorage;
    double const volumetric_strain_rate =
        dt_safe > 0.0
            ? (local_context.volumetric_strain -
               local_context.volumetric_strain_prev) /
                  dt_safe
            : 0.0;
    double const n_l_ceiling =
        use_microstate_storage_mode
            ? boundedMicroWaterContentCeiling(local_context, n_l_floor)
            : std::max(n_l_floor, 1.0);

    auto eval_at = [&](double const n_l)
    {
        auto const micro_potential =
            computeActiveMicroPotential(n_l, rho_LR, local_context, potential_exchange_params);
        double const mu_LR_active = macro_potential.mu_LR;
        double const mu_lR_active = micro_potential.mu_lR;
        auto const exchange = computePotentialDrivenMassExchange(
            alpha_M_effective, mu_LR_active, mu_lR_active);
        auto const micro_liquid_density =
            use_mass_storage
                ? std::optional<ReducedMicroLiquidDensityData>{
                      computeActiveMicroLiquidDensity(
                          n_l, rho_LR, local_context, potential_exchange_params)}
                : std::nullopt;
        return std::tuple{micro_potential, exchange, micro_liquid_density};
    };

    auto const prev_micro_liquid_density =
        use_mass_storage
            ? std::optional<ReducedMicroLiquidDensityData>{
                  computePreviousMicroLiquidDensity(n_l_prev, rho_LR,
                                                      local_context, potential_exchange_params)}
            : std::nullopt;
    // REV-scale previous liquid apparent density: phi_m_prev * rho_lR_prev.
    // local_context.phi_m_prev = (1-phi_M_prev)*n_l_prev (hierarchical split).
    double const rho_l_prev =
        prev_micro_liquid_density
            ? local_context.phi_m_prev * prev_micro_liquid_density->rho_lR
            : 0.0;

    ImplicitMicroWaterContentUpdateData out;
    if (dt_safe <= 0.0)
    {
        out.n_l = std::clamp(n_l_prev, n_l_floor, n_l_ceiling);
        auto const [micro_potential, exchange, micro_liquid_density] =
            eval_at(out.n_l);
        (void)micro_liquid_density;
        out.micro_potential = micro_potential;
        out.exchange = exchange;
        return out;
    }

    if (use_mass_storage)
    {
        auto const coupled_update = solveReferenceMassStorageCoupledState(
            n_l_prev, rho_l_prev,
            prev_micro_liquid_density ? prev_micro_liquid_density->rho_lR
                                      : rho_LR,
            dt_safe, rho_LR, alpha_bar, mu, macro_potential,
            local_context, potential_exchange_params);
        out.n_l = coupled_update.n_l;
        out.rho_lR = coupled_update.rho_lR;
        out.micro_potential = coupled_update.micro_potential;
        out.exchange = coupled_update.exchange;
        out.converged = coupled_update.converged;
        return out;
    }

    double n_l = std::clamp(n_l_prev, n_l_floor, n_l_ceiling);
    constexpr int max_iterations = 25;
    constexpr double residual_tolerance = 1e-12;
    constexpr double increment_tolerance = 1e-12;
    bool converged = false;

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        auto const [micro_potential, exchange, micro_liquid_density] =
            eval_at(n_l);
        double residual = 0.0;
        double jacobian = 0.0;
        double const drho_l_hat_dn_l =
            exchange.drho_l_hat_dmu_lR * micro_potential.dmu_lR_dnl;
        if (use_mass_storage)
        {
            double const rho_l =
                n_l * micro_liquid_density->rho_lR;
            residual =
                rho_l - rho_l_prev - dt_safe * exchange.rho_l_hat;
            // (Unreachable: use_mass_storage returns above; V2 keeps the sign
            // consistent with the live mass-storage solves.)
            double const strain_sign_dead =
                microMassStrainTermSign(potential_exchange_params);  // [-]
            residual -= strain_sign_dead * dt_safe * rho_l *
                        volumetric_strain_rate;

            jacobian = micro_liquid_density->drho_l_dn_l -
                       dt_safe * drho_l_hat_dn_l;
            jacobian -=
                strain_sign_dead * dt_safe *
                micro_liquid_density->drho_l_dn_l * volumetric_strain_rate;
        }
        else
        {
            residual =
                n_l - n_l_prev - dt_safe * exchange.rho_l_hat / rho_LR;
            if (use_microstate_storage_mode)
            {
                residual -= dt_safe * n_l * volumetric_strain_rate;
            }

            jacobian = 1.0 - dt_safe * drho_l_hat_dn_l / rho_LR;
            if (use_microstate_storage_mode)
            {
                jacobian -= dt_safe * volumetric_strain_rate;
            }
        }

        if (std::abs(residual) <=
            residual_tolerance * std::max(1.0, std::abs(n_l_prev)))
        {
            converged = true;
            out.n_l = n_l;
            out.micro_potential = micro_potential;
            out.exchange = exchange;
            break;
        }

        if (!(std::isfinite(jacobian) && std::abs(jacobian) > 1e-20))
        {
            break;
        }

        double const delta_n_l = -residual / jacobian;
        double const n_l_candidate =
            std::clamp(n_l + delta_n_l, n_l_floor, n_l_ceiling);
        if (std::abs(n_l_candidate - n_l) <=
            increment_tolerance * std::max(1.0, std::abs(n_l)))
        {
            auto const [micro_potential_candidate, exchange_candidate,
                        micro_density_candidate] =
                eval_at(n_l_candidate);
            (void)micro_density_candidate;
            out.n_l = n_l_candidate;
            out.micro_potential = micro_potential_candidate;
            out.exchange = exchange_candidate;
            converged = true;
            break;
        }

        n_l = n_l_candidate;
    }

    if (!converged)
    {
        // Fallback to explicit update if local scalar Newton does not converge.
        auto const [micro_potential_prev, exchange_prev, micro_density_prev] =
            eval_at(std::clamp(n_l_prev, n_l_floor, n_l_ceiling));
        double explicit_increment = 0.0;
        if (use_mass_storage)
        {
            double const rho_l_prev_fallback =
                std::clamp(n_l_prev, n_l_floor, n_l_ceiling) *
                micro_density_prev->rho_lR;
            explicit_increment =
                dt_safe * exchange_prev.rho_l_hat /
                std::max(1e-16, micro_density_prev->rho_lR);
            explicit_increment +=
                microMassStrainTermSign(potential_exchange_params) * dt_safe *
                rho_l_prev_fallback * volumetric_strain_rate /
                std::max(1e-16, micro_density_prev->rho_lR);
        }
        else
        {
            explicit_increment =
                dt_safe * exchange_prev.rho_l_hat / rho_LR;
            if (use_microstate_storage_mode)
            {
                explicit_increment +=
                    dt_safe * std::clamp(n_l_prev, n_l_floor, n_l_ceiling) *
                    volumetric_strain_rate;
            }
        }
        out.n_l = std::clamp(n_l_prev + explicit_increment, n_l_floor,
                             n_l_ceiling);
        auto const [micro_potential_fallback, exchange_fallback,
                    micro_density_fallback] =
            eval_at(out.n_l);
        (void)micro_density_fallback;
        out.micro_potential = micro_potential_fallback;
        out.exchange = exchange_fallback;
        out.converged = false;

        static std::once_flag once;
        std::call_once(once, []
        {
            WARN(
                "DSM: local implicit n_l solve did not converge at least once; falling back to explicit n_l update for robustness.");
        });
        return out;
    }

    out.converged = true;
    return out;
}

inline double computeImplicitNlDpL(
    double const n_l_prev, double const p_L_ip, double const dt,
    double const rho_LR, double const drho_LR_dpL,
    double const alpha_bar, double const mu,
    YoungLaplaceMacroPotentialData const& macro_potential,
    VanDerWaalsMicroPotentialData const& micro_potential,
    PotentialDrivenMassExchangeData const& exchange,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params,
    double const n_l_converged = std::numeric_limits<double>::quiet_NaN(),
    double const rho_lR_micro = std::numeric_limits<double>::quiet_NaN())
{
    requirePositiveViscosity("computeImplicitNlDpL", mu);
    double const dt_safe = std::isfinite(dt) && dt > 0.0 ? dt : 0.0;
    if (dt_safe <= 0.0)
    {
        return 0.0;
    }

    // ── F1 (2026-06-06, tangent-only): ScalarReferenceMassStorage REV-mass
    // residual linearization ───────────────────────────────────────────────
    // The previous shared analytic below linearized the n_l-NORMALIZED residual
    // r = (n_l - n_l_prev) - dt*rho_l_hat/rho_LR - dt*eps_v_rate*n_l. But in
    // ScalarReferenceMassStorage mode solveReferenceMassStorageCoupledState
    // actually solves the REV-MASS residual
    //   r1 = rho_l*(1 - dt*eps_v_rate) - rho_l_prev - dt*rho_l_hat,
    //   rho_l = phi_m*rho_lR,  phi_m = (1-phi)/(1-n_l)*n_l,
    // with rho_lR a 2nd unknown slaved to n_l by the density EOS residual r2=0.
    // dn_l/dp_L from the n_l-normalized form is therefore inconsistent with the
    // solved system (it greens the wrong test). Rebuild the tangent on the REV
    // residual with rho_lR eliminated along r2=0 (1x1 reduction; the EOS already
    // returns the slaved drho_lR/dn_l). TANGENT-ONLY: the converged forward
    // solve is untouched; only the Newton Jacobian is made consistent. The
    // converged (n_l, rho_lR) are threaded in by the caller; micro_potential /
    // exchange are the values evaluated at that converged state.
    if (potential_exchange_params.local_nonlinear_solve_mode ==
        LocalNonlinearSolveMode::ScalarReferenceMassStorage)
    {
        double const eps_v_rate =
            (local_context.volumetric_strain -
             local_context.volumetric_strain_prev) /
            dt_safe;
        // V2 (F3): s = +1 shipped, -1 Eulerian (micro_mass_strain_term_eulerian).
        double const time_factor =
            1.0 - microMassStrainTermSign(potential_exchange_params) *
                      dt_safe * eps_v_rate;  // [-]

        // Converged n_l (fall back to n_l_prev only if the caller omitted it).
        double const n_l =
            std::max(1e-16, std::isfinite(n_l_converged) ? n_l_converged
                                                         : n_l_prev);
        // Active micro-solid fraction (LIVE = 1 - n_l in CurrentPorositySplit,
        // constant in Reference) and the EOS at the converged state.
        double const nS = computeActiveMicroSolidVolumeFraction(
            n_l, local_context, potential_exchange_params);
        auto const eos = computeReducedMicroLiquidDensity(
            n_l, rho_LR, nS, potential_exchange_params);
        double const rho_lR = (std::isfinite(rho_lR_micro) && rho_lR_micro > 0.0)
                                  ? rho_lR_micro
                                  : eos.rho_lR;

        // REV macro porosity phi (current step): prefer ctx.phi, else prev sum.
        double const phi = std::isfinite(local_context.phi)
                               ? std::clamp(local_context.phi, 0.0, 1.0 - 1e-12)
                               : std::clamp(local_context.phi_M_prev +
                                                local_context.phi_m_prev,
                                            0.0, 1.0 - 1e-12);
        double const c = 1.0 - phi;
        double const one_minus_n_l = std::max(1e-12, 1.0 - n_l);
        double const f = n_l / one_minus_n_l;
        double const f_prime = 1.0 / (one_minus_n_l * one_minus_n_l);
        double const phi_m = c * f;

        // Total drho_l/dn_l along the EOS-slaved manifold (rho_l = c*f*rho_lR):
        //   drho_l/dn_l = c*( f'*rho_lR + f*drho_lR/dn_l ).
        double const drho_l_dn_l =
            c * (f_prime * rho_lR + f * eos.drho_lR_dnl);

        // Total dmu_lR/dn_l along the manifold (micro_potential carries the
        // partial dmu_lR_dnl and the rho_lR channel dmu_lR_drho_lR).
        double const dmu_lR_dn_l_tot =
            micro_potential.dmu_lR_dnl +
            micro_potential.dmu_lR_drho_lR * eos.drho_lR_dnl;
        double const drho_l_hat_dn_l =
            exchange.drho_l_hat_dmu_lR * dmu_lR_dn_l_tot;

        double const dr_dn_l =
            drho_l_dn_l * time_factor - dt_safe * drho_l_hat_dn_l;
        if (!(std::isfinite(dr_dn_l) && std::abs(dr_dn_l) > 1e-20))
        {
            return 0.0;
        }

        // p_L channel at fixed n_l. rho_lR varies with p_L only through the bulk
        // rho_LR appearing additively in the EOS: g = rho_lR - rho_LR -
        // rho_l0*exp(-a*omega^b) = 0 (omega = n_l*rho_lR/(nS*rho_SR), no rho_LR),
        // so drho_lR/drho_LR|_{fixed n_l} = -dg/drho_LR / dg/drho_lR = 1/dg_drho.
        double drho_lR_dpL_fixed_n = 0.0;
        if (drho_LR_dpL != 0.0)
        {
            double const rho_SR = std::max(
                1e-16, potential_exchange_params.micro_solid_density_reference);
            double const a_rho =
                std::max(1e-16, potential_exchange_params.micro_liquid_density_a);
            double const b_rho =
                std::max(1e-16, potential_exchange_params.micro_liquid_density_b);
            double const denom = std::max(1e-16, nS) * rho_SR;
            double const common = (rho_lR - rho_LR) * a_rho * b_rho *
                                  std::pow(std::max(1e-16, eos.omega_l),
                                           b_rho - 1.0);
            double const dg_drho =
                1.0 + common * (std::max(1e-16, n_l) / denom);
            drho_lR_dpL_fixed_n =
                (std::isfinite(dg_drho) && std::abs(dg_drho) > 1e-20)
                    ? drho_LR_dpL / dg_drho
                    : 0.0;
        }

        double const dalpha_M_dpL = alpha_bar / mu * drho_LR_dpL;
        double const dmu_LR_dpL = macro_potential.dmu_LR_dpLR +
                                  macro_potential.dmu_LR_drho_LR * drho_LR_dpL;
        double const dmu_lR_dpL_fixed_n =
            micro_potential.dmu_lR_drho_lR * drho_lR_dpL_fixed_n;
        double const drho_l_hat_dpL_fixed_n =
            exchange.drho_l_hat_dalpha_M * dalpha_M_dpL +
            exchange.drho_l_hat_dmu_LR * dmu_LR_dpL +
            exchange.drho_l_hat_dmu_lR * dmu_lR_dpL_fixed_n;

        double const drho_l_dpL_fixed_n = phi_m * drho_lR_dpL_fixed_n;
        double const dr_dpL =
            drho_l_dpL_fixed_n * time_factor - dt_safe * drho_l_hat_dpL_fixed_n;
        // L3 (review 2026-06-14) — NOW WIRED (Jacobian-only), supersedes the
        // earlier DOCUMENTED-NOT-WIRED note. In live-K mode the REV-mass
        // residual r also depends on the augmentation prefactor K through
        // micro_potential.mu_lR, and K = K_table(rho_SR*(1-phi)) couples to
        // displacement via phi(eps_v). The converged local n_l therefore carries
        // an implicit strain channel
        //   dn_l/d eps_v |_K = dn_l/dK * dK/dphi * dphi/deps_v,
        //   dn_l/dK = -(dr/dK)/(dr/dn_l)
        //           = (dt * drho_l_hat_dmu_lR * dmu_lR_dK) / dr_dn_l,
        // which THIS dn_l/dpL (a FIXED-n_l-vs-pL sensitivity) does NOT carry.
        // That gap is now closed by the SIBLING helper computeImplicitNlDK
        // (above), consumed at the M2 swelling-eigenstress assembly site to add
        // the implicit-n_l(K) half of the K[u,u]/K[u,p] tangent. The wiring is
        // JACOBIAN-ONLY: this function's return value, the local FORWARD n_l
        // solve, and the residual are all UNCHANGED here — only a new analytic
        // tangent contribution was added at assembly. (The original concern that
        // wiring "risks the converged forward solve" was avoided by NOT touching
        // this return / the solve and adding the sensitivity purely on the
        // Jacobian side; review L3, K_OF_RHO_D_LIVE.md.) PREDICTED (not yet
        // verified by re-run): completes the live-K mass-storage displacement
        // tangent (1b_A form-(a) candidate cure); the converged root is
        // unaffected by construction.
        return -dr_dpL / dr_dn_l;
    }

    // ── ScalarExchange / ScalarReferenceStorage: n_l-normalized residual ─────
    // (unchanged; F1 above leaves these modes bit-for-bit.)
    // P2 fix (2026-06-06): ScalarReferenceStorage previously used a
    // finite-difference dn_l/dp_L here (perturbing the full coupled solve by
    // h ~ 1e-8*|p_L|). At a dry IC the two perturbed coupled solves barely move
    // -> catastrophic cancellation -> random-sign ~3e-13 noise -> corrupts the
    // global pressure-block diagonal (drho_L_hat_dpL_direct, ~line 3824) ->
    // step-1 macro-pressure blow-up. Fall through to the ANALYTIC tangent below.
    // This is TANGENT-ONLY: the converged mass-conserving forward solve is
    // unchanged; only the Newton Jacobian gets a clean, correctly-signed value.

    double const dalpha_M_effective_dpL = alpha_bar / mu * drho_LR_dpL;
    double const dmu_first_dpL_fixed_n =
        macro_potential.dmu_LR_dpLR +
        macro_potential.dmu_LR_drho_LR * drho_LR_dpL;
    double const dmu_second_dpL_fixed_n = micro_potential.dmu_lR_drho_lR *
                                          drho_LR_dpL;

    double const drho_l_hat_dpL_fixed_n =
        exchange.drho_l_hat_dalpha_M * dalpha_M_effective_dpL +
        exchange.drho_l_hat_dmu_LR * dmu_first_dpL_fixed_n +
        exchange.drho_l_hat_dmu_lR * dmu_second_dpL_fixed_n;
    double const drho_l_hat_dn_l =
        exchange.drho_l_hat_dmu_lR * micro_potential.dmu_lR_dnl;

    double dr_dn_l = 1.0 - dt_safe * drho_l_hat_dn_l / rho_LR;
    if (potential_exchange_params.local_nonlinear_solve_mode ==
            LocalNonlinearSolveMode::ScalarReferenceStorage ||
        potential_exchange_params.local_nonlinear_solve_mode ==
            LocalNonlinearSolveMode::ScalarReferenceMassStorage)
    {
        double const volumetric_strain_rate =
            (local_context.volumetric_strain -
             local_context.volumetric_strain_prev) /
            dt_safe;
        dr_dn_l -= dt_safe * volumetric_strain_rate;
    }
    if (!(std::isfinite(dr_dn_l) && std::abs(dr_dn_l) > 1e-20))
    {
        return 0.0;
    }

    double const dr_dp_l =
        -dt_safe * (drho_l_hat_dpL_fixed_n / rho_LR -
                    exchange.rho_l_hat / (rho_LR * rho_LR) * drho_LR_dpL);
    return -dr_dp_l / dr_dn_l;
}

// ── L3 (review 2026-06-14, JACOBIAN-ONLY) ────────────────────────────────────
// Sensitivity of the LOCALLY-SOLVED micro water content n_l to the augmentation
// prefactor K, for ScalarReferenceMassStorage mode. Sibling of
// computeImplicitNlDpL: identical 1x1 REV-mass reduction (rho_lR slaved along
// the density EOS r2=0), differing only in which partial of the residual r is
// taken. K enters r ONLY through the exchange term rho_l_hat (via mu_lR; the
// rho_l = phi_m*rho_lR mass term carries no K), so
//   dr/dK = -dt * drho_l_hat/dK = -dt * exchange.drho_l_hat_dmu_lR
//                                 * micro_potential.dmu_lR_dK,
// and by the implicit-function theorem on r(n_l;K)=0 at the converged state
//   dn_l/dK = -(dr/dK) / (dr/dn_l).                                  [n_l per (J/kg)]
// This is the channel the fixed-n_l dn_l/dpL does NOT carry; in live-K mode it
// closes the implicit n_l(K(phi(eps_v))) strain channel of the swelling
// eigenstress (review L3), wired at the M2 displacement-Jacobian site.
//
// RESIDUAL-SAFE: reads ONLY the already-converged (n_l, rho_lR, micro_potential,
// exchange) the caller threads in -- it never re-solves, never mutates the
// forward state, and is identically 0 outside ScalarReferenceMassStorage (and
// when dt<=0). dr_dn_l is rebuilt here by the SAME expression as
// computeImplicitNlDpL so numerator and denominator are linearized about one
// state. Returns 0 on a singular/non-finite dr_dn_l (matching the dn_l/dpL
// guard), so a degenerate tangent silently drops rather than poisoning K[u,u].
inline double computeImplicitNlDK(
    double const n_l_prev, double const dt, double const rho_LR,
    double const mu,
    VanDerWaalsMicroPotentialData const& micro_potential,
    PotentialDrivenMassExchangeData const& exchange,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params,
    double const n_l_converged = std::numeric_limits<double>::quiet_NaN(),
    double const rho_lR_micro = std::numeric_limits<double>::quiet_NaN())
{
    requirePositiveViscosity("computeImplicitNlDK", mu);
    double const dt_safe = std::isfinite(dt) && dt > 0.0 ? dt : 0.0;
    if (dt_safe <= 0.0)
    {
        return 0.0;
    }
    if (potential_exchange_params.local_nonlinear_solve_mode !=
        LocalNonlinearSolveMode::ScalarReferenceMassStorage)
    {
        // dn_l/dK only defined for the local mass-storage solve; other modes
        // do not solve a K-dependent n_l here (the channel is absent / handled
        // elsewhere), so the L3 chain is exactly zero.
        return 0.0;
    }

    double const eps_v_rate =
        (local_context.volumetric_strain -
         local_context.volumetric_strain_prev) /
        dt_safe;  // 1/s
    double const time_factor =
        1.0 - microMassStrainTermSign(potential_exchange_params) * dt_safe *
                  eps_v_rate;  // [-]  (V2 F3 sign, +1 = shipped)

    // Converged n_l (fall back to n_l_prev only if the caller omitted it) --
    // identical guard to computeImplicitNlDpL so dr_dn_l linearizes about the
    // same state.
    double const n_l =
        std::max(1e-16, std::isfinite(n_l_converged) ? n_l_converged
                                                     : n_l_prev);
    double const nS = computeActiveMicroSolidVolumeFraction(
        n_l, local_context, potential_exchange_params);  // [-]
    auto const eos = computeReducedMicroLiquidDensity(
        n_l, rho_LR, nS, potential_exchange_params);
    double const rho_lR = (std::isfinite(rho_lR_micro) && rho_lR_micro > 0.0)
                              ? rho_lR_micro
                              : eos.rho_lR;  // kg/m^3

    double const phi = std::isfinite(local_context.phi)
                           ? std::clamp(local_context.phi, 0.0, 1.0 - 1e-12)
                           : std::clamp(local_context.phi_M_prev +
                                            local_context.phi_m_prev,
                                        0.0, 1.0 - 1e-12);  // [-]
    double const c = 1.0 - phi;  // [-]
    double const one_minus_n_l = std::max(1e-12, 1.0 - n_l);  // [-]
    double const f = n_l / one_minus_n_l;                     // [-]
    double const f_prime = 1.0 / (one_minus_n_l * one_minus_n_l);  // [1/n_l]

    // dr/dn_l rebuilt EXACTLY as computeImplicitNlDpL (mass-storage branch):
    // r = rho_l*time_factor - rho_l_prev - dt*rho_l_hat, rho_l = c*f*rho_lR.
    double const drho_l_dn_l =
        c * (f_prime * rho_lR + f * eos.drho_lR_dnl);  // kg/m^3 per n_l
    double const dmu_lR_dn_l_tot =
        micro_potential.dmu_lR_dnl +
        micro_potential.dmu_lR_drho_lR * eos.drho_lR_dnl;  // (J/kg) per n_l
    double const drho_l_hat_dn_l =
        exchange.drho_l_hat_dmu_lR * dmu_lR_dn_l_tot;  // (kg/m^3/s) per n_l
    double const dr_dn_l =
        drho_l_dn_l * time_factor - dt_safe * drho_l_hat_dn_l;  // kg/m^3 per n_l
    if (!(std::isfinite(dr_dn_l) && std::abs(dr_dn_l) > 1e-20))
    {
        return 0.0;
    }

    // dr/dK: K enters r ONLY via rho_l_hat = exchange(mu_lR(.;K)); the mass term
    // rho_l = c*f*rho_lR has no K dependence (the EOS omega has no K). So
    //   dr/dK = -dt * drho_l_hat/dmu_lR * dmu_lR/dK.            [kg/m^3 per (J/kg)]
    // micro_potential.dmu_lR_dK is the augmentation channel (linear in K), the
    // SAME field the M2 explicit-K eigenstress chain consumes.
    double const dr_dK =
        -dt_safe * exchange.drho_l_hat_dmu_lR * micro_potential.dmu_lR_dK;
    double const dn_l_dK = -dr_dK / dr_dn_l;  // n_l per (J/kg)
    return std::isfinite(dn_l_dK) ? dn_l_dK : 0.0;
}

template <int DisplacementDim>
inline void updateMicroscaleHydraulicState(
    StatefulData<DisplacementDim>& state_current,
    StatefulDataPrev<DisplacementDim> const& state_previous, double const p_cap_ip,
    double const rho_LR, double const mu, double const dt, double const t,
    MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
    PotentialExchangeLocalSolveContext const& local_context,
    std::optional<MicroPorosityParameters> const& micro_porosity_parameters,
    PotentialExchangeParameters const* const potential_exchange_parameters,
    MicroCeilingTraceTag const* const trace_tag = nullptr)
{
    auto& n_l = std::get<MicroWaterContent>(state_current);
    auto const n_l_prev = std::get<PrevState<MicroWaterContent>>(state_previous);
    auto& rho_lR = std::get<MicroLiquidDensity>(state_current);
    auto const rho_lR_prev = std::get<PrevState<MicroLiquidDensity>>(state_previous);

    double const n_l_prev_value = std::max(1e-16, **n_l_prev);
    *n_l = n_l_prev_value;

    if (!isPotentialExchangeEnabled(potential_exchange_parameters) ||
        !micro_porosity_parameters)
    {
        return;
    }

    auto const& potential_exchange_params = *potential_exchange_parameters;
    if (potential_exchange_params.local_nonlinear_solve_mode ==
        LocalNonlinearSolveMode::ScalarReferenceMassStorage)
    {
        auto const macro_potential = computeYoungLaplaceMacroPotential(
            -p_cap_ip, rho_LR, potential_exchange_params.pressure_tolerance);
        double const rho_lR_prev_value = std::max(1e-16, **rho_lR_prev);
        // REV-scale previous liquid apparent density: phi_m_prev * rho_lR_prev.
        // local_context.phi_m_prev = (1-phi_M_prev)*n_l_prev (hierarchical split).
        double const rho_l_prev = local_context.phi_m_prev * rho_lR_prev_value;
        // KKT micro-water ceiling (branch dsm_mass_conservation_v3_kkt_ceiling_
        // 2026-09-30; NOT adopted): complementarity local solve instead of the
        // projected Newton below. Default treatment = clamp -> this block is
        // skipped and everything below is the shipped code, bitwise.
        // DESIGN.md 3.4.
        if (isKktCeiling(potential_exchange_params))
        {
            auto const kkt = solveReferenceMassStorageKktState(
                n_l_prev_value, rho_l_prev, rho_lR_prev_value, dt, rho_LR,
                micro_porosity_parameters->mass_exchange_coefficient, mu,
                macro_potential, local_context, potential_exchange_params);
            // Status of the previous evaluation at this IP (read before it is
            // overwritten), for the flips counter and the parity counters.
            bool const active_before =
                *std::get<MicroCeilingStatus>(state_current) ==
                static_cast<double>(MicroCeilingKktStatus::Active);
            applyReferenceMassStorageLocalState<DisplacementDim>(
                state_current, state_previous, variables, variables_prev,
                rho_LR, local_context, potential_exchange_params, kkt.local);
            bool const active_now =
                kkt.status == MicroCeilingKktStatus::Active;
            *std::get<MicroCeilingStatus>(state_current) =
                static_cast<double>(kkt.status);
            *std::get<MicroCeilingMultiplier>(state_current) = kkt.multiplier;
            *std::get<MicroExchangeReceived>(state_current) =
                kkt.exchange_received;
            *std::get<MicroCeilingRejectedExchange>(state_current) =
                kkt.rejected_exchange;
            *std::get<MicroCeilingFlips>(state_current) +=
                (active_before != active_now) ? 1.0 : 0.0;
            *std::get<MicroCeilingNonMonotone>(state_current) +=
                kkt.status == MicroCeilingKktStatus::NonMonotone ? 1.0 : 0.0;
            *std::get<MicroCeilingPremise>(state_current) +=
                kkt.status == MicroCeilingKktStatus::PremiseViolated ? 1.0
                                                                     : 0.0;
            // DESIGN.md 3.2 step 7: the KKT rule chooses nothing at status 2 and
            // 3; the point is handled like a clamp point. Counted in the state
            // fields above; printed once per process and status, never hidden.
            if (kkt.status == MicroCeilingKktStatus::NonMonotone)
            {
                static std::once_flag once_non_monotone;
                std::call_once(
                    once_non_monotone,
                    []
                    {
                        WARN(
                            "micro_ceiling_treatment = kkt: at least one "
                            "evaluation has f(n_max) < 0 AND an interior root "
                            "of the micro residual (status NonMonotone): the "
                            "KKT rule chooses nothing there, the base "
                            "(clamp) solve and sink are used at that point. "
                            "Counts: output field micro_ceiling_nonmonotone.");
                    });
            }
            if (kkt.status == MicroCeilingKktStatus::PremiseViolated)
            {
                static std::once_flag once_premise;
                std::call_once(
                    once_premise,
                    []
                    {
                        WARN(
                            "micro_ceiling_treatment = kkt: at least one "
                            "evaluation could not be certified by the scan "
                            "(status PremiseViolated: f(n_floor) >= 0, a "
                            "non-finite value, or alpha_M not positive): the "
                            "base (clamp) solve and sink are used at that "
                            "point. Counts: output field "
                            "micro_ceiling_premise.");
                    });
            }
            updateCeilingIterationDiagnostics(
                *std::get<MicroCeilingEpsSeen>(state_current),
                *std::get<MicroCeilingIncLast>(state_current),
                *std::get<MicroCeilingAttemptT>(state_current),
                *std::get<MicroCeilingIncAlt>(state_current),
                *std::get<MicroCeilingIncSame>(state_current),
                local_context.volumetric_strain, t, active_now, active_before);
            if (trace_tag != nullptr)
            {
                writeMicroCeilingTraceLine(*trace_tag, kkt, local_context);
            }
            return;
        }
        auto const coupled_update = solveReferenceMassStorageCoupledState(
            n_l_prev_value, rho_l_prev, rho_lR_prev_value, dt,
            rho_LR, micro_porosity_parameters->mass_exchange_coefficient, mu,
            macro_potential, local_context, potential_exchange_params);
        applyReferenceMassStorageLocalState<DisplacementDim>(
            state_current, state_previous, variables, variables_prev, rho_LR, local_context, potential_exchange_params,
            coupled_update);
        // DIAGNOSTIC B (mass-strip A/B, 2026-09-30; Vinay R-03; off by default,
        // then this block is skipped and the state is bit-identical to A).
        // At an IP on the micro ceiling n_l = phi the stored MicroExchangeSource
        // is set to the ACTUAL micro storage rate, i.e. the water the micro
        // domain really gained over this step,
        //   rho_l_hat_B = (phi_m*rho_lR - phi_m_prev*rho_lR_prev)/dt
        //                 [kg/(m^3 s)], phi_m*rho_lR = rho_l as in the micro solve
        // (both signs). The macro pressure residual reads this state at clamped
        // IPs (assembleWithJacobian / assemble) instead of alpha_M(mu_LR - mu_lR).
        // Output consequence: in B/B' runs the VTU field micro_exchange_source
        // is the BOOKED sink at clamped IPs, not the potential-driven rho_hat.
        if (potential_exchange_params.ceiling_micro_storage_exchange &&
            dt > 0.0 && microWaterContentIsAtCeiling(coupled_update.n_l,
                                                       local_context.phi))
        {
            double const rho_l_now_B =
                *std::get<MicroPorosity>(state_current) *
                coupled_update.rho_lR;  // kg/m^3
            double const rho_l_prev_B =
                local_context.phi_m_prev * rho_lR_prev_value;  // kg/m^3
            double rho_hat_booked_B =
                (rho_l_now_B - rho_l_prev_B) / dt;  // kg/(m^3 s)
            // V2 (Q1, 2026-09-30): Eulerian micro balance per current bulk
            // volume, d(rho_l)/dt + rho_l*eps_dot = rho_hat (DERIVATION.md):
            // the water the micro actually gained includes the dilution by the
            // volume change. eps_dot = (eps_v - eps_v_prev)/dt [1/s]; off ->
            // bare storage rate (V1, bit-identical).
            if (potential_exchange_params.ceiling_micro_storage_includes_strain)
            {
                rho_hat_booked_B +=
                    rho_l_now_B *
                    (local_context.volumetric_strain -
                     local_context.volumetric_strain_prev) /
                    dt;  // kg/(m^3 s)
            }
            std::get<MicroExchangeSource>(state_current) =
                MicroExchangeSource{rho_hat_booked_B};
        }
        return;
    }

    auto const macro_potential = computeYoungLaplaceMacroPotential(
        -p_cap_ip, rho_LR, potential_exchange_params.pressure_tolerance);
    auto const n_l_update = solveImplicitMicroWaterContent(
        n_l_prev_value, dt, rho_LR,
        micro_porosity_parameters->mass_exchange_coefficient, mu,
        macro_potential, local_context, potential_exchange_params);

    *n_l = n_l_update.n_l;
    // Keep dsm_micromacro-mode rho_lR evolution consistent with the dsm_micromacro bridge:
    // rho_lR is updated from the active reduced micro EOS.
    *rho_lR = computeActiveMicroLiquidDensity(n_l_update.n_l, rho_LR,
                                                local_context, potential_exchange_params)
                  .rho_lR;

    auto& p_L_m = std::get<MicroPressure>(state_current);
    auto& S_L_m = std::get<MicroSaturation>(state_current);
    auto& rho_l_hat = std::get<MicroExchangeSource>(state_current);
    auto const compatibility_output =
        computeCompatibilityMicroHydraulicOutput(
            n_l_update.n_l, rho_LR, local_context, potential_exchange_params);
    *p_L_m = compatibility_output.p_L_m;
    *S_L_m = compatibility_output.S_L_m;
    rho_l_hat = MicroExchangeSource{n_l_update.exchange.rho_l_hat};
}

template <int DisplacementDim>
inline void updatePorositySplitState(
    StatefulData<DisplacementDim>& state_current,
    StatefulDataPrev<DisplacementDim> const& state_previous, double const phi,
    MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
    PotentialExchangeParameters const* const potential_exchange_parameters)
{
    if (!isPotentialExchangeEnabled(potential_exchange_parameters))
    {
        return;
    }

    auto& micro_porosity = std::get<MicroPorosity>(state_current);
    auto& transport_porosity =
        std::get<ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(state_current)
            .phi;
    auto const phi_M_prev = std::get<PrevState<
        ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(state_previous)
                                ->phi;
    auto const n_l = std::max(1e-16, *std::get<MicroWaterContent>(state_current));
    auto const phi_m_prev = **std::get<PrevState<MicroPorosity>>(state_previous);

    auto const transport_porosity_update =
        computeTransportPorosityUpdate(
            phi, phi_M_prev, phi_m_prev, n_l, variables.volumetric_strain,
            variables_prev.volumetric_strain,
            potential_exchange_parameters->macro_porosity_update_mode);

    *micro_porosity = transport_porosity_update.phi_m;
    transport_porosity = transport_porosity_update.phi_M;
    variables.transport_porosity = transport_porosity_update.phi_M;
    variables_prev.transport_porosity = transport_porosity_update.phi_M_prev;
}

template <int DisplacementDim>
inline void updateTotalPorosityState(
    StatefulData<DisplacementDim>& state_current,
    StatefulDataPrev<DisplacementDim> const& state_previous,
    double& phi, MPL::VariableArray& variables,
    MPL::VariableArray& variables_prev,
    PotentialExchangeParameters const* const potential_exchange_parameters)
{
    if (!isPotentialExchangeEnabled(potential_exchange_parameters))
    {
        return;
    }

    if (potential_exchange_parameters->local_nonlinear_solve_mode ==
            LocalNonlinearSolveMode::ScalarReferenceStorage ||
        potential_exchange_parameters->local_nonlinear_solve_mode ==
            LocalNonlinearSolveMode::ScalarReferenceMassStorage)
    {
        // In scalar dsm_micromacro-storage mode, micro porosity is support-state only.
        // Keep the process porosity state on the medium-law carrier.
        return;
    }

    auto const phi_m = *std::get<MicroPorosity>(state_current);
    auto const phi_m_prev = **std::get<PrevState<MicroPorosity>>(state_previous);
    auto const phi_M =
        std::get<ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(state_current)
            .phi;
    auto const phi_M_prev =
        std::get<PrevState<
            ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(
            state_previous)
            ->phi;

    auto& porosity =
        std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(state_current).phi;
    phi = phi_M + phi_m;
    porosity = phi;
    variables.porosity = phi;
    variables_prev.porosity = phi_M_prev + phi_m_prev;
}

template <int DisplacementDim>
inline MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
computeReferenceMicroPorositySwellingStressIncrement(
    double const n_l_prev, double const n_l,
    double const n_S, double const rho_lR, double const rho_lR_prev,
    double const rho_LR,
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C_el,
    PotentialExchangeParameters const& potential_exchange_params,
    double const biot_coefficient = 1.0,
    double const p_conf = std::numeric_limits<double>::quiet_NaN(),
    double const eps_v = std::numeric_limits<double>::quiet_NaN(),
    double const eps_v_prev = std::numeric_limits<double>::quiet_NaN(),
    // TOTAL porosity phi for live K(rho_d) (K_OF_RHO_D_LIVE.md); NaN sentinel
    // (callers without porosity in scope) -> parse-time scalar K.
    double const total_porosity = std::numeric_limits<double>::quiet_NaN())
{
    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    // Live K(rho_d): rho_d = rho_SR*(1-phi) [kg/m^3]; one K for BOTH the prev
    // and curr Pi evaluations of this increment (phi is the current state —
    // mirrors the held-fixed p_conf telescoping convention).
    double const K_aug_sw = effectiveAugmentationPrefactor(
        potential_exchange_params, total_porosity);  // K [J/kg]
    auto const& params = potential_exchange_params;

    // C_el is unused on the OFF and operational film branches (transmitted-
    // pressure form, no drained K needed); the EXACT route (H1) DOES use it for
    // the drained-line eigenstress half. Kept in the signature for call-site
    // stability.

    KV delta_sigma_sw = KV::Zero();
    double const delta_n_l = n_l - n_l_prev;
    if (!(std::isfinite(delta_n_l) &&
          std::abs(delta_n_l) > std::numeric_limits<double>::epsilon()))
    {
        return delta_sigma_sw;
    }

    // ── H1 (review 2026-06-14; RESIDUAL-CHANGING, Vinay-authorized) ──────────
    // Under film_energy_route = Exact, source the eigenstress half from the SAME
    // one-Psi functional whose mu_mech half is folded into mu_lR
    // (applyFilmPressureMicroPotential exact branch, L703), so the assembled
    // Maxwell pair dsigma_sw/dn_l == nS*rho_lR*dmu_mech/deps_v holds in the
    // residual (the §9a operational defect |W|/scale=0.93 the exact route exists
    // to cure is NOT cured if the eigenstress half stays operational). The pair
    // gives the drained-line LEVEL sigma_sw_m = -nS*n_l*(Pi(w_eff) + b*K_d*eps);
    // telescope to the step increment as the operational branch does:
    //   delta_sigma_sw = (sigma_sw_m_curr - sigma_sw_m_prev)*I.
    // PHYSICS TRADEOFF (predicted, §5): the pair's eigenstress is on the DRAINED
    // LINE p_conf = -K_d*eps_v, whereas the operational branch uses the ACTUAL
    // GP p_conf (held fixed across the step). On the drained line the two agree;
    // off it (e.g. fully confined, eps_v~0 with p_conf growing) they differ by
    // the off-line p_conf excursion. This is the deliberate one-Psi consistency
    // choice. Exact route requires Kinematic coupling (create-time validated,
    // mirrored at the mu fold L707).
    if (potential_exchange_params.film_pressure_coupling &&
        potential_exchange_params.film_strain_coupling !=
            FilmStrainCouplingMode::Off &&
        potential_exchange_params.film_energy_route == FilmEnergyRoute::Exact &&
        std::isfinite(eps_v))
    {
        auto const& params_h1 = potential_exchange_params;
        if (!(params_h1.hamaker_constant > 0.0) ||
            !(params_h1.specific_surface > 0.0) ||
            !(params_h1.micro_solid_density_reference > 0.0))
        {
            OGS_FATAL(
                "The exact-route DSM swelling stress requires positive vdW "
                "parameters: hamaker_constant > 0 (got {:g}), specific_surface "
                "> 0 (got {:g}) and micro_solid_density_reference > 0 (got "
                "{:g}).",
                params_h1.hamaker_constant, params_h1.specific_surface,
                params_h1.micro_solid_density_reference);
        }
        auto const& identity2_h1 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::kelvin_vector_dimensions(
                DisplacementDim)>::identity2;
        double const sign_h1 =
            microPotentialSignFactorFromParameters(params_h1);
        double const K_drained_h1 =
            drainedBulkModulusFromStiffness<DisplacementDim>(C_el);  // Pa
        double const eps_v_prev_h1 =
            std::isfinite(eps_v_prev) ? eps_v_prev : eps_v;
        double const rho_pi_prev_h1 =
            params_h1.use_micro_liquid_density_for_micro_pressure ? rho_lR_prev
                                                                  : rho_LR;
        double const rho_pi_curr_h1 =
            params_h1.use_micro_liquid_density_for_micro_pressure ? rho_lR
                                                                  : rho_LR;
        double const active_nS_prev_h1 = computeActiveMicroSolidVolumeFraction(
            n_l_prev, PotentialExchangeLocalSolveContext{}, params_h1);
        double const active_nS_curr_h1 = computeActiveMicroSolidVolumeFraction(
            n_l, PotentialExchangeLocalSolveContext{}, params_h1);
        double const kappa_prev_h1 =
            params_h1.film_strain_kappa == FilmStrainKappaMode::Aggregate
                ? active_nS_prev_h1
                : 1.0;
        double const kappa_curr_h1 =
            params_h1.film_strain_kappa == FilmStrainKappaMode::Aggregate
                ? active_nS_curr_h1
                : 1.0;
        double const sigma_sw_m_prev_h1 =
            computeStrainedFilmEnergyPair(
                n_l_prev, eps_v_prev_h1, kappa_prev_h1, biot_coefficient,
                K_drained_h1, true /*include_S, route R3*/, rho_pi_prev_h1,
                active_nS_prev_h1, params_h1.micro_solid_density_reference,
                params_h1.hamaker_constant, params_h1.specific_surface, sign_h1,
                K_aug_sw, params_h1.potential_augmentation_exponent,
                0.0 /*dnS_dnl: frozen nS (B1)*/,
                params_h1.micro_water_content_floor)
                .sigma_sw_m;  // Pa
        double const sigma_sw_m_curr_h1 =
            computeStrainedFilmEnergyPair(
                n_l, eps_v, kappa_curr_h1, biot_coefficient, K_drained_h1,
                true /*include_S, route R3*/, rho_pi_curr_h1, active_nS_curr_h1,
                params_h1.micro_solid_density_reference,
                params_h1.hamaker_constant, params_h1.specific_surface, sign_h1,
                K_aug_sw, params_h1.potential_augmentation_exponent,
                0.0 /*dnS_dnl: frozen nS (B1)*/,
                params_h1.micro_water_content_floor)
                .sigma_sw_m;  // Pa
        delta_sigma_sw.noalias() +=
            (sigma_sw_m_curr_h1 - sigma_sw_m_prev_h1) * identity2_h1;  // Pa
        return delta_sigma_sw;
    }

    // ── Film-pressure swelling stress (maxwell sec.5, flag ON) ──────────────
    // CORRECTION (Vinay, 2026-06-06): the micro swelling stress is a transmitted
    // PRESSURE over the contact fraction, weighted by the MICRO porosity, NOT an
    // elastic eigenstress. The previous K_sw*b*eps_sw^m form routed a pressure
    // through an elastic strain (the bulk modulus K_sw does NOT belong in a
    // residual stress); the chosen form is
    //   sigma_sw(n) = -phi_m * p_film,
    //   phi_m       = (1 - phi_M)*n_l = n_S*n_l   (MICRO porosity),
    //   p_film      = Pi(n_l) - b*p_conf,
    // tension-positive. Here Pi is the BARE van-der-Waals disjoining pressure at
    // n_l (BEFORE the film delta),
    //   Pi(n_l) = -rho_lR_used * mu_lR_vdw(n_l)   (FULL p^disj, density mirrored to
    //   the hydraulic p_L_m choice exactly as the OFF branch below),
    // b = biot_coefficient and p_conf = -tr(sigma_eff)/3 (confining pressure,
    // supplied by the caller as local_context.confining_pressure_p_conf).
    //
    // Increment over the step (telescoped EXACTLY as the OFF branch telescopes
    // -phi_m*Pi, but with p_film in place of Pi):
    //   delta_sigma_sw = n_S*(n_l_prev*p_film_prev - n_l*p_film_curr)*identity2.
    // p_conf_prev is NOT threaded here (only the current p_conf is), so p_conf is
    // HELD FIXED across the step (current p_conf used for both p_film terms): the
    // p_conf->stress coupling then enters only through the outer (Newton)
    // iteration as p_conf updates step to step. With the n_l weighting the
    // -b*p_conf term does NOT cancel (the weights n_l_prev, n_l differ):
    //   delta_sigma_sw = n_S*(n_l_prev*Pi_prev - n_l*Pi_curr)*identity2
    //                  - n_S*b*p_conf*(n_l_prev - n_l)*identity2.
    // When p_conf is the NaN sentinel (default callers / GP tests) the drain is
    // dropped (p_conf->0), leaving the pure micro-weighted Pi increment.
    // K_sw is GONE from this residual entirely. This makes sigma_sw compressive,
    // ~0 at dry (n_l->0), growing in magnitude with wetting, and proportional to
    // Pi (density-correct).
    if (potential_exchange_params.film_pressure_coupling)
    {
        // The bare-vdW disjoining law needs physical vdW parameters (mirrors the
        // OFF-branch up-front check so the message is swelling-law-specific).
        if (!(params.hamaker_constant > 0.0) ||
            !(params.specific_surface > 0.0) ||
            !(params.micro_solid_density_reference > 0.0))
        {
            OGS_FATAL(
                "The film-pressure DSM swelling stress requires positive vdW "
                "parameters: hamaker_constant > 0 (got {:g}), specific_surface > "
                "0 (got {:g}) and micro_solid_density_reference > 0 (got {:g}).",
                params.hamaker_constant, params.specific_surface,
                params.micro_solid_density_reference);
        }

        auto const& identity2_film = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::kelvin_vector_dimensions(
                DisplacementDim)>::identity2;

        // Bare van-der-Waals micro potential at n_l and n_l_prev, computed EXACTLY
        // as the OFF branch does (same active_nS, sign factor and augmentation
        // args) -> the BARE Pi, BEFORE any film delta.
        double const active_nS_prev_film = computeActiveMicroSolidVolumeFraction(
            n_l_prev, PotentialExchangeLocalSolveContext{}, params);
        double const active_nS_curr_film = computeActiveMicroSolidVolumeFraction(
            n_l, PotentialExchangeLocalSolveContext{}, params);
        double const sign_factor_film =
            microPotentialSignFactorFromParameters(params);

        // ── Strained-film modes (DSM/STRAINED_FILM_IMPLEMENTATION.md) ──────
        // Evaluate Pi at the SAME strained state w_eff the micro-potential
        // fold point uses, so both halves of Psi_film see one film thickness
        // (one-Psi consistency). The density fed to the strained state mirrors
        // the hydraulic p_L_m choice exactly like the Pi evaluation below.
        // eps_v NaN sentinel (callers without strain in scope) or mode Off ->
        // w_eval = n_l, bit-for-bit the frozen-geometry path. p_conf is HELD
        // FIXED across the step for BOTH states (mirrors the telescoping
        // convention for the -b*p_conf drain).
        double w_eval_prev = n_l_prev;
        double w_eval_curr = n_l;
        if (params.film_strain_coupling != FilmStrainCouplingMode::Off &&
            std::isfinite(eps_v))
        {
            double const rho_pi_prev =
                params.use_micro_liquid_density_for_micro_pressure ? rho_lR_prev
                                                                   : rho_LR;
            double const rho_pi_curr =
                params.use_micro_liquid_density_for_micro_pressure ? rho_lR
                                                                   : rho_LR;
            double const eps_v_prev_used =
                std::isfinite(eps_v_prev) ? eps_v_prev : eps_v;
            w_eval_prev =
                computeStrainedFilmState(
                    params.film_strain_coupling, params.film_strain_kappa,
                    n_l_prev, active_nS_prev_film, eps_v_prev_used, p_conf,
                    rho_lR_prev, params.micro_solid_density_reference,
                    params.hamaker_constant, params.specific_surface,
                    sign_factor_film, K_aug_sw /*live K(rho_d), J/kg*/,
                    params.potential_augmentation_exponent,
                    params.micro_water_content_floor, rho_pi_prev)
                    .w_eff;
            w_eval_curr =
                computeStrainedFilmState(
                    params.film_strain_coupling, params.film_strain_kappa, n_l,
                    active_nS_curr_film, eps_v, p_conf, rho_lR,
                    params.micro_solid_density_reference,
                    params.hamaker_constant, params.specific_surface,
                    sign_factor_film, K_aug_sw /*live K(rho_d), J/kg*/,
                    params.potential_augmentation_exponent,
                    params.micro_water_content_floor, rho_pi_curr)
                    .w_eff;
        }
        double const mu_lR_prev_film =
            computeVanDerWaalsMicroPotential(
                w_eval_prev, rho_lR_prev, active_nS_prev_film,
                params.micro_solid_density_reference, params.hamaker_constant,
                params.specific_surface, sign_factor_film,
                K_aug_sw /*live K(rho_d), J/kg*/,
                params.potential_augmentation_exponent, 0.0 /*dnS_dnl*/,
                params.micro_water_content_floor)
                .mu_lR;
        double const mu_lR_curr_film =
            computeVanDerWaalsMicroPotential(
                w_eval_curr, rho_lR, active_nS_curr_film,
                params.micro_solid_density_reference, params.hamaker_constant,
                params.specific_surface, sign_factor_film,
                K_aug_sw /*live K(rho_d), J/kg*/,
                params.potential_augmentation_exponent, 0.0 /*dnS_dnl*/,
                params.micro_water_content_floor)
                .mu_lR;

        // Density mirrors the hydraulic p_L_m choice EXACTLY (see OFF branch and
        // computeCompatibilityMicroHydraulicOutput): confined micro-liquid
        // density when enabled, bulk otherwise.
        double const p_L_m_density_prev_film =
            params.use_micro_liquid_density_for_micro_pressure ? rho_lR_prev
                                                               : rho_LR;
        double const p_L_m_density_curr_film =
            params.use_micro_liquid_density_for_micro_pressure ? rho_lR : rho_LR;

        double const Pi_prev_film = -p_L_m_density_prev_film * mu_lR_prev_film;
        double const Pi_curr_film = -p_L_m_density_curr_film * mu_lR_curr_film;

        // CORRECTION (Vinay, 2026-06-06): the previous form
        //   delta_sigma_sw = -(1 - phi_M)*(Pi_curr - Pi_prev)
        // was the increment of sigma_sw = -(1 - phi_M)*Pi: it (i) dropped the
        // MICRO-porosity weight n_l (so sigma_sw was the wrong magnitude / not
        // density-correct) and (ii) accumulated TENSILE during wetting (Pi_curr >
        // Pi_prev as it wets => -(Pi_curr - Pi_prev) < 0 increment, but starting
        // from a non-zero dry offset the running sum drifted tensile). Vinay's
        // chosen form is
        //   sigma_sw = -phi_m * p_film,   phi_m = (1 - phi_M)*n_l = n_S*n_l,
        //   p_film   = Pi - b*p_conf,
        // i.e. the film/disjoining pressure weighted by the MICRO porosity (NOT
        // just n_S), carrying the p_conf drain, NO bulk modulus. Telescope the
        // increment EXACTLY as the OFF branch telescopes -phi_m*Pi, but with
        // p_film in place of Pi:
        //   delta_sigma_sw = n_S*(n_l_prev*p_film_prev - n_l*p_film_curr)*I.
        // p_conf is HELD FIXED across the step (the current p_conf is used for
        // BOTH p_film_prev and p_film_curr): p_conf_prev is not threaded here, and
        // this is exactly the prior code's intent (the p_conf->stress coupling
        // then enters only through the outer Newton iteration as p_conf updates
        // step to step). Unlike the old (unweighted) form the -b*p_conf term does
        // NOT cancel here, because the n_l weights differ:
        //   delta_sigma_sw = n_S*(n_l_prev*Pi_prev - n_l*Pi_curr)*I
        //                  - n_S*b*p_conf*(n_l_prev - n_l)*I.
        // When p_conf is the NaN sentinel (default callers / GP tests with the
        // film flag but no confining stress) the p_conf drain is dropped so the
        // increment reduces to the pure-Pi micro-weighted form (finite).
        double const b_film = biot_coefficient;
        double const p_conf_film = std::isfinite(p_conf) ? p_conf : 0.0;
        double const p_film_prev = Pi_prev_film - b_film * p_conf_film;
        double const p_film_curr = Pi_curr_film - b_film * p_conf_film;

        // n_S passed in is the REV macro-solid fraction (1 - phi_M); phi_m = n_S*n_l.
        delta_sigma_sw.noalias() +=
            n_S * (n_l_prev * p_film_prev - n_l * p_film_curr) * identity2_film;
        return delta_sigma_sw;
    }

    // Fail loud: the full p^disj swelling law has no fallback branch, so the
    // vdW micro-potential parameters MUST be physical. (computeVanDerWaals-
    // MicroPotential would itself OGS_FATAL on these, but we check up front to
    // emit a swelling-law-specific message and to cover the n_S-reference
    // mode before the helper is ever reached.)
    if (!(params.hamaker_constant > 0.0) || !(params.specific_surface > 0.0) ||
        !(params.micro_solid_density_reference > 0.0))
    {
        OGS_FATAL(
            "The full-p^disj DSM swelling law requires positive vdW "
            "parameters: hamaker_constant > 0 (got {:g}), specific_surface > 0 "
            "(got {:g}) and micro_solid_density_reference > 0 (got {:g}).",
            params.hamaker_constant, params.specific_surface,
            params.micro_solid_density_reference);
    }
    if (params.micro_solid_volume_fraction_mode ==
            MicroSolidVolumeFractionMode::Reference &&
        !(params.micro_solid_volume_fraction_reference > 0.0))
    {
        OGS_FATAL(
            "The full-p^disj DSM swelling law with "
            "micro_solid_volume_fraction_mode='reference' requires "
            "micro_solid_volume_fraction_reference > 0, got {:g}.",
            params.micro_solid_volume_fraction_reference);
    }

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    // ONE unconditional swelling law (full micro disjoining pressure):
    //
    //   sigma_sw eigenstress increment
    //       = n_S * (n_l_prev * Pi_prev - n_l * Pi_curr) * identity2
    //
    // with, for tau in {prev, curr},
    //   Pi_tau           = -p_L_m_density_tau * mu_lR_tau   (FULL p^disj)
    //   mu_lR_tau        = computeVanDerWaalsMicroPotential(...).mu_lR
    //   p_L_m_density_tau = use_micro_liquid_density_for_micro_pressure
    //                           ? rho_lR_tau : rho_LR       (MIRRORS hydraulic)
    //
    // The vdW potential AND its exponential augmentation are BOTH carried
    // through unconditionally (adsorption potential, NOT a plate-plate term),
    // so both are swelling-promoting.
    //
    // Sign (SETTLED — do not change): negative_attractive => mu_lR < 0 =>
    //   Pi = -density * mu_lR > 0 => sigma_sw = -phi_m * Pi compressive
    //   (tension-positive convention), i.e. swelling.
    //   phi_m = (1 - phi_M) * n_l = n_S * n_l in the hierarchical split, so the
    //   eigenstress increment carries the explicit factor n_S * n_l. During
    //   hydration n_l*Pi increases, so (n_l_prev*Pi_prev - n_l*Pi_curr) < 0,
    //   giving a compressive (swelling) increment.
    //
    // NAMING NOTE: the n_S prefactor in scope here is the REV-scale solid
    // fraction (1 - phi_M), passed in from the caller. It is DISTINCT from the
    // aggregate-scale active_nS = computeActiveMicroSolidVolumeFraction(...)
    // used inside the vdW potential (the omega_l denominator). The
    // identification phi_m = n_S * n_l holds only for the REV-scale n_S here.

    // active_nS feeds the vdW potential's nS argument. The helper IGNORES the
    // PotentialExchangeLocalSolveContext entirely (Reference mode returns the
    // reference fraction; CurrentPorositySplit mode uses only n_l), so a
    // default-constructed context is the correct, intentional argument here.
    double const active_nS_prev = computeActiveMicroSolidVolumeFraction(
        n_l_prev, PotentialExchangeLocalSolveContext{}, params);
    double const active_nS_curr = computeActiveMicroSolidVolumeFraction(
        n_l, PotentialExchangeLocalSolveContext{}, params);

    double const sign_factor = microPotentialSignFactorFromParameters(params);

    double const mu_lR_prev =
        computeVanDerWaalsMicroPotential(
            n_l_prev, rho_lR_prev, active_nS_prev,
            params.micro_solid_density_reference, params.hamaker_constant,
            params.specific_surface, sign_factor,
            K_aug_sw /*live K(rho_d), J/kg*/,
            params.potential_augmentation_exponent, 0.0 /*dnS_dnl*/,
            params.micro_water_content_floor)
            .mu_lR;
    double const mu_lR_curr =
        computeVanDerWaalsMicroPotential(
            n_l, rho_lR, active_nS_curr,
            params.micro_solid_density_reference, params.hamaker_constant,
            params.specific_surface, sign_factor,
            K_aug_sw /*live K(rho_d), J/kg*/,
            params.potential_augmentation_exponent, 0.0 /*dnS_dnl*/,
            params.micro_water_content_floor)
            .mu_lR;

    // MIRROR the hydraulic p_L_m density choice exactly (see
    // computeCompatibilityMicroHydraulicOutput): confined micro-liquid density
    // when enabled, bulk density otherwise.
    double const p_L_m_density_prev =
        params.use_micro_liquid_density_for_micro_pressure ? rho_lR_prev
                                                           : rho_LR;
    double const p_L_m_density_curr =
        params.use_micro_liquid_density_for_micro_pressure ? rho_lR : rho_LR;

    double const Pi_prev = -p_L_m_density_prev * mu_lR_prev;
    double const Pi_curr = -p_L_m_density_curr * mu_lR_curr;

    delta_sigma_sw.noalias() +=
        n_S * (n_l_prev * Pi_prev - n_l * Pi_curr) * identity2;
    return delta_sigma_sw;
}

template <int DisplacementDim>
inline MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
computeSwellingStressIncrement(
    double const n_l_prev, double const n_l,
    double const n_S, double const rho_lR,
    double const rho_lR_prev, double const rho_LR,
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C_el,
    PotentialExchangeParameters const& potential_exchange_params,
    double const biot_coefficient = 1.0,
    double const p_conf = std::numeric_limits<double>::quiet_NaN(),
    double const eps_v = std::numeric_limits<double>::quiet_NaN(),
    double const eps_v_prev = std::numeric_limits<double>::quiet_NaN(),
    // TOTAL porosity phi for live K(rho_d) (K_OF_RHO_D_LIVE.md); NaN sentinel
    // -> parse-time scalar K.
    double const total_porosity = std::numeric_limits<double>::quiet_NaN())
{
    return computeReferenceMicroPorositySwellingStressIncrement<DisplacementDim>(
        n_l_prev, n_l, n_S, rho_lR, rho_lR_prev, rho_LR, C_el,
        potential_exchange_params, biot_coefficient, p_conf, eps_v,
        eps_v_prev, total_porosity);
}

template <int DisplacementDim>
inline void updateSwellingState(
    MaterialPropertyLib::Phase const& solid_phase,
    double const rho_LR,
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C_el,
    StatefulData<DisplacementDim>& state_current,
    StatefulDataPrev<DisplacementDim> const& state_previous,
    MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
    ParameterLib::SpatialPosition const& x_position, double const t,
    double const dt,
    PotentialExchangeParameters const* const potential_exchange_parameters,
    double const biot_coefficient = 1.0)
{
    if (!isPotentialExchangeEnabled(potential_exchange_parameters))
    {
        return;
    }

    auto const& potential_exchange_params = *potential_exchange_parameters;
    (void)solid_phase;
    (void)x_position;
    (void)t;
    (void)dt;

    auto const n_l_prev = **std::get<PrevState<MicroWaterContent>>(state_previous);
    auto const n_l = *std::get<MicroWaterContent>(state_current);
    auto const phi_M =
        std::get<ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(
            state_current).phi;
    double const n_S = std::max(1e-16, 1.0 - phi_M);
    double const rho_lR = *std::get<MicroLiquidDensity>(state_current);
    double const rho_lR_prev =
        **std::get<PrevState<MicroLiquidDensity>>(state_previous);

    auto& sigma_sw =
        std::get<ProcessLib::ThermoRichardsMechanics::
                     ConstitutiveStress_StrainTemperature::
                         SwellingDataStateful<DisplacementDim>>(state_current);
    auto const& sigma_sw_prev = std::get<
        PrevState<ProcessLib::ThermoRichardsMechanics::
                      ConstitutiveStress_StrainTemperature::
                          SwellingDataStateful<DisplacementDim>>>(state_previous);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    // Film-pressure coupling (maxwell sec.5): confining pressure
    // p_conf = -tr(sigma_eff)/3 (>0 in compression, OGS tension-positive) from the
    // CURRENT effective stress, threaded into the swelling stress p_film term.
    // Only when the flag is ON; otherwise the NaN sentinel keeps the term out and
    // the OFF (disjoining-eigenstress) branch is bit-for-bit unchanged. Mirrors
    // the p_conf_micro_solve computed at the production call sites.
    double const p_conf_swelling =
        isFilmPressureCouplingEnabled(potential_exchange_parameters)
            ? -std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                   DisplacementDim>>(state_current)
                   .sigma_eff.dot(identity2) /
                  3.0
            : std::numeric_limits<double>::quiet_NaN();

    sigma_sw = *sigma_sw_prev;
    sigma_sw.sigma_sw +=
        computeSwellingStressIncrement<DisplacementDim>(
            n_l_prev, n_l, n_S, rho_lR, rho_lR_prev, rho_LR, C_el,
            potential_exchange_params, biot_coefficient, p_conf_swelling,
            variables.volumetric_strain, variables_prev.volumetric_strain,
            // TOTAL porosity (live K(rho_d); rho_d = rho_SR*(1-phi)).
            std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                state_current)
                .phi);

    auto const C_el_inverse = C_el.inverse().eval();

    variables.volumetric_mechanical_strain =
        variables.volumetric_strain +
        identity2.transpose() * C_el_inverse * sigma_sw.sigma_sw;
    variables_prev.volumetric_mechanical_strain =
        variables_prev.volumetric_strain +
        identity2.transpose() * C_el_inverse * sigma_sw_prev->sigma_sw;
}

template <int DisplacementDim>
void updateSwellingStressAndVolumetricStrain(
    MaterialPropertyLib::Medium const& medium,
    MaterialPropertyLib::Phase const& solid_phase,
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C_el,
    double const rho_LR, double const mu,
    std::optional<MicroPorosityParameters> micro_porosity_parameters,
    PotentialExchangeParameters const* const potential_exchange_parameters,
    double const alpha, double const phi, double const p_cap_ip,
    MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
    ParameterLib::SpatialPosition const& x_position, double const t,
    double const dt,
    ProcessLib::ThermoRichardsMechanics::ConstitutiveStress_StrainTemperature::
        SwellingDataStateful<DisplacementDim>& sigma_sw,
    PrevState<ProcessLib::ThermoRichardsMechanics::
                  ConstitutiveStress_StrainTemperature::SwellingDataStateful<
                      DisplacementDim>> const& sigma_sw_prev,
    PrevState<ProcessLib::ThermoRichardsMechanics::TransportPorosityData> const
        phi_M_prev,
    PrevState<ProcessLib::ThermoRichardsMechanics::PorosityData> const phi_prev,
    ProcessLib::ThermoRichardsMechanics::TransportPorosityData& phi_M,
    PrevState<MicroPressure> const p_L_m_prev,
    PrevState<MicroSaturation> const S_L_m_prev, MicroPressure& p_L_m,
    MicroSaturation& S_L_m)
{
    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;
    bool const potential_exchange_enabled =
        isPotentialExchangeEnabled(potential_exchange_parameters);

    if (!medium.hasProperty(MPL::PropertyType::saturation_micro))
    {
        if (potential_exchange_enabled)
        {
            sigma_sw = *sigma_sw_prev;
            variables.volumetric_mechanical_strain =
                variables.volumetric_strain +
                identity2.transpose() * C_el.inverse() * sigma_sw.sigma_sw;
            variables_prev.volumetric_mechanical_strain =
                variables_prev.volumetric_strain + identity2.transpose() *
                                                       C_el.inverse() *
                                                       sigma_sw_prev->sigma_sw;
            return;
        }

        // If there is swelling, compute it. Update volumetric strain rate,
        // s.t. it corresponds to the mechanical part only.
        sigma_sw = *sigma_sw_prev;
        if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
        {
            auto const sigma_sw_dot =
                MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                    MPL::formEigenTensor<3>(
                        solid_phase[MPL::PropertyType::swelling_stress_rate]
                            .value(variables, variables_prev, x_position, t,
                                   dt)));
            sigma_sw.sigma_sw += sigma_sw_dot * dt;

            variables.volumetric_mechanical_strain =
                variables.volumetric_strain +
                identity2.transpose() * C_el.inverse() * sigma_sw.sigma_sw;
            variables_prev.volumetric_mechanical_strain =
                variables_prev.volumetric_strain + identity2.transpose() *
                                                       C_el.inverse() *
                                                       sigma_sw_prev->sigma_sw;
        }
        else
        {
            variables.volumetric_mechanical_strain =
                variables.volumetric_strain;
            variables_prev.volumetric_mechanical_strain =
                variables_prev.volumetric_strain;
        }
    }

    // TODO (naumov) saturation_micro must be always defined together with
    // the micro_porosity_parameters.
    if (medium.hasProperty(MPL::PropertyType::saturation_micro))
    {
        if (potential_exchange_enabled)
        {
            phi_M.phi = phi_M_prev->phi;
            // Prevent propagation of a non-physical negative macro porosity
            // from previous-step state into current assembly/output.
            phi_M.phi = std::max(0.0, phi_M.phi);
            variables_prev.transport_porosity = phi_M_prev->phi;
            variables.transport_porosity = phi_M.phi;

            *p_L_m = **p_L_m_prev;
            *S_L_m = **S_L_m_prev;
            sigma_sw = *sigma_sw_prev;

            variables.volumetric_mechanical_strain =
                variables.volumetric_strain +
                identity2.transpose() * C_el.inverse() * sigma_sw.sigma_sw;
            variables_prev.volumetric_mechanical_strain =
                variables_prev.volumetric_strain +
                identity2.transpose() * C_el.inverse() *
                    sigma_sw_prev->sigma_sw;
            return;
        }

        double const phi_m_prev = phi_prev->phi - phi_M_prev->phi;

        auto const [delta_phi_m, delta_e_sw, delta_p_L_m, delta_sigma_sw] =
            computeMicroPorosity<DisplacementDim>(
                identity2.transpose() * C_el.inverse(), rho_LR, mu,
                *micro_porosity_parameters, alpha, phi, -p_cap_ip, **p_L_m_prev,
                variables_prev, **S_L_m_prev, phi_m_prev, x_position, t, dt,
                medium.property(MPL::PropertyType::saturation_micro),
                solid_phase.property(MPL::PropertyType::swelling_stress_rate));

        phi_M.phi = phi - (phi_m_prev + delta_phi_m);
        variables_prev.transport_porosity = phi_M_prev->phi;
        variables.transport_porosity = phi_M.phi;

        *p_L_m = **p_L_m_prev + delta_p_L_m;
        {  // Update micro saturation.
            MPL::VariableArray variables_prev;
            variables_prev.capillary_pressure = -**p_L_m_prev;
            MPL::VariableArray variables;
            variables.capillary_pressure = -*p_L_m;

            *S_L_m = medium.property(MPL::PropertyType::saturation_micro)
                         .template value<double>(variables, x_position, t, dt);
        }
        sigma_sw.sigma_sw = sigma_sw_prev->sigma_sw + delta_sigma_sw;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                ShapeFunctionPressure, DisplacementDim>::
    RichardsMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        RichardsMechanicsProcessData<DisplacementDim>& process_data)
    : LocalAssemblerInterface<DisplacementDim>{
          e, integration_method, is_axially_symmetric, process_data}
{
    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    ip_data_.resize(n_integration_points);
    secondary_data_.N_u.resize(n_integration_points);

    auto const shape_matrices_u =
        NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                  ShapeMatricesTypeDisplacement,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   this->integration_method_);

    auto const shape_matrices_p =
        NumLib::initShapeMatrices<ShapeFunctionPressure,
                                  ShapeMatricesTypePressure, DisplacementDim>(
            e, is_axially_symmetric, this->integration_method_);

    auto const& medium =
        this->process_data_.media_map.getMedium(this->element_.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = ip_data_[ip];
        auto const& sm_u = shape_matrices_u[ip];
        ip_data_[ip].integration_weight =
            this->integration_method_.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        ParameterLib::SpatialPosition x_position = {
            std::nullopt, this->element_.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, ip_data.N_u))};

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

        // Initial porosity. Could be read from integration point data or mesh.
        auto& porosity =
            std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                this->current_states_[ip])
                .phi;
        porosity = medium->property(MPL::porosity)
                       .template initialValue<double>(
                           x_position,
                           std::numeric_limits<
                               double>::quiet_NaN() /* t independent */);

        auto& transport_porosity =
            std::get<
                ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(
                this->current_states_[ip])
                .phi;
        transport_porosity = porosity;
        if (medium->hasProperty(MPL::PropertyType::transport_porosity))
        {
            transport_porosity =
                medium->property(MPL::transport_porosity)
                    .template initialValue<double>(
                        x_position,
                        std::numeric_limits<
                            double>::quiet_NaN() /* t independent */);
        }

        secondary_data_.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    setInitialConditionsConcrete(Eigen::VectorXd const local_x,
                                 double const t,
                                 int const /*process_id*/)
{
    assert(local_x.size() == pressure_size + displacement_size);

    auto const [p_L, u] = localDOF(local_x);

    constexpr double dt = std::numeric_limits<double>::quiet_NaN();
    auto const& medium =
        this->process_data_.media_map.getMedium(this->element_.getID());
    MPL::VariableArray variables;

    auto const& solid_phase = medium->phase(MaterialPropertyLib::PhaseName::Solid);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N_p = ip_data_[ip].N_p;

        ParameterLib::SpatialPosition x_position = {
            std::nullopt, this->element_.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionPressure,
                                               ShapeMatricesTypePressure>(
                    this->element_, N_p))};

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;

        {
            auto& p_L_m = std::get<MicroPressure>(this->current_states_[ip]);
            auto& p_L_m_prev =
                std::get<PrevState<MicroPressure>>(this->prev_states_[ip]);
            **p_L_m_prev = -p_cap_ip;
            *p_L_m = -p_cap_ip;
        }

        auto const temperature =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(variables, x_position, t, dt);
        variables.temperature = temperature;

        auto& S_L_prev =
            std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::SaturationData>>(
                this->prev_states_[ip])
                ->S_L;
        S_L_prev = medium->property(MPL::PropertyType::saturation)
                       .template value<double>(variables, x_position, t, dt);

        if (this->process_data_.initial_stress.isTotalStress())
        {
            auto const alpha_b =
                medium->property(MPL::PropertyType::biot_coefficient)
                    .template value<double>(variables, x_position, t, dt);

            variables.liquid_saturation = S_L_prev;
            double const chi_S_L =
                medium->property(MPL::PropertyType::bishops_effective_stress)
                    .template value<double>(variables, x_position, t, dt);

            // Initial stresses are total stress, which were assigned to
            // sigma_eff in
            // RichardsMechanicsLocalAssembler::initializeConcrete().
            auto& sigma_eff =
                std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                    DisplacementDim>>(this->current_states_[ip]);

            auto& sigma_eff_prev =
                std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                       EffectiveStressData<DisplacementDim>>>(
                    this->prev_states_[ip]);

            // Reset sigma_eff to effective stress
            sigma_eff.sigma_eff.noalias() +=
                chi_S_L * alpha_b * (-p_cap_ip) * identity2;
            sigma_eff_prev->sigma_eff = sigma_eff.sigma_eff;
        }

        if (medium->hasProperty(MPL::PropertyType::saturation_micro))
        {
            MPL::VariableArray vars;
            vars.capillary_pressure = p_cap_ip;

            auto& S_L_m = std::get<MicroSaturation>(this->current_states_[ip]);
            auto& S_L_m_prev =
                std::get<PrevState<MicroSaturation>>(this->prev_states_[ip]);

            *S_L_m = medium->property(MPL::PropertyType::saturation_micro)
                         .template value<double>(vars, x_position, t, dt);
            *S_L_m_prev = S_L_m;
        }

        {
            auto& n_l = std::get<MicroWaterContent>(this->current_states_[ip]);
            auto& n_l_prev =
        std::get<PrevState<MicroWaterContent>>(this->prev_states_[ip]);
    auto& rho_lR =
        std::get<MicroLiquidDensity>(this->current_states_[ip]);
    auto& rho_lR_prev =
        std::get<PrevState<MicroLiquidDensity>>(this->prev_states_[ip]);
    auto& phi_m = std::get<MicroPorosity>(this->current_states_[ip]);
    auto& phi_m_prev =
        std::get<PrevState<MicroPorosity>>(this->prev_states_[ip]);

    // Default fallback keeps state positive for vdW algebra.
    double n_l_initial = 1e-6;
    double rho_lR_initial = 1.0;
            if (medium->hasProperty(MPL::PropertyType::saturation_micro))
            {
                auto const S_L_m_init =
                    *std::get<MicroSaturation>(this->current_states_[ip]);
                n_l_initial = std::max(1e-12, S_L_m_init);
            }

            auto const* const potential_exchange_params_ptr =
                this->getPotentialExchangeParameters();

            if (isPotentialExchangeEnabled(potential_exchange_params_ptr))
            {
                auto const porosity =
                    std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                        this->current_states_[ip])
                        .phi;
                auto const transport_porosity =
                    std::get<ProcessLib::ThermoRichardsMechanics::
                                 TransportPorosityData>(this->current_states_[ip])
                        .phi;
                double const porosity_safe = std::clamp(
                    std::max(0.0, porosity), 0.0, 1.0 - 1e-12);
                double const transport_porosity_safe = std::clamp(
                    std::max(0.0, transport_porosity), 0.0, porosity_safe);
                double const one_minus_phi_M =
                    std::max(1e-12, 1.0 - transport_porosity_safe);
                n_l_initial =
                    std::clamp((porosity_safe - transport_porosity_safe) /
                                   one_minus_phi_M,
                               1e-12, porosity_safe);

                n_l_initial =
                    std::clamp(
                        potential_exchange_params_ptr->initial_micro_water_content
                            .value_or(n_l_initial),
                        1e-12, porosity_safe);

                // Without a transport_porosity property the initial macro
                // porosity phi_M0 is the total porosity (constructor:
                // transport_porosity = porosity), so the previous micro
                // porosity stored below would be 0 and the first solve would
                // rebuild the micro water from empty pores. Take phi_M0 from
                // the hierarchical split of the declared
                // initial_micro_water_content instead (same formula as
                // computeTransportPorosityUpdate): phi_M0 = (phi - n_l0) /
                // (1 - n_l0), phi_m0 = (1 - phi_M0) * n_l0, phi_M0 + phi_m0 =
                // phi. Decks with a transport_porosity property, or without
                // initial_micro_water_content, are unchanged.
                if (potential_exchange_params_ptr->initial_micro_water_content
                        .has_value() &&
                    !medium->hasProperty(MPL::PropertyType::transport_porosity))
                {
                    std::get<ProcessLib::ThermoRichardsMechanics::
                                 TransportPorosityData>(
                        this->current_states_[ip])
                        .phi = std::clamp(
                        (porosity_safe - n_l_initial) /
                            std::max(1e-12, 1.0 - n_l_initial),
                        0.0, porosity_safe);
                }

                rho_lR_initial = std::max(
                    1e-16,
                    potential_exchange_params_ptr
                        ->micro_liquid_density_reference);
            }

            double phi_m_initial = std::max(1e-16, n_l_initial);
            if (isPotentialExchangeEnabled(potential_exchange_params_ptr))
            {
                auto const porosity_init_for_split =
                    std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                        this->current_states_[ip])
                        .phi;
                auto const transport_porosity_init_for_split =
                    std::get<ProcessLib::ThermoRichardsMechanics::
                                 TransportPorosityData>(this->current_states_[ip])
                        .phi;
                double const porosity_safe_for_split = std::clamp(
                    std::max(0.0, porosity_init_for_split), 0.0, 1.0 - 1e-12);
                double const transport_safe_for_split = std::clamp(
                    std::max(0.0, transport_porosity_init_for_split), 0.0,
                    porosity_safe_for_split);
                phi_m_initial = std::clamp(
                    (1.0 - transport_safe_for_split) * n_l_initial, 1e-16,
                    porosity_safe_for_split);
            }
            *n_l = n_l_initial;
            **n_l_prev = n_l_initial;
            *rho_lR = rho_lR_initial;
            **rho_lR_prev = rho_lR_initial;
            *phi_m = phi_m_initial;
            **phi_m_prev = phi_m_initial;

            if (isPotentialExchangeEnabled(potential_exchange_params_ptr))
            {
                auto const porosity =
                    std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                        this->current_states_[ip])
                        .phi;
                auto const transport_porosity_init =
                    std::get<ProcessLib::ThermoRichardsMechanics::
                                 TransportPorosityData>(this->current_states_[ip])
                        .phi;
                auto const rho_LR_initial =
                    medium->phase(MaterialPropertyLib::PhaseName::AqueousLiquid)
                        .property(MPL::PropertyType::density)
                        .template value<double>(variables, x_position, t, dt);
                auto const compatibility_output =
                    computeCompatibilityMicroHydraulicOutput(
                        n_l_initial, rho_LR_initial,
                        {.phi = porosity,
                         .phi_M_prev = transport_porosity_init,
                         .phi_m_prev = phi_m_initial,
                         .volumetric_strain = 0.0,
                         .volumetric_strain_prev = 0.0},
                        *potential_exchange_params_ptr);
                auto& p_L_m =
                    std::get<MicroPressure>(this->current_states_[ip]);
                auto& p_L_m_prev =
                    std::get<PrevState<MicroPressure>>(this->prev_states_[ip]);
                auto& S_L_m =
                    std::get<MicroSaturation>(this->current_states_[ip]);
                auto& S_L_m_prev =
                    std::get<PrevState<MicroSaturation>>(
                        this->prev_states_[ip]);
                *p_L_m = compatibility_output.p_L_m;
                **p_L_m_prev = compatibility_output.p_L_m;
                *S_L_m = compatibility_output.S_L_m;
                **S_L_m_prev = compatibility_output.S_L_m;

                auto& transport_porosity =
                    std::get<ProcessLib::ThermoRichardsMechanics::
                                 TransportPorosityData>(
                        this->current_states_[ip])
                        .phi;
                auto& transport_porosity_prev = std::get<PrevState<
                    ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(
                    this->prev_states_[ip]);
                auto const transport_porosity_update =
                    computeTransportPorosityUpdate(
                        porosity, transport_porosity, phi_m_initial, n_l_initial,
                        /*volumetric_strain=*/0.0,
                        /*volumetric_strain_prev=*/0.0,
                        potential_exchange_params_ptr
                            ->macro_porosity_update_mode);
                *phi_m = transport_porosity_update.phi_m;
                **phi_m_prev = transport_porosity_update.phi_m_prev;
                transport_porosity = transport_porosity_update.phi_M;
                transport_porosity_prev->phi =
                    transport_porosity_update.phi_M_prev;

                // Correct the micro liquid density initial state.
                // micro_liquid_density_reference (used above, line ~2015) is a
                // trivial EOS placeholder (e.g. 1e-6 kg/m³), NOT the physical
                // initial density.  In the first time step the exchange solve
                // updates rho_lR to ~rho_LR (~1000 kg/m³), so rho_lR_prev = 1e-6
                // while rho_lR = 1000.  When micro density enters the Pi-path
                // swelling stress (Pi = rho_lR * mu_lR), the density mismatch
                // Pi_prev = 1e-6 * K * exp(-xi_prev) ≈ 0 while
                // Pi_curr = 1000 * K * exp(-xi_curr) >> 0 produces a ~10^6×
                // tensile sigma_sw spike that permanently corrupts the accumulation.
                // Fix: initialise rho_lR and rho_lR_prev from the actual EOS at
                // the initial state so the first-step Pi difference is physical.
                {
                    PotentialExchangeLocalSolveContext const local_ctx_rho{
                        .phi = porosity,
                        .phi_M_prev = transport_porosity_update.phi_M_prev,
                        .phi_m_prev = transport_porosity_update.phi_m_prev,
                        .volumetric_strain = 0.0,
                        .volumetric_strain_prev = 0.0};
                    auto const rho_lR_data = computeActiveMicroLiquidDensity(
                        n_l_initial, rho_LR_initial, local_ctx_rho,
                        *potential_exchange_params_ptr);
                    double const rho_lR_corrected =
                        std::max(1e-16, rho_lR_data.rho_lR);
                    *rho_lR = rho_lR_corrected;
                    **rho_lR_prev = rho_lR_corrected;
                }
            }
        }

        // Set eps_m_prev from potentially non-zero eps and sigma_sw from
        // restart.
        auto& state_current = this->current_states_[ip];
        variables.stress =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                DisplacementDim>>(state_current)
                .sigma_eff;

        auto const& N_u = ip_data_[ip].N_u;
        auto const& dNdx_u = ip_data_[ip].dNdx_u;
        auto const x_coord =
            x_position.getCoordinates().value()[0];  // r for axisymetric
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, this->is_axially_symmetric_);
        auto& eps =
            std::get<StrainData<DisplacementDim>>(this->current_states_[ip])
                .eps;
        eps.noalias() = B * u;

        // Set mechanical strain temporary to compute tangent stiffness.
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);

        auto const C_el = ip_data_[ip].computeElasticTangentStiffness(
            variables, t, x_position, dt, this->solid_material_,
            *this->material_states_[ip].material_state_variables);

        auto const& sigma_sw =
            std::get<ProcessLib::ThermoRichardsMechanics::
                         ConstitutiveStress_StrainTemperature::
                             SwellingDataStateful<DisplacementDim>>(
                this->current_states_[ip])
                .sigma_sw;
        auto& eps_m_prev =
            std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                   MechanicalStrainData<DisplacementDim>>>(
                this->prev_states_[ip])
                ->eps_m;

        bool const swelling_stress_active =
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate) ||
            isPotentialExchangeEnabled(this->getPotentialExchangeParameters());
        eps_m_prev.noalias() =
            swelling_stress_active ? eps + C_el.inverse() * sigma_sw : eps;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::assemble(double const t, double const dt,
                               std::vector<double> const& local_x,
                               std::vector<double> const& local_x_prev,
                               std::vector<double>& local_M_data,
                               std::vector<double>& local_K_data,
                               std::vector<double>& local_rhs_data)
{
    assert(local_x.size() == pressure_size + displacement_size);

    auto const [p_L, u] = localDOF(local_x);
    auto const [p_L_prev, u_prev] = localDOF(local_x_prev);

    auto K = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_K_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto M = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_M_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size + pressure_size>>(
        local_rhs_data, displacement_size + pressure_size);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    auto const& medium =
        this->process_data_.media_map.getMedium(this->element_.getID());
    auto const& liquid_phase = medium->phase(MaterialPropertyLib::PhaseName::AqueousLiquid);
    auto const& solid_phase = medium->phase(MaterialPropertyLib::PhaseName::Solid);
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(this->element_.getID());

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = ip_data_[ip].integration_weight;

        auto const& N_u = ip_data_[ip].N_u;
        auto const& dNdx_u = ip_data_[ip].dNdx_u;

        auto const& N_p = ip_data_[ip].N_p;
        auto const& dNdx_p = ip_data_[ip].dNdx_p;

        x_position = {
            std::nullopt, this->element_.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, N_u))};
        auto const x_coord = x_position.getCoordinates().value()[0];

        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, this->is_axially_symmetric_);

        auto& eps =
            std::get<StrainData<DisplacementDim>>(this->current_states_[ip]);
        eps.eps.noalias() = B * u;

        auto& S_L =
            std::get<ProcessLib::ThermoRichardsMechanics::SaturationData>(
                this->current_states_[ip])
                .S_L;
        auto const S_L_prev =
            std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::SaturationData>>(
                this->prev_states_[ip])
                ->S_L;

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        double p_cap_prev_ip;
        NumLib::shapeFunctionInterpolate(-p_L_prev, N_p, p_cap_prev_ip);

        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;

        auto const temperature =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(variables, x_position, t, dt);
        variables.temperature = temperature;

        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);
        auto& state_current = this->current_states_[ip];
        variables.stress =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                DisplacementDim>>(state_current)
                .sigma_eff;
        // Set mechanical strain temporary to compute tangent stiffness.
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps.eps);
        auto const C_el = ip_data_[ip].computeElasticTangentStiffness(
            variables, t, x_position, dt, this->solid_material_,
            *this->material_states_[ip].material_state_variables);

        auto const beta_SR = (1 - alpha) / this->solid_material_.getBulkModulus(
                                               t, x_position, &C_el);
        variables.grain_compressibility = beta_SR;

        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        variables.density = rho_LR;
        auto const& b = this->process_data_.specific_body_force;

        S_L = medium->property(MPL::PropertyType::saturation)
                  .template value<double>(variables, x_position, t, dt);
        variables.liquid_saturation = S_L;
        variables_prev.liquid_saturation = S_L_prev;

        // tangent derivative for Jacobian
        double const dS_L_dp_cap =
            medium->property(MPL::PropertyType::saturation)
                .template dValue<double>(variables,
                                         MPL::Variable::capillary_pressure,
                                         x_position, t, dt);
        // secant derivative from time discretization for storage
        // use tangent, if secant is not available
        double const DeltaS_L_Deltap_cap =
            (p_cap_ip == p_cap_prev_ip)
                ? dS_L_dp_cap
                : (S_L - S_L_prev) / (p_cap_ip - p_cap_prev_ip);

        auto const chi = [medium, x_position, t, dt](double const S_L)
        {
            MPL::VariableArray vs;
            vs.liquid_saturation = S_L;
            return medium->property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(vs, x_position, t, dt);
        };
        double const chi_S_L = chi(S_L);
        double const chi_S_L_prev = chi(S_L_prev);

        double const p_FR = -chi_S_L * p_cap_ip;
        variables.effective_pore_pressure = p_FR;
        variables_prev.effective_pore_pressure = -chi_S_L_prev * p_cap_prev_ip;

        // Set volumetric strain rate for the general case without swelling.
        variables.volumetric_strain = Invariants::trace(eps.eps);
        variables_prev.volumetric_strain = Invariants::trace(B * u_prev);

        auto& phi = std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                        this->current_states_[ip])
                        .phi;
        {  // Porosity update
            auto const phi_prev = std::get<PrevState<
                ProcessLib::ThermoRichardsMechanics::PorosityData>>(
                                      this->prev_states_[ip])
                                      ->phi;
            variables_prev.porosity = phi_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables.porosity = phi;
        }

        if (alpha < phi)
        {
            OGS_FATAL(
                "RichardsMechanics: Biot-coefficient {} is smaller than "
                "porosity {} in element/integration point {}/{}.",
                alpha, phi, this->element_.getID(), ip);
        }

        // Swelling and possibly volumetric strain rate update.
        {
            auto& sigma_sw =
                std::get<ProcessLib::ThermoRichardsMechanics::
                             ConstitutiveStress_StrainTemperature::
                                 SwellingDataStateful<DisplacementDim>>(
                    this->current_states_[ip])
                    .sigma_sw;
            auto const& sigma_sw_prev = std::get<PrevState<
                ProcessLib::ThermoRichardsMechanics::
                    ConstitutiveStress_StrainTemperature::SwellingDataStateful<
                        DisplacementDim>>>(this->prev_states_[ip])
                                            ->sigma_sw;

            // If there is swelling, compute it. Update volumetric strain rate,
            // s.t. it corresponds to the mechanical part only.
            sigma_sw = sigma_sw_prev;
            if (solid_phase.hasProperty(
                    MPL::PropertyType::swelling_stress_rate))
            {
                auto const sigma_sw_dot =
                    MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                        MPL::formEigenTensor<3>(
                            solid_phase[MPL::PropertyType::swelling_stress_rate]
                                .value(variables, variables_prev, x_position, t,
                                       dt)));
                sigma_sw += sigma_sw_dot * dt;

                variables.volumetric_mechanical_strain =
                    variables.volumetric_strain +
                    identity2.transpose() * C_el.inverse() * sigma_sw;
                variables_prev.volumetric_mechanical_strain =
                    variables_prev.volumetric_strain +
                    identity2.transpose() * C_el.inverse() * sigma_sw_prev;
            }
            else
            {
                variables.volumetric_mechanical_strain =
                    variables.volumetric_strain;
                variables_prev.volumetric_mechanical_strain =
                    variables_prev.volumetric_strain;
            }

            if (medium->hasProperty(MPL::PropertyType::transport_porosity))
            {
                auto& transport_porosity =
                    std::get<ProcessLib::ThermoRichardsMechanics::
                                 TransportPorosityData>(
                        this->current_states_[ip])
                        .phi;
                auto const transport_porosity_prev =
                    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                                           TransportPorosityData>>(
                        this->prev_states_[ip])
                        ->phi;
                variables_prev.transport_porosity = transport_porosity_prev;

                transport_porosity =
                    medium->property(MPL::PropertyType::transport_porosity)
                        .template value<double>(variables, variables_prev,
                                                x_position, t, dt);
                variables.transport_porosity = transport_porosity;
            }
            else
            {
                variables.transport_porosity = phi;
            }
        }

        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);
        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);

        auto const& sigma_sw =
            std::get<ProcessLib::ThermoRichardsMechanics::
                         ConstitutiveStress_StrainTemperature::
                             SwellingDataStateful<DisplacementDim>>(
                this->current_states_[ip])
                .sigma_sw;
        auto const& sigma_eff =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                DisplacementDim>>(this->current_states_[ip])
                .sigma_eff;

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        {
            auto const sigma_total =
                (sigma_eff - alpha * p_FR * identity2).eval();

            // For stress dependent permeability.
            variables.total_stress.emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    sigma_total));
        }

        variables.equivalent_plastic_strain =
            this->material_states_[ip]
                .material_state_variables->getEquivalentPlasticStrain();

        auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        GlobalDimMatrixType const rho_K_over_mu =
            K_intrinsic * rho_LR * k_rel / mu;

        //
        // displacement equation, displacement part
        //
        {
            auto& eps_m = std::get<ProcessLib::ConstitutiveRelations::
                                       MechanicalStrainData<DisplacementDim>>(
                              this->current_states_[ip])
                              .eps_m;
            bool const swelling_stress_active =
                solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate) ||
                isPotentialExchangeEnabled(
                    this->getPotentialExchangeParameters());
            eps_m.noalias() = swelling_stress_active
                                  ? eps.eps + C_el.inverse() * sigma_sw
                                  : eps.eps;
            variables.mechanical_strain.emplace<
                MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);
        }

        {
            auto& state_current = this->current_states_[ip];
            auto const& state_previous = this->prev_states_[ip];
            auto& sigma_eff =
                std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                    DisplacementDim>>(state_current);
            auto const& sigma_eff_prev =
                std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                       EffectiveStressData<DisplacementDim>>>(
                    state_previous);
            auto const& eps_m =
                std::get<ProcessLib::ConstitutiveRelations::
                             MechanicalStrainData<DisplacementDim>>(state_current);
            auto& eps_m_prev =
                std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                       MechanicalStrainData<DisplacementDim>>>(
                    state_previous);

            auto const C = ip_data_[ip].updateConstitutiveRelation(
                variables, t, x_position, dt, temperature, sigma_eff,
                sigma_eff_prev, eps_m, eps_m_prev, this->solid_material_,
                this->material_states_[ip].material_state_variables);

            if (this->process_data_.use_numerical_jacobian)
            {
                K.template block<displacement_size, displacement_size>(
                     displacement_index, displacement_index)
                    .noalias() += B.transpose() * C * B * w;
            }
        }

        // p_SR
        variables.solid_grain_pressure =
            p_FR - sigma_eff.dot(identity2) / (3 * (1 - phi));
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);

        //
        // displacement equation, displacement part
        //
        double const rho = rho_SR * (1 - phi) + S_L * phi * rho_LR;
        rhs.template segment<displacement_size>(displacement_index).noalias() -=
            (B.transpose() * sigma_eff - N_u_op(N_u).transpose() * rho * b) * w;

        //
        // pressure equation, pressure part.
        //
        auto const beta_LR =
            1 / rho_LR *
            liquid_phase.property(MPL::PropertyType::density)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_phase_pressure,
                                         x_position, t, dt);

        double const a0 = S_L * (alpha - phi) * beta_SR;
        // DIAGNOSTIC B' (2026-09-30, off by default -> phi_storage == phi,
        // bit-identical): pore space of the macro water storage, phi_M instead
        // of the total porosity phi.
        double const phi_storage =
            (this->getPotentialExchangeParameters() &&
             this->getPotentialExchangeParameters()
                 ->macro_storage_uses_macro_porosity)
                ? std::get<ProcessLib::ThermoRichardsMechanics::
                               TransportPorosityData>(this->current_states_[ip])
                      .phi
                : phi;  // [-]
        // Volumetric average specific storage of the solid and fluid phases.
        double const specific_storage =
            DeltaS_L_Deltap_cap * (p_cap_ip * a0 - phi_storage) +
            S_L * (phi_storage * beta_LR + a0);
        M.template block<pressure_size, pressure_size>(pressure_index,
                                                       pressure_index)
            .noalias() += N_p.transpose() * rho_LR * specific_storage * N_p * w;

        K.template block<pressure_size, pressure_size>(pressure_index,
                                                       pressure_index)
            .noalias() += dNdx_p.transpose() * rho_K_over_mu * dNdx_p * w;

        rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_LR * rho_K_over_mu * b * w;

        auto const* const potential_exchange_params_ptr =
            this->getPotentialExchangeParameters();
        bool const potential_exchange_enabled =
            isPotentialExchangeEnabled(potential_exchange_params_ptr);
        if ((medium->hasProperty(MPL::PropertyType::saturation_micro) ||
             potential_exchange_enabled) &&
            this->process_data_.micro_porosity_parameters)
        {
            double const alpha_bar =
                this->process_data_.micro_porosity_parameters
                    ->mass_exchange_coefficient;
            auto const p_L_m =
                *std::get<MicroPressure>(this->current_states_[ip]);
            double const p_L_ip = -p_cap_ip;
            double const pressure_tolerance =
                getPotentialPressureTolerance(
                    potential_exchange_params_ptr);

            bool use_vdw_micro_potential_for_active_exchange = false;
            double mu_lR_vdw = 0.0;
            double dmu_lR_vdw_drho_lR = 0.0;
            double rho_lR_exchange_input =
                std::numeric_limits<double>::quiet_NaN();
            double drho_lR_exchange_input_dpL =
                std::numeric_limits<double>::quiet_NaN();
            // DIAGNOSTIC B (2026-09-30, off by default): residual-only twin of
            // the assembleWithJacobian block (see there).
            bool ceiling_B_active = false;
            double ceiling_B_rho_L_hat = 0.0;  // kg/(m^3 s), macro source
            // KKT micro-water ceiling (DESIGN.md 3.5, residual-only twin): reads
            // the state left by the last Newton-type evaluation (this path does
            // not run the micro update chain, exactly as variant B; nothing is
            // claimed for Picard, risk R4).
            bool kkt_active = false;
            double kkt_rho_L_hat = 0.0;  // kg/(m^3 s), macro source

            if (potential_exchange_enabled)
            {
                auto const n_l =
                    std::max(1e-16,
                             *std::get<MicroWaterContent>(
                                 this->current_states_[ip]));
                auto const transport_porosity_prev =
                    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                                           TransportPorosityData>>(
                        this->prev_states_[ip])
                        ->phi;
                auto const phi_m_prev =
                    **std::get<PrevState<MicroPorosity>>(this->prev_states_[ip]);
                // Confining pressure p_conf = -tr(sigma_eff)/3 (>0 in
                // compression). Threaded into the context ONLY when film-pressure
                // coupling is ON, so computeActiveMicroPotential folds the SAME
                // mu_lR(p_film) used by the n_l solve. Flag OFF -> NaN -> no film
                // term -> bit-for-bit identical to the pre-film code.
                bool const film_pressure_coupling =
                    potential_exchange_params_ptr->film_pressure_coupling;
                double const p_conf_assembly =
                    film_pressure_coupling
                        ? -std::get<ProcessLib::ConstitutiveRelations::
                                        EffectiveStressData<DisplacementDim>>(
                               this->current_states_[ip])
                               .sigma_eff.dot(identity2) /
                              3.0
                        : std::numeric_limits<double>::quiet_NaN();
                // Drained bulk modulus for the INTEGRABLE Maxwell partner
                // (mu_lR_mech), threaded so the macro-exchange mu_lR equals the
                // micro-solve mu_lR (equipresence). Only under film coupling; NaN
                // sentinel otherwise -> partner inert, flag-off bit-for-bit. C_el
                // here is the elastic stiffness evaluated above in this assemble()
                // ip loop (line ~2833).
                double const K_drained_assembly =
                    film_pressure_coupling
                        ? drainedBulkModulusFromStiffness<DisplacementDim>(C_el)
                        : std::numeric_limits<double>::quiet_NaN();
                PotentialExchangeLocalSolveContext const local_solve_context{
                    .phi = phi,
                    .phi_M_prev = transport_porosity_prev,
                    .phi_m_prev = phi_m_prev,
                    .volumetric_strain = variables.volumetric_strain,
                    .volumetric_strain_prev = variables_prev.volumetric_strain,
                    .confining_pressure_p_conf = p_conf_assembly,
                    .biot_coefficient = alpha,
                    .drained_bulk_modulus = K_drained_assembly};
                // DIAGNOSTIC B: rho_L_hat_B = -(booked micro gain) [kg/(m^3 s)]
                // at a clamped IP (see assembleWithJacobian / the micro state
                // update for the definition and the tangent).
                ceiling_B_active =
                    potential_exchange_params_ptr
                        ->ceiling_micro_storage_exchange &&
                    microWaterContentIsAtCeiling(n_l, phi);
                if (ceiling_B_active)
                {
                    ceiling_B_rho_L_hat = -*std::get<MicroExchangeSource>(
                        this->current_states_[ip]);  // kg/(m^3 s)
                }
                if (isKktCeiling(*potential_exchange_params_ptr))
                {
                    kkt_active =
                        *std::get<MicroCeilingStatus>(
                            this->current_states_[ip]) ==
                        static_cast<double>(MicroCeilingKktStatus::Active);
                    if (kkt_active)
                    {
                        kkt_rho_L_hat = -*std::get<MicroExchangeReceived>(
                            this->current_states_[ip]);  // kg/(m^3 s)
                    }
                }
                auto const micro_potential = computeActiveMicroPotential(
                    n_l, rho_LR, local_solve_context,
                    *potential_exchange_params_ptr);
                if (potential_exchange_params_ptr
                        ->use_micro_liquid_density_for_micro_pressure)
                {
                    auto const rho_lR_state =
                        *std::get<MicroLiquidDensity>(this->current_states_[ip]);
                    if (std::isfinite(rho_lR_state) && rho_lR_state > 0.0)
                    {
                        rho_lR_exchange_input = rho_lR_state;
                        drho_lR_exchange_input_dpL = rho_lR_state * beta_LR;
                    }
                }
                use_vdw_micro_potential_for_active_exchange = true;
                mu_lR_vdw = micro_potential.mu_lR;
                dmu_lR_vdw_drho_lR = micro_potential.dmu_lR_drho_lR;
                // ── DSM Maxwell-conjugate term (B1) ──────────────────────────
                // Restore the partner of the swelling eigenstress so (sigma,
                // mu_lR) come from one Psi. Sharp gate p'>=phi_m*Pi (opt.1); freeze phi (B1).
                // See ProcessLib/RichardsMechanics/DSM/
                // MAXWELL_CONJUGATE_IMPLEMENTATION.md. EXACTLY zero below the
                // gate, so gate-closed runs are unchanged bit-for-bit.
                //
                // Film-pressure coupling (maxwell sec.5): when the flag is ON,
                // computeActiveMicroPotential ALREADY folded the film term into
                // mu_lR above (the consolidated mu_lR(p_film) supersedes this
                // strain-view S1*eps_v block), so SKIP this inline block to
                // avoid double-counting. Flag OFF -> this block runs unchanged.
                if (!film_pressure_coupling)
                {
                    double const Pi_mc = p_L_m;  // Pa, disjoining = -rho*mu_lR > 0
                    // The INTEGRABLE partner takes eps_v + K_drained (not p_conf)
                    // and carries NO explicit (1-phi_M) factor (the n_S referencing
                    // cancels in the specific potential -- equipresence note in
                    // PotentialExchange.h), so p_conf / n_S / one-minus-n_l are no
                    // longer formed here; mirrors the L717 micro-solve fold.
                    double const rho_mc =
                        (std::isfinite(rho_lR_exchange_input) &&
                         rho_lR_exchange_input > 0.0)
                            ? rho_lR_exchange_input
                            : rho_LR;
                    // dPi/dn_l = (Pi/mu_lR)*dmu_lR_dnl (Pi=-rho*mu_lR; density-agnostic)
                    double const dPi_dnl_mc =
                        (std::abs(mu_lR_vdw) > 1e-300)
                            ? (Pi_mc / mu_lR_vdw) * micro_potential.dmu_lR_dnl
                            : 0.0;
                    // Pi'' = -rho*mu_lR''; density-agnostic (Pi=-rho*mu_lR) ->
                    // Pi''/Pi = mu_lR''/mu_lR, matching the dPi_dnl_mc convention.
                    double const d2Pi_dnl2_mc =
                        (std::abs(mu_lR_vdw) > 1e-300)
                            ? (Pi_mc / mu_lR_vdw) * micro_potential.d2mu_lR_dnl2
                            : 0.0;
                    // INTEGRABLE Maxwell mechanical partner (Vinay's Option-B,
                    // 2026-06-08): the SINGLE mu_lR_mech form, UNGATED (no Pi
                    // Heaviside) -- mirrors the L717 micro-solve fold so the
                    // macro-exchange mu_lR equals the micro-solve mu_lR
                    // (equipresence). PARAMETER-FREE: Pi/Pi'/Pi'' from the vdW
                    // potential, eps_v + b (=alpha) from the assembly, K_drained
                    // from the elastic stiffness C_el evaluated in this ip loop.
                    double const K_drained_mc =
                        drainedBulkModulusFromStiffness<DisplacementDim>(C_el);
                    auto const mc = computeIntegrableMechanicalMicroPotential(
                        Pi_mc, dPi_dnl_mc, d2Pi_dnl2_mc, n_l,
                        variables.volumetric_strain, /*biot_b=*/alpha,
                        K_drained_mc, rho_mc);
                    mu_lR_vdw += mc.mu_lR_mech;  // J/kg, additive (ungated)
                    dmu_lR_vdw_drho_lR += mc.dmu_lR_mech_drho_lR;
                }
            }

            auto const potential_exchange_result = computePotentialExchangeUpdate(
                alpha_bar, mu, p_L_ip, p_L_m, rho_LR, beta_LR,
                rho_lR_exchange_input, drho_lR_exchange_input_dpL,
                pressure_tolerance, potential_exchange_enabled,
                use_vdw_micro_potential_for_active_exchange, mu_lR_vdw,
                dmu_lR_vdw_drho_lR,
                /*use_custom_dmu_lR_vdw_dpL=*/false, /*dmu_lR_vdw_dpL=*/0.0,
                /*use_fd_jacobian_for_direct_macro_derivative=*/false,
                /*fd_jacobian_perturbation=*/1e-8);
            rhs.template segment<pressure_size>(pressure_index).noalias() +=
                N_p.transpose() *
                (ceiling_B_active
                     ? ceiling_B_rho_L_hat
                     : kkt_active
                           ? kkt_rho_L_hat
                           : potential_exchange_result.exchange.rho_L_hat) *
                w;
        }

        //
        // displacement equation, pressure part
        //
        K.template block<displacement_size, pressure_size>(displacement_index,
                                                           pressure_index)
            .noalias() -= B.transpose() * alpha * chi_S_L * identity2 * N_p * w;

        //
        // pressure equation, displacement part.
        //
        M.template block<pressure_size, displacement_size>(pressure_index,
                                                           displacement_index)
            .noalias() += N_p.transpose() * S_L * rho_LR * alpha *
                          identity2.transpose() * B * w;
    }

    if (this->process_data_.apply_mass_lumping)
    {
        auto pressure_mass_block_diag = M.template block<pressure_size, pressure_size>(
            pressure_index, pressure_index);
        pressure_mass_block_diag = pressure_mass_block_diag.colwise().sum().eval().asDiagonal();
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianEvalConstitutiveSetting(
        double const t, double const dt,
        ParameterLib::SpatialPosition const& x_position,
        RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure,
                                        DisplacementDim>::IpData& ip_data,
        MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
        MPL::Medium const* const medium, TemperatureData const T_data,
        CapillaryPressureData<DisplacementDim> const& p_cap_data,
        ConstitutiveData<DisplacementDim>& constitutive_data,
        StatefulData<DisplacementDim>& state_current,
        StatefulDataPrev<DisplacementDim> const& state_previous,
        OutputData<DisplacementDim>& OD,
        std::optional<MicroPorosityParameters> const& micro_porosity_parameters,
        PotentialExchangeParameters const* const
            potential_exchange_parameters,
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material,
        ProcessLib::ThermoRichardsMechanics::MaterialStateData<DisplacementDim>&
            material_state_data,
        MicroCeilingTraceTag const* const trace_tag)
{
    auto const& liquid_phase = medium->phase(MaterialPropertyLib::PhaseName::AqueousLiquid);
    auto const& solid_phase = medium->phase(MaterialPropertyLib::PhaseName::Solid);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    double const temperature = T_data();
    double const p_cap_ip = p_cap_data.p_cap;
    double const p_cap_prev_ip = p_cap_data.p_cap_prev;

    auto const& eps = std::get<StrainData<DisplacementDim>>(state_current);
    auto& S_L =
        std::get<ProcessLib::ThermoRichardsMechanics::SaturationData>(state_current).S_L;
    auto const S_L_prev =
        std::get<
            PrevState<ProcessLib::ThermoRichardsMechanics::SaturationData>>(
            state_previous)
            ->S_L;
    auto const alpha =
        medium->property(MPL::PropertyType::biot_coefficient)
            .template value<double>(variables, x_position, t, dt);
    *std::get<ProcessLib::ThermoRichardsMechanics::BiotData>(constitutive_data) = alpha;

    variables.stress =
        std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
            DisplacementDim>>(state_current)
            .sigma_eff;
    // Set mechanical strain temporary to compute tangent stiffness.
    variables.mechanical_strain
        .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
            eps.eps);
    auto const C_el = ip_data.computeElasticTangentStiffness(
        variables, t, x_position, dt, solid_material,
        *material_state_data.material_state_variables);

    auto const beta_SR =
        (1 - alpha) / solid_material.getBulkModulus(t, x_position, &C_el);
    variables.grain_compressibility = beta_SR;
    std::get<ProcessLib::ThermoRichardsMechanics::SolidCompressibilityData>(constitutive_data)
        .beta_SR = beta_SR;

    auto const rho_LR =
        liquid_phase.property(MPL::PropertyType::density)
            .template value<double>(variables, x_position, t, dt);
    variables.density = rho_LR;
    *std::get<LiquidDensity>(constitutive_data) = rho_LR;

    S_L = medium->property(MPL::PropertyType::saturation)
              .template value<double>(variables, x_position, t, dt);
    variables.liquid_saturation = S_L;
    variables_prev.liquid_saturation = S_L_prev;

    // tangent derivative for Jacobian
    double const dS_L_dp_cap =
        medium->property(MPL::PropertyType::saturation)
            .template dValue<double>(variables,
                                     MPL::Variable::capillary_pressure,
                                     x_position, t, dt);
    std::get<ProcessLib::ThermoRichardsMechanics::SaturationDataDeriv>(constitutive_data)
        .dS_L_dp_cap = dS_L_dp_cap;
    // secant derivative from time discretization for storage
    // use tangent, if secant is not available
    double const DeltaS_L_Deltap_cap =
        (p_cap_ip == p_cap_prev_ip)
            ? dS_L_dp_cap
            : (S_L - S_L_prev) / (p_cap_ip - p_cap_prev_ip);
    std::get<SaturationSecantDerivative>(constitutive_data).DeltaS_L_Deltap_cap =
        DeltaS_L_Deltap_cap;

    auto const chi = [medium, x_position, t, dt](double const S_L)
    {
        MPL::VariableArray vs;
        vs.liquid_saturation = S_L;
        return medium->property(MPL::PropertyType::bishops_effective_stress)
            .template value<double>(vs, x_position, t, dt);
    };
    double const chi_S_L = chi(S_L);
    std::get<ProcessLib::ThermoRichardsMechanics::BishopsData>(constitutive_data).chi_S_L =
        chi_S_L;
    double const chi_S_L_prev = chi(S_L_prev);
    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::BishopsData>>(constitutive_data)
        ->chi_S_L = chi_S_L_prev;

    auto const dchi_dS_L =
        medium->property(MPL::PropertyType::bishops_effective_stress)
            .template dValue<double>(
                variables, MPL::Variable::liquid_saturation, x_position, t, dt);
    std::get<ProcessLib::ThermoRichardsMechanics::BishopsData>(constitutive_data).dchi_dS_L =
        dchi_dS_L;

    double const p_FR = -chi_S_L * p_cap_ip;
    variables.effective_pore_pressure = p_FR;
    variables_prev.effective_pore_pressure = -chi_S_L_prev * p_cap_prev_ip;

    // Set volumetric strain rate for the general case without swelling.
    variables.volumetric_strain = Invariants::trace(eps.eps);
    // TODO (CL) changed that, using eps_prev for the moment, not B * u_prev
    // variables_prev.volumetric_strain = Invariants::trace(B * u_prev);
    variables_prev.volumetric_strain = Invariants::trace(
        std::get<PrevState<StrainData<DisplacementDim>>>(state_previous)->eps);

    auto& phi =
        std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(state_current).phi;
    {  // Porosity update
        auto const phi_prev =
            std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::PorosityData>>(
                state_previous)
                ->phi;
        variables_prev.porosity = phi_prev;
        phi = medium->property(MPL::PropertyType::porosity)
                  .template value<double>(variables, variables_prev, x_position,
                                          t, dt);
        variables.porosity = phi;
    }
    std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(constitutive_data).phi = phi;

    if (alpha < phi)
    {
        auto const eid =
            x_position.getElementID()
                ? static_cast<std::ptrdiff_t>(*x_position.getElementID())
                : static_cast<std::ptrdiff_t>(-1);
        OGS_FATAL(
            "RichardsMechanics: Biot-coefficient {} is smaller than porosity "
            "{} in element {}.",
            alpha, phi, eid);
    }

    auto const mu = liquid_phase.property(MPL::PropertyType::viscosity)
                        .template value<double>(variables, x_position, t, dt);
    *std::get<ProcessLib::ThermoRichardsMechanics::LiquidViscosityData>(constitutive_data) =
        mu;

    {
        // Swelling and possibly volumetric strain rate update.
        auto& sigma_sw =
            std::get<ProcessLib::ThermoRichardsMechanics::
                         ConstitutiveStress_StrainTemperature::
                             SwellingDataStateful<DisplacementDim>>(state_current);
        auto const& sigma_sw_prev =
            std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                                   ConstitutiveStress_StrainTemperature::
                                       SwellingDataStateful<DisplacementDim>>>(
                state_previous);
        auto const transport_porosity_prev = std::get<PrevState<
            ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(
            state_previous);
        auto const phi_prev = std::get<
            PrevState<ProcessLib::ThermoRichardsMechanics::PorosityData>>(
            state_previous);
        auto& transport_porosity = std::get<
            ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(state_current);
        auto& p_L_m = std::get<MicroPressure>(state_current);
        auto const p_L_m_prev = std::get<PrevState<MicroPressure>>(state_previous);
        auto& S_L_m = std::get<MicroSaturation>(state_current);
        auto const S_L_m_prev = std::get<PrevState<MicroSaturation>>(state_previous);

        updateSwellingStressAndVolumetricStrain<DisplacementDim>(
            *medium, solid_phase, C_el, rho_LR, mu, micro_porosity_parameters,
            potential_exchange_parameters, alpha, phi, p_cap_ip, variables,
            variables_prev, x_position, t, dt, sigma_sw, sigma_sw_prev,
            transport_porosity_prev, phi_prev, transport_porosity, p_L_m_prev,
            S_L_m_prev, p_L_m, S_L_m);
    }

    auto const transport_porosity_prev_value = std::get<PrevState<
        ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(state_previous)
                                                    ->phi;
    auto const phi_m_prev_value =
        **std::get<PrevState<MicroPorosity>>(state_previous);

    // Film-pressure coupling: supply the confining pressure p_conf =
    // -tr(sigma_eff)/3 so the n_l local solve sees the film term. Threaded only
    // when the flag is ON; otherwise the NaN sentinel keeps the term disabled.
    double const p_conf_micro_solve =
        isFilmPressureCouplingEnabled(potential_exchange_parameters)
            ? -std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                   DisplacementDim>>(state_current)
                   .sigma_eff.dot(identity2) /
                  3.0
            : std::numeric_limits<double>::quiet_NaN();
    // Drained bulk modulus for the INTEGRABLE Maxwell partner (mu_lR_mech). Only
    // under film coupling; NaN sentinel otherwise -> the partner is inert and the
    // flag-off micro solve is bit-for-bit unchanged.
    double const K_drained_micro_solve =
        isFilmPressureCouplingEnabled(potential_exchange_parameters)
            ? drainedBulkModulusFromStiffness<DisplacementDim>(C_el)
            : std::numeric_limits<double>::quiet_NaN();
    updateMicroscaleHydraulicState<DisplacementDim>(
        state_current, state_previous, p_cap_ip, rho_LR, mu, dt, t, variables, variables_prev,
        {.phi = phi,
         .phi_M_prev = transport_porosity_prev_value,
         .phi_m_prev = phi_m_prev_value,
         .volumetric_strain = variables.volumetric_strain,
         .volumetric_strain_prev = variables_prev.volumetric_strain,
         .confining_pressure_p_conf = p_conf_micro_solve,
         .biot_coefficient = alpha,
         .drained_bulk_modulus = K_drained_micro_solve},
        micro_porosity_parameters, potential_exchange_parameters, trace_tag);
    updatePorositySplitState<DisplacementDim>(
        state_current, state_previous, phi, variables, variables_prev,
        potential_exchange_parameters);
    updateTotalPorosityState<DisplacementDim>(
        state_current, state_previous, phi, variables, variables_prev,
        potential_exchange_parameters);
    std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(constitutive_data).phi =
        std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(state_current).phi;
    updateSwellingState<DisplacementDim>(
        solid_phase, rho_LR, C_el, state_current, state_previous, variables,
        variables_prev, x_position, t, dt,
        potential_exchange_parameters, /*biot_coefficient=*/alpha);

    // Gate 1/2 fix for DSM micro-porosity mode: enforce phi_m <= phi_total and
    // phi_M = phi_total - phi_m >= 0 in constitutive state.
    if (micro_porosity_parameters.has_value())
    {
        auto& phi_M_cs =
            std::get<ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(
                state_current)
                .phi;
        auto& phi_m_cs = *std::get<MicroPorosity>(state_current);
        double const phi_total_cs =
            std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                state_current)
                .phi;
        phi_m_cs = std::min(phi_m_cs, phi_total_cs);
        phi_M_cs = phi_total_cs - phi_m_cs;  // >= 0
        variables.transport_porosity = phi_M_cs;
        variables_prev.transport_porosity =
            std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                                   TransportPorosityData>>(state_previous)
                ->phi;
    }

    if (medium->hasProperty(MPL::PropertyType::transport_porosity))
    {
        if (!medium->hasProperty(MPL::PropertyType::saturation_micro) &&
            !isPotentialExchangeEnabled(potential_exchange_parameters))
        {
            auto& transport_porosity =
                std::get<
                    ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(
                    state_current)
                    .phi;
            auto const transport_porosity_prev = std::get<PrevState<
                ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(
                                                     state_previous)
                                                     ->phi;
            variables_prev.transport_porosity = transport_porosity_prev;

            transport_porosity =
                medium->property(MPL::PropertyType::transport_porosity)
                    .template value<double>(variables, variables_prev,
                                            x_position, t, dt);
            variables.transport_porosity = transport_porosity;
        }
        // Pi-path and legacy-DSM modes: phi_M already set by updatePorositySplitState /
        // updateSwellingStressAndVolumetricStrain; no property evaluation needed.
    }
    else
    {
        // No transport_porosity medium property. In Pi-path / DSM mode
        // variables.transport_porosity is already phi_M from updatePorositySplitState.
        // Only fall back to total porosity for plain RM (no micro-porosity split).
        if (!isPotentialExchangeEnabled(potential_exchange_parameters) &&
            !medium->hasProperty(MPL::PropertyType::saturation_micro))
        {
            variables.transport_porosity = phi;
        }
    }

    // Set mechanical variables for the intrinsic permeability model
    // For stress dependent permeability.
    {
        // TODO mechanical constitutive relation will be evaluated afterwards
        auto const sigma_total =
            (std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                 DisplacementDim>>(state_current)
                 .sigma_eff +
             alpha * p_FR * identity2)
                .eval();
        // For stress dependent permeability.
        variables.total_stress.emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_total));
    }

    variables.equivalent_plastic_strain =
        material_state_data.material_state_variables
            ->getEquivalentPlasticStrain();

    double const k_rel =
        medium->property(MPL::PropertyType::relative_permeability)
            .template value<double>(variables, x_position, t, dt);

    auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
        medium->property(MPL::PropertyType::permeability)
            .value(variables, x_position, t, dt));

    std::get<
        ProcessLib::ThermoRichardsMechanics::PermeabilityData<DisplacementDim>>(
        OD)
        .k_rel = k_rel;
    std::get<
        ProcessLib::ThermoRichardsMechanics::PermeabilityData<DisplacementDim>>(
        OD)
        .Ki = K_intrinsic;

    //
    // displacement equation, displacement part
    //

    {
        auto& sigma_sw =
            std::get<ProcessLib::ThermoRichardsMechanics::
                         ConstitutiveStress_StrainTemperature::
                             SwellingDataStateful<DisplacementDim>>(state_current)
                .sigma_sw;

        auto& eps_m =
            std::get<ProcessLib::ConstitutiveRelations::MechanicalStrainData<
                DisplacementDim>>(state_current)
                .eps_m;
        bool const swelling_stress_active =
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate) ||
            isPotentialExchangeEnabled(potential_exchange_parameters);
        eps_m.noalias() =
            swelling_stress_active ? eps.eps + C_el.inverse() * sigma_sw
                                   : eps.eps;
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);
    }

    {
        auto& sigma_eff =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                DisplacementDim>>(state_current);
        auto const& sigma_eff_prev =
            std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                   EffectiveStressData<DisplacementDim>>>(
                state_previous);
        auto const& eps_m =
            std::get<ProcessLib::ConstitutiveRelations::MechanicalStrainData<
                DisplacementDim>>(state_current);
        auto& eps_m_prev =
            std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                   MechanicalStrainData<DisplacementDim>>>(
                state_previous);

        auto C = ip_data.updateConstitutiveRelation(
            variables, t, x_position, dt, temperature, sigma_eff,
            sigma_eff_prev, eps_m, eps_m_prev, solid_material,
            material_state_data.material_state_variables);

        *std::get<StiffnessTensor<DisplacementDim>>(constitutive_data) = std::move(C);
    }

    // p_SR
    variables.solid_grain_pressure =
        p_FR - std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                   DisplacementDim>>(state_current)
                       .sigma_eff.dot(identity2) /
                   (3 * (1 - phi));
    auto const rho_SR =
        solid_phase.property(MPL::PropertyType::density)
            .template value<double>(variables, x_position, t, dt);

    double const rho = rho_SR * (1 - phi) + S_L * phi * rho_LR;
    *std::get<Density>(constitutive_data) = rho;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_x_prev,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    assert(local_x.size() == pressure_size + displacement_size);

    auto const [p_L, u] = localDOF(local_x);
    auto const [p_L_prev, u_prev] = localDOF(local_x_prev);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_Jac_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size + pressure_size>>(
        local_rhs_data, displacement_size + pressure_size);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    typename ShapeMatricesTypePressure::NodalMatrixType laplace_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p_a_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p_a_S_Jpp =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p_a_S =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, pressure_size>
        Kup = ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, pressure_size>::Zero(displacement_size,
                                                    pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        pressure_size, displacement_size>
        Kpu = ShapeMatricesTypeDisplacement::template MatrixType<
            pressure_size, displacement_size>::Zero(pressure_size,
                                                    displacement_size);

    // KKT micro-water ceiling (DESIGN.md 3.5): consistent exchange p-u entries of
    // the KKT-active IPs (and, with micro_ceiling_pu_tangent = all_exchange, the
    // accumulation target for the added Kpu/dt). Kept apart from local_Jac
    // because the assembly line `local_Jac.pu = Kpu/dt` at the end of this
    // function erases anything accumulated there (Q9). Zero unless a KKT-active
    // IP exists; used only when micro_ceiling_pu_tangent != overwritten.
    typename ShapeMatricesTypeDisplacement::template MatrixType<
        pressure_size, displacement_size>
        kkt_Kpu_exchange = ShapeMatricesTypeDisplacement::template MatrixType<
            pressure_size, displacement_size>::Zero(pressure_size,
                                                    displacement_size);

    auto const& medium =
        this->process_data_.media_map.getMedium(this->element_.getID());
    auto const& liquid_phase = medium->phase(MaterialPropertyLib::PhaseName::AqueousLiquid);
    auto const& solid_phase = medium->phase(MaterialPropertyLib::PhaseName::Solid);
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        ConstitutiveData<DisplacementDim> constitutive_data;
        auto& state_current = this->current_states_[ip];
        auto const& state_previous = this->prev_states_[ip];
        [[maybe_unused]] auto models = createConstitutiveModels(
            this->process_data_, this->solid_material_);

        auto const& w = ip_data_[ip].integration_weight;

        auto const& N_u = ip_data_[ip].N_u;
        auto const& dNdx_u = ip_data_[ip].dNdx_u;

        auto const& N_p = ip_data_[ip].N_p;
        auto const& dNdx_p = ip_data_[ip].dNdx_p;

        ParameterLib::SpatialPosition x_position = {
            std::nullopt, this->element_.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, N_u))};
        auto const x_coord = x_position.getCoordinates().value()[0];

        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, this->is_axially_symmetric_);

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        double p_cap_prev_ip;
        NumLib::shapeFunctionInterpolate(-p_L_prev, N_p, p_cap_prev_ip);

        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;

        auto const temperature =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(variables, x_position, t, dt);
        variables.temperature = temperature;

        std::get<StrainData<DisplacementDim>>(state_current).eps.noalias() = B * u;

        MicroCeilingTraceTag trace_tag_newton;
        assembleWithJacobianEvalConstitutiveSetting(
            t, dt, x_position, ip_data_[ip], variables, variables_prev, medium,
            TemperatureData{temperature},
            CapillaryPressureData<DisplacementDim>{
                p_cap_ip, p_cap_prev_ip,
                Eigen::Vector<double, DisplacementDim>::Zero()},
            constitutive_data, state_current, state_previous,
            this->output_data_[ip],
            this->process_data_.micro_porosity_parameters,
            this->getPotentialExchangeParameters(),
            this->solid_material_, this->material_states_[ip],
            // KKT iteration trace (DESIGN.md 3.8): tag only for the listed
            // elements.
            [&]() -> MicroCeilingTraceTag const*
            {
                auto const* const pep_trace =
                    this->getPotentialExchangeParameters();
                if (pep_trace == nullptr ||
                    pep_trace->micro_ceiling_trace_elements.empty() ||
                    std::find(pep_trace->micro_ceiling_trace_elements.begin(),
                              pep_trace->micro_ceiling_trace_elements.end(),
                              static_cast<std::size_t>(
                                  this->element_.getID())) ==
                        pep_trace->micro_ceiling_trace_elements.end())
                {
                    return nullptr;
                }
                trace_tag_newton = {
                    .element_id =
                        static_cast<std::size_t>(this->element_.getID()),
                    .integration_point = static_cast<std::size_t>(ip),
                    .t = t,
                    .dt = dt,
                    .output_reevaluation = false};
                return &trace_tag_newton;
            }());

        {
            auto const& C = *std::get<StiffnessTensor<DisplacementDim>>(constitutive_data);
            local_Jac
                .template block<displacement_size, displacement_size>(
                    displacement_index, displacement_index)
                .noalias() += B.transpose() * C * B * w;
        }

        auto const& b = this->process_data_.specific_body_force;

        {
            auto const& sigma_eff =
                std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                    DisplacementDim>>(this->current_states_[ip])
                    .sigma_eff;
            double const rho = *std::get<Density>(constitutive_data);
            local_rhs.template segment<displacement_size>(displacement_index)
                .noalias() -= (B.transpose() * sigma_eff -
                               N_u_op(N_u).transpose() * rho * b) *
                              w;
        }

        //
        // displacement equation, pressure part
        //

        double const alpha =
            *std::get<ProcessLib::ThermoRichardsMechanics::BiotData>(constitutive_data);
        double const dS_L_dp_cap =
            std::get<ProcessLib::ThermoRichardsMechanics::SaturationDataDeriv>(
                constitutive_data)
                .dS_L_dp_cap;

        {
            double const chi_S_L =
                std::get<ProcessLib::ThermoRichardsMechanics::BishopsData>(constitutive_data)
                    .chi_S_L;
            Kup.noalias() +=
                B.transpose() * alpha * chi_S_L * identity2 * N_p * w;
            double const dchi_dS_L =
                std::get<ProcessLib::ThermoRichardsMechanics::BishopsData>(constitutive_data)
                    .dchi_dS_L;

            local_Jac
                .template block<displacement_size, pressure_size>(
                    displacement_index, pressure_index)
                .noalias() -= B.transpose() * alpha *
                              (chi_S_L + dchi_dS_L * p_cap_ip * dS_L_dp_cap) *
                              identity2 * N_p * w;
        }

        double const phi =
            std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(constitutive_data).phi;
        double const rho_LR = *std::get<LiquidDensity>(constitutive_data);
        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() +=
            N_u_op(N_u).transpose() * phi * rho_LR * dS_L_dp_cap * b * N_p * w;

        // For the swelling stress with double structure model the corresponding
        // Jacobian u-p entry would be required, but it does not improve
        // convergence and sometimes worsens it:
        // if (medium->hasProperty(MPL::PropertyType::saturation_micro))
        // {
        //     -B.transpose() *
        //         dsigma_sw_dS_L_m* dS_L_m_dp_cap_m*(p_L_m - p_L_m_prev) /
        //         (p_cap_ip - p_cap_prev_ip) * N_p* w;
        // }
        if (!medium->hasProperty(MPL::PropertyType::saturation_micro) &&
            !isPotentialExchangeEnabled(
                this->getPotentialExchangeParameters()) &&
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
        {
            using DimMatrix = Eigen::Matrix<double, 3, 3>;
            auto const dsigma_sw_dS_L =
                MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                    solid_phase
                        .property(MPL::PropertyType::swelling_stress_rate)
                        .template dValue<DimMatrix>(
                            variables, variables_prev,
                            MPL::Variable::liquid_saturation, x_position, t,
                            dt));
            local_Jac
                .template block<displacement_size, pressure_size>(
                    displacement_index, pressure_index)
                .noalias() +=
                B.transpose() * dsigma_sw_dS_L * dS_L_dp_cap * N_p * w;
        }
        //
        // pressure equation, displacement part.
        //
        double const S_L =
            std::get<ProcessLib::ThermoRichardsMechanics::SaturationData>(
                this->current_states_[ip])
                .S_L;
        // OPEN (code review 2026-09-30, NOT a change; Vinay's ruling needed, see
        // ~/ogs-models/scratch/2026-09-30_massfix_V2/DERIVATION.md section 6):
        // Kpu below is the Biot volume-change storage S_L*rho_LR*alpha*div(u_dot)
        // of the WHOLE pore space. With alpha = 1 and beta_SR = 0 the porosity law
        // PorosityFromMassBalance gives alpha*dEps = dphi + phi*dEps, and
        // phi = phi_M + phi_m, so the term contains the micro part
        // S_L*rho_LR*(dphi_m + phi_m*dEps)/dt. macro_storage_uses_macro_porosity
        // (B') moves only a_p and a_S to phi_M and leaves this term unchanged:
        // the variant is 'phi_M in a_p/a_S only', not 'macro storage = macro pores
        // only'. Not removed here (formulation decision, CLAUDE.md section 9).
        if (this->process_data_.explicit_hm_coupling_in_unsaturated_zone)
        {
            double const chi_S_L_prev = std::get<PrevState<
                ProcessLib::ThermoRichardsMechanics::BishopsData>>(constitutive_data)
                                            ->chi_S_L;
            Kpu.noalias() += N_p.transpose() * chi_S_L_prev * rho_LR * alpha *
                             identity2.transpose() * B * w;
        }
        else
        {
            Kpu.noalias() += N_p.transpose() * S_L * rho_LR * alpha *
                             identity2.transpose() * B * w;
        }

        //
        // pressure equation, pressure part.
        //

        double const k_rel =
            std::get<ProcessLib::ThermoRichardsMechanics::PermeabilityData<
                DisplacementDim>>(this->output_data_[ip])
                .k_rel;
        auto const& K_intrinsic =
            std::get<ProcessLib::ThermoRichardsMechanics::PermeabilityData<
                DisplacementDim>>(this->output_data_[ip])
                .Ki;
        double const mu =
            *std::get<ProcessLib::ThermoRichardsMechanics::LiquidViscosityData>(
                constitutive_data);

        GlobalDimMatrixType const rho_Ki_over_mu = K_intrinsic * rho_LR / mu;

        laplace_p.noalias() +=
            dNdx_p.transpose() * k_rel * rho_Ki_over_mu * dNdx_p * w;

        auto const beta_LR =
            1 / rho_LR *
            liquid_phase.property(MPL::PropertyType::density)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_phase_pressure,
                                         x_position, t, dt);

        double const beta_SR =
            std::get<
                ProcessLib::ThermoRichardsMechanics::SolidCompressibilityData>(
                constitutive_data)
                .beta_SR;
        double const a0 = (alpha - phi) * beta_SR;
        // DIAGNOSTIC B' (2026-09-30, off by default -> phi_storage == phi,
        // bit-identical): pore space of the macro water storage, phi_M instead
        // of the total porosity phi, in the pore-fluid coefficients a_p, a_S
        // (and dS_L-derivative of a_p). a0 (skeleton/grain, alpha - phi) untouched.
        double const phi_storage =
            (this->getPotentialExchangeParameters() &&
             this->getPotentialExchangeParameters()
                 ->macro_storage_uses_macro_porosity)
                ? std::get<ProcessLib::ThermoRichardsMechanics::
                               TransportPorosityData>(this->current_states_[ip])
                      .phi
                : phi;  // [-]
        double const specific_storage_a_p =
            S_L * (phi_storage * beta_LR + S_L * a0);
        double const specific_storage_a_S = phi_storage - p_cap_ip * S_L * a0;

        double const dspecific_storage_a_p_dp_cap =
            dS_L_dp_cap * (phi_storage * beta_LR + 2 * S_L * a0);
        double const dspecific_storage_a_S_dp_cap =
            -a0 * (S_L + p_cap_ip * dS_L_dp_cap);

        storage_p_a_p.noalias() +=
            N_p.transpose() * rho_LR * specific_storage_a_p * N_p * w;

        double const DeltaS_L_Deltap_cap =
            std::get<SaturationSecantDerivative>(constitutive_data).DeltaS_L_Deltap_cap;
        storage_p_a_S.noalias() -= N_p.transpose() * rho_LR *
                                   specific_storage_a_S * DeltaS_L_Deltap_cap *
                                   N_p * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += N_p.transpose() * (p_cap_ip - p_cap_prev_ip) / dt *
                          rho_LR * dspecific_storage_a_p_dp_cap * N_p * w;

        double const S_L_prev =
            std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::SaturationData>>(
                this->prev_states_[ip])
                ->S_L;
        storage_p_a_S_Jpp.noalias() -=
            N_p.transpose() * rho_LR *
            ((S_L - S_L_prev) * dspecific_storage_a_S_dp_cap +
             specific_storage_a_S * dS_L_dp_cap) /
            dt * N_p * w;

        if (!this->process_data_.explicit_hm_coupling_in_unsaturated_zone)
        {
            local_Jac
                .template block<pressure_size, pressure_size>(pressure_index,
                                                              pressure_index)
                .noalias() -= N_p.transpose() * rho_LR * dS_L_dp_cap * alpha *
                              identity2.transpose() * B * (u - u_prev) / dt *
                              N_p * w;
        }

        double const dk_rel_dS_l =
            medium->property(MPL::PropertyType::relative_permeability)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_saturation,
                                         x_position, t, dt);
        typename ShapeMatricesTypeDisplacement::GlobalDimVectorType const
            grad_p_cap = -dNdx_p * p_L;
        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += dNdx_p.transpose() * rho_Ki_over_mu * grad_p_cap *
                          dk_rel_dS_l * dS_L_dp_cap * N_p * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += dNdx_p.transpose() * rho_LR * rho_Ki_over_mu * b *
                          dk_rel_dS_l * dS_L_dp_cap * N_p * w;

        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_LR * k_rel * rho_Ki_over_mu * b * w;

        auto const* const potential_exchange_params_ptr =
            this->getPotentialExchangeParameters();
        bool const potential_exchange_enabled =
            isPotentialExchangeEnabled(potential_exchange_params_ptr);
        if ((medium->hasProperty(MPL::PropertyType::saturation_micro) ||
             potential_exchange_enabled) &&
            this->process_data_.micro_porosity_parameters)
        {
            double const alpha_bar =
                this->process_data_.micro_porosity_parameters
                    ->mass_exchange_coefficient;
            auto const p_L_m =
                *std::get<MicroPressure>(this->current_states_[ip]);
            double const p_L_ip = -p_cap_ip;
            double const pressure_tolerance =
                getPotentialPressureTolerance(
                    potential_exchange_params_ptr);

            bool use_vdw_micro_potential_for_active_exchange = false;
            double mu_lR_vdw = 0.0;
            double dmu_lR_vdw_drho_lR = 0.0;
            double rho_lR_exchange_input =
                std::numeric_limits<double>::quiet_NaN();
            double drho_lR_exchange_input_dpL =
                std::numeric_limits<double>::quiet_NaN();
            bool use_custom_dmu_lR_vdw_dpL = false;
            double dmu_lR_vdw_dpL = 0.0;
            bool use_fd_jacobian_for_direct_macro_derivative = false;
            double fd_jacobian_perturbation = 1e-8;
            // DIAGNOSTIC B (mass-strip A/B, 2026-09-30; Vinay R-03; off by
            // default). At IPs on the micro ceiling n_l = phi the macro residual
            // books the ACTUAL micro storage rate instead of rho_hat.
            bool ceiling_B_active = false;
            double ceiling_B_rho_L_hat = 0.0;       // kg/(m^3 s), macro source
            double ceiling_B_drho_L_hat_dpL = 0.0;  // kg/(m^3 s)/Pa
            double ceiling_B_drho_L_hat_deps_v = 0.0;  // kg/(m^3 s) per unit eps_v
            // KKT micro-water ceiling (DESIGN.md 3.5): active flag read from the
            // state MicroCeilingStatus (no detector, no tolerance), macro
            // source convention rho_L_hat = -rhohat.
            bool kkt_active = false;
            double kkt_rho_L_hat = 0.0;            // kg/(m^3 s), macro source
            double kkt_drho_L_hat_dpL = 0.0;       // kg/(m^3 s)/Pa
            double kkt_drho_L_hat_deps_v = 0.0;    // kg/(m^3 s) per unit eps_v
            if (potential_exchange_enabled)
            {
                auto const n_l =
                    std::max(1e-16,
                             *std::get<MicroWaterContent>(
                                 this->current_states_[ip]));
                auto const transport_porosity_prev =
                    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                                           TransportPorosityData>>(
                        this->prev_states_[ip])
                        ->phi;
                auto const n_l_prev = **std::get<PrevState<MicroWaterContent>>(
                    this->prev_states_[ip]);
                auto const phi_m_prev =
                    **std::get<PrevState<MicroPorosity>>(this->prev_states_[ip]);
                // Film-pressure coupling (maxwell sec.5): thread p_conf into the
                // context so computeActiveMicroPotential folds the SAME
                // mu_lR(p_film) the residual path uses. Flag OFF -> NaN -> no
                // film term -> bit-for-bit identical Jacobian.
                bool const film_pressure_coupling =
                    potential_exchange_params_ptr->film_pressure_coupling;
                double const p_conf_assembly =
                    film_pressure_coupling
                        ? -std::get<ProcessLib::ConstitutiveRelations::
                                        EffectiveStressData<DisplacementDim>>(
                               this->current_states_[ip])
                               .sigma_eff.dot(identity2) /
                              3.0
                        : std::numeric_limits<double>::quiet_NaN();
                // Drained bulk modulus for the INTEGRABLE Maxwell partner
                // (mu_lR_mech), threaded so the macro-exchange mu_lR equals the
                // micro-solve mu_lR (equipresence). Only under film coupling; NaN
                // sentinel otherwise -> partner inert, flag-off bit-for-bit.
                double const K_drained_assembly =
                    film_pressure_coupling
                        ? drainedBulkModulusFromStiffness<DisplacementDim>(
                              ip_data_[ip].computeElasticTangentStiffness(
                                  variables, t, x_position, dt,
                                  this->solid_material_,
                                  *this->material_states_[ip]
                                       .material_state_variables))
                        : std::numeric_limits<double>::quiet_NaN();
                PotentialExchangeLocalSolveContext const local_solve_context{
                    .phi = phi,
                    .phi_M_prev = transport_porosity_prev,
                    .phi_m_prev = phi_m_prev,
                    .volumetric_strain = variables.volumetric_strain,
                    .volumetric_strain_prev = variables_prev.volumetric_strain,
                    .confining_pressure_p_conf = p_conf_assembly,
                    .biot_coefficient = alpha,
                    .drained_bulk_modulus = K_drained_assembly};
                // ── DIAGNOSTIC B (2026-09-30) ────────────────────────────────
                // Actual micro storage rate at a clamped IP, with its tangent.
                //   rho_l = phi_m * rho_lR   [kg/m^3 REV], same definition as the
                //   micro solve (rho_l_prev = phi_m_prev * rho_lR_prev).
                //   rho_L_hat_B = -(rho_l - rho_l_prev)/dt   [kg/(m^3 s)], read
                //   from the MicroExchangeSource state that the micro update
                //   booked at this IP (macro source = minus the micro gain; both
                //   directions).
                // Tangent at the clamped state (n_l = phi, phi_m = phi):
                //   rho_l(eps_v, p_L) = phi * rho_lR(n_l = phi, rho_LR(p_L)),
                //   d rho_l/d eps_v = dphi/deps_v * (rho_lR + phi*drho_lR/dn_l),
                //   d rho_l/d p_L   = phi * (d rho_lR/d rho_LR)*rho_LR*beta_LR,
                //   dphi/deps_v as PorosityFromMassBalance (clamp -> 0), the same
                //   law as the live-K p-u block below. The dphi/dp_eff*beta_SR
                //   chain is omitted (beta_SR is checked to be 0 below).
                ceiling_B_active =
                    potential_exchange_params_ptr
                        ->ceiling_micro_storage_exchange &&
                    microWaterContentIsAtCeiling(n_l, phi);
                if (ceiling_B_active)
                {
                    if (beta_SR != 0.0)
                    {
                        OGS_FATAL(
                            "ceiling_micro_storage_exchange: beta_SR = {:g} != 0; "
                            "the dphi/dp_eff chain of the clamped-IP tangent is "
                            "not implemented.",
                            beta_SR);
                    }
                    double const rho_lR_now =
                        *std::get<MicroLiquidDensity>(this->current_states_[ip]);
                    double const phi_m_now =
                        *std::get<MicroPorosity>(this->current_states_[ip]);
                    // Booked micro gain, set by the micro state update at this
                    // clamped IP (see updateMicroscaleHydraulicState); the macro
                    // source is its negative.
                    ceiling_B_rho_L_hat = -*std::get<MicroExchangeSource>(
                        this->current_states_[ip]);  // kg/(m^3 s)

                    auto const rho_lR_eos = computeActiveMicroLiquidDensity(
                        n_l, rho_LR, local_solve_context,
                        *potential_exchange_params_ptr);
                    double const drho_LR_dpL_B = rho_LR * beta_LR;  // kg/m^3/Pa
                    // d rho_lR/d rho_LR = 1 (rho_lR = rho_LR + rho_l0*exp(-a*omega),
                    // a*omega ~ 1e-16 -> the exp factor is constant to round-off).
                    // V2 (Q1): with the strain term S_B = [(rho_l - rho_l_prev)
                    // + rho_l*d_eps]/dt, dS_B/d(.) = [d rho_l/d(.)*(1 + d_eps)
                    // (+ rho_l for eps_v)]/dt; d_eps = eps_v - eps_v_prev [-].
                    double const d_eps_B =
                        potential_exchange_params_ptr
                            ->ceiling_micro_storage_includes_strain
                            ? (variables.volumetric_strain -
                               variables_prev.volumetric_strain)
                            : 0.0;  // [-]
                    ceiling_B_drho_L_hat_dpL =
                        -phi_m_now * drho_LR_dpL_B * (1.0 + d_eps_B) /
                        dt;  // kg/(m^3 s)/Pa

                    double dphi_deps_v_B = 0.0;  // [-]
                    if (dynamic_cast<MPL::PorosityFromMassBalance const*>(
                            &medium->property(MPL::PropertyType::porosity)))
                    {
                        double const w_B =
                            (variables.volumetric_strain -
                             variables_prev.volumetric_strain) +
                            (variables.effective_pore_pressure -
                             variables_prev.effective_pore_pressure) *
                                beta_SR;  // [-]
                        double const phi_unclamped_B =
                            (variables_prev.porosity + alpha * w_B) /
                            (1.0 + w_B);  // [-]
                        bool const clamp_active_B =
                            std::abs(phi_unclamped_B - phi) >
                            1e-12 * std::max(1.0, std::abs(phi));
                        dphi_deps_v_B =
                            clamp_active_B ? 0.0
                                           : (alpha - phi) / (1.0 + w_B);  // [-]
                    }
                    double const drho_l_deps_v_B =
                        dphi_deps_v_B *
                        (rho_lR_now + phi_m_now * rho_lR_eos.drho_lR_dnl);
                    ceiling_B_drho_L_hat_deps_v =
                        -(drho_l_deps_v_B * (1.0 + d_eps_B) +
                          (potential_exchange_params_ptr
                               ->ceiling_micro_storage_includes_strain
                               ? phi_m_now * rho_lR_now
                               : 0.0)) /
                        dt;  // kg/(m^3 s) per unit eps_v
                }
                // ── KKT micro-water ceiling (DESIGN.md 3.5; NOT adopted) ─────
                // At a KKT-active IP (status written by the micro update of this
                // iterate) the macro sink is the booked exchange rhohat = S_s
                // (MicroExchangeReceived), its tangents are those of S_s (D-5.2):
                //   d rhohat/d p_L, d rhohat/d eps_v  (computeCeilingKktActive-
                // ExchangeTangents). status 2, 3 (NonMonotone, PremiseViolated)
                // are handled as clamp points: base sink, base tangent.
                if (isKktCeiling(*potential_exchange_params_ptr))
                {
                    if (beta_SR != 0.0)
                    {
                        // The branch decision uses n_max(eps_v) only; with
                        // beta_SR != 0 it would depend on p_L (W-6).
                        OGS_FATAL(
                            "micro_ceiling_treatment = kkt: beta_SR = {:g} != "
                            "0; the dphi/dp_eff chain of the active tangent is "
                            "not implemented.",
                            beta_SR);
                    }
                    kkt_active =
                        *std::get<MicroCeilingStatus>(
                            this->current_states_[ip]) ==
                        static_cast<double>(MicroCeilingKktStatus::Active);
                    if (kkt_active)
                    {
                        double const s_kkt = microMassStrainTermSign(
                            *potential_exchange_params_ptr);  // [-]
                        kkt_rho_L_hat = -*std::get<MicroExchangeReceived>(
                            this->current_states_[ip]);  // kg/(m^3 s)
                        double const rho_lR_kkt = *std::get<MicroLiquidDensity>(
                            this->current_states_[ip]);  // kg/m^3
                        auto const eos_kkt = computeActiveMicroLiquidDensity(
                            n_l, rho_LR, local_solve_context,
                            *potential_exchange_params_ptr);
                        double const d_eps_kkt =
                            variables.volumetric_strain -
                            variables_prev.volumetric_strain;  // [-]
                        double const dphi_deps_v_kkt =
                            dynamic_cast<MPL::PorosityFromMassBalance const*>(
                                &medium->property(
                                    MPL::PropertyType::porosity))
                                ? porosityDerivativeWrtVolumetricStrain(
                                      alpha, phi, variables_prev.porosity,
                                      d_eps_kkt)
                                : 0.0;  // [-]
                        auto const tangents_kkt =
                            computeCeilingKktActiveExchangeTangents(
                                s_kkt, d_eps_kkt, phi, rho_lR_kkt,
                                eos_kkt.drho_lR_dnl, eos_kkt.drho_lR_drho_LR,
                                rho_LR * beta_LR, dphi_deps_v_kkt, dt);
                        kkt_drho_L_hat_dpL = -tangents_kkt.drhohat_dpL;
                        kkt_drho_L_hat_deps_v = -tangents_kkt.drhohat_deps_v;
                    }
                }
                auto const micro_potential = computeActiveMicroPotential(
                    n_l, rho_LR, local_solve_context,
                    *potential_exchange_params_ptr);
                if (potential_exchange_params_ptr
                        ->use_micro_liquid_density_for_micro_pressure)
                {
                    auto const rho_lR_state =
                        *std::get<MicroLiquidDensity>(this->current_states_[ip]);
                    if (std::isfinite(rho_lR_state) && rho_lR_state > 0.0)
                    {
                        rho_lR_exchange_input = rho_lR_state;
                        drho_lR_exchange_input_dpL = rho_lR_state * beta_LR;
                    }
                }
                use_vdw_micro_potential_for_active_exchange = true;
                mu_lR_vdw = micro_potential.mu_lR;
                dmu_lR_vdw_drho_lR = micro_potential.dmu_lR_drho_lR;
                // ── DSM Maxwell-conjugate term (B1) ──────────────────────────
                // Maxwell partner of the swelling eigenstress (one Psi). Sharp
                // gate p'>=Pi; freeze phi (B1). ==0 below the gate (gate-closed
                // runs unchanged bit-for-bit). Mirrors the assemble() path.
                //
                // Film-pressure coupling (maxwell sec.5): when ON, mu_lR ALREADY
                // carries the film term (folded in computeActiveMicroPotential),
                // so SKIP this strain-view block to avoid double-counting; the
                // film p-u tangent is added in the dedicated ON branch below.
                // Flag OFF -> this block runs unchanged (bit-for-bit Jacobian).
                if (!film_pressure_coupling)
                {
                    double const Pi_mc = p_L_m;  // Pa, disjoining = -rho*mu_lR > 0
                    // The INTEGRABLE partner takes eps_v + K_drained (not p_conf)
                    // and carries NO explicit (1-phi_M) factor (the n_S referencing
                    // cancels in the specific potential -- equipresence note in
                    // PotentialExchange.h), so p_conf / n_S / one-minus-n_l are no
                    // longer formed here; mirrors the L717 micro-solve fold.
                    double const rho_mc =
                        (std::isfinite(rho_lR_exchange_input) &&
                         rho_lR_exchange_input > 0.0)
                            ? rho_lR_exchange_input
                            : rho_LR;
                    // dPi/dn_l = (Pi/mu_lR)*dmu_lR_dnl (Pi=-rho*mu_lR; density-agnostic)
                    double const dPi_dnl_mc =
                        (std::abs(mu_lR_vdw) > 1e-300)
                            ? (Pi_mc / mu_lR_vdw) * micro_potential.dmu_lR_dnl
                            : 0.0;
                    // Pi'' = -rho*mu_lR''; density-agnostic (Pi=-rho*mu_lR) ->
                    // Pi''/Pi = mu_lR''/mu_lR, matching the dPi_dnl_mc convention.
                    double const d2Pi_dnl2_mc =
                        (std::abs(mu_lR_vdw) > 1e-300)
                            ? (Pi_mc / mu_lR_vdw) * micro_potential.d2mu_lR_dnl2
                            : 0.0;
                    // Drained bulk modulus from the elastic tangent stiffness, the
                    // SAME mechanical K the L717 fold uses (here reconstructed via
                    // computeElasticTangentStiffness, mirroring the ON p-u block
                    // below). dp_conf/deps_v = -dsigma'_m/deps_v.
                    double const K_drained_mc =
                        drainedBulkModulusFromStiffness<DisplacementDim>(
                            ip_data_[ip].computeElasticTangentStiffness(
                                variables, t, x_position, dt,
                                this->solid_material_,
                                *this->material_states_[ip]
                                     .material_state_variables));
                    // INTEGRABLE Maxwell mechanical partner (Vinay's Option-B,
                    // 2026-06-08): the SINGLE mu_lR_mech form, UNGATED (no Pi
                    // Heaviside) -- mirrors the L717 micro-solve fold so the
                    // macro-exchange mu_lR equals the micro-solve mu_lR
                    // (equipresence). PARAMETER-FREE.
                    auto const mc = computeIntegrableMechanicalMicroPotential(
                        Pi_mc, dPi_dnl_mc, d2Pi_dnl2_mc, n_l,
                        variables.volumetric_strain, /*biot_b=*/alpha,
                        K_drained_mc, rho_mc);
                    mu_lR_vdw += mc.mu_lR_mech;  // J/kg, additive (ungated)
                    dmu_lR_vdw_drho_lR += mc.dmu_lR_mech_drho_lR;
                    // DSM Maxwell-conjugate: exchange<->displacement tangent --
                    // the transpose partner of the swelling-eigenstress block.
                    // mu_lR_mech makes the exchange depend on eps_v, so the
                    // pressure residual rho_L_hat = alpha_M*(mu_lR - mu_LR) has
                    //   d rho_L_hat/d u = alpha_M * (dmu_lR_mech/deps_v) * m^T B,
                    // m = identity2 (eps_v = m^T B u). The integrable partner is
                    // SMOOTH in eps_v (NO Heaviside gate), so this fires whenever
                    // the exchange is active (no gate_open guard). Without this
                    // block the tangent is inconsistent and Newton diverges
                    // (Task-13 Pi blow-up). Sign mirrors the K[p,p] exchange (-=).
                    if (mc.dmu_lR_mech_deps_v != 0.0 && mu > 0.0)
                    {
                        double const alpha_M_eff_mc = alpha_bar * rho_LR / mu;
                        local_Jac
                            .template block<pressure_size, displacement_size>(
                                pressure_index, displacement_index)
                            .noalias() -=
                            N_p.transpose() *
                            (alpha_M_eff_mc * mc.dmu_lR_mech_deps_v) *
                            identity2.transpose() * B * w;
                    }
                }
                // ── INTEGRABLE Maxwell p-u tangent (spec item 2/3) ────────────
                // The exchange depends on eps_v through the integrable partner
                //   mu_lR_mech = -[ (Pi + n_l*Pi')*eps_v
                //                   + 0.5*b*K_drained*eps_v^2 ] / rho_lR,
                // so d mu_lR/d eps_v = -[ (Pi + n_l*Pi') + b*K_drained*eps_v ]
                // / rho_lR (computeIntegrableMechanicalMicroPotential.
                // dmu_lR_mech_deps_v). Then d rho_L_hat/d u =
                //   alpha_M*(d mu_lR/d eps_v)*identity2^T B, assembled with the
                // SAME sign (-=) as the K[p,p] exchange block. This is the
                // one-Psi transpose partner of the swelling-stress u-eps block
                // d sigma_sw/d eps_v = +(1-phi_M)*n_l*b*K (residual sigma_sw =
                // -phi_m*p_film). The partner is SMOOTH in eps_v (NO Heaviside
                // gate), so the block fires whenever film coupling is ON; it is
                // still EXACTLY zero when the flag is OFF (no block reached).
                // K_drained = dp_conf/deps_v (mechanical drained bulk modulus),
                // the SAME K the residual fold used (threaded via the context),
                // so the residual and tangent are consistent. PARAMETER-FREE.
                // DIAGNOSTIC B: at a clamped IP the exchange source is the
                // micro storage rate (not mu_lR-driven), so this mu_lR p-u
                // tangent is replaced by ceiling_B_drho_L_hat_deps_v below.
                // KKT active IP: the exchange is S_s, not mu_lR-driven, so the
                // Maxwell / live-K p-u entries of mu_lR do not apply (D-5.2);
                // effective only at micro_ceiling_pu_tangent = all_exchange
                // (otherwise erased by the Kpu/dt assignment below).
                if (film_pressure_coupling && mu > 0.0 && !ceiling_B_active &&
                    !kkt_active)
                {
                    double const rho_film =
                        (std::isfinite(rho_lR_exchange_input) &&
                         rho_lR_exchange_input > 0.0)
                            ? rho_lR_exchange_input
                            : rho_LR;
                    // BARE vdW Pi, Pi', Pi'' at this state (density mirrored to
                    // rho_film exactly as the residual fold), recomputed so the
                    // tangent is built from the SAME disjoining law. active_nS is
                    // context-independent (Reference -> ref fraction;
                    // CurrentPorositySplit -> 1 - n_l). Thread dnS_dnl exactly as
                    // computeActiveMicroPotential does (so Pi' here equals the
                    // residual's Pi' including the F2 live-nS chain).
                    double const active_nS_pu =
                        computeActiveMicroSolidVolumeFraction(
                            n_l, PotentialExchangeLocalSolveContext{},
                            *potential_exchange_params_ptr);
                    double const dnS_dnl_pu =
                        potential_exchange_params_ptr
                                    ->micro_solid_volume_fraction_mode ==
                                MicroSolidVolumeFractionMode::CurrentPorositySplit
                            ? -1.0
                            : 0.0;
                    auto const vdw_pu = computeVanDerWaalsMicroPotential(
                        n_l, rho_film, active_nS_pu,
                        potential_exchange_params_ptr
                            ->micro_solid_density_reference,
                        potential_exchange_params_ptr->hamaker_constant,
                        potential_exchange_params_ptr->specific_surface,
                        microPotentialSignFactorFromParameters(
                            *potential_exchange_params_ptr),
                        // Live K(rho_d): rho_d = rho_SR*(1-phi) (total
                        // porosity in scope); off -> parse scalar.
                        effectiveAugmentationPrefactor(
                            *potential_exchange_params_ptr, phi),  // K [J/kg]
                        potential_exchange_params_ptr
                            ->potential_augmentation_exponent,
                        dnS_dnl_pu,
                        potential_exchange_params_ptr
                            ->micro_water_content_floor);
                    double const Pi_pu = -rho_film * vdw_pu.mu_lR;          // Pa
                    double const dPi_dnl_pu = -rho_film * vdw_pu.dmu_lR_dnl;  // Pa/n_l
                    double const d2Pi_dnl2_pu =
                        -rho_film * vdw_pu.d2mu_lR_dnl2;  // Pa/n_l^2
                    // Drained bulk modulus = dp_conf/deps_v, the SAME mechanical
                    // K threaded into the residual fold. (Reconstruct C_el as the
                    // enable_dsm_swelling_up_jacobian block does below.)
                    auto const C_el_film =
                        ip_data_[ip].computeElasticTangentStiffness(
                            variables, t, x_position, dt,
                            this->solid_material_,
                            *this->material_states_[ip]
                                 .material_state_variables);
                    double const K_drained =
                        drainedBulkModulusFromStiffness<DisplacementDim>(
                            C_el_film);
                    auto const mech_pu =
                        computeIntegrableMechanicalMicroPotential(
                            Pi_pu, dPi_dnl_pu, d2Pi_dnl2_pu, n_l,
                            variables.volumetric_strain, /*biot_b=*/alpha,
                            K_drained, rho_film);
                    double const alpha_M_eff_film = alpha_bar * rho_LR / mu;

                    // ── M1 route dispatch (review fix 2026-06-14) ────────────
                    // The residual mu_lR's eps_v dependence differs by route:
                    //   Off (no strain coupling): the integrable partner
                    //     mu_lR_mech (folded at L844), d/deps_v = mech_pu
                    //     .dmu_lR_mech_deps_v -- the historical default below.
                    //   Operational-strained (film_strain_coupling != Off,
                    //     route != Exact, L759-817): mu_lR = bare(w_eff(eps_v))
                    //     + b*p_conf(eps_v)/rho, so
                    //       dmu_lR/deps_v = bare'(w_eff)*dw_eff/deps_v
                    //                       + b*(dp_conf/deps_v)/rho,
                    //     with bare'(w_eff) = the law's d/d(arg) (its dmu_lR_dnl
                    //     BEFORE the dw_eff_dnl chain at L796) and dp_conf/deps_v
                    //     = -K_drained (drainedBulkModulusFromStiffness = dsigma'
                    //     _m/deps_v = -dp_conf/deps_v; L61).
                    //   Exact (film_energy_route == Exact, L703-748): mu_lR +=
                    //     g_cut*pair.mu_mech with the bare evaluated at the TRUE
                    //     n_l (no eps_v), so dmu_lR/deps_v = g_cut*pair
                    //     .dmu_mech_deps_v.
                    // JACOBIAN-ONLY: residual untouched. Off-mode reaches the
                    // default branch unchanged (bit-for-bit).
                    double dmu_lR_deps_v_film =
                        mech_pu.dmu_lR_mech_deps_v;  // J/kg per strain (Off)
                    auto const& pep_m1 = *potential_exchange_params_ptr;
                    if (pep_m1.film_strain_coupling !=
                        FilmStrainCouplingMode::Off)
                    {
                        double const sign_m1 =
                            microPotentialSignFactorFromParameters(pep_m1);
                        double const K_aug_m1 = effectiveAugmentationPrefactor(
                            pep_m1, phi);  // K [J/kg]
                        if (pep_m1.film_energy_route == FilmEnergyRoute::Exact)
                        {
                            double const kappa_m1 =
                                pep_m1.film_strain_kappa ==
                                        FilmStrainKappaMode::Aggregate
                                    ? active_nS_pu
                                    : 1.0;  // [-]
                            auto const pair_m1 = computeStrainedFilmEnergyPair(
                                n_l, variables.volumetric_strain, kappa_m1,
                                alpha, K_drained, true /*include_S, R3*/,
                                rho_film, active_nS_pu,
                                pep_m1.micro_solid_density_reference,
                                pep_m1.hamaker_constant, pep_m1.specific_surface,
                                sign_m1, K_aug_m1,
                                pep_m1.potential_augmentation_exponent,
                                0.0 /*dnS_dnl: frozen nS (B1)*/,
                                pep_m1.micro_water_content_floor);
                            // g_cut = mu_lR(post macro-floor cutoff)/mu_bare_pre,
                            // matching the residual fold (L738). The residual
                            // mu_lR_vdw already carries the cutoff; recover g via
                            // bare_pre. |mu_bare_pre| > 0 by the bare law FATALs.
                            double const g_cut_m1 =
                                vdw_pu.mu_lR / pair_m1.mu_bare_pre;  // [-]
                            dmu_lR_deps_v_film =
                                g_cut_m1 * pair_m1.dmu_mech_deps_v;  // J/kg/strain
                        }
                        else
                        {
                            // Operational-strained route.
                            auto const film_state_m1 = computeStrainedFilmState(
                                pep_m1.film_strain_coupling,
                                pep_m1.film_strain_kappa, n_l, active_nS_pu,
                                variables.volumetric_strain, p_conf_assembly,
                                rho_film,
                                pep_m1.micro_solid_density_reference,
                                pep_m1.hamaker_constant, pep_m1.specific_surface,
                                sign_m1, K_aug_m1,
                                pep_m1.potential_augmentation_exponent,
                                pep_m1.micro_water_content_floor,
                                rho_film /*rho_pi*/);
                            // Bare law at w_eff; its dmu_lR_dnl is d(bare)/d(arg)
                            // (the arg plays the role of n_l), so multiplying by
                            // dw_eff/deps_v gives d(bare(w_eff))/deps_v.
                            auto const bare_weff_m1 =
                                computeVanDerWaalsMicroPotential(
                                    film_state_m1.w_eff, rho_film, active_nS_pu,
                                    pep_m1.micro_solid_density_reference,
                                    pep_m1.hamaker_constant,
                                    pep_m1.specific_surface, sign_m1, K_aug_m1,
                                    pep_m1.potential_augmentation_exponent,
                                    0.0 /*dnS_dnl: frozen nS (B1)*/,
                                    pep_m1.micro_water_content_floor);
                            double const dbare_deps_v =
                                bare_weff_m1.dmu_lR_dnl *
                                film_state_m1.dw_eff_deps_v;  // J/kg per strain
                            // d(mu_load)/deps_v = b*(dp_conf/deps_v)/rho,
                            // dp_conf/deps_v = -K_drained (L61 sign).
                            double const dmuload_deps_v =
                                alpha * (-K_drained) / rho_film;  // J/kg/strain
                            dmu_lR_deps_v_film =
                                dbare_deps_v + dmuload_deps_v;  // J/kg per strain
                        }
                    }
                    local_Jac
                        .template block<pressure_size, displacement_size>(
                            pressure_index, displacement_index)
                        .noalias() -=
                        N_p.transpose() *
                        (alpha_M_eff_film * dmu_lR_deps_v_film) *
                        identity2.transpose() * B * w;

                    // ── Live K(rho_d) analytic tangent (K_OF_RHO_D_LIVE.md;
                    // Vinay 2026-06-12 approved Jacobian completion) ─────────
                    // JACOBIAN-ONLY: the residual's live K is untouched. In
                    // live mode K = K_table(rho_SR*(1-phi)) makes mu_lR depend
                    // on eps_v through phi, so the exchange p-u tangent gains
                    // the product-rule chain
                    //   d mu_lR/d eps_v |_K = dmu_lR/dK * dK/dphi * dphi/deps_v
                    // with dK/dphi = -rho_SR*(table segment slope) (exact,
                    // clamped-edge slope 0; PotentialExchangeParameters.h) and
                    // dmu_lR/dK covering the two LINEAR K-channels the residual
                    // mu_lR carries: the direct mu_aug term and, via Pi/Pi' in
                    // the integrable mechanical partner mu_lR_mech,
                    //   d mu_lR_mech/dK = (dmu_lR/dK + n_l*d(dmu_lR/dnl)/dK)
                    //                     * eps_v
                    // (computeIntegrableMechanicalMicroPotential form; Pi =
                    // -rho*mu_lR makes the rho factors cancel). Off-mode /
                    // frozen table / clamped edge -> dK/dphi == 0 -> block
                    // skipped, bit-for-bit identical Jacobian.
                    double const dK_dphi_pu =
                        effectiveAugmentationPrefactorPhiDerivative(
                            *potential_exchange_params_ptr,
                            phi);  // J/kg per unit phi
                    if (dK_dphi_pu != 0.0)
                    {
                        // dphi/deps_v of the porosity law the residual actually
                        // evaluated. PorosityFromMassBalance (the MS33 law):
                        //   phi = (phi_prev + alpha*w)/(1 + w),
                        //   w = delta_eps_v + delta_p_eff*beta_SR
                        // => dphi/deps_v = (alpha - phi)/(1 + w)  (analytic,
                        // derived here from that law's value()). Any other
                        // porosity law (e.g. Constant: exact) is treated as
                        // strain-independent -> chain 0.
                        double dphi_deps_v_pu = 0.0;  // [-]
                        if (dynamic_cast<
                                MPL::PorosityFromMassBalance const*>(
                                &medium->property(MPL::PropertyType::porosity)))
                        {
                            double const w_phi_pu =
                                (variables.volumetric_strain -
                                 variables_prev.volumetric_strain) +
                                (variables.effective_pore_pressure -
                                 variables_prev.effective_pore_pressure) *
                                    beta_SR;  // [-]
                            // N1 fix (2026-06-14): PorosityFromMassBalance
                            // CLAMPS phi to [phi_min, phi_max] (its value(),
                            // PorosityFromMassBalance.cpp L56). When the clamp
                            // is active phi is flat in eps_v -> dphi/deps_v = 0;
                            // the analytic (alpha-phi)/(1+w) would be a spurious
                            // nonzero tangent. The bounds are private, so detect
                            // the clamp by comparing the unclamped law value to
                            // the stored (clamped) phi: a mismatch beyond a
                            // relative floor means a bound is active.
                            double const phi_unclamped_pu =
                                (variables_prev.porosity + alpha * w_phi_pu) /
                                (1.0 + w_phi_pu);  // [-]
                            bool const clamp_active_pu =
                                std::abs(phi_unclamped_pu - phi) >
                                1e-12 * std::max(1.0, std::abs(phi));  // [-]
                            dphi_deps_v_pu =
                                clamp_active_pu
                                    ? 0.0
                                    : (alpha - phi) / (1.0 + w_phi_pu);  // [-]
                        }
                        if (dphi_deps_v_pu != 0.0)
                        {
                            double const dmu_lR_dK_tot_pu =
                                vdw_pu.dmu_lR_dK +
                                (vdw_pu.dmu_lR_dK +
                                 n_l * vdw_pu.ddmu_lR_dnl_dK) *
                                    variables.volumetric_strain;  // [-]
                            // Same sign/shape as the dmu_lR_mech_deps_v block
                            // above: an additional contribution to
                            // d mu_lR/d eps_v in d rho_L_hat/d u.
                            local_Jac
                                .template block<pressure_size,
                                                displacement_size>(
                                    pressure_index, displacement_index)
                                .noalias() -=
                                N_p.transpose() *
                                (alpha_M_eff_film * dmu_lR_dK_tot_pu *
                                 dK_dphi_pu * dphi_deps_v_pu) *  // J/kg per eps_v
                                identity2.transpose() * B * w;
                        }
                    }
                }
                use_fd_jacobian_for_direct_macro_derivative =
                    potential_exchange_params_ptr
                        ->use_fd_jacobian_for_exchange;
                fd_jacobian_perturbation =
                    potential_exchange_params_ptr->fd_jacobian_perturbation;

                // In analytic mode, include implicit n_l(p_L) chain coupling
                // in dmu_lR/dp_L for the active exchange Jacobian term.
                if (!use_fd_jacobian_for_direct_macro_derivative)
                {
                    requirePositiveViscosity(
                        "RichardsMechanics local exchange Jacobian assembly",
                        mu);
                    double const drho_LR_dpL = rho_LR * beta_LR;
                    auto const macro_potential = computeYoungLaplaceMacroPotential(
                        p_L_ip, rho_LR, pressure_tolerance);
                    double const alpha_M_effective = alpha_bar * rho_LR / mu;
                    auto const exchange = computePotentialDrivenMassExchange(
                        alpha_M_effective, macro_potential.mu_LR,
                        micro_potential.mu_lR);
                    // F1: thread the converged (n_l, micro rho_lR) so the
                    // ScalarReferenceMassStorage branch linearizes the REV-mass
                    // residual at the converged state. rho_lR_exchange_input is
                    // the converged micro liquid density (set above when
                    // use_micro_liquid_density_for_micro_pressure, i.e. always in
                    // mass-storage mode); NaN -> the helper recomputes via EOS.
                    // KKT active IP: n_l = n_max(eps_v), dn_l/dp_L = 0 (D-5.3.1);
                    // the unconstrained root's derivative would be wrong here.
                    double const dn_l_dpL = kkt_active ? 0.0 : computeImplicitNlDpL(
                        n_l_prev, p_L_ip, dt, rho_LR, drho_LR_dpL, alpha_bar, mu,
                        macro_potential, micro_potential, exchange,
                        local_solve_context,
                        *potential_exchange_params_ptr,
                        /*n_l_converged=*/n_l,
                        /*rho_lR_micro=*/rho_lR_exchange_input);

                    // Full total derivative of the vdW micro potential w.r.t.
                    // pL. NOTE (on-disk): dmu_lR_drho_lR is NON-zero here
                    // (= -mu_lR/rho_lR; PotentialExchange.h line 181, "non-zero
                    // after /rho_lR fix"), despite the stale struct comment at
                    // line 64 ("exactly zero in the reduced algebraic form").
                    // It is paired with the BULK drho_LR_dpL, matching both
                    // computeImplicitNlDpL (line ~1327) and the
                    // computePotentialExchangeUpdate fallback (line ~217). The
                    // dominant contribution is the implicit n_l(p_L) chain
                    // dmu_lR_dnl * dn_l_dpL.
                    // Film-pressure coupling (maxwell sec.5, increment D-ii):
                    // micro_potential.dmu_lR_dnl and .dmu_lR_drho_lR ALREADY
                    // carry the film contributions (folded in
                    // computeActiveMicroPotential, B3), so this single line
                    // captures BOTH film p_L channels: (i) the implicit
                    // n_l(p_L) chain through the gate (dmu_lR_film_dnl*dn_l_dpL,
                    // since dmu_lR_dnl includes the gate's dPi_gate/dn_l term)
                    // and (ii) the rho_lR(p_L) channel
                    // dmu_lR_film_drho_lR*drho_lR_dpL (here rho_lR == bulk rho_LR
                    // in ScalarReferenceStorage mode, so the bulk drho_LR_dpL is
                    // the right pairing, matching the vdW convention above). NOT
                    // included (by design): the direct p_conf(p_L) channel via
                    // Bishop effective stress -- that lagged u-p coupling is left
                    // to the consistent-tangent block (enable_dsm_swelling_up_
                    // jacobian, default OFF) per the spec.
                    dmu_lR_vdw_dpL = micro_potential.dmu_lR_dnl * dn_l_dpL +
                                     micro_potential.dmu_lR_drho_lR *
                                         drho_LR_dpL;
                    use_custom_dmu_lR_vdw_dpL = true;

                    // ── M2+L2: displacement-side live-K swelling-eigenstress
                    // tangent (review fix 2026-06-14; the 1b compliant-top cure)
                    // ───────────────────────────────────────────────────────
                    // The residual swelling eigenstress (computeReferenceMicro-
                    // PorositySwellingStressIncrement, L1983/L2090-2169) sources
                    // its K from K_aug_sw = effectiveAugmentationPrefactor(params,
                    // total_porosity=phi), so in live-K mode delta_sigma_sw
                    // depends on eps_v and p_L through phi. That dependence was
                    // present in the residual but ABSENT from K[u,u]/K[u,p] (the
                    // only swelling strain tangent was the constexpr-false dead
                    // block below, which carries no dK/dphi chain). Wire ONLY the
                    // live-K chain here -- NOT the pre-existing swelling u-p/u-u
                    // term Vinay set OFF on 2026-06-01 (enable_dsm_swelling_up_
                    // jacobian stays at its default; the two are independent).
                    //   d(delta_sigma_sw)/d(.) = d(delta_sigma_sw)/dK
                    //                            * dK/dphi * dphi/d(.),
                    //   delta_sigma_sw = n_S*(n_l_prev*p_film_prev
                    //                          - n_l*p_film_curr)*I,
                    //   p_film = Pi - b*p_conf,  Pi = -rho*mu_lR(w_eff; K),
                    //   => d(delta_sigma_sw)/dK
                    //      = -n_S*( n_l_prev*rho_prev*dmu_lR_prev/dK
                    //               - n_l*rho_curr*dmu_lR_curr/dK )*I,
                    // mapped to R_u via dsigma'/d(.) = C*C_el^{-1}
                    // *d(delta_sigma_sw)/d(.). dphi/deps_v = (alpha-phi)/(1+w)
                    // (PorosityFromMassBalance; w = delta_eps_v
                    // + delta_p_eff*beta_SR) and dphi/dp = dphi/deps_v*beta_SR
                    // (dw/dp_eff = beta_SR). Gated on live-K mode AND dK/dphi !=
                    // 0 -> fires ONLY when live K is active inside the table
                    // interior; off / frozen / clamped edge -> skipped ->
                    // Jacobian bit-for-bit. The residual eigenstress uses live K
                    // in BOTH its film-ON branch and its OFF (disjoining) branch
                    // (computeReferenceMicroPorositySwellingStressIncrement
                    // L2380/L2389 use K_aug_sw), so this gate is the live-K flag,
                    // NOT film_pressure_coupling -- otherwise the 1b family (live
                    // K, film coupling OFF) keeps the residual dependence without
                    // the matching tangent (the 1b_A step-1 divergence).
                    // JACOBIAN-ONLY: residual untouched.
                    if (potential_exchange_params_ptr
                            ->potential_augmentation_prefactor_live_dry_density)
                    {
                        double const dK_dphi_sw =
                            effectiveAugmentationPrefactorPhiDerivative(
                                *potential_exchange_params_ptr,
                                phi);  // J/kg per unit phi
                        bool const is_pfmb =
                            dynamic_cast<MPL::PorosityFromMassBalance const*>(
                                &medium->property(
                                    MPL::PropertyType::porosity)) != nullptr;
                        if (dK_dphi_sw != 0.0 && is_pfmb)
                        {
                            auto const& pep_sw = *potential_exchange_params_ptr;
                            // Residual-identical state inputs (mirror
                            // updateSwellingStateWithMicroPorosity, L2325-2371).
                            double const n_l_prev_sw =
                                **std::get<PrevState<MicroWaterContent>>(
                                    this->prev_states_[ip]);
                            double const phi_M_sw =
                                std::get<ProcessLib::ThermoRichardsMechanics::
                                             TransportPorosityData>(
                                    this->current_states_[ip])
                                    .phi;
                            double const n_S_sw =
                                std::max(1e-16, 1.0 - phi_M_sw);  // [-]
                            double const rho_lR_curr_micro_sw =
                                *std::get<MicroLiquidDensity>(
                                    this->current_states_[ip]);
                            double const rho_lR_prev_micro_sw =
                                **std::get<PrevState<MicroLiquidDensity>>(
                                    this->prev_states_[ip]);
                            // Density mirrors the residual's p_L_m choice
                            // (micro rho_lR when enabled, bulk otherwise).
                            double const rho_pi_prev_sw =
                                pep_sw
                                        .use_micro_liquid_density_for_micro_pressure
                                    ? rho_lR_prev_micro_sw
                                    : rho_LR;
                            double const rho_pi_curr_sw =
                                pep_sw
                                        .use_micro_liquid_density_for_micro_pressure
                                    ? rho_lR_curr_micro_sw
                                    : rho_LR;
                            double const eps_v_sw = variables.volumetric_strain;
                            double const eps_v_prev_sw =
                                variables_prev.volumetric_strain;
                            double const sign_sw =
                                microPotentialSignFactorFromParameters(pep_sw);
                            double const K_aug_sw_j =
                                effectiveAugmentationPrefactor(pep_sw,
                                                               phi);  // K [J/kg]
                            double const active_nS_prev_sw =
                                computeActiveMicroSolidVolumeFraction(
                                    n_l_prev_sw,
                                    PotentialExchangeLocalSolveContext{},
                                    pep_sw);  // [-]
                            double const active_nS_curr_sw =
                                computeActiveMicroSolidVolumeFraction(
                                    n_l, PotentialExchangeLocalSolveContext{},
                                    pep_sw);  // [-]
                            // w_eff at prev/curr EXACTLY as the residual: the
                            // strained w_eval is taken ONLY inside the residual's
                            // film_pressure_coupling block (L2073); the OFF
                            // (disjoining) branch uses w_eval = n_l. So gate the
                            // strained state on film_pressure_coupling here too --
                            // NOT on film_strain_coupling -- otherwise 1b_B
                            // (kinematic coupling, film OFF) would strain w_eval
                            // while its residual eigenstress does not.
                            double w_eval_prev_sw = n_l_prev_sw;
                            double w_eval_curr_sw = n_l;
                            // L3 also needs d(w_eval_curr)/dn_l of the residual's
                            // eval argument: 1 on the OFF branch (w_eval = n_l),
                            // dw_eff/dn_l on the strained branch.
                            double dw_eval_curr_dnl_sw = 1.0;  // [-]
                            if (film_pressure_coupling &&
                                pep_sw.film_strain_coupling !=
                                    FilmStrainCouplingMode::Off &&
                                std::isfinite(eps_v_sw))
                            {
                                double const eps_v_prev_used_sw =
                                    std::isfinite(eps_v_prev_sw) ? eps_v_prev_sw
                                                                 : eps_v_sw;
                                w_eval_prev_sw =
                                    computeStrainedFilmState(
                                        pep_sw.film_strain_coupling,
                                        pep_sw.film_strain_kappa, n_l_prev_sw,
                                        active_nS_prev_sw, eps_v_prev_used_sw,
                                        p_conf_assembly, rho_lR_prev_micro_sw,
                                        pep_sw.micro_solid_density_reference,
                                        pep_sw.hamaker_constant,
                                        pep_sw.specific_surface, sign_sw,
                                        K_aug_sw_j,
                                        pep_sw.potential_augmentation_exponent,
                                        pep_sw.micro_water_content_floor,
                                        rho_pi_prev_sw)
                                        .w_eff;
                                auto const film_state_curr_sw =
                                    computeStrainedFilmState(
                                        pep_sw.film_strain_coupling,
                                        pep_sw.film_strain_kappa, n_l,
                                        active_nS_curr_sw, eps_v_sw,
                                        p_conf_assembly, rho_lR_curr_micro_sw,
                                        pep_sw.micro_solid_density_reference,
                                        pep_sw.hamaker_constant,
                                        pep_sw.specific_surface, sign_sw,
                                        K_aug_sw_j,
                                        pep_sw.potential_augmentation_exponent,
                                        pep_sw.micro_water_content_floor,
                                        rho_pi_curr_sw);
                                w_eval_curr_sw = film_state_curr_sw.w_eff;
                                dw_eval_curr_dnl_sw =
                                    film_state_curr_sw.dw_eff_dnl;  // [-]
                            }
                            // dmu_lR/dK of the bare law at the two states (the
                            // augmentation channel is linear in K -> dmu_lR_dK).
                            double const dmu_lR_prev_dK_sw =
                                computeVanDerWaalsMicroPotential(
                                    w_eval_prev_sw, rho_lR_prev_micro_sw,
                                    active_nS_prev_sw,
                                    pep_sw.micro_solid_density_reference,
                                    pep_sw.hamaker_constant,
                                    pep_sw.specific_surface, sign_sw, K_aug_sw_j,
                                    pep_sw.potential_augmentation_exponent,
                                    0.0 /*dnS_dnl*/,
                                    pep_sw.micro_water_content_floor)
                                    .dmu_lR_dK;  // [-] (J/kg per J/kg)
                            // Full curr-state vdW struct (L3 also reads .mu_lR
                            // and .dmu_lR_dnl from it; dnS_dnl frozen to 0
                            // matches the residual eigenstress, L2390).
                            auto const vdw_curr_sw =
                                computeVanDerWaalsMicroPotential(
                                    w_eval_curr_sw, rho_lR_curr_micro_sw,
                                    active_nS_curr_sw,
                                    pep_sw.micro_solid_density_reference,
                                    pep_sw.hamaker_constant,
                                    pep_sw.specific_surface, sign_sw, K_aug_sw_j,
                                    pep_sw.potential_augmentation_exponent,
                                    0.0 /*dnS_dnl*/,
                                    pep_sw.micro_water_content_floor);
                            double const dmu_lR_curr_dK_sw =
                                vdw_curr_sw.dmu_lR_dK;  // [-]
                            // d(delta_sigma_sw)/dK scalar (on identity2):
                            // -n_S*( n_l_prev*rho_prev*dmu_prev/dK
                            //        - n_l*rho_curr*dmu_curr/dK ).  [Pa per J/kg]
                            double const d_delta_sigma_sw_dK_scalar =
                                -n_S_sw *
                                (n_l_prev_sw * rho_pi_prev_sw *
                                     dmu_lR_prev_dK_sw -
                                 n_l * rho_pi_curr_sw * dmu_lR_curr_dK_sw);
                            // dphi/deps_v and dphi/dp_eff (PorosityFromMassBalance).
                            double const w_phi_sw =
                                (variables.volumetric_strain -
                                 variables_prev.volumetric_strain) +
                                (variables.effective_pore_pressure -
                                 variables_prev.effective_pore_pressure) *
                                    beta_SR;  // [-]
                            // N1 clamp guard (2026-06-14): zero the chain when
                            // PorosityFromMassBalance clamps phi to a bound
                            // (flat in eps_v/p there); detect via unclamped vs
                            // stored phi (bounds are private). Same logic as the
                            // live-K p-u block below.
                            double const phi_unclamped_sw =
                                (variables_prev.porosity + alpha * w_phi_sw) /
                                (1.0 + w_phi_sw);  // [-]
                            bool const clamp_active_sw =
                                std::abs(phi_unclamped_sw - phi) >
                                1e-12 * std::max(1.0, std::abs(phi));  // [-]
                            double const dphi_deps_v_sw =
                                clamp_active_sw
                                    ? 0.0
                                    : (alpha - phi) / (1.0 + w_phi_sw);  // [-]
                            double const dphi_dp_sw =
                                dphi_deps_v_sw * beta_SR;  // 1/Pa
                            // d(delta_sigma_sw)/deps_v and /dp via the live-K
                            // chain. This is the EXPLICIT-K channel (M2): K(phi)
                            // varies with strain at FIXED converged n_l.
                            double dsig_sw_deps_v_scalar =
                                d_delta_sigma_sw_dK_scalar * dK_dphi_sw *
                                dphi_deps_v_sw;  // Pa per unit strain
                            double dsig_sw_dp_scalar =
                                d_delta_sigma_sw_dK_scalar * dK_dphi_sw *
                                dphi_dp_sw;  // Pa/Pa

                            // ── L3 (review 2026-06-14, JACOBIAN-ONLY) ─────────
                            // IMPLICIT n_l(K) channel of the SAME eigenstress.
                            // In ScalarReferenceMassStorage mode the local solve
                            // returns n_l satisfying the REV-mass residual
                            // r(n_l;K)=0, and K=K(phi(eps_v)), so the converged
                            // n_l ALSO moves with strain through K -- a channel
                            // the explicit-K term above (n_l held) omits and the
                            // fixed-n_l dn_l/dpL does not carry (review L3). By
                            // the chain rule on delta_sigma_sw(n_l(K(phi)),K):
                            //   d(delta_sigma_sw)/d(.) |_implicit
                            //     = d(delta_sigma_sw)/dn_l * dn_l/dK
                            //       * dK/dphi * dphi/d(.),
                            // with dn_l/dK from computeImplicitNlDK (same 1x1
                            // REV-mass reduction as dn_l/dpL; returns 0 outside
                            // mass-storage mode -> this whole channel vanishes
                            // for ScalarExchange/ReferenceStorage). The explicit
                            // partial d(delta_sigma_sw)/dn_l is taken at FIXED K
                            // and FIXED prev state, from the residual increment
                            // delta_sigma_sw = n_S*(n_l_prev*Pi_prev
                            //                        - n_l*Pi_curr) (L2406/L2301),
                            // Pi_curr = -rho_curr*mu_lR_curr(n_l;K):
                            //   d(delta_sigma_sw)/dn_l
                            //     = -n_S*( Pi_curr + n_l*dPi_curr/dn_l )
                            //     = -n_S*( Pi_curr - n_l*rho_curr
                            //              * dmu_lR_curr/dw_eval
                            //              * dw_eval/dn_l ),
                            // where dmu_lR_curr/dw_eval = vdw_curr_sw.dmu_lR_dnl
                            // (the law's derivative w.r.t. its eval argument) and
                            // dw_eval/dn_l = 1 on the OFF branch (w_eval=n_l) or
                            // film_state.dw_eff_dnl on the strained branch --
                            // EXACTLY the residual's argument-chain (L2384/L2073).
                            // Mirrors the residual's dnS_dnl=0 (active_nS frozen)
                            // and held rho_curr. JACOBIAN-ONLY: the local FORWARD
                            // n_l solve and the residual are untouched; this only
                            // completes the assembled displacement tangent.
                            // Folded into the SAME scalars so it rides the
                            // identical C*C_el^{-1}*identity2 map below. Off /
                            // frozen / clamped edge / non-mass-storage -> the
                            // factors are 0 -> bit-for-bit.
                            // Exact route (H1) sources the residual eigenstress
                            // from the one-Psi pair, NOT this telescoped form, so
                            // the telescoped d(delta_sigma_sw)/dn_l here is not
                            // the assembled residual's partial there. Restrict
                            // the implicit chain to the telescoped-residual
                            // regimes (OFF + film-ON-operational) so it stays a
                            // true partial of the assembled eigenstress and is a
                            // bit-for-bit no-op on the exact-route 1b_B (which is
                            // already converging; its tangent is left untouched).
                            bool const exact_route_l3 =
                                film_pressure_coupling &&
                                pep_sw.film_strain_coupling !=
                                    FilmStrainCouplingMode::Off &&
                                pep_sw.film_energy_route ==
                                    FilmEnergyRoute::Exact &&
                                std::isfinite(eps_v_sw);
                            // KKT active IP: dn_l/dK = 0 (D-5.3.1, n_l = n_max).
                            double const dn_l_dK_sw =
                                (exact_route_l3 || kkt_active)
                                    ? 0.0
                                    : computeImplicitNlDK(
                                          n_l_prev_sw, dt, rho_LR, mu,
                                          micro_potential, exchange,
                                          local_solve_context,
                                          *potential_exchange_params_ptr,
                                          /*n_l_converged=*/n_l,
                                          /*rho_lR_micro=*/
                                          rho_lR_exchange_input);  // n_l/(J/kg)
                            if (dn_l_dK_sw != 0.0)
                            {
                                double const Pi_curr_sw =
                                    -rho_pi_curr_sw * vdw_curr_sw.mu_lR;  // Pa
                                // d(delta_sigma_sw)/dn_l at fixed K, prev
                                // (argument chain dw_eval/dn_l included):
                                double const d_delta_sigma_sw_dnl_sw =
                                    -n_S_sw *
                                    (Pi_curr_sw -
                                     n_l * rho_pi_curr_sw *
                                         vdw_curr_sw.dmu_lR_dnl *
                                         dw_eval_curr_dnl_sw);  // Pa per n_l
                                // [Pa per n_l]*[n_l/(J/kg)]*[(J/kg)/phi] = Pa/phi.
                                double const d_delta_sigma_sw_dphi_implicit_sw =
                                    d_delta_sigma_sw_dnl_sw * dn_l_dK_sw *
                                    dK_dphi_sw;  // Pa per unit phi
                                dsig_sw_deps_v_scalar +=
                                    d_delta_sigma_sw_dphi_implicit_sw *
                                    dphi_deps_v_sw;  // Pa per unit strain
                                dsig_sw_dp_scalar +=
                                    d_delta_sigma_sw_dphi_implicit_sw *
                                    dphi_dp_sw;  // Pa/Pa
                            }
                            // Map to R_u: dsigma'/d(.) = C*C_el^{-1}
                            // *d(delta_sigma_sw)/d(.) (the swelling eigenstress
                            // enters eps_m = eps + C_el^{-1}:sigma_sw).
                            auto const& C_consistent_sw =
                                *std::get<StiffnessTensor<DisplacementDim>>(
                                    constitutive_data);
                            auto const C_el_sw =
                                ip_data_[ip].computeElasticTangentStiffness(
                                    variables, t, x_position, dt,
                                    this->solid_material_,
                                    *this->material_states_[ip]
                                         .material_state_variables);
                            auto const C_el_inv_sw = C_el_sw.inverse().eval();
                            MathLib::KelvinVector::KelvinVectorType<
                                DisplacementDim> const dsig_sw_deps_v =
                                dsig_sw_deps_v_scalar * identity2;  // Pa
                            MathLib::KelvinVector::KelvinVectorType<
                                DisplacementDim> const dsig_sw_dp =
                                dsig_sw_dp_scalar * identity2;  // Pa
                            // K[u,u]: eps_v = identity2^T B u.
                            local_Jac
                                .template block<displacement_size,
                                                displacement_size>(
                                    displacement_index, displacement_index)
                                .noalias() += B.transpose() * C_consistent_sw *
                                              C_el_inv_sw * dsig_sw_deps_v *
                                              identity2.transpose() * B * w;
                            // K[u,p].
                            local_Jac
                                .template block<displacement_size,
                                                pressure_size>(
                                    displacement_index, pressure_index)
                                .noalias() += B.transpose() * C_consistent_sw *
                                              C_el_inv_sw * dsig_sw_dp * N_p * w;
                        }
                    }

                    // --- DSM swelling-eigenstress u-p Jacobian (full p^disj) -
                    // Consistent-tangent completeness term for the swelling
                    // eigenstress that enters R_u through the mechanical strain
                    // (eps_m = eps + C_el^{-1} : sigma_sw, see line ~3080 and
                    // the swelling-state update at line ~1688). Differentiating
                    // the DSM eigenstress w.r.t. pL propagates as
                    //   dsigma'/dpL = C * C_el^{-1} * d(delta_sigma_sw)/dpL,
                    // where delta_sigma_sw =
                    //   n_S*(n_l_prev*Pi_prev - n_l*Pi_curr)*identity2 (the
                    // *_prev terms are frozen) and -n_l*Pi_curr =
                    //   +n_l * rho_d * mu_lR, with rho_d = micro liquid density
                    // (when use_micro_liquid_density_for_micro_pressure) else
                    // bulk rho_LR (mirrors the residual, line ~1618). It reuses
                    // the analytic dn_l_dpL just computed, hence its placement
                    // inside this !use_fd_jacobian guard.
                    //
                    // CAVEAT (carried forward from the upstream authors, see
                    // the commented block at line ~3315): for the classical
                    // saturation_micro swelling path this u-p coupling "does
                    // not improve convergence and sometimes worsens it". It is
                    // included here only for consistent-tangent completeness on
                    // the DSM potential-exchange path and is ISOLATED behind
                    // the opt-in flag below so it can be compiled out by
                    // flipping one line without disturbing the residual or any
                    // other Jacobian block.
                    // ENABLED 2026-06-09 (dsm_maxwell_jac_parallel): the u-side
                    // swelling tangent completes the consistent Maxwell Jacobian.
                    // The film-ON block below (K[u,p] + K[u,u]) implements the
                    // Maxwell-identity pair: d sigma_sw/d eps_v =
                    // +(1-phi_M)*n_l*b*K_drained and d sigma_sw/d n_l =
                    // -(1-phi_M)*(p_film + n_l*Pi') routed through dn_l/dpL — the
                    // exact transpose of the eigenstress block (THM_DSM_Richards_
                    // maxwell_web.wl L54, L135).
                    //
                    // u-side swelling Jacobian blocks (K[u,p]/K[u,u]) kept in code
                    // but OFF by default — enabling them singularizes the assembled
                    // tangent on stiff/dense MS33 cases (dd1800, ModelIII gap2mm:
                    // SparseLU compute() failure) per the at-scale comparison
                    // 2026-06-09; flip to true to opt in. Analytic micro 2x2 (a)
                    // is unchanged — the strict win.
                    constexpr bool enable_dsm_swelling_up_jacobian = false;
                    // Scope: already inside `if (potential_exchange_enabled)`,
                    // which IS the p^disj (Pi-path) DSM path. Do NOT additionally
                    // gate on saturation_micro: the Pi-path models REMOVE that
                    // MPL property as vestigial (see e.g. ms33_modelI_dd1400.prj
                    // line ~205), so a hasProperty(saturation_micro) gate would
                    // make this term silently never fire on the real models.
                    if (enable_dsm_swelling_up_jacobian && film_pressure_coupling)
                    {
                        // ── Film-pressure swelling-stress tangent (maxwell sec.5,
                        // 2026-06-06 MICRO-WEIGHTED PRESSURE form) ──────────────
                        // Residual sigma_sw = -(1 - phi_M)*n_l*(Pi - b*p_conf)
                        //                   = -(1 - phi_M)*n_l*p_film, so by the
                        // product rule (p_film = Pi - b*p_conf):
                        //   d sigma_sw/d n_l   = -(1 - phi_M)*( p_film + n_l*dPi/dn_l ),
                        //   d sigma_sw/d eps_v = -(1 - phi_M)*n_l*b*(dp_conf/deps_v)
                        //                      = +(1 - phi_M)*n_l*b*K_drained
                        // (dp_conf/deps_v = -K_drained, the MECHANICAL drained bulk
                        // modulus; Pi and n_l do not depend on eps_v at fixed n_l).
                        // K_sw is GONE here -- it no longer routes a pressure
                        // through an elastic strain. Pi is the BARE vdW disjoining
                        // pressure (recomputed below, NOT micro_potential.mu_lR
                        // which already carries the film delta). Both pieces map to
                        // R_u through dsigma'/d(.) = C * C_el^{-1} * d sigma_sw/d(.).
                        double const phi_M_swj =
                            std::get<ProcessLib::ThermoRichardsMechanics::
                                         TransportPorosityData>(
                                this->current_states_[ip])
                                .phi;
                        double const n_S_swj = std::max(1e-16, 1.0 - phi_M_swj);
                        double const b_swj = alpha;  // Biot b == poroelastic alpha

                        // BARE vdW Pi(n_l) and dPi/dn_l, density mirrored to the
                        // residual's p_L_m choice (micro rho_lR when enabled, bulk
                        // otherwise) and treated density-agnostically in dPi/dn_l
                        // (matching the LIVE Maxwell-conjugate block: dPi/dn_l =
                        // -rho * dmu_lR_bare/dn_l).
                        double const active_nS_bare_swj =
                            computeActiveMicroSolidVolumeFraction(
                                n_l, PotentialExchangeLocalSolveContext{},
                                *potential_exchange_params_ptr);
                        auto const vdw_bare_swj = computeVanDerWaalsMicroPotential(
                            n_l, rho_LR, active_nS_bare_swj,
                            potential_exchange_params_ptr
                                ->micro_solid_density_reference,
                            potential_exchange_params_ptr->hamaker_constant,
                            potential_exchange_params_ptr->specific_surface,
                            microPotentialSignFactorFromParameters(
                                *potential_exchange_params_ptr),
                            // Live K(rho_d): rho_d = rho_SR*(1-phi) (total
                            // porosity in scope); off -> parse scalar.
                            // NOTE (2026-06-12): the live-K dK/dphi chain is
                            // NOT added in this default-OFF block (dead code,
                            // enable_dsm_swelling_up_jacobian=false); wire it
                            // per the live p-u block above if ever enabled.
                            effectiveAugmentationPrefactor(
                                *potential_exchange_params_ptr,
                                phi),  // K [J/kg]
                            potential_exchange_params_ptr
                                ->potential_augmentation_exponent,
                            0.0 /*dnS_dnl*/,
                            potential_exchange_params_ptr
                                ->micro_water_content_floor);
                        double rho_d_film_swj = rho_LR;
                        if (potential_exchange_params_ptr
                                ->use_micro_liquid_density_for_micro_pressure)
                        {
                            double const rho_lR_state_swj =
                                *std::get<MicroLiquidDensity>(
                                    this->current_states_[ip]);
                            if (std::isfinite(rho_lR_state_swj) &&
                                rho_lR_state_swj > 0.0)
                            {
                                rho_d_film_swj = rho_lR_state_swj;
                            }
                        }
                        double const Pi_bare_swj =
                            -rho_d_film_swj * vdw_bare_swj.mu_lR;
                        double const dPi_dnl_swj =
                            -rho_d_film_swj * vdw_bare_swj.dmu_lR_dnl;

                        // p_film = Pi - b*p_conf, with p_conf = p_conf_assembly
                        // (= -tr(sigma_eff)/3, current GP effective stress, the
                        // SAME value threaded into the residual via
                        // local_context.confining_pressure_p_conf above); NaN
                        // sentinel -> drop the drain.
                        double const p_conf_drain_swj =
                            std::isfinite(p_conf_assembly) ? p_conf_assembly : 0.0;
                        double const p_film_swj =
                            Pi_bare_swj - b_swj * p_conf_drain_swj;

                        // u-p (via n_l): d sigma_sw/d n_l * dn_l/dpL on identity2,
                        // with d sigma_sw/d n_l = -(1-phi_M)*( p_film + n_l*dPi/dn_l ).
                        double const dsigma_sw_dnl_scalar_swj =
                            -n_S_swj * (p_film_swj + n_l * dPi_dnl_swj);
                        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
                            const d_delta_sigma_sw_dpL =
                                (dsigma_sw_dnl_scalar_swj * dn_l_dpL) * identity2;

                        auto const& C_consistent_swj =
                            *std::get<StiffnessTensor<DisplacementDim>>(
                                constitutive_data);
                        auto const C_el_swj =
                            ip_data_[ip].computeElasticTangentStiffness(
                                variables, t, x_position, dt,
                                this->solid_material_,
                                *this->material_states_[ip]
                                     .material_state_variables);
                        auto const C_el_inv_swj = C_el_swj.inverse().eval();

                        local_Jac
                            .template block<displacement_size, pressure_size>(
                                displacement_index, pressure_index)
                            .noalias() += B.transpose() * C_consistent_swj *
                                          C_el_inv_swj * d_delta_sigma_sw_dpL *
                                          N_p * w;

                        // u-eps (K[u,u]): d sigma_sw/d eps_v =
                        //   +(1-phi_M)*n_l*b*K_drained.
                        // eps_v = identity2^T B u, so the block is
                        //   B^T C C_el^{-1} (dsigma_sw/deps_v * identity2)
                        //        identity2^T B w.
                        double const K_drained_swj =
                            drainedBulkModulusFromStiffness<DisplacementDim>(
                                C_el_swj);
                        double const dsigma_sw_deps_v_scalar_swj =
                            n_S_swj * n_l * b_swj * K_drained_swj;
                        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
                            const dsigma_sw_deps_v =
                                dsigma_sw_deps_v_scalar_swj * identity2;
                        local_Jac
                            .template block<displacement_size, displacement_size>(
                                displacement_index, displacement_index)
                            .noalias() += B.transpose() * C_consistent_swj *
                                          C_el_inv_swj * dsigma_sw_deps_v *
                                          identity2.transpose() * B * w;
                    }
                    else if (enable_dsm_swelling_up_jacobian)
                    {
                        // Current transport porosity phi_M -> REV solid
                        // fraction n_S = 1 - phi_M, matching the residual caller
                        // (updateSwellingStateWithMicroPorosity, line ~1671).
                        double const phi_M_swj =
                            std::get<ProcessLib::ThermoRichardsMechanics::
                                         TransportPorosityData>(
                                this->current_states_[ip])
                                .phi;
                        double const n_S_swj = std::max(1e-16, 1.0 - phi_M_swj);

                        // rho_d and its pL-derivative: mirror the residual's
                        // p_L_m_density choice (micro liquid density when
                        // enabled, bulk otherwise; line ~1618-1622). Inside
                        // mu_lR the density argument is the MICRO rho_lR, so
                        // dmu_lR_drho_lR is paired with the MICRO drho_lR/dpL
                        // here (distinct from the bulk pairing used for the
                        // exchange equation above).
                        double rho_d_swj = rho_LR;
                        double drho_d_dpL_swj = drho_LR_dpL;
                        double dmu_lR_dpL_swj =
                            micro_potential.dmu_lR_dnl * dn_l_dpL +
                            micro_potential.dmu_lR_drho_lR * drho_LR_dpL;
                        if (potential_exchange_params_ptr
                                ->use_micro_liquid_density_for_micro_pressure)
                        {
                            double const rho_lR_state_swj =
                                *std::get<MicroLiquidDensity>(
                                    this->current_states_[ip]);
                            double const drho_lR_dpL_swj =
                                rho_lR_state_swj * beta_LR;
                            rho_d_swj = rho_lR_state_swj;
                            drho_d_dpL_swj = drho_lR_dpL_swj;
                            dmu_lR_dpL_swj =
                                micro_potential.dmu_lR_dnl * dn_l_dpL +
                                micro_potential.dmu_lR_drho_lR *
                                    drho_lR_dpL_swj;
                        }

                        // d(delta_sigma_sw)/dpL scalar on identity2 (product
                        // rule on the three current-pL-dependent factors n_l,
                        // rho_d, mu_lR; -n_l*Pi_curr = +n_l*rho_d*mu_lR).
                        double const d_delta_sigma_sw_dpL_scalar =
                            n_S_swj *
                            (dn_l_dpL * rho_d_swj * micro_potential.mu_lR +
                             n_l * drho_d_dpL_swj * micro_potential.mu_lR +
                             n_l * rho_d_swj * dmu_lR_dpL_swj);

                        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
                            const d_delta_sigma_sw_dpL =
                                d_delta_sigma_sw_dpL_scalar * identity2;

                        // Consistent tangent C (as fetched at line ~3256) and
                        // the elastic tangent C_el (reconstructed exactly as at
                        // line ~3773; not a local in this function). For a
                        // linear-elastic solid C == C_el so C*C_el^{-1} ==
                        // Identity and this block reduces to B^T *
                        // d(delta_sigma_sw)/dpL * N_p * w; the remap factor only
                        // matters for nonlinear tangents.
                        auto const& C_consistent_swj =
                            *std::get<StiffnessTensor<DisplacementDim>>(
                                constitutive_data);
                        auto const C_el_swj =
                            ip_data_[ip].computeElasticTangentStiffness(
                                variables, t, x_position, dt,
                                this->solid_material_,
                                *this->material_states_[ip]
                                     .material_state_variables);

                        local_Jac
                            .template block<displacement_size, pressure_size>(
                                displacement_index, pressure_index)
                            .noalias() += B.transpose() * C_consistent_swj *
                                          C_el_swj.inverse() *
                                          d_delta_sigma_sw_dpL * N_p * w;
                    }
                }
            }

            auto const potential_exchange_result = computePotentialExchangeUpdate(
                alpha_bar, mu, p_L_ip, p_L_m, rho_LR, beta_LR,
                rho_lR_exchange_input, drho_lR_exchange_input_dpL,
                pressure_tolerance, potential_exchange_enabled,
                use_vdw_micro_potential_for_active_exchange, mu_lR_vdw,
                dmu_lR_vdw_drho_lR, use_custom_dmu_lR_vdw_dpL,
                dmu_lR_vdw_dpL,
                use_fd_jacobian_for_direct_macro_derivative,
                fd_jacobian_perturbation);
            // DIAGNOSTIC B: variant A values are passed through unchanged.
            double const rho_L_hat_residual =
                ceiling_B_active
                    ? ceiling_B_rho_L_hat
                    : kkt_active
                          ? kkt_rho_L_hat
                          : potential_exchange_result.exchange.rho_L_hat;
            double const drho_L_hat_dpL_jacobian =
                ceiling_B_active
                    ? ceiling_B_drho_L_hat_dpL
                    : kkt_active
                          ? kkt_drho_L_hat_dpL
                          : potential_exchange_result.drho_L_hat_dpL_direct;
            local_rhs.template segment<pressure_size>(pressure_index)
                .noalias() += N_p.transpose() * rho_L_hat_residual * w;

            // Direct macro Jacobian term for the exchange source. In analytic
            // mode this includes the implicit n_l(p_L) chain contribution.
            local_Jac
                .template block<pressure_size, pressure_size>(pressure_index,
                                                              pressure_index)
                .noalias() -= N_p.transpose() * drho_L_hat_dpL_jacobian * N_p * w;
            if (kkt_active && kkt_drho_L_hat_deps_v != 0.0)
            {
                // KKT active p-u entry, same sign/shape as the B entry below
                // (-=), but into the separate matrix: the assembly line
                // `local_Jac.pu = Kpu/dt` erases local_Jac.pu (Q9). Added after
                // that line at micro_ceiling_pu_tangent = kkt_active /
                // all_exchange; unused (overwritten) otherwise.
                kkt_Kpu_exchange.noalias() -= N_p.transpose() *
                                              kkt_drho_L_hat_deps_v *
                                              identity2.transpose() * B * w;
            }
            if (ceiling_B_active && ceiling_B_drho_L_hat_deps_v != 0.0)
            {
                // Same sign/shape as the exchange p-u blocks above (-=).
                local_Jac
                    .template block<pressure_size, displacement_size>(
                        pressure_index, displacement_index)
                    .noalias() -= N_p.transpose() *
                                  ceiling_B_drho_L_hat_deps_v *
                                  identity2.transpose() * B * w;
            }

            // Keep the microscale pressure-state sensitivity lagged via the
            // secant term only in the placeholder microscale path. In the
            // vdW+ n_l opt-in path this term is intentionally omitted because
            // the active microscale potential is no longer p_L_m/rho_LR.
            if (!use_vdw_micro_potential_for_active_exchange &&
                p_cap_ip != p_cap_prev_ip)
            {
                requirePositiveViscosity(
                    "RichardsMechanics local secant exchange Jacobian assembly",
                    mu);
                auto const p_L_m_prev = **std::get<PrevState<MicroPressure>>(
                    this->prev_states_[ip]);
                local_Jac
                    .template block<pressure_size, pressure_size>(
                        pressure_index, pressure_index)
                    .noalias() += N_p.transpose() * alpha_bar / mu *
                                  (p_L_m - p_L_m_prev) /
                                  (p_cap_ip - p_cap_prev_ip) * N_p * w;
            }
        }
    }

    if (this->process_data_.apply_mass_lumping)
    {
        storage_p_a_p = storage_p_a_p.colwise().sum().eval().asDiagonal();
        storage_p_a_S = storage_p_a_S.colwise().sum().eval().asDiagonal();
        storage_p_a_S_Jpp =
            storage_p_a_S_Jpp.colwise().sum().eval().asDiagonal();
    }

    // pressure equation, pressure part.
    local_Jac
        .template block<pressure_size, pressure_size>(pressure_index,
                                                      pressure_index)
        .noalias() += laplace_p + storage_p_a_p / dt + storage_p_a_S_Jpp;

    // pressure equation, displacement part.
    // Q9 (DESIGN.md D1): micro_ceiling_pu_tangent. overwritten (default, the
    // shipped line, bitwise): every exchange p-u entry accumulated above is
    // erased. all_exchange: Kpu/dt is ADDED to the accumulated exchange entries
    // (repo rule 4.1, accumulators use +=), which also revives the Maxwell,
    // film and live-K p-u entries at inactive points. kkt_active /
    // all_exchange: the KKT-active entries are added after the assignment.
    {
        auto const* const pep_q9 = this->getPotentialExchangeParameters();
        bool const q9_all_exchange =
            pep_q9 != nullptr &&
            pep_q9->micro_ceiling_pu_tangent ==
                MicroCeilingPuTangent::AllExchange;
        bool const q9_add_kkt =
            pep_q9 != nullptr &&
            pep_q9->micro_ceiling_pu_tangent != MicroCeilingPuTangent::Overwritten;
        if (q9_all_exchange)
        {
            local_Jac
                .template block<pressure_size, displacement_size>(
                    pressure_index, displacement_index)
                .noalias() += Kpu / dt;
        }
        else
        {
            local_Jac
                .template block<pressure_size, displacement_size>(
                    pressure_index, displacement_index)
                .noalias() = Kpu / dt;
        }
        if (q9_add_kkt)
        {
            local_Jac
                .template block<pressure_size, displacement_size>(
                    pressure_index, displacement_index)
                .noalias() += kkt_Kpu_exchange;
        }
    }

    // pressure equation
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p_L +
        (storage_p_a_p + storage_p_a_S) * (p_L - p_L_prev) / dt +
        Kpu * (u - u_prev) / dt;

    // displacement equation
    local_rhs.template segment<displacement_size>(displacement_index)
        .noalias() += Kup * p_L;

    // ── Route-B debug flag micro_ceiling_fd_check (DESIGN.md 3.7; NOT adopted) ──
    // Central-difference check of the ASSEMBLED element Jacobian against the
    // analytic one: the only measurement that can say whether the assembled p-u
    // block is consistent (W-1) and, at micro_ceiling_pu_tangent = overwritten
    // against kkt_active, that Q9 changes exactly the block the derivation says.
    // Not for production runs: 2*(n_p + n_u) extra element assemblies per call.
    // Step (DESIGN.md 4.2): h_j = eps_mach^(1/3) * max(|x_j|, x_scale), x_scale =
    // the largest |value| of the element's own dofs of that block (1 SI unit only
    // when all of them are exactly zero; flagged in the log line). Residual
    // convention: local_rhs = -R, so J_fd = -d(local_rhs)/dx. The perturbed
    // assemblies run the same code with the flag suppressed (thread-local); the
    // element's integration-point states and output data are restored after every
    // perturbation. History-dependent solids carry material state that is not
    // snapshotted: only LinearElasticIsotropic is accepted (FATAL otherwise).
    // The comparison blames nothing on the KKT entries for the known inexact
    // terms (W-4, W-5: phi_s tangents in J_pp, J_pu; compiled-out swelling u-side
    // entries; Bishop and gravity terms): the log states per block, it makes no
    // pass/fail claim.
    {
        auto const* const pep_fd = this->getPotentialExchangeParameters();
        static thread_local bool fd_check_running = false;
        if (pep_fd != nullptr && pep_fd->micro_ceiling_fd_check &&
            !fd_check_running)
        {
            if (dynamic_cast<MaterialLib::Solids::LinearElasticIsotropic<
                    DisplacementDim> const*>(&this->solid_material_) == nullptr)
            {
                OGS_FATAL(
                    "micro_ceiling_fd_check supports LinearElasticIsotropic "
                    "solids only: MFront and history-dependent materials carry "
                    "material state that the in-assembler FD check does not "
                    "snapshot.");
            }
            struct FlagReset
            {
                bool& flag;
                ~FlagReset() { flag = false; }
            };
            fd_check_running = true;
            FlagReset const flag_reset{fd_check_running};

            constexpr int n_dof =
                static_cast<int>(pressure_size + displacement_size);
            Eigen::MatrixXd const J_analytic = local_Jac;
            auto const states_backup = this->current_states_;
            auto const output_backup = this->output_data_;

            // Element class by the status written by the base evaluation above.
            unsigned n_active_ip = 0;
            for (unsigned ip = 0; ip < n_integration_points; ++ip)
            {
                n_active_ip +=
                    *std::get<MicroCeilingStatus>(this->current_states_[ip]) ==
                            static_cast<double>(MicroCeilingKktStatus::Active)
                        ? 1
                        : 0;
            }
            char const* const element_class =
                n_active_ip == 0 ? "none_active"
                                 : (n_active_ip == n_integration_points
                                        ? "all_active"
                                        : "mixed");

            auto const block_scale = [&](bool const pressure_block,
                                         bool& used_fallback)
            {
                double m = 0.0;
                unsigned const begin = pressure_block ? 0 : pressure_size;
                unsigned const end =
                    pressure_block ? pressure_size : pressure_size + displacement_size;
                for (unsigned i = begin; i < end; ++i)
                {
                    m = std::max(m, std::abs(local_x[i]));
                }
                used_fallback = !(m > 0.0);
                return used_fallback ? 1.0 : m;
            };
            bool fallback_p = false;
            bool fallback_u = false;
            double const scale_p = block_scale(true, fallback_p);  // [Pa]
            double const scale_u = block_scale(false, fallback_u);  // [m]
            double const h_rel =
                std::cbrt(std::numeric_limits<double>::epsilon());

            Eigen::MatrixXd J_fd = Eigen::MatrixXd::Zero(n_dof, n_dof);
            for (int j = 0; j < n_dof; ++j)
            {
                double const x_scale =
                    static_cast<unsigned>(j) < pressure_size ? scale_p : scale_u;
                double const h = h_rel * std::max(std::abs(local_x[j]), x_scale);
                auto x_plus = local_x;
                auto x_minus = local_x;
                x_plus[j] += h;
                x_minus[j] -= h;
                // the representable steps actually taken
                double const h_total = (x_plus[j] - local_x[j]) +
                                       (local_x[j] - x_minus[j]);
                std::vector<double> rhs_plus;
                std::vector<double> rhs_minus;
                std::vector<double> jac_scratch;
                this->assembleWithJacobian(t, dt, x_plus, local_x_prev,
                                           rhs_plus, jac_scratch);
                this->current_states_ = states_backup;
                this->output_data_ = output_backup;
                this->assembleWithJacobian(t, dt, x_minus, local_x_prev,
                                           rhs_minus, jac_scratch);
                this->current_states_ = states_backup;
                this->output_data_ = output_backup;
                for (int i = 0; i < n_dof; ++i)
                {
                    J_fd(i, j) = -(rhs_plus[i] - rhs_minus[i]) / h_total;
                }
            }

            auto const report = [&](char const* const name, int const r0,
                                    int const c0, int const nr, int const nc)
            {
                auto const a = J_analytic.block(r0, c0, nr, nc);
                auto const f = J_fd.block(r0, c0, nr, nc);
                double const max_abs_dev = (a - f).cwiseAbs().maxCoeff();
                double const max_abs_an = a.cwiseAbs().maxCoeff();
                double const max_abs_fd = f.cwiseAbs().maxCoeff();
                double const denominator = std::max(max_abs_an, max_abs_fd);
                double const max_rel_dev =
                    denominator > 0.0 ? max_abs_dev / denominator : 0.0;
                INFO(
                    "KKT-FD t={:.10g} dt={:.6g} elem={} class={} "
                    "active_ip={}/{} block={} max_abs_dev={:.6e} "
                    "max_rel_dev={:.6e} (relative to the larger of the two block "
                    "maxima) max_abs_analytic={:.6e} max_abs_fd={:.6e} "
                    "h_rel={:.3e} x_scale_p={:.3e}{} x_scale_u={:.3e}{} "
                    "pu_tangent={}",
                    t, dt, this->element_.getID(), element_class, n_active_ip,
                    n_integration_points, name, max_abs_dev, max_rel_dev,
                    max_abs_an, max_abs_fd, h_rel, scale_p,
                    fallback_p ? "(fallback 1)" : "", scale_u,
                    fallback_u ? "(fallback 1)" : "",
                    toString(pep_fd->micro_ceiling_pu_tangent));
            };
            constexpr int np = static_cast<int>(pressure_size);
            constexpr int nu = static_cast<int>(displacement_size);
            report("pp", 0, 0, np, np);
            report("pu", 0, np, np, nu);
            report("up", np, 0, nu, np);
            report("uu", np, np, nu, nu);
        }
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianForPressureEquations(
        const double /*t*/, double const /*dt*/,
        Eigen::VectorXd const& /*local_x*/,
        Eigen::VectorXd const& /*local_x_prev*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/)
{
    OGS_FATAL("RichardsMechanics; The staggered scheme is not implemented.");
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        const double /*t*/, double const /*dt*/,
        Eigen::VectorXd const& /*local_x*/,
        Eigen::VectorXd const& /*local_x_prev*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/)
{
    OGS_FATAL("RichardsMechanics; The staggered scheme is not implemented.");
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(double const t, double const dt,
                                           Eigen::VectorXd const& local_x,
                                           Eigen::VectorXd const& local_x_prev,
                                           int const process_id,
                                           std::vector<double>& local_b_data,
                                           std::vector<double>& local_Jac_data)
{
    // For the equations with pressure
    if (process_id == 0)
    {
        assembleWithJacobianForPressureEquations(t, dt, local_x, local_x_prev,
                                                 local_b_data, local_Jac_data);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(t, dt, local_x, local_x_prev,
                                                local_b_data, local_Jac_data);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& local_x_prev)
{
    auto const [p_L, u] = localDOF(local_x);
    auto const [p_L_prev, u_prev] = localDOF(local_x_prev);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    auto const& medium =
        this->process_data_.media_map.getMedium(this->element_.getID());
    auto const& liquid_phase = medium->phase(MaterialPropertyLib::PhaseName::AqueousLiquid);
    auto const& solid_phase = medium->phase(MaterialPropertyLib::PhaseName::Solid);
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    double saturation_avg = 0;
    double porosity_avg = 0;

    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    KV sigma_avg = KV::Zero();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N_p = ip_data_[ip].N_p;
        auto const& N_u = ip_data_[ip].N_u;
        auto const& dNdx_u = ip_data_[ip].dNdx_u;

        ParameterLib::SpatialPosition x_position = {
            std::nullopt, this->element_.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, N_u))};
        auto const x_coord = x_position.getCoordinates().value()[0];

        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, this->is_axially_symmetric_);

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        double p_cap_prev_ip;
        NumLib::shapeFunctionInterpolate(-p_L_prev, N_p, p_cap_prev_ip);

        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;

        auto const temperature =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(variables, x_position, t, dt);
        variables.temperature = temperature;

        auto& eps =
            std::get<StrainData<DisplacementDim>>(this->current_states_[ip])
                .eps;
        eps.noalias() = B * u;
        auto& S_L =
            std::get<ProcessLib::ThermoRichardsMechanics::SaturationData>(
                this->current_states_[ip])
                .S_L;
        auto const S_L_prev =
            std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::SaturationData>>(
                this->prev_states_[ip])
                ->S_L;
        S_L = medium->property(MPL::PropertyType::saturation)
                  .template value<double>(variables, x_position, t, dt);
        variables.liquid_saturation = S_L;
        variables_prev.liquid_saturation = S_L_prev;

        auto const chi = [medium, x_position, t, dt](double const S_L)
        {
            MPL::VariableArray vs;
            vs.liquid_saturation = S_L;
            return medium->property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(vs, x_position, t, dt);
        };
        double const chi_S_L = chi(S_L);
        double const chi_S_L_prev = chi(S_L_prev);

        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);
        auto& state_current = this->current_states_[ip];
        variables.stress =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                DisplacementDim>>(state_current)
                .sigma_eff;
        // Set mechanical strain temporary to compute tangent stiffness.
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);
        auto const C_el = ip_data_[ip].computeElasticTangentStiffness(
            variables, t, x_position, dt, this->solid_material_,
            *this->material_states_[ip].material_state_variables);

        auto const beta_SR = (1 - alpha) / this->solid_material_.getBulkModulus(
                                               t, x_position, &C_el);
        variables.grain_compressibility = beta_SR;

        variables.effective_pore_pressure = -chi_S_L * p_cap_ip;
        variables_prev.effective_pore_pressure = -chi_S_L_prev * p_cap_prev_ip;

        // Set volumetric strain rate for the general case without swelling.
        variables.volumetric_strain = Invariants::trace(eps);
        variables_prev.volumetric_strain = Invariants::trace(B * u_prev);

        auto& phi = std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                        this->current_states_[ip])
                        .phi;
        {  // Porosity update
            auto const phi_prev = std::get<PrevState<
                ProcessLib::ThermoRichardsMechanics::PorosityData>>(
                                      this->prev_states_[ip])
                                      ->phi;
            variables_prev.porosity = phi_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables.porosity = phi;
        }

        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        variables.density = rho_LR;
        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);

        {
            // Swelling and possibly volumetric strain rate update.
            auto& sigma_sw =
                std::get<ProcessLib::ThermoRichardsMechanics::
                             ConstitutiveStress_StrainTemperature::
                                 SwellingDataStateful<DisplacementDim>>(
                    this->current_states_[ip]);
            auto const& sigma_sw_prev = std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::
                              ConstitutiveStress_StrainTemperature::
                                  SwellingDataStateful<DisplacementDim>>>(
                this->prev_states_[ip]);
            auto const transport_porosity_prev = std::get<PrevState<
                ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(
                this->prev_states_[ip]);
            auto const phi_prev = std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::PorosityData>>(
                this->prev_states_[ip]);
            auto& transport_porosity = std::get<
                ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(
                this->current_states_[ip]);
            auto& p_L_m = std::get<MicroPressure>(this->current_states_[ip]);
            auto const p_L_m_prev =
                std::get<PrevState<MicroPressure>>(this->prev_states_[ip]);
            auto& S_L_m = std::get<MicroSaturation>(this->current_states_[ip]);
            auto const S_L_m_prev =
                std::get<PrevState<MicroSaturation>>(this->prev_states_[ip]);

            updateSwellingStressAndVolumetricStrain<DisplacementDim>(
                *medium, solid_phase, C_el, rho_LR, mu,
                this->process_data_.micro_porosity_parameters,
                this->getPotentialExchangeParameters(), alpha, phi, p_cap_ip,
                variables, variables_prev, x_position, t, dt, sigma_sw,
                sigma_sw_prev, transport_porosity_prev, phi_prev,
                transport_porosity, p_L_m_prev, S_L_m_prev, p_L_m, S_L_m);
        }

        auto const transport_porosity_prev_value = std::get<PrevState<
            ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(
            this->prev_states_[ip])
                                                        ->phi;
        auto const phi_m_prev_value =
            **std::get<PrevState<MicroPorosity>>(
                this->prev_states_[ip]);

        // Film-pressure coupling: supply p_conf = -tr(sigma_eff)/3 to the n_l
        // local solve (assemble path). NaN sentinel keeps the term off when the
        // flag is OFF.
        double const p_conf_micro_solve =
            isFilmPressureCouplingEnabled(
                this->getPotentialExchangeParameters())
                ? -std::get<ProcessLib::ConstitutiveRelations::
                                EffectiveStressData<DisplacementDim>>(
                       this->current_states_[ip])
                       .sigma_eff.dot(identity2) /
                      3.0
                : std::numeric_limits<double>::quiet_NaN();
        // Drained bulk modulus for the INTEGRABLE Maxwell partner (mu_lR_mech),
        // assembleWithJacobian micro-solve path. Only under film coupling; NaN
        // sentinel otherwise -> partner inert, flag-off bit-for-bit. C_el is the
        // elastic stiffness evaluated above in this ip loop (line ~4761).
        double const K_drained_micro_solve =
            isFilmPressureCouplingEnabled(
                this->getPotentialExchangeParameters())
                ? drainedBulkModulusFromStiffness<DisplacementDim>(C_el)
                : std::numeric_limits<double>::quiet_NaN();
        // KKT iteration trace (DESIGN.md 3.8), output re-evaluation flagged.
        MicroCeilingTraceTag trace_tag_storage;
        MicroCeilingTraceTag const* trace_tag = nullptr;
        if (auto const* const pep_trace = this->getPotentialExchangeParameters();
            pep_trace != nullptr &&
            !pep_trace->micro_ceiling_trace_elements.empty() &&
            std::find(pep_trace->micro_ceiling_trace_elements.begin(),
                      pep_trace->micro_ceiling_trace_elements.end(),
                      static_cast<std::size_t>(this->element_.getID())) !=
                pep_trace->micro_ceiling_trace_elements.end())
        {
            trace_tag_storage = {
                .element_id = static_cast<std::size_t>(this->element_.getID()),
                .integration_point = static_cast<std::size_t>(ip),
                .t = t,
                .dt = dt,
                .output_reevaluation = true};
            trace_tag = &trace_tag_storage;
        }
        updateMicroscaleHydraulicState<DisplacementDim>(
            this->current_states_[ip], this->prev_states_[ip], p_cap_ip,
            rho_LR, mu, dt, t, variables, variables_prev,
            {.phi = phi,
             .phi_M_prev = transport_porosity_prev_value,
             .phi_m_prev = phi_m_prev_value,
             .volumetric_strain = variables.volumetric_strain,
             .volumetric_strain_prev = variables_prev.volumetric_strain,
             .confining_pressure_p_conf = p_conf_micro_solve,
             .biot_coefficient = alpha,
             .drained_bulk_modulus = K_drained_micro_solve},
            this->process_data_.micro_porosity_parameters,
            this->getPotentialExchangeParameters(), trace_tag);
        updatePorositySplitState<DisplacementDim>(
            this->current_states_[ip], this->prev_states_[ip], phi, variables,
            variables_prev, this->getPotentialExchangeParameters());
        updateTotalPorosityState<DisplacementDim>(
            this->current_states_[ip], this->prev_states_[ip], phi, variables,
            variables_prev, this->getPotentialExchangeParameters());
        updateSwellingState<DisplacementDim>(
            solid_phase, rho_LR, C_el, this->current_states_[ip],
            this->prev_states_[ip], variables, variables_prev, x_position, t,
            dt, this->getPotentialExchangeParameters(),
            /*biot_coefficient=*/alpha);

        // Gate 1/2 fix for DSM micro-porosity mode: enforce phi_m <= phi_total
        // and phi_M = phi_total - phi_m >= 0. When the
        // micro water content n_l approaches the total porosity phi under
        // confinement, the hierarchical split can produce phi_M < 0. Cap
        // micro porosity at the total and set phi_M = phi - phi_m >= 0.
        // This is the output-field clamp; the constitutive root cause (missing
        // micro-swelling saturation) is tracked separately.
        if (this->process_data_.micro_porosity_parameters.has_value())
        {
            auto& phi_M_out =
                std::get<ProcessLib::ThermoRichardsMechanics::
                             TransportPorosityData>(
                    this->current_states_[ip])
                    .phi;
            auto& phi_m_out =
                *std::get<MicroPorosity>(this->current_states_[ip]);
            double const phi_total_out =
                std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                    this->current_states_[ip])
                    .phi;
            phi_m_out = std::min(phi_m_out, phi_total_out);
            phi_M_out = phi_total_out - phi_m_out;  // >= 0
            variables.transport_porosity = phi_M_out;
        }

        if (medium->hasProperty(MPL::PropertyType::transport_porosity))
        {
            if (!medium->hasProperty(MPL::PropertyType::saturation_micro) &&
                !isPotentialExchangeEnabled(
                    this->getPotentialExchangeParameters()))
            {
                auto& transport_porosity =
                    std::get<ProcessLib::ThermoRichardsMechanics::
                                 TransportPorosityData>(
                        this->current_states_[ip])
                        .phi;
                auto const transport_porosity_prev =
                    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                                           TransportPorosityData>>(
                        this->prev_states_[ip])
                        ->phi;

                variables_prev.transport_porosity = transport_porosity_prev;

                transport_porosity =
                    medium->property(MPL::PropertyType::transport_porosity)
                        .template value<double>(variables, variables_prev,
                                                x_position, t, dt);
                variables.transport_porosity = transport_porosity;
            }
            // Pi-path and legacy-DSM modes: phi_M already set by the porosity
            // split / Gate fix above; no property evaluation needed.
        }
        else
        {
            // No transport_porosity medium property. In Pi-path / DSM mode
            // variables.transport_porosity is already phi_M from the Gate fix above.
            // Only fall back to total porosity for plain RM (no micro-porosity split).
            if (!isPotentialExchangeEnabled(
                    this->getPotentialExchangeParameters()) &&
                !medium->hasProperty(MPL::PropertyType::saturation_micro))
            {
                variables.transport_porosity = phi;
            }
        }

        auto const& sigma_eff =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                DisplacementDim>>(this->current_states_[ip])
                .sigma_eff;

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        {
            auto const sigma_total =
                (sigma_eff + alpha * chi_S_L * identity2 * p_cap_ip).eval();
            // For stress dependent permeability.
            variables.total_stress.emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    sigma_total));
        }

        variables.equivalent_plastic_strain =
            this->material_states_[ip]
                .material_state_variables->getEquivalentPlasticStrain();

        auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);

        std::get<
            ProcessLib::ThermoRichardsMechanics::PermeabilityData<DisplacementDim>>(
            this->output_data_[ip])
            .Ki = K_intrinsic;
        std::get<
            ProcessLib::ThermoRichardsMechanics::PermeabilityData<DisplacementDim>>(
            this->output_data_[ip])
            .k_rel = k_rel;

        GlobalDimMatrixType const K_over_mu = k_rel * K_intrinsic / mu;

        double const p_FR = -chi_S_L * p_cap_ip;
        // p_SR
        variables.solid_grain_pressure =
            p_FR - sigma_eff.dot(identity2) / (3 * (1 - phi));
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        *std::get<DrySolidDensity>(this->output_data_[ip]) = (1 - phi) * rho_SR;

        {
            auto& state_current = this->current_states_[ip];
            auto const& sigma_sw =
                std::get<ProcessLib::ThermoRichardsMechanics::
                             ConstitutiveStress_StrainTemperature::
                                 SwellingDataStateful<DisplacementDim>>(state_current)
                    .sigma_sw;
            auto& eps_m =
                std::get<ProcessLib::ConstitutiveRelations::
                             MechanicalStrainData<DisplacementDim>>(state_current)
                    .eps_m;
            bool const swelling_stress_active =
                solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate) ||
                isPotentialExchangeEnabled(
                    this->getPotentialExchangeParameters());
            eps_m.noalias() = swelling_stress_active
                                  ? eps + C_el.inverse() * sigma_sw
                                  : eps;
            variables.mechanical_strain.emplace<
                MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);
        }

        {
            auto& state_current = this->current_states_[ip];
            auto const& state_previous = this->prev_states_[ip];
            auto& sigma_eff =
                std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                    DisplacementDim>>(state_current);
            auto const& sigma_eff_prev =
                std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                       EffectiveStressData<DisplacementDim>>>(
                    state_previous);
            auto const& eps_m =
                std::get<ProcessLib::ConstitutiveRelations::
                             MechanicalStrainData<DisplacementDim>>(state_current);
            auto const& eps_m_prev =
                std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                       MechanicalStrainData<DisplacementDim>>>(
                    state_previous);

            ip_data_[ip].updateConstitutiveRelation(
                variables, t, x_position, dt, temperature, sigma_eff,
                sigma_eff_prev, eps_m, eps_m_prev, this->solid_material_,
                this->material_states_[ip].material_state_variables);
        }

        auto const& b = this->process_data_.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = ip_data_[ip].dNdx_p;
        std::get<
            ProcessLib::ThermoRichardsMechanics::DarcyLawData<DisplacementDim>>(
            this->output_data_[ip])
            ->noalias() = -K_over_mu * dNdx_p * p_L + rho_LR * K_over_mu * b;

        saturation_avg += S_L;
        porosity_avg += phi;
        sigma_avg += sigma_eff;
    }
    saturation_avg /= n_integration_points;
    porosity_avg /= n_integration_points;
    sigma_avg /= n_integration_points;

    (*this->process_data_.element_saturation)[this->element_.getID()] =
        saturation_avg;
    (*this->process_data_.element_porosity)[this->element_.getID()] =
        porosity_avg;

    Eigen::Map<KV>(
        &(*this->process_data_.element_stresses)[this->element_.getID() *
                                                 KV::RowsAtCompileTime]) =
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_avg);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(this->element_, this->is_axially_symmetric_, p_L,
                         *this->process_data_.pressure_interpolated);
}
}  // namespace RichardsMechanics
}  // namespace ProcessLib
