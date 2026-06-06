// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <cmath>
#include <Eigen/LU>
#include <cassert>
#include <limits>
#include <mutex>

#include "BaseLib/Logging.h"
#include "ComputeMicroPorosity.h"
#include "ConstitutiveRelations/ConstitutiveModels.h"
#include "ConstitutiveRelations/PotentialExchange.h"
#include "IntegrationPointData.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
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

    return {
        .rho_lR = rho_lR,
        .omega_l = omega_l,
        .drho_lR_dnl = drho_lR_dnl,
        .drho_l_dn_l = rho_lR + n_l_safe * drho_lR_dnl,
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
    // Pi = -rho_lR_used * mu_lR_vdw > 0 (mirrors p_L_m and the eigenstress; the
    // SAME density used inside the vdW formula at this call site).
    double const Pi = -rho_lR_used * out.mu_lR;
    double const one_minus_nl = std::max(1e-12, 1.0 - n_l);
    double const n_S_rev =
        std::isfinite(local_context.phi)
            ? std::clamp((1.0 - local_context.phi) / one_minus_nl, 0.0, 1.0)
            : active_nS;  // = 1 - phi_M (REV macro-solid fraction)
    auto const film = computeFilmPressureMicroPotential(
        Pi, local_context.confining_pressure_p_conf, rho_lR_used,
        std::max(1e-16, n_S_rev), n_l, out.dmu_lR_dnl,
        potential_exchange_params.film_pressure_gate_width,
        local_context.biot_coefficient);
    out.mu_lR += film.mu_lR_film_delta;
    out.dmu_lR_dnl += film.dmu_lR_film_dnl;
    out.dmu_lR_drho_lR += film.dmu_lR_film_drho_lR;
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
            potential_exchange_params.potential_augmentation_prefactor,
            potential_exchange_params.potential_augmentation_exponent);
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
                                dt_safe * rho_l * volumetric_strain_rate;
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
                                dt_safe * drho_l_REV_dn_l *
                                    volumetric_strain_rate;
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

inline MicroMacroMassStorageCoupledSolveData
solveReferenceMassStorageCoupledState(
    double const n_l_prev, double const rho_l_prev, double const rho_lR_prev,
    double const dt, double const rho_LR, double const alpha_bar,
    double const mu, YoungLaplaceMacroPotentialData const& macro_potential,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    requirePositiveViscosity("solveReferenceMassStorageCoupledState", mu);
    constexpr double n_l_floor = 1e-16;
    constexpr double rho_floor = 1e-16;
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
            potential_exchange_params.potential_augmentation_prefactor,
            potential_exchange_params.potential_augmentation_exponent);
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
                                     dt_safe * rho_l * volumetric_strain_rate;
        auto const density = computeReducedMicroLiquidDensity(
            n_l, rho_LR, active_nS, potential_exchange_params);
        double const density_residual = rho_lR - density.rho_lR;
        return std::tuple{mass_residual, density_residual, micro_potential,
                          exchange};
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

        double const J11 = (r1_n_plus - r1_n_minus) / denom_n;
        double const J21 = (r2_n_plus - r2_n_minus) / denom_n;
        double const J12 = (r1_rho_plus - r1_rho_minus) / denom_rho;
        double const J22 = (r2_rho_plus - r2_rho_minus) / denom_rho;

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
            potential_exchange_params.potential_augmentation_prefactor,
            potential_exchange_params.potential_augmentation_exponent, dnS_dnl);

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
            residual -= dt_safe * rho_l * volumetric_strain_rate;

            jacobian = micro_liquid_density->drho_l_dn_l -
                       dt_safe * drho_l_hat_dn_l;
            jacobian -=
                dt_safe * micro_liquid_density->drho_l_dn_l *
                volumetric_strain_rate;
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
                dt_safe * rho_l_prev_fallback * volumetric_strain_rate /
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
        double const time_factor = 1.0 - dt_safe * eps_v_rate;

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

template <int DisplacementDim>
inline void updateMicroscaleHydraulicState(
    StatefulData<DisplacementDim>& state_current,
    StatefulDataPrev<DisplacementDim> const& state_previous, double const p_cap_ip,
    double const rho_LR, double const mu, double const dt,
    MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
    PotentialExchangeLocalSolveContext const& local_context,
    std::optional<MicroPorosityParameters> const& micro_porosity_parameters,
    PotentialExchangeParameters const* const potential_exchange_parameters)
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
        auto const coupled_update = solveReferenceMassStorageCoupledState(
            n_l_prev_value, rho_l_prev, rho_lR_prev_value, dt,
            rho_LR, micro_porosity_parameters->mass_exchange_coefficient, mu,
            macro_potential, local_context, potential_exchange_params);
        applyReferenceMassStorageLocalState<DisplacementDim>(
            state_current, state_previous, variables, variables_prev, rho_LR, local_context, potential_exchange_params,
            coupled_update);
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
    double const p_conf = std::numeric_limits<double>::quiet_NaN())
{
    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    auto const& params = potential_exchange_params;

    // C_el is unused on BOTH branches now (the film-ON branch is a transmitted
    // PRESSURE over the contact fraction, no longer an elastic eigenstress, so it
    // no longer needs a drained K here; the OFF disjoining-eigenstress path never
    // touched it). Kept in the signature for call-site stability.
    (void)C_el;

    KV delta_sigma_sw = KV::Zero();
    double const delta_n_l = n_l - n_l_prev;
    if (!(std::isfinite(delta_n_l) &&
          std::abs(delta_n_l) > std::numeric_limits<double>::epsilon()))
    {
        return delta_sigma_sw;
    }

    // ── Film-pressure swelling stress (maxwell sec.5, flag ON) ──────────────
    // CORRECTION (Vinay, 2026-06-06): the micro swelling stress is a transmitted
    // PRESSURE over the solid/aggregate contact fraction, NOT an elastic
    // eigenstress. The previous K_sw*b*eps_sw^m form routed a pressure through an
    // elastic strain (the bulk modulus K_sw does NOT belong in a residual
    // stress); it is replaced by
    //   sigma_sw(n) = -(1 - phi_M)*p_film,   p_film = Pi(n_l) - b*p_conf,
    // tension-positive, with the (1 - phi_M) = n_S contact-area factor. Here Pi is
    // the BARE van-der-Waals disjoining pressure at n_l (BEFORE the film delta),
    //   Pi(n_l) = -rho_lR_used * mu_lR_vdw(n_l)   (FULL p^disj, density mirrored to
    //   the hydraulic p_L_m choice exactly as the OFF branch below),
    // b = biot_coefficient and p_conf = -tr(sigma_eff)/3 (confining pressure,
    // supplied by the caller as local_context.confining_pressure_p_conf).
    //
    // Increment over the step:
    //   delta_sigma_sw = sigma_sw(n) - sigma_sw_prev,
    //   sigma_sw_prev  = -(1 - phi_M)*(Pi(n_l_prev) - b*p_conf_prev).
    // p_conf_prev is NOT available in this function (only the current p_conf is
    // threaded). We therefore HOLD p_conf FIXED across the step (use the current
    // p_conf for both terms): the -b*p_conf contribution then cancels EXACTLY in
    // the increment, so the Pi part is exact and the p_conf->stress coupling
    // enters only through the assembly's outer (Newton) iteration as p_conf
    // updates step to step. Consequently the increment is independent of p_conf
    // (and well-defined even when p_conf is the NaN sentinel for default callers /
    // GP tests):
    //   delta_sigma_sw = -(1 - phi_M)*(Pi(n_l) - Pi(n_l_prev))*identity2.
    // K_sw is GONE from this residual entirely.
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
        double const mu_lR_prev_film =
            computeVanDerWaalsMicroPotential(
                n_l_prev, rho_lR_prev, active_nS_prev_film,
                params.micro_solid_density_reference, params.hamaker_constant,
                params.specific_surface, sign_factor_film,
                params.potential_augmentation_prefactor,
                params.potential_augmentation_exponent)
                .mu_lR;
        double const mu_lR_curr_film =
            computeVanDerWaalsMicroPotential(
                n_l, rho_lR, active_nS_curr_film,
                params.micro_solid_density_reference, params.hamaker_constant,
                params.specific_surface, sign_factor_film,
                params.potential_augmentation_prefactor,
                params.potential_augmentation_exponent)
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

        (void)biot_coefficient;  // b enters only via the p_conf term, which
        (void)p_conf;            // cancels in the (p_conf-fixed) increment.

        // n_S passed in is the REV macro-solid fraction (1 - phi_M).
        // delta_sigma_sw = -(1 - phi_M)*(Pi_curr - Pi_prev)*identity2.
        delta_sigma_sw.noalias() +=
            -n_S * (Pi_curr_film - Pi_prev_film) * identity2_film;
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
            params.potential_augmentation_prefactor,
            params.potential_augmentation_exponent)
            .mu_lR;
    double const mu_lR_curr =
        computeVanDerWaalsMicroPotential(
            n_l, rho_lR, active_nS_curr,
            params.micro_solid_density_reference, params.hamaker_constant,
            params.specific_surface, sign_factor,
            params.potential_augmentation_prefactor,
            params.potential_augmentation_exponent)
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
    double const p_conf = std::numeric_limits<double>::quiet_NaN())
{
    return computeReferenceMicroPorositySwellingStressIncrement<DisplacementDim>(
        n_l_prev, n_l, n_S, rho_lR, rho_lR_prev, rho_LR, C_el,
        potential_exchange_params, biot_coefficient, p_conf);
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
            potential_exchange_params, biot_coefficient, p_conf_swelling);

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
        // Volumetric average specific storage of the solid and fluid phases.
        double const specific_storage =
            DeltaS_L_Deltap_cap * (p_cap_ip * a0 - phi) +
            S_L * (phi * beta_LR + a0);
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
                PotentialExchangeLocalSolveContext const local_solve_context{
                    .phi = phi,
                    .phi_M_prev = transport_porosity_prev,
                    .phi_m_prev = phi_m_prev,
                    .volumetric_strain = variables.volumetric_strain,
                    .volumetric_strain_prev = variables_prev.volumetric_strain,
                    .confining_pressure_p_conf = p_conf_assembly,
                    .biot_coefficient = alpha};
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
                    double const p_conf_mc =
                        -std::get<ProcessLib::ConstitutiveRelations::
                                      EffectiveStressData<DisplacementDim>>(
                             this->current_states_[ip])
                             .sigma_eff.dot(identity2) /
                        3.0;  // Pa, confining pressure = -tr(sigma')/3
                    double const one_minus_nl_mc = std::max(1e-12, 1.0 - n_l);
                    double const n_S_mc = std::clamp(
                        (1.0 - phi) / one_minus_nl_mc, 0.0, 1.0);  // = 1 - phi_M
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
                    double const S1_mc =
                        -n_S_mc * (Pi_mc + n_l * dPi_dnl_mc);  // Pa
                    // Gate-scale = option 1 (decision #4, 2026-06-02): the gate
                    // threshold is the REV-consistent PARTIAL stress phi_m*Pi =
                    // n_S*n_l*Pi, NOT the intrinsic micro Pi. p'(REV) reaches
                    // phi_m*Pi but never Pi (confined probe 7.15<16.6; EPFL
                    // path2 10.4>9.2). Trigger only -- S1_mc keeps the full Pi
                    // for the term magnitude. See tex/maxwell_gate_scale.
                    double const Pi_gate_mc = n_S_mc * n_l * Pi_mc;  // = phi_m*Pi
                    auto const mc = computeMaxwellConjugateMicroPotential(
                        S1_mc, /*dS1_dnl=*/0.0, variables.volumetric_strain,
                        p_conf_mc, Pi_gate_mc, rho_mc, /*n_S=*/n_S_mc);
                    mu_lR_vdw += mc.mu_lR_mech;  // J/kg, additive; ==0 below gate
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
                N_p.transpose() * potential_exchange_result.exchange.rho_L_hat * w;
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
            material_state_data)
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
    updateMicroscaleHydraulicState<DisplacementDim>(
        state_current, state_previous, p_cap_ip, rho_LR, mu, dt, variables, variables_prev,
        {.phi = phi,
         .phi_M_prev = transport_porosity_prev_value,
         .phi_m_prev = phi_m_prev_value,
         .volumetric_strain = variables.volumetric_strain,
         .volumetric_strain_prev = variables_prev.volumetric_strain,
         .confining_pressure_p_conf = p_conf_micro_solve,
         .biot_coefficient = alpha},
        micro_porosity_parameters, potential_exchange_parameters);
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
            this->solid_material_, this->material_states_[ip]);

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
        double const specific_storage_a_p = S_L * (phi * beta_LR + S_L * a0);
        double const specific_storage_a_S = phi - p_cap_ip * S_L * a0;

        double const dspecific_storage_a_p_dp_cap =
            dS_L_dp_cap * (phi * beta_LR + 2 * S_L * a0);
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
                PotentialExchangeLocalSolveContext const local_solve_context{
                    .phi = phi,
                    .phi_M_prev = transport_porosity_prev,
                    .phi_m_prev = phi_m_prev,
                    .volumetric_strain = variables.volumetric_strain,
                    .volumetric_strain_prev = variables_prev.volumetric_strain,
                    .confining_pressure_p_conf = p_conf_assembly,
                    .biot_coefficient = alpha};
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
                    double const p_conf_mc =
                        -std::get<ProcessLib::ConstitutiveRelations::
                                      EffectiveStressData<DisplacementDim>>(
                             this->current_states_[ip])
                             .sigma_eff.dot(identity2) /
                        3.0;  // Pa, confining pressure = -tr(sigma')/3
                    double const one_minus_nl_mc = std::max(1e-12, 1.0 - n_l);
                    double const n_S_mc = std::clamp(
                        (1.0 - phi) / one_minus_nl_mc, 0.0, 1.0);  // = 1 - phi_M
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
                    double const S1_mc =
                        -n_S_mc * (Pi_mc + n_l * dPi_dnl_mc);  // Pa
                    // Gate-scale = option 1 (decision #4, 2026-06-02): the gate
                    // threshold is the REV-consistent PARTIAL stress phi_m*Pi =
                    // n_S*n_l*Pi, NOT the intrinsic micro Pi. p'(REV) reaches
                    // phi_m*Pi but never Pi (confined probe 7.15<16.6; EPFL
                    // path2 10.4>9.2). Trigger only -- S1_mc keeps the full Pi
                    // for the term magnitude. See tex/maxwell_gate_scale.
                    double const Pi_gate_mc = n_S_mc * n_l * Pi_mc;  // = phi_m*Pi
                    auto const mc = computeMaxwellConjugateMicroPotential(
                        S1_mc, /*dS1_dnl=*/0.0, variables.volumetric_strain,
                        p_conf_mc, Pi_gate_mc, rho_mc, /*n_S=*/n_S_mc);
                    mu_lR_vdw += mc.mu_lR_mech;  // J/kg, additive; ==0 below gate
                    dmu_lR_vdw_drho_lR += mc.dmu_lR_mech_drho_lR;
                    // DSM Maxwell-conjugate: exchange<->displacement tangent --
                    // the transpose partner of the swelling-eigenstress block.
                    // mu_lR_mech makes the exchange depend on eps_v, so the
                    // pressure residual rho_L_hat = alpha_M*(mu_lR - mu_LR) has
                    //   d rho_L_hat/d u = alpha_M * (dmu_lR_mech/deps_v) * m^T B,
                    // m = identity2 (eps_v = m^T B u). ==0 below the gate, so
                    // gate-closed runs stay bit-for-bit. Without this block the
                    // tangent is inconsistent and Newton diverges (Task-13 Pi
                    // blow-up). Sign mirrors the K[p,p] exchange block (-=).
                    if (mc.gate_open && mc.dmu_lR_mech_deps_v != 0.0 && mu > 0.0)
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
                // ── Film-pressure p-u tangent (maxwell sec.5, increment D-i) ──
                // The drain block K_{n_l u}: with mu_lR = -(Pi - b*p_conf)/rho_lR,
                //   dmu_lR/deps_v = (dmu_lR/dp_conf)*(dp_conf/deps_v)
                //                 = (g*b/rho_lR)*(-K) = -(K/rho_lR)*g*b,
                // since p_conf = -sigma'_m and dsigma'_m/deps_v = K (drained bulk
                // modulus). Then d rho_L_hat/d u = alpha_M*(dmu_lR/deps_v)*m^T B,
                // assembled with the SAME sign (-=) as the K[p,p] exchange block.
                // This is the one-Psi transpose partner of the NEW swelling-stress
                // u-eps block d sigma_sw/d eps_v = +(1-phi_M)*b*K (2026-06-06
                // PRESSURE form; the old -(1-phi_M)*K_sw*b eigenstress is retired,
                // K_sw no longer appears anywhere). Here K is the MECHANICAL
                // drained bulk modulus = dp_conf/deps_v -- NOT K_sw. ==0 when the
                // smooth gate is closed (g==0), so a run that never opens the gate
                // keeps a film-free Jacobian.
                if (film_pressure_coupling && mu > 0.0)
                {
                    double const Pi_film = p_L_m;  // disjoining = -rho*mu_lR > 0
                    double const one_minus_nl_film =
                        std::max(1e-12, 1.0 - n_l);
                    double const n_S_film = std::clamp(
                        (1.0 - phi) / one_minus_nl_film, 0.0, 1.0);
                    double const rho_film =
                        (std::isfinite(rho_lR_exchange_input) &&
                         rho_lR_exchange_input > 0.0)
                            ? rho_lR_exchange_input
                            : rho_LR;
                    // Recover the gate value g at this state (same evaluator as
                    // computeActiveMicroPotential used for mu_lR above).
                    auto const film = computeFilmPressureMicroPotential(
                        Pi_film, p_conf_assembly, rho_film,
                        std::max(1e-16, n_S_film), n_l,
                        micro_potential.dmu_lR_dnl,
                        potential_exchange_params_ptr->film_pressure_gate_width,
                        /*biot_b=*/alpha);
                    if (film.gate_open)
                    {
                        // Drained (elastic) bulk modulus = dp_conf/deps_v, the
                        // MECHANICAL K that the NEW swelling-stress u-eps block
                        // (+(1-phi_M)*b*K) also uses, so the p-u and u-eps blocks
                        // are a true one-Psi transpose pair. (Reconstruct C_el
                        // exactly as the enable_dsm_swelling_up_jacobian block
                        // does below.)
                        auto const C_el_film =
                            ip_data_[ip].computeElasticTangentStiffness(
                                variables, t, x_position, dt,
                                this->solid_material_,
                                *this->material_states_[ip]
                                     .material_state_variables);
                        double const K_drained =
                            drainedBulkModulusFromStiffness<DisplacementDim>(
                                C_el_film);
                        // Eigenstrain Biot b == poroelastic Biot alpha (threaded
                        // from the medium biot_coefficient; one-Psi consistency).
                        double const b_film = alpha;
                        // dmu_lR/deps_v = -(K/rho_lR)*g*b.
                        double const dmu_lR_film_deps_v =
                            -(K_drained / rho_film) * film.gate_value * b_film;
                        double const alpha_M_eff_film =
                            alpha_bar * rho_LR / mu;
                        local_Jac
                            .template block<pressure_size, displacement_size>(
                                pressure_index, displacement_index)
                            .noalias() -=
                            N_p.transpose() *
                            (alpha_M_eff_film * dmu_lR_film_deps_v) *
                            identity2.transpose() * B * w;
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
                    double const dn_l_dpL = computeImplicitNlDpL(
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
                    // DEFAULT OFF (verified 2026-06-01): with the corrected
                    // scope (below) this term DOES fire on the Pi-path models,
                    // but on dd1400 it is residual-neutral (P_sw 4.91637 MPa
                    // unchanged) AND convergence-neutral (824 Newton iters
                    // unchanged) because the micro-macro coupling here is tiny
                    // (dn_l_dpL ~ 3e-13, d(sigma_sw)/dpL ~ 4e-5, negligible vs the
                    // O(1) Biot term in K_up). It costs an elastic-tangent
                    // reconstruction per ip and the upstream authors disabled the
                    // analogous classical term (line ~3315). Flip to true to
                    // enable the consistent-tangent contribution.
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
                        // 2026-06-06 PRESSURE form) ─────────────────────────────
                        // New residual sigma_sw = -(1 - phi_M)*(Pi - b*p_conf), so
                        //   d sigma_sw/d n_l   = -(1 - phi_M)*dPi/dn_l,
                        //   d sigma_sw/d eps_v = -(1 - phi_M)*b*(dp_conf/deps_v)
                        //                      = +(1 - phi_M)*b*K_drained
                        // (dp_conf/deps_v = -K_drained, the MECHANICAL drained bulk
                        // modulus). K_sw is GONE here -- it no longer routes a
                        // pressure through an elastic strain. Pi is the BARE vdW
                        // disjoining pressure (recomputed below, NOT
                        // micro_potential.mu_lR which already carries the film
                        // delta). Both pieces map to R_u through
                        //   dsigma'/d(.) = C * C_el^{-1} * d sigma_sw/d(.).
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
                            potential_exchange_params_ptr
                                ->potential_augmentation_prefactor,
                            potential_exchange_params_ptr
                                ->potential_augmentation_exponent);
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
                        double const dPi_dnl_swj =
                            -rho_d_film_swj * vdw_bare_swj.dmu_lR_dnl;

                        // u-p (via n_l): d sigma_sw/d n_l * dn_l/dpL on identity2.
                        double const dsigma_sw_dnl_scalar_swj =
                            -n_S_swj * dPi_dnl_swj;
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

                        // u-eps (K[u,u]): d sigma_sw/d eps_v = +(1-phi_M)*b*K_drained.
                        // eps_v = identity2^T B u, so the block is
                        //   B^T C C_el^{-1} (dsigma_sw/deps_v * identity2)
                        //        identity2^T B w.
                        double const K_drained_swj =
                            drainedBulkModulusFromStiffness<DisplacementDim>(
                                C_el_swj);
                        double const dsigma_sw_deps_v_scalar_swj =
                            n_S_swj * b_swj * K_drained_swj;
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
            local_rhs.template segment<pressure_size>(pressure_index)
                .noalias() += N_p.transpose() * potential_exchange_result.exchange.rho_L_hat *
                              w;

            // Direct macro Jacobian term for the exchange source. In analytic
            // mode this includes the implicit n_l(p_L) chain contribution.
            local_Jac
                .template block<pressure_size, pressure_size>(pressure_index,
                                                              pressure_index)
                .noalias() -= N_p.transpose() *
                              potential_exchange_result.drho_L_hat_dpL_direct * N_p * w;

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
    local_Jac
        .template block<pressure_size, displacement_size>(pressure_index,
                                                          displacement_index)
        .noalias() = Kpu / dt;

    // pressure equation
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p_L +
        (storage_p_a_p + storage_p_a_S) * (p_L - p_L_prev) / dt +
        Kpu * (u - u_prev) / dt;

    // displacement equation
    local_rhs.template segment<displacement_size>(displacement_index)
        .noalias() += Kup * p_L;
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
        updateMicroscaleHydraulicState<DisplacementDim>(
            this->current_states_[ip], this->prev_states_[ip], p_cap_ip,
            rho_LR, mu, dt, variables, variables_prev,
            {.phi = phi,
             .phi_M_prev = transport_porosity_prev_value,
             .phi_m_prev = phi_m_prev_value,
             .volumetric_strain = variables.volumetric_strain,
             .volumetric_strain_prev = variables_prev.volumetric_strain,
             .confining_pressure_p_conf = p_conf_micro_solve,
             .biot_coefficient = alpha},
            this->process_data_.micro_porosity_parameters,
            this->getPotentialExchangeParameters());
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
