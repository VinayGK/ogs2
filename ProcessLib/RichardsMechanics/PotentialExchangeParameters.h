// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <memory>
#include <optional>
#include <vector>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace ProcessLib::RichardsMechanics
{
// Piecewise-linear K(rho_d) table with an EXACT per-segment slope accessor.
// MathLib::PiecewiseLinearInterpolation::getDerivative blends the two
// adjacent segment slopes (quadratic smoothing, see its .cpp), which is NOT
// the derivative of getValue's clamped piecewise-linear evaluation. The
// live-K(rho_d) Jacobian (K_OF_RHO_D_LIVE.md) needs the slope of the VALUE
// actually fed into the residual, so this thin subclass exposes the exact
// segment slope via the protected knot vectors.
class AugmentationPrefactorTable final
    : public MathLib::PiecewiseLinearInterpolation
{
public:
    using MathLib::PiecewiseLinearInterpolation::PiecewiseLinearInterpolation;

    // Exact d(getValue)/dx of the clamped piecewise-linear evaluation.
    // Convention (documented choice, mirrors getValue's branch structure):
    //  - x <= x_min or x >= x_max: 0 (getValue holds the endpoint value, so
    //    the clamped evaluation is FLAT there; this is the one-sided outward
    //    slope AT the edge knots as well).
    //  - interior knots: the LEFT segment slope (one-sided), consistent with
    //    getValue's lower_bound interval selection (idx = lower_bound - 1).
    double getSegmentSlope(double const x) const
    {
        if (x <= supp_pnts_.front() || supp_pnts_.back() <= x)
        {
            return 0.0;
        }
        auto const it =
            std::lower_bound(supp_pnts_.begin(), supp_pnts_.end(), x);
        std::size_t const i = std::distance(supp_pnts_.begin(), it) - 1;
        return (values_at_supp_pnts_[i + 1] - values_at_supp_pnts_[i]) /
               (supp_pnts_[i + 1] - supp_pnts_[i]);
    }
};

enum class MicroPotentialConvention
{
    PositiveReduced,
    NegativeAttractive
};

enum class LocalNonlinearSolveMode
{
    ScalarExchange,
    ScalarReferenceStorage,
    ScalarReferenceMassStorage
};

enum class MacroPorosityUpdateMode
{
    AlgebraicSplit,
    ReferenceAdditiveRate
};

enum class MicroSolidVolumeFractionMode
{
    Reference,
    CurrentPorositySplit
};

// ── Strained-film disjoining law h(w_m, eps_v) (DSM/STRAINED_FILM_IMPLEMENTATION.md) ──
// Off:         film geometry frozen (current behavior, bit-for-bit).
// Kinematic:   variant A — spacing follows the volumetric strain,
//              h = h0(n_l)*(1 + kappa*eps_v)  <=>  evaluate the bare law at
//              w_eff = n_l*(1 + kappa*eps_v).
// Equilibrium: variant B — spacing tracks the film force balance once the load
//              can compress the film: w_eff solves Pi(w_eff) = p_conf on the
//              loaded branch (p_conf > Pi(n_l)), else w_eff = n_l (emergent
//              branch point; no bolted-on gate).
enum class FilmStrainCouplingMode
{
    Off,
    Kinematic,
    Equilibrium
};

// Spacing-strain weighting kappa in dh/deps_v = kappa*h0 (design doc §3, D1):
// Aggregate: kappa = (1 - phi_M) (active_nS at the GP) — the integrable
//            completion of the existing eigenstress scale (recommended).
// Unity:     kappa = 1 — naive geometric reading (spacing follows REV strain
//            one-to-one); kept PRJ-selectable for discrimination (Vinay,
//            2026-06-09).
enum class FilmStrainKappaMode
{
    Aggregate,
    Unity
};

// ── Film energy route (DSM/PI_OF_NL_EV_IMPLEMENTATION.md, Vinay 2026-06-11) ──
// Operational: the shipped Derjaguin cut — bare law evaluated at w_eff plus the
//              hand-added load term +b*p_conf/rho_lR (NOT Maxwell-exact; defect
//              O(Pi*eps_v), strained-film design doc §9a). Default, bit-for-bit.
// Exact:       the one-Psi energy route — Psi_film(n_l, eps_v) with closed-form
//              strain integrals of the disjoining law along the kinematic
//              h-law; mu_mech = (1/(nS*rho_lR)) dPsi/dn_l. Maxwell holds
//              identically; kappa->0 reduces EXACTLY to the shipped integrable
//              partner. Requires film_strain_coupling == Kinematic (the closed
//              forms are for the kinematic h-law).
enum class FilmEnergyRoute
{
    Operational,
    Exact
};

// Create-time admissibility of the (film_strain_coupling, film_energy_route)
// combination (PI_OF_NL_EV_IMPLEMENTATION.md §3 mode matrix). Pure predicate so
// it is unit-testable; the OGS_FATAL lives at the parse site.
inline constexpr bool isValidFilmEnergyRouteCombination(
    FilmStrainCouplingMode const mode, FilmEnergyRoute const route)
{
    return route == FilmEnergyRoute::Operational ||
           mode == FilmStrainCouplingMode::Kinematic;
}

inline constexpr char const* toString(
    MicroPotentialConvention const convention)
{
    switch (convention)
    {
        case MicroPotentialConvention::PositiveReduced:
            return "positive_reduced";
        case MicroPotentialConvention::NegativeAttractive:
            return "negative_attractive";
    }
    return "unknown";
}

inline constexpr double microPotentialSignFactor(
    MicroPotentialConvention const convention)
{
    return convention == MicroPotentialConvention::NegativeAttractive ? -1.0
                                                                       : 1.0;
}

inline constexpr char const* toString(LocalNonlinearSolveMode const mode)
{
    switch (mode)
    {
        case LocalNonlinearSolveMode::ScalarExchange:
            return "scalar_exchange";
        case LocalNonlinearSolveMode::ScalarReferenceStorage:
            return "scalar_microstate_storage_mode";
        case LocalNonlinearSolveMode::ScalarReferenceMassStorage:
            return "scalar_micro_macro_mass_storage_mode";
    }
    return "unknown";
}

inline constexpr char const* toString(MacroPorosityUpdateMode const mode)
{
    switch (mode)
    {
        case MacroPorosityUpdateMode::AlgebraicSplit:
            return "algebraic_split";
        case MacroPorosityUpdateMode::ReferenceAdditiveRate:
            return "additive_macro_porosity_rate_mode";
    }
    return "unknown";
}

inline constexpr char const* toString(
    MicroSolidVolumeFractionMode const mode)
{
    switch (mode)
    {
        case MicroSolidVolumeFractionMode::Reference:
            return "reference";
        case MicroSolidVolumeFractionMode::CurrentPorositySplit:
            return "current_porosity_split";
    }
    return "unknown";
}

inline constexpr char const* toString(FilmStrainCouplingMode const mode)
{
    switch (mode)
    {
        case FilmStrainCouplingMode::Off:
            return "off";
        case FilmStrainCouplingMode::Kinematic:
            return "kinematic";
        case FilmStrainCouplingMode::Equilibrium:
            return "equilibrium";
    }
    return "unknown";
}

inline constexpr char const* toString(FilmStrainKappaMode const mode)
{
    switch (mode)
    {
        case FilmStrainKappaMode::Aggregate:
            return "aggregate";
        case FilmStrainKappaMode::Unity:
            return "unity";
    }
    return "unknown";
}

inline constexpr char const* toString(FilmEnergyRoute const route)
{
    switch (route)
    {
        case FilmEnergyRoute::Operational:
            return "operational";
        case FilmEnergyRoute::Exact:
            return "exact";
    }
    return "unknown";
}

struct PotentialExchangeParameters
{
    bool enabled = false;

    // Young-Laplace macro potential branch tolerance.
    double pressure_tolerance = 0.0;

    // vdW microscale potential parameters / reference state constants.
    double hamaker_constant = 0.0;
    double specific_surface = 0.0;
    double micro_solid_density_reference = 0.0;          // rho_SR
    double micro_solid_volume_fraction_reference = 0.0;  // n_S
    double micro_liquid_density_reference = 0.0;         // rho_l0
    double micro_liquid_density_a = 0.0;                 // a_rho
    double micro_liquid_density_b = 0.0;                 // b_rho
    MicroPotentialConvention micro_potential_convention =
        MicroPotentialConvention::PositiveReduced;
    LocalNonlinearSolveMode local_nonlinear_solve_mode =
        LocalNonlinearSolveMode::ScalarExchange;
    MacroPorosityUpdateMode macro_porosity_update_mode =
        MacroPorosityUpdateMode::AlgebraicSplit;
    MicroSolidVolumeFractionMode micro_solid_volume_fraction_mode =
        MicroSolidVolumeFractionMode::Reference;

    // Optional GP-local n_l initialization (future full 2C path).
    std::optional<double> initial_micro_water_content;

    // Optional Jacobian approximation for DSM exchange contribution only.
    // If true, drho_L_hat/dp_L is computed by finite difference in the local
    // helper path.
    bool use_fd_jacobian_for_exchange = false;
    double fd_jacobian_perturbation = 1e-8;

    // Finite-difference step for the implicit n_l(p_L) chain-rule derivative
    // used in ScalarReferenceMassStorage mode.
    double local_jacobian_perturbation = 1e-8;

    // Lumped exponential force augmentation to the vdW micro-potential.
    // h = n_l / (nS * rho_SR * Sa)   [mean water film thickness, m]
    // mu_lR_aug = sign * K * exp(-h / lambda)
    // Zero prefactor (default) disables augmentation and preserves
    // existing behaviour.
    double potential_augmentation_prefactor = 0.0;     // K      [J/kg], must be >= 0
    double potential_augmentation_exponent = 0.0;  // lambda [m],    must be > 0 if K > 0

    // ── Disjoining-pressure FLOOR via a micro-water-content lower bound ───────
    // Optional lower bound n_l,min [-] on the water content USED IN THE vdW
    // DISJOINING LAW ONLY (Pi ~ 1/n_l^3). When > 0, the law is evaluated at
    // max(n_l, micro_water_content_floor), so Pi is CAPPED at Pi(floor) instead
    // of diverging as n_l -> 0. This is local to the disjoining evaluation: it
    // does NOT change the global n_l, the exchange, or the porosity. Below the
    // floor the clamped Pi is FLAT in n_l, so its n_l-derivatives are 0 there.
    // 0.0 (default) -> no floor -> evaluation is byte-identical to before.
    // Value source: PRJ-supplied (Vinay's call), not defaulted in code.
    double micro_water_content_floor = 0.0;  // n_l,min [-], must be >= 0

    // Optional consistency switches for the hierarchical DSM branch.
    // Default micro-pressure density is the confined micro-liquid density.
    bool use_micro_liquid_density_for_micro_pressure = true;

    // ── Film-pressure coupling (maxwell beamer sec.5) ──────────────────────
    // Default ON (2026-06-08, Vinay): the model is CONSOLIDATED on the film
    // coupling. mu_lR carries the effective-stress (film) term mu_lR(p_film =
    // p_disj + sigma') in ALL local solves and the macro exchange, the swelling
    // stress is the eigenstrain form (S1 < 0 -> compression drains), biot=alpha
    // (incompressible grains), and the sharp gate is a C1 activation of width
    // film_pressure_gate_width. The bare-Pi OFF formulation is RETIRED: it is
    // forced true at parse (CreateRichardsMechanicsProcess), so OFF is unrunnable;
    // the residual OFF code branches are dead and pending physical removal.
    bool film_pressure_coupling = true;
    // NOTE: the eigenstrain Biot b is NO LONGER a separate film parameter. It is
    // unified with the poroelastic biot_coefficient MPL medium property (same
    // solid-fluid volume partitioning; one-Psi consistency) and threaded into the
    // local solve via PotentialExchangeLocalSolveContext::biot_coefficient.
    double film_pressure_gate_width = 0.0;        // smooth-gate width w [Pa]; 0 -> sharp fallback  [Vinay's call]
    // DEPRECATED 2026-06-06: swelling stress is now (1-phi_M)*p_film; this modulus is unused.
    double film_pressure_swelling_modulus = 0.0;  // eigenstrain modulus K_sw [Pa]; 0 -> drained K  [Vinay's call]

    // ── Macro-porosity floor (Vinay 2026-06-06) ────────────────────────────
    // phi_M,min (REV macro porosity). Prevents the macro pore from collapsing
    // into the interlayer: the interlayer water n_l is capped at
    // n_l_cap = (phi - macro_porosity_floor)/(1 - macro_porosity_floor), so the
    // hierarchical split phi_M = (phi - n_l)/(1 - n_l) >= macro_porosity_floor.
    // Beyond the cap the film is saturated and further water stays bulk (macro):
    // porosity- and water-conserving (phi = phi_M + phi_m held; the capped micro
    // uptake remains in the macro mass balance). Value source: EPFL MIP bimodal
    // pore structure (Seiphoori 2014 / Acta 2022) [Vinay's call]. 0 (default) ->
    // no floor -> bit-for-bit unchanged.
    double macro_porosity_floor = 0.0;
    double macro_floor_cutoff_width = 0.0;  // film-to-bulk cutoff width in n_l [-]; 0 -> default 5% of n_l_cap [Vinay's call]

    // ── Strained-film disjoining law (DSM/STRAINED_FILM_IMPLEMENTATION.md) ──
    // When != Off, the bare disjoining law is evaluated at the strained film
    // state w_eff and mu_lR gains the load term +b*p_conf/rho_lR; the shipped
    // integrable mechanical partner is REPLACED (it is the frozen-h, O(eps_v)
    // truncation of the same physics — running both double-counts; D3
    // provisional, demonstrated by the shipped-limit unit test). Off (default)
    // is bit-for-bit the current behavior.
    FilmStrainCouplingMode film_strain_coupling = FilmStrainCouplingMode::Off;
    FilmStrainKappaMode film_strain_kappa = FilmStrainKappaMode::Aggregate;

    // ── Film energy route (DSM/PI_OF_NL_EV_IMPLEMENTATION.md) ───────────────
    // Operational (default): shipped Derjaguin cut, bit-for-bit. Exact: the
    // one-Psi pair — REPLACES the operational mu assembly when ON (kinematic
    // only; create-time validated). The eigenstress half is identical in both
    // routes (Pi at w_eff with the actual p_conf), so only the fold-point mu
    // assembly differs.
    FilmEnergyRoute film_energy_route = FilmEnergyRoute::Operational;

    // ── K(rho_d): augmentation prefactor as a function of dry density ──────
    // Optional piecewise-linear table K = K(rho_d) [J/kg vs kg/m^3]. When set
    // together with `dry_density`, the augmentation prefactor above is
    // RESOLVED at parse time to K(dry_density) and stored back into
    // `potential_augmentation_prefactor` — i.e. K is the *initial/target*
    // dry-density value, a per-material constant in time (Vinay 2026-06-08).
    // Because resolution is parse-time and time-constant, the downstream
    // potential/exchange tangent is unchanged (no dK/drho_d term). The table
    // and dry density are carried here only so a per-<medium id> override can
    // inherit the shared table from the global block as its default.
    // getValue() clamps outside [rho_d_min, rho_d_max] (endpoint hold).
    std::shared_ptr<AugmentationPrefactorTable const>
        potential_augmentation_prefactor_vs_dry_density = nullptr;
    std::optional<double> dry_density;  // rho_d [kg/m^3], initial/target

    // ── LIVE K(rho_d) (K_OF_RHO_D_LIVE.md; Vinay 2026-06-10 "K(rho_d) try
    // it") ──. When true, the table above is NOT frozen at parse time;
    // instead K is re-evaluated at the EVOLVING dry density rho_d =
    // rho_SR*(1-phi) at every evaluation site that has the current total
    // porosity phi in scope (see effectiveAugmentationPrefactor below).
    // Sites without phi fall back to the scalar `potential_augmentation_
    // prefactor`. The analytic dK/dphi = -rho_SR*(table segment slope)
    // tangent is wired into the live p-u augmentation Jacobian block since
    // 2026-06-12 (Vinay's approved completion; see
    // effectiveAugmentationPrefactorPhiDerivative below and
    // K_OF_RHO_D_LIVE.md) — the first cut's omission note is historical.
    // false (default) -> parse-time freeze, bit-for-bit the existing
    // behavior.
    bool potential_augmentation_prefactor_live_dry_density = false;
};

// Effective augmentation prefactor K [J/kg] at the current state.
// Live mode + table + finite phi -> K(rho_d) with rho_d = rho_SR*(1-phi)
// [kg/m^3] (rho_SR = micro_solid_density_reference; phi = current TOTAL
// porosity). PiecewiseLinearInterpolation::getValue holds the endpoint
// values outside [rho_d_min, rho_d_max] (verified: MathLib/
// InterpolationAlgorithms/PiecewiseLinearInterpolation.cpp, getValue),
// so K is clamped at the table range ends. Any other case (mode off, no
// table, phi sentinel/NaN) -> the parse-time scalar, bit-for-bit.
inline double effectiveAugmentationPrefactor(
    PotentialExchangeParameters const& params, double const phi)
{
    if (params.potential_augmentation_prefactor_live_dry_density &&
        params.potential_augmentation_prefactor_vs_dry_density &&
        std::isfinite(phi))
    {
        // rho_d = rho_SR * (1 - phi)  [kg/m^3]
        return params.potential_augmentation_prefactor_vs_dry_density
            ->getValue(params.micro_solid_density_reference *
                       (1.0 - phi));  // K [J/kg]
    }
    return params.potential_augmentation_prefactor;  // K [J/kg]
}

// d K_eff/d phi of effectiveAugmentationPrefactor above, at the same state.
// Chain (analytic derivation, this file; Vinay 2026-06-12 approved Jacobian
// completion of live K(rho_d)): rho_d = rho_SR*(1-phi) [kg/m^3], so
//   dK/dphi = (dK/drho_d) * (drho_d/dphi) = (table segment slope) * (-rho_SR).
// Returns 0 in EVERY case where effectiveAugmentationPrefactor returns the
// parse-time scalar (mode off, no table, phi sentinel/NaN) and at/outside the
// clamped table edges (where the clamped value is flat in rho_d) — exactly
// the one-sided/zero-slope convention documented on getSegmentSlope. The
// RESIDUAL is untouched by this helper; it feeds the Jacobian only.
inline double effectiveAugmentationPrefactorPhiDerivative(
    PotentialExchangeParameters const& params, double const phi)
{
    if (params.potential_augmentation_prefactor_live_dry_density &&
        params.potential_augmentation_prefactor_vs_dry_density &&
        std::isfinite(phi))
    {
        double const rho_SR = params.micro_solid_density_reference;  // kg/m^3
        return -rho_SR *
               params.potential_augmentation_prefactor_vs_dry_density
                   ->getSegmentSlope(
                       rho_SR * (1.0 - phi));  // dK/dphi [J/kg per unit phi]:
                                               // [kg/m^3]*[J/kg / (kg/m^3)]
    }
    return 0.0;  // J/kg per unit phi
}
}  // namespace ProcessLib::RichardsMechanics
