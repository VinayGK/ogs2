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

    // ── Log-linear interpolant (production for the LIVE K(rho_d) path per
    // Vinay's interpolation-scheme decision, 2026-08-26) ──────────────────
    // ln(K) linear in rho_d between knots instead of the K-linear chord
    // used by getValue() above (Dixon 2023's own exponential
    // swelling-pressure-vs-density law is the physical motivation;
    // getValue()/getSegmentSlope() remain in use by the parse-time
    // frozen-K path). Node-preserving at the knots by construction (t=0 ->
    // K_l, t=1 -> K_r); flat-clamped outside [x_min, x_max] exactly like
    // getValue() -- never slope- or exp-extends past the table edges (the
    // extrapolation scheme is a separate question, out of scope here).
    // Mirrors getValue's <=/>= clamp branches and lower_bound interval
    // selection exactly, so this and getValue agree on which segment owns
    // an exact interior-knot argument.
    double getValueLogLinear(double const x) const
    {
        if (x <= supp_pnts_.front())
        {
            return values_at_supp_pnts_.front();
        }
        if (supp_pnts_.back() <= x)
        {
            return values_at_supp_pnts_.back();
        }
        auto const it =
            std::lower_bound(supp_pnts_.begin(), supp_pnts_.end(), x);
        std::size_t const i = std::distance(supp_pnts_.begin(), it) - 1;
        double const x_l = supp_pnts_[i];
        double const x_r = supp_pnts_[i + 1];
        double const K_l = values_at_supp_pnts_[i];
        double const K_r = values_at_supp_pnts_[i + 1];
        double const t = (x - x_l) / (x_r - x_l);
        // K(x) = K_l * exp(t * ln(K_r/K_l)) = K_l * (K_r/K_l)^t; exact at
        // t=0 -> K_l and t=1 -> K_r (node-preserving), matching getValue()
        // at every knot. [J/kg]
        return K_l * std::exp(t * std::log(K_r / K_l));
    }

    // d(getValueLogLinear)/dx, the EXACT companion tangent (chain rule of
    // d/dx[K_l * exp(t*ln(K_r/K_l))], t=(x-x_l)/(x_r-x_l)):
    //   dK/dx = K(x) * ln(K_r/K_l) / (x_r - x_l)   [(J/kg)/(kg/m^3)]
    // i.e. PROPORTIONAL TO THE LOCAL K VALUE, not a per-segment constant
    // as getSegmentSlope() is (that mismatch would give an inconsistent
    // residual/Jacobian tangent under this scheme). Same clamp/one-sided
    // convention as getSegmentSlope: 0 at/outside the boundary knots,
    // LEFT-segment value at an interior knot (idx = lower_bound - 1),
    // matching getValueLogLinear's and getValue's interval selection.
    double getSegmentSlopeLogLinear(double const x) const
    {
        if (x <= supp_pnts_.front() || supp_pnts_.back() <= x)
        {
            return 0.0;
        }
        auto const it =
            std::lower_bound(supp_pnts_.begin(), supp_pnts_.end(), x);
        std::size_t const i = std::distance(supp_pnts_.begin(), it) - 1;
        double const x_l = supp_pnts_[i];
        double const x_r = supp_pnts_[i + 1];
        double const K_l = values_at_supp_pnts_[i];
        double const K_r = values_at_supp_pnts_[i + 1];
        double const K_x = getValueLogLinear(x);  // J/kg
        return K_x * std::log(K_r / K_l) / (x_r - x_l);  // (J/kg)/(kg/m^3)
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

// ── KKT micro-water ceiling (branch dsm_mass_conservation_v3_kkt_ceiling_2026-09-30) ──
// Treatment of the micro water content ceiling n_l <= n_max(eps_v) = phi in the
// scalar_micro_macro_mass_storage_mode local solve.
//   Clamp (default): the shipped path, bitwise (projected Newton, the exchange
//     that the micro cannot take is deleted from the macro).
//   Kkt: complementarity problem, multiplier lambda >= 0, active set
//     f(n_max) < 0 and no interior root (DERIVATION.md 2.3-2.5, 4.4 of the
//     record folder ~/ogs-models/scratch/2026-09-30_kkt_ceiling_impl/, D-n).
//     NOT adopted; Vinay's ruling is open.
enum class MicroCeilingTreatment
{
    Clamp,
    Kkt
};

// Q9 (DESIGN.md D1): does the consistent exchange p-u entry reach Newton?
//   Overwritten (default): the assembly line local_Jac.pu = Kpu/dt stays as is
//     and erases every exchange p-u entry (shipped behaviour, bitwise).
//   KktActive: the KKT-active entries (D-5.2) are accumulated in a separate
//     matrix and added after that line.
//   AllExchange: the line becomes Kpu/dt plus the accumulated exchange entries
//     (the Maxwell, film and live-K p-u entries at inactive points as well).
enum class MicroCeilingPuTangent
{
    Overwritten,
    KktActive,
    AllExchange
};

// Model IV tangent term (DESIGN_FIXES.md part A, 2026-10-01): does the swelling
// eigenstress delta_sigma_sw reach Newton through the strain derivative that it
// has on the KKT-active branch (n_l = n_max(eps_v) = phi(eps_v))?
//   Overwritten (default): the K_uu swelling block stays as shipped (explicit
//     live-K chain only, read with variables_prev.porosity) and is bitwise the
//     tree without the switch.
//   KktActive: at the KKT-active integration points K_uu gets the total
//     derivative d(delta_sigma_sw)/d eps_v = [dK-chain + (partial in n_l) +
//     (partial in rho_lR) d rho_lR/d n_l] dphi/d eps_v, with the previous
//     porosity of the porosity law in the dphi/d eps_v of the live-K chain.
//     The residual is NOT changed (tangent only).
enum class MicroCeilingSwTangent
{
    Overwritten,
    KktActive
};

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

inline constexpr char const* toString(MicroCeilingTreatment const t)
{
    switch (t)
    {
        case MicroCeilingTreatment::Clamp:
            return "clamp";
        case MicroCeilingTreatment::Kkt:
            return "kkt";
    }
    return "unknown";
}

inline constexpr char const* toString(MicroCeilingSwTangent const t)
{
    switch (t)
    {
        case MicroCeilingSwTangent::Overwritten:
            return "overwritten";
        case MicroCeilingSwTangent::KktActive:
            return "kkt_active";
    }
    return "unknown";
}

inline constexpr char const* toString(MicroCeilingPuTangent const t)
{
    switch (t)
    {
        case MicroCeilingPuTangent::Overwritten:
            return "overwritten";
        case MicroCeilingPuTangent::KktActive:
            return "kkt_active";
        case MicroCeilingPuTangent::AllExchange:
            return "all_exchange";
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
    // 0.0 -> no floor -> evaluation is byte-identical to before.
    // Value source: PRJ-supplied (Vinay's call), not defaulted in code.
    // MANDATORY in the PRJ (Vinay 2026-06-17): the top-level <potential_exchange>
    // MUST declare <micro_water_content_floor>; the parser no longer defaults it
    // (see parsePotentialExchangeParameters). The 0.0 here is only the in-struct
    // fallback for medium-override inheritance, never a parse default.
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

    // -- DIAGNOSTIC switches, mass-strip A/B test (2026-09-30) ---------------
    // Implemented as specified by Vinay (R-03, 2026-09-30: "disclose and run
    // the test in parallel"; test text: mass_audit README l.52-57). NOTHING is
    // adopted: both default false -> the assembly is bit-identical to bed3e395
    // (variant A). See DSM/AGENTS.md worklog entry 2026-09-30.
    //
    // B: at IPs where the micro water content sits on its ceiling n_l = phi
    // (boundedMicroWaterContentCeiling), the macro pressure residual books the
    // ACTUAL micro storage rate (rho_l - rho_l_prev)/dt, rho_l = phi_m*rho_lR,
    // instead of the potential-driven rho_hat = alpha_M (mu_LR - mu_lR), in both
    // directions (drain and return), with the matching Jacobian. The micro
    // state update itself is untouched.
    bool ceiling_micro_storage_exchange = false;
    // B': the macro storage term uses phi_M (macro pore space) in place of the
    // total porosity phi in the two pore-fluid storage coefficients (a_p, a_S).
    bool macro_storage_uses_macro_porosity = false;

    // -- V2 switches, mass-fix trees (2026-09-30; Vinay: "do both ... in two
    // different trees"; Q1 = book the volume-change term or not). Both default
    // false -> bit-identical to the tree without them (V1 / diag B').
    // Derivation: DERIVATION.md in the record folder of the V2 tree
    // (~/ogs-models/scratch/2026-09-30_massfix_V2/DERIVATION.md). FORMULATION
    // CHANGE, flagged for Vinay's ruling; nothing adopted.
    //
    // F3 sign: the micro mass residual of the scalar_micro_macro_mass_storage_mode
    // carries  - dt*rho_l*eps_dot  (as shipped at bed3e395); the Eulerian micro
    // balance per current bulk volume, d(rho_l)/dt + rho_l*eps_dot = rho_hat,
    // gives  + dt*rho_l*eps_dot. true -> use the Eulerian sign in the residual
    // and the matching tangents of that mode (n_l-normalised scalar modes are
    // left unchanged).
    bool micro_mass_strain_term_eulerian = false;
    // Q1: at the IPs booked by ceiling_micro_storage_exchange the booked rate
    // includes the volume-change term: rho_hat_booked = (rho_l - rho_l_prev)/dt
    // + rho_l*eps_dot (the same Eulerian balance as above) instead of the bare
    // storage rate. Requires ceiling_micro_storage_exchange.
    bool ceiling_micro_storage_includes_strain = false;

    // -- KKT micro-water ceiling (branch dsm_mass_conservation_v3_kkt_ceiling_
    // 2026-09-30; Vinay 2026-09-30 "derivation, report, beamer, design docs,
    // implementation, weak forms, unit tests and then the ms33 suite"). ALL
    // default to the shipped behaviour; micro_ceiling_treatment = clamp is
    // bitwise the tree without them. NOT adopted. Design: DESIGN.md of the record
    // folder ~/ogs-models/scratch/2026-09-30_kkt_ceiling_impl/ (D4: enum tag).
    MicroCeilingTreatment micro_ceiling_treatment =
        MicroCeilingTreatment::Clamp;
    // Q9 (DESIGN.md D1); only legal with micro_ceiling_treatment = kkt.
    MicroCeilingPuTangent micro_ceiling_pu_tangent =
        MicroCeilingPuTangent::Overwritten;
    // Route-B debug flag (DESIGN.md 3.7): in-assembler central-difference check
    // of the assembled element Jacobian, log lines prefixed KKT-FD.
    bool micro_ceiling_fd_check = false;
    // N_dec of the bracketed scan (DERIVATION.md 2.5; DESIGN.md D7, a PROPOSAL
    // of 8 that needs Vinay's approval under the repo rule 1.2). >= 2.
    int micro_ceiling_scan_nodes_per_decade = 8;
    // Element ids whose integration points write the iteration trace
    // (DESIGN.md 3.8). Empty = off (no cost).
    std::vector<std::size_t> micro_ceiling_trace_elements;
    // Model IV tangent term (DESIGN_FIXES.md part A); only legal with
    // micro_ceiling_treatment = kkt. Default = shipped behaviour, bitwise.
    MicroCeilingSwTangent micro_ceiling_sw_tangent =
        MicroCeilingSwTangent::Overwritten;
};

// True when the KKT treatment of the micro-water ceiling is selected.
inline bool isKktCeiling(PotentialExchangeParameters const& p)
{
    return p.micro_ceiling_treatment == MicroCeilingTreatment::Kkt;
}

// Sign s of the micro mass residual's volume-change term in the mass-storage
// mode, written as  residual -= s*dt*rho_l*eps_dot  and in the tangents as
// (1 - s*dt*eps_dot):  s = +1 is the shipped form (F3),  s = -1 the Eulerian
// balance (micro_mass_strain_term_eulerian). Dimensionless.
inline double microMassStrainTermSign(PotentialExchangeParameters const& p)
{
    return p.micro_mass_strain_term_eulerian ? -1.0 : 1.0;  // [-]
}

// Effective augmentation prefactor K [J/kg] at the current state.
//
// LOG-LINEAR SCHEME (production per Vinay's K(rho_d) interpolation-scheme
// decision, 2026-08-26): evaluates the table via getValueLogLinear() --
// ln(K) linear in rho_d between knots -- instead of the K-linear
// getValue() (which remains in use by the parse-time frozen-K path in
// CreateRichardsMechanicsProcess.cpp). Node values are unchanged (both
// schemes are node-preserving); only the INTERIOR chord shape differs.
// Outside [rho_d_min, rho_d_max] the flat endpoint-hold clamp is
// preserved exactly as in getValue() (verified: getValueLogLinear mirrors
// getValue's <=/>= branches) -- this change does NOT touch the
// extrapolation scheme, only the interior interpolant.
// Live mode + table + finite phi -> K(rho_d) with rho_d = rho_SR*(1-phi)
// [kg/m^3] (rho_SR = micro_solid_density_reference; phi = current TOTAL
// porosity). Any other case (mode off, no table, phi sentinel/NaN) -> the
// parse-time scalar, bit-for-bit (unchanged).
inline double effectiveAugmentationPrefactor(
    PotentialExchangeParameters const& params, double const phi)
{
    if (params.potential_augmentation_prefactor_live_dry_density &&
        params.potential_augmentation_prefactor_vs_dry_density &&
        std::isfinite(phi))
    {
        // rho_d = rho_SR * (1 - phi)  [kg/m^3]
        return params.potential_augmentation_prefactor_vs_dry_density
            ->getValueLogLinear(params.micro_solid_density_reference *
                                (1.0 - phi));  // K [J/kg]
    }
    return params.potential_augmentation_prefactor;  // K [J/kg]
}

// d K_eff/d phi of effectiveAugmentationPrefactor above, at the same state.
//
// LOG-LINEAR SCHEME (companion to the log-linear value above, production
// per Vinay's decision 2026-08-26 -- REQUIRED so the residual and its
// Jacobian stay tangent-consistent; see getSegmentSlopeLogLinear
// doc for the chain-rule derivation). Chain (analytic derivation, this
// file): rho_d = rho_SR*(1-phi) [kg/m^3], so
//   dK/dphi = (dK/drho_d) * (drho_d/dphi)
//           = [K(rho_d) * ln(K_r/K_l)/(x_r-x_l)] * (-rho_SR)
// i.e. proportional to the LOCAL K value (getSegmentSlopeLogLinear), not the
// old per-segment constant (getSegmentSlope). Returns 0 in EVERY case where
// effectiveAugmentationPrefactor returns the parse-time scalar (mode off, no
// table, phi sentinel/NaN) and at/outside the clamped table edges (where the
// clamped value is flat in rho_d) -- exactly the one-sided/zero-slope
// convention documented on getSegmentSlopeLogLinear, unchanged from the
// standing getSegmentSlope convention. The RESIDUAL is untouched by this
// helper; it feeds the Jacobian only.
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
                   ->getSegmentSlopeLogLinear(
                       rho_SR * (1.0 - phi));  // dK/dphi [J/kg per unit phi]:
                                               // [kg/m^3]*[J/kg / (kg/m^3)]
    }
    return 0.0;  // J/kg per unit phi
}
}  // namespace ProcessLib::RichardsMechanics
