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
// K(rho_d) table carrying TWO value/slope pairs, each slope the EXACT
// derivative of the value it belongs to:
//   getValue()          / getSegmentSlope()          -- K     linear in rho_d
//   getValueLogLinear() / getSegmentSlopeLogLinear() -- ln(K) linear in rho_d
// The LIVE K(rho_d) path (K_OF_RHO_D_LIVE.md) uses the LOG-LINEAR pair, per
// Vinay's interpolation-scheme decision of 2026-08-26 (commit 1bb414ac05):
// effectiveAugmentationPrefactor() calls getValueLogLinear(), and its Jacobian
// companion effectiveAugmentationPrefactorPhiDerivative() calls
// getSegmentSlopeLogLinear(), so residual and tangent share one interpolant.
// The K-linear pair is retained, not superseded: getValue() still resolves the
// parse-time frozen-K scalar (CreateRichardsMechanicsProcess.cpp), and
// getSegmentSlope() is kept as that pair's exact companion -- unit-test
// covered, with no production caller since the live Jacobian moved to the
// log-linear slope in 1bb414ac05.
// HISTORICAL: the class was introduced as the K-linear slope accessor alone,
// because MathLib::PiecewiseLinearInterpolation::getDerivative blends the two
// adjacent segment slopes (quadratic smoothing, see its .cpp), which is NOT
// the derivative of getValue's clamped piecewise-linear evaluation; the live
// Jacobian needs the slope of the VALUE actually fed into the residual, so
// this thin subclass exposes exact segment slopes via the protected knot
// vectors. That reason is unchanged -- it now covers both pairs.
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
    // Both branches come from locateSegment() below, which is the single
    // copy of that clamp + interval selection shared by every accessor in
    // this class (getValue() itself lives in the MathLib base and is
    // deliberately neither touched nor shadowed; locateSegment reproduces
    // its branches).
    double getSegmentSlope(double const x) const
    {
        auto const s = locateSegment(x);
        if (!s)
        {
            return 0.0;
        }
        // Unchanged chord slope (K_{i+1} - K_i) / (x_{i+1} - x_i).
        return (s->K_r - s->K_l) / (s->x_r - s->x_l);  // (J/kg)/(kg/m^3)
    }

    // ── Log-linear interpolant (production for the LIVE K(rho_d) path per
    // Vinay's interpolation-scheme decision, 2026-08-26) ──────────────────
    // ln(K) linear in rho_d between knots instead of the K-linear chord
    // used by getValue() above (Dixon 2023's own exponential
    // swelling-pressure-vs-density law is the physical motivation). The
    // K-linear pair is not retired, but its two halves are not equally
    // live: getValue() has exactly ONE production caller -- the parse-time
    // frozen-K resolution in CreateRichardsMechanicsProcess.cpp -- while
    // getSegmentSlope() has NONE. That is structural, not an oversight:
    // the frozen-K path resolves K to a scalar at parse time, so it
    // introduces no Jacobian term and can never want a slope.
    // getSegmentSlope() is exercised only by
    // Tests/ProcessLib/RichardsMechanics/StrainedFilmPotential.cpp (call
    // sites re-grepped over ProcessLib/ and Tests/ProcessLib/ on
    // 2026-08-31). Flat-clamped outside [x_min, x_max] exactly like
    // getValue() -- never slope- or exp-extends past the table edges (the
    // extrapolation scheme is a separate question, out of scope here).
    // Mirrors getValue's <=/>= clamp branches and lower_bound interval
    // selection exactly -- both go through locateSegment() below -- so this
    // and getValue agree on which segment owns an exact interior-knot
    // argument.
    //
    // NODE PRESERVATION, and why it needs a branch: a BOUNDARY knot returns
    // the stored endpoint through the clamp, but an INTERIOR knot lands in
    // the LEFT segment at t == 1, where K_l*exp(1*ln(K_r/K_l)) round-trips
    // through log/exp and can land ~1 ULP off K_r instead of on it
    // (measured 2026-08-31 on the SUPERSEDED 900-knot table
    // K(900)=4367.2277: the 1400 knot came out low by 1 ULP, rel -1.6e-16,
    // the miss recorded in DSM/AGENTS.md; the table shipped at 7ec39ecf4c
    // happens to hit all four of its knots bit-exactly, so this is a latent
    // defect, not a live one). logLinearValueOnSegment() below returns the
    // STORED knot value at t == 1 (and at t == 0), which is what makes the
    // node-preservation claim true of the code rather than nearly true.
    //
    // PRECONDITION: strictly positive table values; <prefactors> is
    // validated > 0 at parse time in CreateRichardsMechanicsProcess.cpp.
    // The non-positive cases do NOT share one failure mode: probed
    // 2026-08-31 on a standalone transcription of logLinearValueOnSegment
    // (IEEE-754 double), they are three distinct chains, and only one of
    // them starts at the logarithm:
    //   K_l == 0: K_r/K_l = +inf and std::log(+inf) = +inf -- still no
    //             NaN. The NaN is born one step later, in
    //             K_l*exp(t*r) = 0*inf. At an interior knot (t == 1) the
    //             stored-knot return even yields a clean K_r, while the
    //             companion slope goes to +inf.
    //   K_l <  0: the ratio is negative and std::log of it IS NaN, which
    //             then propagates through value and slope alike. This is
    //             the only case the logarithm itself catches.
    //   K_r == 0: the ratio is 0, std::log(0) = -inf, and the VALUE comes
    //             out a clean, plausible K_l*exp(-inf) = 0 that nothing
    //             downstream can distinguish from a legitimately small K,
    //             while the slope 0*(-inf)/(x_r-x_l) is NaN. The silent
    //             case is the dangerous one: the poison enters through
    //             the Jacobian, not the value.
    // Two knots of the SAME negative sign produce no NaN anywhere (the
    // ratio is positive, the logarithm finite) and would carry a negative
    // K through the whole chain untouched. The K-linear getValue() has
    // none of these failure modes, which is why the precondition belongs
    // to the log-linear pair.
    double getValueLogLinear(double const x) const
    {
        auto const s = locateSegment(x);
        if (!s)
        {
            // Endpoint hold: the same two branches, in the same order, as
            // getValue().
            return x <= supp_pnts_.front() ? values_at_supp_pnts_.front()
                                           : values_at_supp_pnts_.back();
        }
        double const t = (x - s->x_l) / (s->x_r - s->x_l);
        // K(x) = K_l * exp(t * ln(K_r/K_l)) = K_l * (K_r/K_l)^t. [J/kg]
        return logLinearValueOnSegment(*s, t,
                                       std::log(s->K_r / s->K_l));  // J/kg
    }

    // d(getValueLogLinear)/dx, the EXACT companion tangent (chain rule of
    // d/dx[K_l * exp(t*ln(K_r/K_l))], t=(x-x_l)/(x_r-x_l)):
    //   dK/dx = K(x) * ln(K_r/K_l) / (x_r - x_l)   [(J/kg)/(kg/m^3)]
    // i.e. PROPORTIONAL TO THE LOCAL K VALUE, not a per-segment constant
    // as getSegmentSlope() is (that mismatch would give an inconsistent
    // residual/Jacobian tangent under this scheme). Same clamp/one-sided
    // convention as getSegmentSlope: 0 at/outside the boundary knots,
    // LEFT-segment value at an interior knot (idx = lower_bound - 1) --
    // both go through locateSegment(), so that interval selection is one
    // object, not a copy kept in step by hand.
    // The segment is evaluated ONCE here (one interval lookup, one
    // std::log, one std::exp) instead of re-entering getValueLogLinear,
    // which would repeat the clamp, the binary search and the logarithm:
    // it sits on the per-integration-point path of assembleWithJacobian,
    // reached from BOTH live-K tangent sites. Those are two DIFFERENT
    // blocks of that integration-point loop, not one (line numbers at
    // 7ec39ecf4c, RichardsMechanicsFEM-impl.h): the p-u augmentation
    // exchange tangent at 5052, and the displacement-side
    // swelling-eigenstress tangent at 5223.
    double getSegmentSlopeLogLinear(double const x) const
    {
        auto const s = locateSegment(x);
        if (!s)
        {
            return 0.0;
        }
        double const t = (x - s->x_l) / (s->x_r - s->x_l);
        double const r = std::log(s->K_r / s->K_l);  // ln(K_r/K_l) [-]
        return logLinearValueOnSegment(*s, t, r) * r /
               (s->x_r - s->x_l);  // (J/kg)/(kg/m^3)
    }

private:
    // One located table segment [x_l, x_r] with its two knot values.
    struct Segment
    {
        double x_l;  // kg/m^3
        double x_r;  // kg/m^3
        double K_l;  // J/kg
        double K_r;  // J/kg
    };

    // THE single copy of the clamp + interval-selection logic that the three
    // accessors above share (getValue() itself lives in the MathLib base
    // class and is non-virtual, so it is deliberately neither touched nor
    // shadowed; this reproduces its branches). Returns nullopt exactly when
    // getValue() would return a clamped endpoint -- x <= x_min or
    // x_max <= x, the same <=/<= tests in the same order. Otherwise the
    // interval is the one getValue() uses, idx = lower_bound(x) - 1, so an
    // exact interior knot belongs to the LEFT segment and sits at t == 1
    // there.
    std::optional<Segment> locateSegment(double const x) const
    {
        if (x <= supp_pnts_.front() || supp_pnts_.back() <= x)
        {
            return std::nullopt;
        }
        auto const it =
            std::lower_bound(supp_pnts_.begin(), supp_pnts_.end(), x);
        std::size_t const i = std::distance(supp_pnts_.begin(), it) - 1;
        return Segment{supp_pnts_[i], supp_pnts_[i + 1],
                       values_at_supp_pnts_[i], values_at_supp_pnts_[i + 1]};
    }

    // Log-linear value on an ALREADY-LOCATED segment, at the segment
    // coordinate t = (x - x_l)/(x_r - x_l), with r = ln(K_r/K_l) for that
    // segment passed in so getSegmentSlopeLogLinear can reuse the one
    // logarithm it needs anyway. t == 1 and t == 0 return the STORED knot
    // values -- that is what makes node preservation bit-exact instead of
    // 1-ULP-approximate. (t == 0 is unreachable through locateSegment,
    // whose lower_bound places an exact knot at t == 1; it is kept so the
    // helper is total.)
    static double logLinearValueOnSegment(Segment const& s, double const t,
                                          double const r)
    {
        if (t == 0.0)
        {
            return s.K_l;  // J/kg
        }
        if (t == 1.0)
        {
            return s.K_r;  // J/kg
        }
        return s.K_l * std::exp(t * r);  // K_l*(K_r/K_l)^t  [J/kg]
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
//     Kkt carried forward with the F3 ruling 2026-10-01; on in MS33 cand. 2a/2b.
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

// Latched saturation gate for the Bishop factor chi and the relative
// permeability k_rel at KKT-active, previously saturated integration points
// (branch dsm_mass_conservation_v3_kkt_vii_gate_2026-10-01; design part B.4 of
// ~/ogs-models/scratch/2026-10-01_kkt_iv_vii_fixes/DESIGN_FIXES.md; ruled by
// Vinay 2026-10-01 ~15:15 CEST, "yes to both, keep full weight and k_rel = 1").
//   Off (default): the shipped rules, bitwise.
//   BishopRelperm: at a latched Active point chi = chi_deck(S = 1), dchi/dS = 0,
//     chi_prev = chi_deck(S = 1), k_rel = k_rel(S = 1), dk_rel/dS = 0.
//   Bishop: the same for chi only (a labelled PROBE, never a fix).
// S_L itself, the retention law, storage, exchange, the Biot term and the
// output saturation are not changed. Requires micro_ceiling_treatment = kkt.
enum class MicroCeilingSaturationGate
{
    Off,
    BishopRelperm,
    Bishop
};

inline constexpr char const* toString(MicroCeilingSaturationGate const gate)
{
    switch (gate)
    {
        case MicroCeilingSaturationGate::Off:
            return "off";
        case MicroCeilingSaturationGate::BishopRelperm:
            return "bishop_relperm";
        case MicroCeilingSaturationGate::Bishop:
            return "bishop";
    }
    return "unknown";
}

// ── v4 switches (branch dsm_mass_conservation_v4_tm_krel_2026-10-02; Vinay's
// ruling 2026-10-02 "(go with L + drop T_m) x (1a, 1b separate)"; design
// ~/ogs-models/scratch/2026-10-02_kkt_v4_tm_krel/DESIGN_V4.md). All defaults =
// the AB code (tip 35fbd4149b), bitwise.
//
// darcy_relative_permeability_mobility (variant 1a, DESIGN_V4.md 2.3):
//   GaussPoint (default): k_rel(S_L(p_c,ip)) at every integration point, the
//     shipped Galerkin mobility, bitwise.
//   KirchhoffElementMean: one mobility per element, the mean of the deck law
//     k_rel(S_L(p_c)) over [min, max] of the element's nodal p_c, evaluated
//     through the Kirchhoff potential of a piecewise-linear table of the deck
//     law (KirchhoffMobility.h). The range-mean over the nodal [min, max] is
//     the DESIGN's reading of the ruled "element-mean (Kirchhoff) mobility"
//     for 2D quads (other readings exist: edge-wise / two-point Kirchhoff, a
//     mean along the gradient; DESIGN_V4.md 2.3.1, Q6), chosen because it is
//     monotone in the boundary-layer mode and has an exact tangent.
enum class DarcyRelativePermeabilityMobility
{
    GaussPoint,
    KirchhoffElementMean
};

inline constexpr char const* toString(
    DarcyRelativePermeabilityMobility const mobility)
{
    switch (mobility)
    {
        case DarcyRelativePermeabilityMobility::GaussPoint:
            return "gauss_point";
        case DarcyRelativePermeabilityMobility::KirchhoffElementMean:
            return "kirchhoff_element_mean";
    }
    return "unknown";
}

// micro_ceiling_closed_macro_gate (variant 1b, DESIGN_V4.md 2.4; offered to
// Vinay as "Treat closed macro pores as gas-free (k_rel = 1)", ruled
// 2026-10-02): at a KKT-active integration point whose macro pores are closed
// (phi_M = 0) the macro pore space holds no gas, so k_rel = k_deck(S = 1).
//   Off (default): the AB rules, bitwise.
//   Relperm (the run level): k_rel = k_deck(S = 1), dk_rel/dS = 0 where the
//     predicate holds; chi, chi_prev, dchi/dS and p_FR are NOT touched.
//   BishopRelperm (BUILT, NOT RUN; Vinay's call): additionally chi =
//     chi_deck(S = 1), dchi/dS = 0, gated p_FR, and chi_prev = chi_deck(S = 1)
//     only if the gate acted at the end of the previous step
//     (MicroClosedMacroGateActed), else the deck/Fix-B value.
// The predicate (closedMacroGateActs below) mirrors the Fix B latch reading:
// the previous converged step's status Active AND its phi_M == 0 AND this
// iterate's status Active. The one-step lag (prev_active) is the design's
// reading of the ruling text "k_rel = 1 at KKT-active points with phi_M = 0",
// mirroring the Fix B latch; it is NOT in the ruling text.
enum class MicroCeilingClosedMacroGate
{
    Off,
    Relperm,
    BishopRelperm
};

inline constexpr char const* toString(MicroCeilingClosedMacroGate const gate)
{
    switch (gate)
    {
        case MicroCeilingClosedMacroGate::Off:
            return "off";
        case MicroCeilingClosedMacroGate::Relperm:
            return "relperm";
        case MicroCeilingClosedMacroGate::BishopRelperm:
            return "bishop_relperm";
    }
    return "unknown";
}

// ── Swelling-stress form (DIAGNOSTIC, NOT FOR PRODUCTION; fix (b) of
// ~/ogs-models/scratch/2026-10-04_swelling_stress_fixes_abc/). ────────────────
// step (default): the shipped telescoped step rule, bitwise 2a.
//   d sigma_sw = n_S (n_l_prev p_film_prev - n_l p_film_curr) I,
//   p_film = Pi - b p_conf with p_conf HELD FIXED (current iterate) in both
//   terms, n_S current in both terms, K(rho_d) of the current porosity in both.
// level: the same eigenstress as a function of the state,
//   L(n_l, n_S, K, sigma') = -n_S n_l [Pi(n_l; K) + b sigma'_mean],
//   d sigma_sw = L(curr) - L(prev) (I),
//   with L(prev) at the previous accepted state (n_l_prev, n_S_prev,
//   K(rho_d,prev), sigma'_mean_prev) and L(curr) at the current iterate. The
//   drain sigma'_mean is the effective stress of the previous Newton evaluation
//   (state_current sigma_eff) carried to the current strain by the ELASTIC
//   response, m_hat = m_lag - s_lag + K_d (eps_v - eps_v_lag), and the level
//   equation s = F - c (m_hat + s), c = n_S n_l b, is solved in closed form
//   (exact for a linear-elastic skeleton; at convergence m_hat + s = m, the
//   current sigma'_mean, for any skeleton). Implies a per-level K(rho_d) (fix
//   (a)).
enum class SwellingStressForm
{
    Step,
    Level
};

inline constexpr char const* toString(SwellingStressForm const form)
{
    switch (form)
    {
        case SwellingStressForm::Step:
            return "step";
        case SwellingStressForm::Level:
            return "level";
    }
    return "unknown";
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
    // prefactor`. The analytic dK/dphi tangent is wired in since 2026-06-12
    // (Vinay's approved completion), and into TWO Jacobian blocks, not one:
    // the p-u augmentation exchange tangent and the displacement-side
    // swelling-eigenstress tangent (RichardsMechanicsFEM-impl.h lines 5052
    // and 5223 at 7ec39ecf4c). Under the log-linear scheme in force since
    // 2026-08-26 (commit 1bb414ac05) that tangent is
    //   dK/dphi = -rho_SR * K(rho_d) * ln(K_r/K_l)/(x_r - x_l),
    // i.e. PROPORTIONAL TO THE LOCAL K value -- NOT the
    // -rho_SR*(K-linear segment slope) of the pre-2026-08-26 wording, which
    // is a per-segment constant and is superseded here (see
    // effectiveAugmentationPrefactorPhiDerivative below,
    // getSegmentSlopeLogLinear above and K_OF_RHO_D_LIVE.md) — the first
    // cut's omission note is historical.
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
    // Latched saturation gate of chi and k_rel at KKT-active points (design
    // part B.4 of DESIGN_FIXES.md; Vinay 2026-10-01 ~15:15 CEST). Appended LAST
    // so that the aggregate initialisation order of the earlier members is
    // unchanged. Default off = the shipped rules, bitwise.
    MicroCeilingSaturationGate micro_ceiling_saturation_gate =
        MicroCeilingSaturationGate::Off;
    // ── v4 switches (DESIGN_V4.md 2.1; ruling 2026-10-02). Appended LAST, in
    // this order, so the aggregate initialisation order of every earlier
    // member is unchanged. Defaults = the AB code, bitwise.
    // Drop T_m (DESIGN_V4.md 2.2): the micro part of the Biot volume-change
    // term, T_m = S_L rho_LR [(phi_m - phi_m,prev) + phi_m Delta eps_v]/dt, is
    // subtracted from the macro mass balance (extends Vinay's Q2 ruling "macro
    // storage is only macropores" to the volume-change term).
    bool macro_balance_drops_micro_biot_term = false;
    // 1a (DESIGN_V4.md 2.3).
    DarcyRelativePermeabilityMobility darcy_relative_permeability_mobility =
        DarcyRelativePermeabilityMobility::GaussPoint;
    // 1b (DESIGN_V4.md 2.4).
    MicroCeilingClosedMacroGate micro_ceiling_closed_macro_gate =
        MicroCeilingClosedMacroGate::Off;
    // Cells per decade of the log-uniform p_c grid of the 1a Kirchhoff table
    // (DESIGN_V4.md 2.3.2 item 4). A NUMERICAL choice, not physics: 2048 gives
    // a k interpolation error of 5.1e-6 and a Kirchhoff-potential error of
    // 1.9e-7 on the AB deck pair (MEASURED, REC/design_facts/
    // out_table_resolution.txt); pending Vinay's approval (DESIGN_V4.md Q5),
    // like micro_ceiling_scan_nodes_per_decade = 8. A new deck name beyond the
    // three briefed switches. >= 2; only meaningful with 1a on.
    int darcy_kirchhoff_cells_per_decade = 2048;
    // ── v5 probe switch (branch dsm_mass_conservation_v5_P_exact_2026-10-02,
    // ~/ogs-models/scratch/2026-10-02_kkt_v5_P_exact/DESIGN_V5.md 2.1-2.5).
    // On in MS33 cand. 2a/2b (main-loop reading, open). Appended LAST so that the aggregate initialisation
    // order of every earlier member is unchanged. Default false = the v4 code,
    // bitwise. When true, the a_S coefficient of the macro storage uses phi_M
    // of the previous converged step (PrevState<TransportPorosityData>)
    // instead of phi_M of the iterate, so that, together with the Biot term
    // (S_L at the new level) after the T_m drop, the discrete macro
    // accumulation is the exact difference Delta(rho_LR S_L phi_M) plus the
    // Biot strain term: the product term P = rho_LR Delta S_L Delta phi_M of
    // the L-books is removed. The uniqueness of phi_M,prev is CONDITIONAL on
    // keeping v4's S_L-new Biot/T_m form (an implementation choice, not a
    // ruling; review must-fix M2). Legal only with
    // macro_balance_drops_micro_biot_term = true and
    // macro_storage_uses_macro_porosity = true; FATAL at runtime on a
    // non-zero beta_LR or a0 (constant liquid density only, DESIGN_V5.md Q-A).
    bool macro_storage_exact_time_levels = false;
    // ── DIAGNOSTIC switches, swelling-stress fixes (a) and (b) (2026-10-04;
    // Vinay: "do all three on my mbp and compare"; scope in
    // ~/ogs-models/scratch/2026-10-04_swelling_stress_fixes_abc/STEP0.md).
    // DIAGNOSTIC, NOT FOR PRODUCTION; nothing adopted. Both default off ->
    // bit-identical to candidate 2a (fc14f19fb9). Appended LAST so the
    // aggregate initialisation order of every earlier member is unchanged.
    // (a) swelling_stress_K_level: in the step rule the Pi of the PREVIOUS level
    //     is evaluated at its own K(rho_d,prev) (previous accepted porosity)
    //     instead of the K of the current porosity. Acts only with the live
    //     K(rho_d) table; with a frozen K it is a no-op.
    bool swelling_stress_K_level = false;
    // (b) swelling_stress_form: step | level, see SwellingStressForm.
    SwellingStressForm swelling_stress_form = SwellingStressForm::Step;
};

// True when the previous level of the swelling stress is evaluated at its own
// K(rho_d): switch (a), or implied by the level form (b).
inline bool swellingStressPerLevelK(PotentialExchangeParameters const& p)
{
    return p.swelling_stress_K_level ||
           p.swelling_stress_form == SwellingStressForm::Level;
}

inline bool isSwellingStressLevelForm(PotentialExchangeParameters const* p)
{
    return p != nullptr && p->swelling_stress_form == SwellingStressForm::Level;
}

// ── Latched saturation gate: pure logic (design part B.4) ──────────────────
//
// Persistent per-integration-point latch L. At every evaluation of an
// integration point of the step that starts from the converged state with the
// latch L_old (the previous-step value, a constant of the step):
//
//   L_new = (status == Active) AND (L_old OR chi_deck(S_L) == 1)
//
// L_new is stored in the current state and becomes L_old when the step is
// accepted (the previous-state copy of the other history). chi_deck is the
// deck's own bishops_effective_stress property at the macro saturation of the
// evaluation; the trigger is the exact comparison chi_deck == 1.0 (the deck law
// BishopsSaturationCutoff returns exactly 0 or 1; no tolerance literal).
// L has no derivative: it is a constant of the step.
inline bool nextSaturatedLatch(bool const status_active, bool const latch_old,
                               double const chi_deck_at_S_L)
{
    return status_active && (latch_old || chi_deck_at_S_L == 1.0);
}

// The gate acts in an iterate when the point was latched at the end of the
// previous step AND the status of THIS iterate is Active. Otherwise (L_old
// false, or a released point) the shipped rules apply, bitwise.
inline bool saturationGateActs(MicroCeilingSaturationGate const level,
                               bool const latch_old, bool const status_active)
{
    return level != MicroCeilingSaturationGate::Off && latch_old &&
           status_active;
}

// Values of the Bishop factor at one integration point of one iterate.
struct BishopFactorValues
{
    double chi;         // [-] chi(S_L) of this iterate
    double chi_prev;    // [-] chi(S_L_prev), the previous converged step
    double dchi_dS_L;   // [-] dchi/dS_L of this iterate
};

// The gated Bishop factors of design part B.4: where the gate acts (previous
// latch true AND status Active in this iterate) chi = chi_deck(S = 1), the
// previous-step factor is the same (L_old true: the point was Active with the
// factor chi_deck(S = 1) at the end of the previous step), and dchi/dS = 0
// because S is a constant there. Everywhere else the deck's values pass
// through unchanged (bitwise). chi_deck is the deck's own
// bishops_effective_stress property as a function of S_L.
template <typename ChiDeck>
inline BishopFactorValues saturationGatedBishopFactors(
    MicroCeilingSaturationGate const level, bool const latch_old,
    bool const status_active, BishopFactorValues const& deck_values,
    ChiDeck&& chi_deck)
{
    if (!saturationGateActs(level, latch_old, status_active))
    {
        return deck_values;
    }
    double const chi_unit_saturation = chi_deck(1.0);  // [-]
    return {chi_unit_saturation, chi_unit_saturation, 0.0};
}

// k_rel is gated at level BishopRelperm only.
inline bool saturationGateActsOnRelativePermeability(
    MicroCeilingSaturationGate const level, bool const latch_old,
    bool const status_active)
{
    return level == MicroCeilingSaturationGate::BishopRelperm && latch_old &&
           status_active;
}

// ── 1b closed-macro gate: pure logic (DESIGN_V4.md 2.4.1) ───────────────────
// The gate acts in an iterate when the point was KKT-active at the end of the
// previous converged step with closed macro pores (phi_M,prev == 0, an exact
// comparison: on the KKT Active branch n_l = phi_s and phi_M = 0 exactly,
// DESIGN_V4.md 1.2; no tolerance literal) AND the status of THIS iterate is
// Active. prev_active and prev_phiM_zero are constants of the step; the
// iterate's status is not (release semantics, DESIGN_V4.md 2.4.5).
inline bool closedMacroGateActs(MicroCeilingClosedMacroGate const level,
                                bool const prev_active,
                                bool const prev_phiM_zero,
                                bool const status_active)
{
    return level != MicroCeilingClosedMacroGate::Off && prev_active &&
           prev_phiM_zero && status_active;
}

// Bishop factors at a point where the 1b gate acts at level bishop_relperm
// (DESIGN_V4.md 2.4.4): chi = chi_deck(1), dchi/dS = 0, and chi_prev =
// chi_deck(1) only if the gate acted at the end of the previous step
// (acted_prev); otherwise chi_prev is the incoming value (the deck's
// chi(S_L,prev), or Fix B's chi_deck(1) where Fix B acted). Elsewhere, and at
// level relperm, the incoming values pass through unchanged (bitwise).
template <typename ChiDeck>
inline BishopFactorValues closedMacroGatedBishopFactors(
    MicroCeilingClosedMacroGate const level, bool const gate_acts,
    bool const acted_prev, BishopFactorValues const& incoming,
    ChiDeck&& chi_deck)
{
    if (level != MicroCeilingClosedMacroGate::BishopRelperm || !gate_acts)
    {
        return incoming;
    }
    double const chi_unit_saturation = chi_deck(1.0);  // [-]
    return {chi_unit_saturation,
            acted_prev ? chi_unit_saturation : incoming.chi_prev, 0.0};
}

// True when any of the three v4 switches is on (Picard guard, labels).
inline bool anyV4SwitchOn(PotentialExchangeParameters const& p)
{
    return p.macro_balance_drops_micro_biot_term ||
           p.darcy_relative_permeability_mobility !=
               DarcyRelativePermeabilityMobility::GaussPoint ||
           p.micro_ceiling_closed_macro_gate != MicroCeilingClosedMacroGate::Off;
}

// ── v5 probe: exact time levels of the macro storage (DESIGN_V5.md 2.5) ────
// The a_S coefficient: phi_M of the previous converged step when the switch
// is on, else the incoming phi_M (pass-through, the same double, bitwise).
// No floating-point operation either way.
inline double macroStorageCoefficient(bool const exact_time_levels,
                                      double const phi_M,
                                      double const phi_M_prev)
{
    return exact_time_levels ? phi_M_prev : phi_M;
}

// Runtime admissibility of the switch at one integration point (DESIGN_V5.md
// 2.4): the exact-difference statement is derived for a constant liquid
// density (beta_LR == 0) and a0 == 0 (beta_SR == 0, enforced by the KKT
// FATAL as well). Exact comparisons, no tolerance literal: a Constant density
// gives dValue == 0 exactly.
inline bool macroStorageExactTimeLevelsAdmissible(double const beta_LR,
                                                  double const a0)
{
    return beta_LR == 0.0 && a0 == 0.0;
}

inline bool isMacroStorageExactTimeLevels(PotentialExchangeParameters const* p)
{
    return p != nullptr && p->macro_storage_exact_time_levels;
}

inline bool isKirchhoffElementMeanMobility(PotentialExchangeParameters const* p)
{
    return p != nullptr && p->darcy_relative_permeability_mobility ==
                               DarcyRelativePermeabilityMobility::
                                   KirchhoffElementMean;
}

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
