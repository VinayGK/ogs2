// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause
//
// DIAGNOSTIC, NOT FOR PRODUCTION. Unit tests of the swelling-stress fixes (a)
// swelling_stress_K_level and (b) swelling_stress_form = level
// (branch dsm_sw_fix_ab_2026-10-04; scope and equations in
// ~/ogs-models/scratch/2026-10-04_swelling_stress_fixes_abc/).
//
// Physics anchors (no expected value of a material parameter is asserted):
//   - bit identity of the default (both switches off) against the shipped
//     call signature;
//   - algebraic identities of the coded forms against an independent
//     reconstruction from the public vdW helper (computeVanDerWaalsMicro-
//     Potential), i.e. Pi at each level and its own K(rho_d);
//   - the STATE-FUNCTION property of the level form: the increments over a
//     path telescope to L(end) - L(start) and are path independent, which the
//     shipped step rule is not once K, n_S or p_conf change between levels;
//   - finite-difference vs analytic identity of the drain-feedback tangent on
//     a linear model of the closed loop.
// Sample-state values mirror the prior approved unit tests of this directory
// (hamaker 6.0e-20 J, Sa 1000 m^2/kg, rho_SR 2650 kg/m^3), as in
// StrainedFilmPotential.cpp.

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include "ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h"

using namespace ProcessLib::RichardsMechanics;

namespace
{
using SwKv = MathLib::KelvinVector::KelvinVectorType<2>;
using SwKm = MathLib::KelvinVector::KelvinMatrixType<2>;

PotentialExchangeParameters sampleParams()
{
    PotentialExchangeParameters p;
    p.enabled = true;
    p.film_pressure_coupling = true;
    p.film_strain_coupling = FilmStrainCouplingMode::Off;
    p.film_energy_route = FilmEnergyRoute::Operational;
    p.micro_potential_convention = MicroPotentialConvention::NegativeAttractive;
    p.hamaker_constant = 6.0e-20;
    p.specific_surface = 1000.0;
    p.micro_solid_density_reference = 2650.0;
    p.micro_solid_volume_fraction_mode = MicroSolidVolumeFractionMode::Reference;
    p.micro_solid_volume_fraction_reference = 0.6;  // aggregate nS of the sample state
    p.potential_augmentation_exponent = 0.30 / (0.70 * 2650.0 * 1000.0);  // m
    p.potential_augmentation_prefactor = 20.0;  // J/kg, scalar fallback
    p.use_micro_liquid_density_for_micro_pressure = true;
    p.potential_augmentation_prefactor_live_dry_density = true;
    p.potential_augmentation_prefactor_vs_dry_density =
        std::make_shared<AugmentationPrefactorTable const>(
            std::vector<double>{1000.0, 2000.0},
            std::vector<double>{10.0, 50.0});  // J/kg vs kg/m^3
    return p;
}

// Pi at the level (n_l, rho_lR, phi) with its own K(rho_d), independent
// reconstruction from the public helper (film branch, strain coupling off).
double Pi_level(PotentialExchangeParameters const& p, double const n_l,
                double const rho_lR, double const phi)
{
    double const sign = microPotentialSignFactorFromParameters(p);
    double const nS_act = computeActiveMicroSolidVolumeFraction(
        n_l, PotentialExchangeLocalSolveContext{}, p);
    double const K = effectiveAugmentationPrefactor(p, phi);
    double const mu =
        computeVanDerWaalsMicroPotential(
            n_l, rho_lR, nS_act, p.micro_solid_density_reference,
            p.hamaker_constant, p.specific_surface, sign, K,
            p.potential_augmentation_exponent, 0.0,
            p.micro_water_content_floor)
            .mu_lR;
    return -rho_lR * mu;  // Pa
}

struct Level
{
    double n_l, rho_lR, phi, n_S, p_conf;
};

// L = -n_S n_l [Pi - b p_conf]
double L_of(PotentialExchangeParameters const& p, Level const& s, double b)
{
    return -s.n_S * s.n_l * (Pi_level(p, s.n_l, s.rho_lR, s.phi) - b * s.p_conf);
}

double increment(PotentialExchangeParameters const& p, Level const& prev,
                 Level const& curr, double b,
                 double const level_prev_used =
                     std::numeric_limits<double>::quiet_NaN(),
                 double* const level_curr_out = nullptr)
{
    SwKm const C_el = SwKm::Identity() * 1.5e8;
    double const rho_LR = 1000.0;
    SwKv const inc = computeSwellingStressIncrement<2>(
        prev.n_l, curr.n_l, curr.n_S, curr.rho_lR, prev.rho_lR, rho_LR, C_el, p,
        b, curr.p_conf, /*eps_v*/ 0.0, /*eps_v_prev*/ 0.0, curr.phi, prev.phi,
        prev.n_S, prev.p_conf, level_prev_used, level_curr_out);
    auto const& I2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(2)>::identity2;
    return inc.dot(I2) / I2.dot(I2);  // scalar on identity2 [Pa]
}
}  // namespace

// Default (both switches off): the new trailing arguments are never read; the
// increment is bitwise the one of the shipped signature.
TEST(SwellingStressLevelForm, DefaultOffIsBitIdenticalToShippedSignature)
{
    auto const p = sampleParams();
    SwKm const C_el = SwKm::Identity() * 1.5e8;
    Level const prev{0.20, 1100.0, 0.40, 0.55, 2.0e6};
    Level const curr{0.25, 1100.0, 0.39, 0.57, 3.0e6};
    SwKv const shipped = computeSwellingStressIncrement<2>(
        prev.n_l, curr.n_l, curr.n_S, curr.rho_lR, prev.rho_lR, 1000.0, C_el, p,
        1.0, curr.p_conf, 0.0, 0.0, curr.phi);
    SwKv const with_prev_inputs = computeSwellingStressIncrement<2>(
        prev.n_l, curr.n_l, curr.n_S, curr.rho_lR, prev.rho_lR, 1000.0, C_el, p,
        1.0, curr.p_conf, 0.0, 0.0, curr.phi, prev.phi, prev.n_S, prev.p_conf);
    for (int i = 0; i < shipped.size(); ++i)
    {
        EXPECT_EQ(shipped[i], with_prev_inputs[i]);  // bitwise
    }
}

// Fix (a): the previous level's Pi at its own K(rho_d,prev); with equal
// porosity at both levels it is a no-op.
TEST(SwellingStressLevelForm, KLevelUsesOwnKAtPreviousLevel)
{
    auto p = sampleParams();
    double const b = 1.0;
    Level const prev{0.20, 1100.0, 0.40, 0.55, 2.0e6};
    Level const curr{0.25, 1100.0, 0.36, 0.55, 2.0e6};  // phi changed: K moves

    // Shipped: K of the CURRENT porosity in both Pi.
    double const K_curr = effectiveAugmentationPrefactor(p, curr.phi);
    ASSERT_NE(K_curr, effectiveAugmentationPrefactor(p, prev.phi));
    double const shipped = increment(p, prev, curr, b);
    {
        double const expected =
            curr.n_S *
            (prev.n_l * (Pi_level(p, prev.n_l, prev.rho_lR, curr.phi) -
                         b * curr.p_conf) -
             curr.n_l * (Pi_level(p, curr.n_l, curr.rho_lR, curr.phi) -
                         b * curr.p_conf));
        EXPECT_NEAR(shipped, expected, 1e-12 * std::abs(expected));
    }

    // (a) on: prev Pi at K(phi_prev), curr Pi at K(phi_curr), p_conf and n_S
    // held at the CURRENT values (step rule kept).
    p.swelling_stress_K_level = true;
    double const with_a = increment(p, prev, curr, b);
    double const expected_a =
        curr.n_S *
        (prev.n_l * (Pi_level(p, prev.n_l, prev.rho_lR, prev.phi) -
                     b * curr.p_conf) -
         curr.n_l * (Pi_level(p, curr.n_l, curr.rho_lR, curr.phi) -
                     b * curr.p_conf));
    EXPECT_NEAR(with_a, expected_a, 1e-12 * std::abs(expected_a));
    EXPECT_NE(with_a, shipped);

    // No-op at equal porosity.
    Level const curr_same_phi{0.25, 1100.0, prev.phi, 0.55, 2.0e6};
    auto p_off = sampleParams();
    EXPECT_EQ(increment(p, prev, curr_same_phi, b),
              increment(p_off, prev, curr_same_phi, b));
}

// Fix (a), must-fix 2026-10-04: a step with UNCHANGED n_l and a changed K books
// the path, -W n_l [Pi(n_l, K) - Pi(n_l, K_prev)] (ruled equation (2.1), the
// loop of (1.5)). The early return at |dn_l| <= eps must not swallow it. The
// 2a step rule (switch off) keeps the return and books nothing; with (a) on and
// the same K at both levels (no change) the increment is exactly zero.
TEST(SwellingStressLevelForm, KLevelBooksChangedKAtUnchangedNl)
{
    auto p = sampleParams();
    double const b = 1.0;
    Level const prev{0.25, 1100.0, 0.40, 0.55, 2.0e6};
    Level const curr{0.25, 1100.0, 0.36, 0.55, 2.0e6};  // n_l same, phi -> K moves
    ASSERT_EQ(prev.n_l, curr.n_l);
    ASSERT_NE(effectiveAugmentationPrefactor(p, curr.phi),
              effectiveAugmentationPrefactor(p, prev.phi));

    // 2a (switch off): early return, nothing booked.
    EXPECT_EQ(increment(p, prev, curr, b), 0.0);

    // (a) on: -W n_l [Pi(n_l, K_curr) - Pi(n_l, K_prev)], W = curr.n_S.
    p.swelling_stress_K_level = true;
    double const expected =
        -curr.n_S * curr.n_l *
        (Pi_level(p, curr.n_l, curr.rho_lR, curr.phi) -
         Pi_level(p, prev.n_l, prev.rho_lR, prev.phi));
    double const got = increment(p, prev, curr, b);
    ASSERT_NE(expected, 0.0);
    EXPECT_NEAR(got, expected, 1e-12 * std::abs(expected));

    // (a) on, same K at both levels and unchanged n_l: exactly zero.
    Level const curr_same_phi{0.25, 1100.0, prev.phi, 0.55, 2.0e6};
    EXPECT_EQ(increment(p, prev, curr_same_phi, b), 0.0);

    // Telescoping: n_l unchanged over two steps with K moving 0.40 -> 0.38 ->
    // 0.36 books the same sum as one step 0.40 -> 0.36 (each Pi level at its
    // own K, the shared W and b p_conf terms cancel at constant n_l).
    Level const mid{0.25, 1100.0, 0.38, 0.55, 2.0e6};
    double const two_steps =
        increment(p, prev, mid, b) + increment(p, mid, curr, b);
    EXPECT_NEAR(two_steps, got, 1e-12 * std::abs(got));
}

// Fix (b), level form: increments telescope to L(end) - L(start) and are path
// independent; the step rule is not (K, n_S, p_conf change along the path).
TEST(SwellingStressLevelForm, LevelFormIsAStateFunction)
{
    auto p = sampleParams();
    double const b = 1.0;
    std::vector<Level> const pathA = {
        {0.050, 1100.0, 0.45, 0.52, 0.0},     {0.100, 1100.0, 0.44, 0.53, 1.0e6},
        {0.160, 1100.0, 0.42, 0.54, 3.0e6},   {0.230, 1100.0, 0.40, 0.56, 6.0e6},
        {0.300, 1100.0, 0.38, 0.58, 1.2e7}};
    // Same endpoints, different intermediate states.
    std::vector<Level> const pathB = {
        pathA.front(), {0.120, 1100.0, 0.43, 0.55, 8.0e6},
        {0.280, 1100.0, 0.39, 0.50, 2.0e6}, pathA.back()};

    auto p_step = sampleParams();
    p.swelling_stress_form = SwellingStressForm::Level;
    auto sum = [&](auto const& pp, std::vector<Level> const& path)
    {
        double s = 0.0;
        for (std::size_t i = 1; i < path.size(); ++i)
        {
            s += increment(pp, path[i - 1], path[i], b);
        }
        return s;
    };
    double const L_end_minus_start =
        L_of(p, pathA.back(), b) - L_of(p, pathA.front(), b);
    double const levelA = sum(p, pathA);
    double const levelB = sum(p, pathB);
    EXPECT_NEAR(levelA, L_end_minus_start, 1e-12 * std::abs(L_end_minus_start));
    EXPECT_NEAR(levelB, L_end_minus_start, 1e-12 * std::abs(L_end_minus_start));

    // The step rule (shipped) is path dependent on the same two paths.
    double const stepA = sum(p_step, pathA);
    double const stepB = sum(p_step, pathB);
    EXPECT_GT(std::abs(stepA - stepB), 1e-3 * std::abs(L_end_minus_start));
}

// Level form = step rule when K, n_S and p_conf are equal at both levels.
TEST(SwellingStressLevelForm, LevelEqualsStepWhenOnlyNlChanges)
{
    auto p_step = sampleParams();
    auto p_level = sampleParams();
    p_level.swelling_stress_form = SwellingStressForm::Level;
    double const b = 1.0;
    Level const prev{0.20, 1100.0, 0.40, 0.55, 2.0e6};
    Level const curr{0.26, 1100.0, 0.40, 0.55, 2.0e6};
    double const s = increment(p_step, prev, curr, b);
    double const l = increment(p_level, prev, curr, b);
    EXPECT_NEAR(l, s, 1e-12 * std::abs(s));
}

// Level form: no early return at dn_l = 0; sigma' and n_S moving at constant
// n_l change the stress by the drain term alone: dL = n_S n_l b dp - dn_S n_l p.
TEST(SwellingStressLevelForm, LevelFormRespondsToDrainAtConstantNl)
{
    auto p_step = sampleParams();
    auto p_level = sampleParams();
    p_level.swelling_stress_form = SwellingStressForm::Level;
    double const b = 1.0;
    Level const prev{0.25, 1100.0, 0.40, 0.55, 2.0e6};
    Level const curr{0.25, 1100.0, 0.40, 0.55, 5.0e6};
    EXPECT_EQ(increment(p_step, prev, curr, b), 0.0);  // shipped early return
    double const expected = 0.55 * 0.25 * b * (5.0e6 - 2.0e6);
    EXPECT_NEAR(increment(p_level, prev, curr, b), expected,
                1e-12 * std::abs(expected));
}

// Drain-feedback derivative ds/d eps of the level form on a linear model of
// the closed loop: sigma' = C_cons (eps + C_el^{-1} I s), s = s0 - c m,
// m = tr(sigma')/3. FD of the solved s(eps) against the helper.
TEST(SwellingStressLevelForm, DrainFeedbackTangentMatchesFiniteDifference)
{
    // Anisotropic-looking but symmetric positive definite C_el, and a softer
    // consistent tangent C_cons (as in a plastic state).
    SwKm C_el = SwKm::Zero();
    C_el.diagonal() << 3.0e8, 2.5e8, 2.8e8, 1.1e8;
    C_el(0, 1) = C_el(1, 0) = 0.9e8;
    C_el(0, 2) = C_el(2, 0) = 0.8e8;
    C_el(1, 2) = C_el(2, 1) = 0.7e8;
    SwKm const C_cons = 0.6 * C_el;
    auto const& I2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(2)>::identity2;
    double const c = 0.23;
    double const s0 = -4.0e6;
    SwKv const eps0 = (SwKv() << 1.0e-3, -2.0e-3, 5.0e-4, 3.0e-4).finished();

    auto solve_s = [&](SwKv const& eps)
    {
        // s = s0 - c m(eps, s); m linear in s: m = a + g s.
        SwKv const sig_at_s0 = C_cons * eps;                     // s = 0 part
        SwKv const dsig_ds = C_cons * (C_el.inverse() * I2);     // per unit s
        double const a = I2.dot(sig_at_s0) / 3.0;
        double const g = I2.dot(dsig_ds) / 3.0;
        return (s0 - c * a) / (1.0 + c * g);
    };
    SwKv const analytic = swellingLevelDrainFeedbackDsDeps<2>(C_cons, C_el, c);
    for (int k = 0; k < 4; ++k)
    {
        double const h = 1e-7;
        SwKv ep = eps0, em = eps0;
        ep[k] += h;
        em[k] -= h;
        double const fd = (solve_s(ep) - solve_s(em)) / (2.0 * h);
        EXPECT_NEAR(fd, analytic[k], 1e-6 * std::abs(analytic[k]) + 1e-3)
            << "component " << k;
    }

    // Linear-elastic isotropic limit: ds/d eps_v = -c K/(1+c) on identity2.
    double const K_d = 1.5e8;
    SwKm const C_iso = SwKm::Identity() * K_d;  // sample; g = 1/3 * tr(I I^T)...
    SwKv const a_iso = swellingLevelDrainFeedbackDsDeps<2>(C_iso, C_iso, c);
    // For C = C_el: g = I^T I / 3 = 1 (3 normal Kelvin components in 2D plane
    // strain: identity2 = (1,1,1,0)) and r = (1/3) C^T I.
    EXPECT_NEAR(a_iso[0], -c * (K_d / 3.0) / (1.0 + c * 1.0),
                1e-12 * K_d);
}

// Level form with the stored level L_prev (state SwellingLevelUsed): the lag of
// sigma'_mean inside a step (p_conf of the LAGGED iterate differs from the
// accepted value) does not accumulate over the steps. sigma_sw(end) =
// L(end, lagged) - L(start) exactly, whatever the lag of the earlier steps.
// Recomputing L_prev from the accepted previous state (no stored level) lets the
// lag of every step enter the sum.
TEST(SwellingStressLevelForm, StoredLevelPreventsLagAccumulation)
{
    auto p = sampleParams();
    p.swelling_stress_form = SwellingStressForm::Level;
    double const b = 1.0;
    // Accepted (true) states along the path and the lagged p_conf used in the
    // evaluation of each step (true value + lag).
    std::vector<Level> const accepted = {
        {0.050, 1100.0, 0.45, 0.52, 0.0},     {0.100, 1100.0, 0.44, 0.53, 1.0e6},
        {0.160, 1100.0, 0.42, 0.54, 3.0e6},   {0.230, 1100.0, 0.40, 0.56, 6.0e6},
        {0.300, 1100.0, 0.38, 0.58, 1.2e7}};
    std::vector<double> const lag = {0.0, 2.0e5, -3.0e5, 4.0e5, -1.0e5};
    auto lagged = [&](std::size_t i)
    {
        Level l = accepted[i];
        l.p_conf += lag[i];
        return l;
    };

    // (i) stored level chain.
    double sum_chain = 0.0;
    double used = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t i = 1; i < accepted.size(); ++i)
    {
        double L_curr = std::numeric_limits<double>::quiet_NaN();
        sum_chain += increment(p, accepted[i - 1], lagged(i), b, used, &L_curr);
        used = L_curr;
    }
    double const L_end_lagged = L_of(p, lagged(accepted.size() - 1), b);
    double const L_start = L_of(p, accepted.front(), b);
    EXPECT_NEAR(sum_chain, L_end_lagged - L_start,
                1e-12 * std::abs(L_end_lagged - L_start));

    // (ii) recomputed from the accepted previous state: lags accumulate.
    double sum_recomputed = 0.0;
    for (std::size_t i = 1; i < accepted.size(); ++i)
    {
        sum_recomputed += increment(p, accepted[i - 1], lagged(i), b);
    }
    double expected_accum = 0.0;  // sum of n_S n_l b lag over the steps 1..4
    for (std::size_t i = 1; i < accepted.size(); ++i)
    {
        expected_accum += accepted[i].n_S * accepted[i].n_l * b * lag[i];
    }
    EXPECT_NEAR(sum_recomputed - (L_of(p, accepted.back(), b) - L_start),
                expected_accum, 1e-9 * std::abs(expected_accum));
}

// The factor on a frozen-sigma' partial dF (the live-K chain of the Jacobian):
// the closed loop s = F(eps) - c m(eps, s), m = tr(sigma')/3,
// sigma' = C_cons (eps + C_el^{-1} I s), has ds/d eps_k = (dF/d eps_k - c r_k) /
// (1 + c g). Finite difference of the solved s(eps) against
// swellingLevelImplicitFactor and swellingLevelDrainFeedbackDsDeps, with F
// carrying a strain dependence (F = F0 + f1 eps_v, the way K(phi(eps_v)) enters).
TEST(SwellingStressLevelForm, ImplicitFactorOnFrozenPartialMatchesFiniteDifference)
{
    SwKm C_el = SwKm::Zero();
    C_el.diagonal() << 3.0e8, 2.5e8, 2.8e8, 1.1e8;
    C_el(0, 1) = C_el(1, 0) = 0.9e8;
    C_el(0, 2) = C_el(2, 0) = 0.8e8;
    C_el(1, 2) = C_el(2, 1) = 0.7e8;
    SwKm const C_cons = 0.6 * C_el;
    auto const& I2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(2)>::identity2;
    double const c = 0.31;
    double const F0 = -6.0e6;
    double const f1 = -2.0e7;  // dF/d eps_v [Pa]
    SwKv const eps0 = (SwKv() << 2.0e-3, -1.0e-3, 4.0e-4, 2.0e-4).finished();

    auto solve_s = [&](SwKv const& eps)
    {
        double const F = F0 + f1 * I2.dot(eps);
        SwKv const dsig_ds = C_cons * (C_el.inverse() * I2);
        double const a = I2.dot(C_cons * eps) / 3.0;
        double const g = I2.dot(dsig_ds) / 3.0;
        return (F - c * a) / (1.0 + c * g);
    };
    double const f = swellingLevelImplicitFactor<2>(
        C_cons, C_el.inverse().eval(), c);
    SwKv const drain = swellingLevelDrainFeedbackDsDeps<2>(C_cons, C_el, c);
    for (int k = 0; k < 4; ++k)
    {
        double const h = 1e-7;
        SwKv ep = eps0, em = eps0;
        ep[k] += h;
        em[k] -= h;
        double const fd = (solve_s(ep) - solve_s(em)) / (2.0 * h);
        // coded: factor on the frozen partial (dF/d eps_k = f1 I2[k]) plus the
        // drain-feedback vector.
        double const coded = f * f1 * I2[k] + drain[k];
        EXPECT_NEAR(fd, coded, 1e-6 * std::abs(coded) + 1e-3) << "component " << k;
    }
    // c = 0: the factor is exactly 1 (no drain, the frozen partial is the whole).
    EXPECT_EQ(swellingLevelImplicitFactor<2>(C_cons, C_el.inverse().eval(), 0.0),
              1.0);
}

// Closed form of the level equation with the elastic prediction m_hat: the
// returned increment satisfies s_new = s_prev + inc and
// s_new = L(n_l, n_S, K, m = m_hat + s_new) - L_prev + s_prev, i.e. the level
// formula evaluated at the CURRENT mean effective stress (no lag), for any m_hat.
TEST(SwellingStressLevelForm, ClosedFormSatisfiesLevelEquationAtCurrentStress)
{
    auto p = sampleParams();
    p.swelling_stress_form = SwellingStressForm::Level;
    double const b = 1.0;
    Level const prev{0.20, 1100.0, 0.40, 0.55, 2.0e6};
    Level const curr{0.27, 1100.0, 0.39, 0.57, 0.0};  // p_conf unused here
    double const L_prev_used = L_of(p, prev, b);
    double const s_prev = -3.0e6;
    double const m_hat = -4.0e6;  // [Pa], tension positive
    SwKm const C_el = SwKm::Identity() * 1.5e8;
    double L_curr = 0.0;
    SwKv const inc = computeSwellingStressIncrement<2>(
        prev.n_l, curr.n_l, curr.n_S, curr.rho_lR, prev.rho_lR, 1000.0, C_el, p,
        b, /*p_conf lagged, unused*/ 1.0e6, 0.0, 0.0, curr.phi, prev.phi,
        prev.n_S, prev.p_conf, L_prev_used, &L_curr, m_hat, s_prev);
    auto const& I2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(2)>::identity2;
    double const s_new = s_prev + inc.dot(I2) / I2.dot(I2);
    // level formula at the current mean effective stress m = m_hat + s_new
    Level at_new = curr;
    at_new.p_conf = -(m_hat + s_new);  // p_conf = -m
    double const L_formula = L_of(p, at_new, b);
    EXPECT_NEAR(L_curr, L_formula, 1e-9 * std::abs(L_formula));
    EXPECT_NEAR(s_new, L_formula - L_prev_used + s_prev,
                1e-9 * std::abs(s_new));
    // c = n_S n_l b -> 0 (b = 0): no drain, closed form = plain level.
    double L0 = 0.0;
    computeSwellingStressIncrement<2>(
        prev.n_l, curr.n_l, curr.n_S, curr.rho_lR, prev.rho_lR, 1000.0, C_el, p,
        0.0, 1.0e6, 0.0, 0.0, curr.phi, prev.phi, prev.n_S, prev.p_conf,
        L_of(p, prev, 0.0), &L0, m_hat, s_prev);
    EXPECT_NEAR(L0, L_of(p, curr, 0.0), 1e-12 * std::abs(L0));
}

namespace
{
// A one-integration-point state for updateSwellingState: previous accepted
// state (n_l, phi, zero stress, no stored level) and a current state.
struct SwState
{
    StatefulData<2> cur;
    StatefulDataPrev<2> prev;
};

SwState makeInitialSwState(double const n_l_prev, double const n_l,
                           double const phi_prev, double const phi)
{
    SwState s;
    auto const& I2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(2)>::identity2;
    (void)I2;
    std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<2>>(s.cur)
        .sigma_eff.setZero();
    std::get<ProcessLib::ThermoRichardsMechanics::
                 ConstitutiveStress_StrainTemperature::SwellingDataStateful<2>>(
        s.cur)
        .sigma_sw.setZero();
    std::get<ProcessLib::ConstitutiveRelations::MechanicalStrainData<2>>(s.cur)
        .eps_m.setZero();
    std::get<StrainData<2>>(s.cur).eps.setZero();
    std::get<PrevState<ProcessLib::ConstitutiveRelations::
                           EffectiveStressData<2>>>(s.prev)
        ->sigma_eff.setZero();
    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                           ConstitutiveStress_StrainTemperature::
                               SwellingDataStateful<2>>>(s.prev)
        ->sigma_sw.setZero();
    **std::get<PrevState<MicroWaterContent>>(s.prev) = n_l_prev;
    *std::get<MicroWaterContent>(s.cur) = n_l;
    *std::get<MicroLiquidDensity>(s.cur) = 1100.0;
    **std::get<PrevState<MicroLiquidDensity>>(s.prev) = 1100.0;
    std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(s.cur).phi = phi;
    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::PorosityData>>(
        s.prev)
        ->phi = phi_prev;
    // transport porosity phi_M = (phi - n_l)/(1 - n_l), the split of the code
    std::get<ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(s.cur)
        .phi = (phi - n_l) / (1.0 - n_l);
    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                           TransportPorosityData>>(s.prev)
        ->phi = (phi_prev - n_l_prev) / (1.0 - n_l_prev);
    *std::get<SwellingLevelUsed>(s.cur) = 0.0;
    **std::get<PrevState<SwellingLevelUsed>>(s.prev) = 0.0;
    *std::get<SwellingLagVolRatio>(s.cur) = 0.0;
    *std::get<SwellingLagStress>(s.cur) = 0.0;
    *std::get<SwellingLevelAssumedStress>(s.cur) = 0.0;
    return s;
}

double meanOf(SwKv const& v)
{
    auto const& I2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(2)>::identity2;
    return v.dot(I2) / 3.0;
}
}  // namespace

// Review must-fix 2026-10-04 (fix (b) tangent search, MEASURED on the Model I
// live-K probe): updateSwellingState read the lagged swelling stress s_lag from
// the state sigma_sw, which updateSwellingStressAndVolumetricStrain has reset to
// the previous accepted value earlier in the same evaluation. m_hat then carried
// the whole lagged s, and the closed form iterated a loop of gain -c/(1+c) at a
// FIXED strain (residual ratio -0.155 to -0.30, Newton eps_v ratio -0.40).
// Here a linear-elastic skeleton is emulated by hand (m = m0 + K_d eps_v + s,
// the incremental form of the code) and the evaluation is repeated at ONE strain
// with the reset of sigma_sw in between: from the second evaluation on the
// swelling stress must not change (no lag loop), and it must satisfy the level
// equation at the actual mean effective stress.
TEST(SwellingStressLevelForm, RepeatedEvaluationAtFixedStrainHasNoLagLoop)
{
    auto p = sampleParams();
    p.swelling_stress_form = SwellingStressForm::Level;
    p.swelling_stress_K_level = true;
    double const n_l_prev = 0.20, n_l = 0.27, phi_prev = 0.40, phi = 0.39;
    SwState st = makeInitialSwState(n_l_prev, n_l, phi_prev, phi);
    SwKm const C_el = SwKm::Identity() * 1.5e8;
    double const K_d = drainedBulkModulusFromStiffness<2>(C_el);
    auto const& I2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(2)>::identity2;
    double const eps_v = -2.0e-3, eps_v_prev = 0.0;
    double const m_prev = -1.0e6;  // stress of the previous accepted state
    std::get<PrevState<ProcessLib::ConstitutiveRelations::
                           EffectiveStressData<2>>>(st.prev)
        ->sigma_eff = m_prev * I2;
    std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<2>>(st.cur)
        .sigma_eff = m_prev * I2;
    MaterialPropertyLib::Phase solid{
        MaterialPropertyLib::PhaseName::Solid, {}, nullptr};
    ParameterLib::SpatialPosition const x_position;
    MPL::VariableArray variables, variables_prev;
    variables.volumetric_strain = eps_v;
    variables_prev.volumetric_strain = eps_v_prev;

    auto const& sigma_sw_prev = *std::get<
        PrevState<ProcessLib::ThermoRichardsMechanics::
                      ConstitutiveStress_StrainTemperature::
                          SwellingDataStateful<2>>>(st.prev);
    auto& sigma_sw_cur =
        std::get<ProcessLib::ThermoRichardsMechanics::
                     ConstitutiveStress_StrainTemperature::
                         SwellingDataStateful<2>>(st.cur);
    std::vector<double> s_hist;
    std::vector<double> L_hist;
    for (int k = 0; k < 5; ++k)
    {
        // what updateSwellingStressAndVolumetricStrain does first
        sigma_sw_cur.sigma_sw = sigma_sw_prev.sigma_sw;
        updateSwellingState<2>(solid, 1000.0, C_el, st.cur, st.prev, variables,
                               variables_prev, x_position, 0.0, 1.0, &p, 1.0);
        double const s = meanOf(sigma_sw_cur.sigma_sw);
        s_hist.push_back(s);
        L_hist.push_back(*std::get<SwellingLevelUsed>(st.cur));
        // linear-elastic stress update (incremental form of the code)
        double const m = m_prev + K_d * (eps_v - eps_v_prev) + s - 0.0;
        std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<2>>(
            st.cur)
            .sigma_eff = m * I2;
    }
    // evaluation 0 uses the lagged form (nothing stored yet); from evaluation 1
    // on the closed form with the stored lag: the iterates must be identical.
    for (int k = 2; k < 5; ++k)
    {
        EXPECT_NEAR(s_hist[k], s_hist[1], 1e-9 * std::abs(s_hist[1]))
            << "evaluation " << k << " (a lag loop of gain -c/(1+c) would give "
            << "a geometric sequence here)";
    }
    // level equation at the actual mean effective stress of the last evaluation
    double const m_last = m_prev + K_d * (eps_v - eps_v_prev) + s_hist[4];
    double const n_S = 1.0 - (phi - n_l) / (1.0 - n_l);
    double const n_S_prev =
        1.0 - (phi_prev - n_l_prev) / (1.0 - n_l_prev);
    // L_prev recomputed at the previous state (no stored level), as the code
    Level const prev_lv{n_l_prev, 1100.0, phi_prev, n_S_prev, -m_prev};
    double const L_prev = L_of(p, prev_lv, 1.0);
    Level const curr_lv{n_l, 1100.0, phi, n_S, -m_last};
    double const L_curr = L_of(p, curr_lv, 1.0);
    EXPECT_NEAR(s_hist[4], 0.0 + L_curr - L_prev, 1e-8 * std::abs(s_hist[4]));
    EXPECT_NEAR(L_hist[4], L_curr, 1e-8 * std::abs(L_curr));
}

// Review must-fix 3: with a skeleton that is softer than the elastic prediction
// (consistent tangent g < 1 on the swelling stress), the loop solves
// phi(s) = (s - s_0) + c (m(s) - m_ev) = 0 by Newton and leaves the state on
// the level equation; for the elastic skeleton (g = 1, m_ev = m) it does not
// run. m(s) = m_a + g (s - s_0) is the emulated stress response.
TEST(SwellingStressLevelForm, LevelCorrectionLoopSolvesTheEquationForASofterSkeleton)
{
    auto const& I2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(2)>::identity2;
    SwKm const C_el = SwKm::Identity() * 1.5e8;
    for (double const g : {1.0, 0.6, 0.3})
    {
        SwState st = makeInitialSwState(0.20, 0.27, 0.40, 0.39);
        auto& sigma_sw =
            std::get<ProcessLib::ThermoRichardsMechanics::
                         ConstitutiveStress_StrainTemperature::
                             SwellingDataStateful<2>>(st.cur);
        auto& sigma_eff =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<2>>(
                st.cur);
        double const s_0 = -5.0e6, m_ev = -3.0e6, L_0 = -6.0e6;
        double const m_a = (g == 1.0) ? m_ev : -2.0e6;  // actual m at s_0
        double const n_S = std::max(
            1e-16,
            1.0 - std::get<ProcessLib::ThermoRichardsMechanics::
                               TransportPorosityData>(st.cur)
                      .phi);
        double const c = n_S * 0.27 * 1.0;
        sigma_sw.sigma_sw = s_0 * I2;
        sigma_eff.sigma_eff = m_a * I2;
        *std::get<SwellingLevelAssumedStress>(st.cur) = m_ev;
        *std::get<SwellingLevelUsed>(st.cur) = L_0;
        *std::get<SwellingLagStress>(st.cur) = s_0;
        SwKm C_cons = g * C_el;
        int n_updates = 0;
        int const passes = iterateSwellingLevelWithStressUpdate<2>(
            st.cur, C_el, C_cons, 1.0,
            [&]()
            {
                ++n_updates;
                double const s = meanOf(sigma_sw.sigma_sw);
                sigma_eff.sigma_eff = (m_a + g * (s - s_0)) * I2;  // emulated
                return SwKm(g * C_el);
            });
        double const s = meanOf(sigma_sw.sigma_sw);
        double const m = meanOf(sigma_eff.sigma_eff);
        double const phi_end = (s - s_0) + c * (m - m_ev);
        EXPECT_NEAR(phi_end, 0.0, 1e-3) << "g = " << g;  // Pa
        EXPECT_EQ(passes, n_updates);
        if (g == 1.0)
        {
            EXPECT_EQ(passes, 0);  // elastic: m_a = m_ev, nothing to correct
        }
        else
        {
            EXPECT_GE(passes, 1);
            EXPECT_LE(passes, 3);  // the emulated response is linear in s
            // closed solution of s = s_0 - c (m_a + g (s - s_0) - m_ev)
            double const s_exact =
                s_0 - c * (m_a - m_ev) / (1.0 + c * g);
            EXPECT_NEAR(s, s_exact, 1e-6 * std::abs(s_exact));
            // stored level and lag stress follow the corrected s
            EXPECT_NEAR(*std::get<SwellingLevelUsed>(st.cur), L_0 + (s - s_0),
                        1e-6);
            EXPECT_NEAR(*std::get<SwellingLagStress>(st.cur), s, 1e-6);
        }
    }
}
