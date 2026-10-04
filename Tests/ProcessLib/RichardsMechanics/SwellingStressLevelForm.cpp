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
