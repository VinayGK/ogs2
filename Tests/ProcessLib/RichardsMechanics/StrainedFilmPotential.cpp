// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause
//
// Strained-film disjoining law h(w_m, eps_v) — unit tests.
// Design: ProcessLib/RichardsMechanics/DSM/STRAINED_FILM_IMPLEMENTATION.md.
// Physics anchors (CLAUDE.md §3): derived identities (FD-vs-analytic chains,
// force-balance inversion residual), analytical limits (zero strain + zero
// load reduction), sign-only physical limits (load raises the potential).
// No fitted expected values; tolerances derive from the FD step / solver
// residual scales.
//
// Sample-state parameter values mirror the prior approved unit tests in
// DSMMicroMacroSingleIntegrationPoint.cpp (hamaker 6.0e-20 J, Sa 1000 m^2/kg,
// rho_SR 2650 kg/m^3, nS 0.6) — citation source: prior user-approved test
// code in this repository.

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h"

using namespace ProcessLib::RichardsMechanics;

namespace
{
struct StrainedFilmSampleState
{
    double n_l = 0.3;
    double rho_lR = 1100.0;  // confined micro-liquid density scale, mirrors
                             // the memory note "micro EOS ~1100 kg/m^3"
    double active_nS = 0.6;
    double rho_SR = 2650.0;
    double hamaker = 6.0e-20;
    double Sa = 1000.0;
    // NegativeAttractive convention (the MS33 PRJ family): mu_lR < 0,
    // Pi = -rho*mu_lR > 0 (repulsive operational disjoining pressure).
    double sign = -1.0;
    double K_aug = 0.0;
    double lambda_aug = 0.0;
    double floor = 0.0;
};

double bareMu(StrainedFilmSampleState const& st, double const w)
{
    return computeVanDerWaalsMicroPotential(
               w, st.rho_lR, st.active_nS, st.rho_SR, st.hamaker, st.Sa,
               st.sign, st.K_aug, st.lambda_aug, 0.0, st.floor)
        .mu_lR;
}

double barePi(StrainedFilmSampleState const& st, double const w)
{
    return -st.rho_lR * bareMu(st, w);
}

StrainedFilmStateData state(StrainedFilmSampleState const& st,
                            FilmStrainCouplingMode const mode,
                            FilmStrainKappaMode const kappa,
                            double const eps_v, double const p_conf)
{
    return computeStrainedFilmState(mode, kappa, st.n_l, st.active_nS, eps_v,
                                    p_conf, st.rho_lR, st.rho_SR, st.hamaker,
                                    st.Sa, st.sign, st.K_aug, st.lambda_aug,
                                    st.floor, st.rho_lR);
}
}  // namespace

// Anchor: approved baseline — the defaults must stay Off/Aggregate so every
// existing PRJ is bit-for-bit unaffected.
TEST(RichardsMechanicsStrainedFilm, DefaultsAreOffAggregate)
{
    PotentialExchangeParameters params;
    EXPECT_EQ(params.film_strain_coupling, FilmStrainCouplingMode::Off);
    EXPECT_EQ(params.film_strain_kappa, FilmStrainKappaMode::Aggregate);
}

// Anchor: derived identity — FD-vs-analytic chains of the kinematic state.
TEST(RichardsMechanicsStrainedFilm, KinematicChainsFDConsistent)
{
    StrainedFilmSampleState st;
    double const eps_v = -0.02;  // compression
    double const p_conf = 1.0e6;
    double const d = 1e-7;

    for (auto const kappa :
         {FilmStrainKappaMode::Aggregate, FilmStrainKappaMode::Unity})
    {
        auto const s0 =
            state(st, FilmStrainCouplingMode::Kinematic, kappa, eps_v, p_conf);
        double const kappa_value =
            kappa == FilmStrainKappaMode::Aggregate ? st.active_nS : 1.0;
        EXPECT_NEAR(s0.w_eff, st.n_l * (1.0 + kappa_value * eps_v),
                    1e-14 * st.n_l);

        // d w_eff / d eps_v by central FD on the helper itself.
        auto const sp = state(st, FilmStrainCouplingMode::Kinematic, kappa,
                              eps_v + d, p_conf);
        auto const sm = state(st, FilmStrainCouplingMode::Kinematic, kappa,
                              eps_v - d, p_conf);
        double const fd_deps = (sp.w_eff - sm.w_eff) / (2.0 * d);
        EXPECT_NEAR(s0.dw_eff_deps_v, fd_deps, 1e-6 * std::abs(fd_deps));

        // d w_eff / d n_l by central FD in n_l.
        StrainedFilmSampleState stp = st;
        stp.n_l += d;
        StrainedFilmSampleState stm = st;
        stm.n_l -= d;
        double const fd_dnl =
            (state(stp, FilmStrainCouplingMode::Kinematic, kappa, eps_v,
                   p_conf)
                 .w_eff -
             state(stm, FilmStrainCouplingMode::Kinematic, kappa, eps_v,
                   p_conf)
                 .w_eff) /
            (2.0 * d);
        EXPECT_NEAR(s0.dw_eff_dnl, fd_dnl, 1e-6 * std::abs(fd_dnl));
    }
}

// Anchor: derived identity — on the loaded branch the inverted state must
// satisfy the film force balance Pi(w_eff) = p_conf to the solver residual
// scale (1e-9 relative; the Newton terminates at 1e-12 relative).
TEST(RichardsMechanicsStrainedFilm, EquilibriumInversionSolvesForceBalance)
{
    StrainedFilmSampleState st;
    double const Pi_unloaded = barePi(st, st.n_l);
    ASSERT_GT(Pi_unloaded, 0.0);

    // Loaded branch: p_conf above the unloaded disjoining pressure.
    double const p_loaded = 2.0 * Pi_unloaded;
    auto const s_loaded =
        state(st, FilmStrainCouplingMode::Equilibrium,
              FilmStrainKappaMode::Aggregate, 0.0, p_loaded);
    EXPECT_TRUE(s_loaded.loaded_branch);
    EXPECT_LT(s_loaded.w_eff, st.n_l);  // squeezed film
    EXPECT_NEAR(barePi(st, s_loaded.w_eff), p_loaded, 1e-9 * p_loaded);

    // Unloaded branch: p_conf below the branch point -> identity state.
    auto const s_unloaded =
        state(st, FilmStrainCouplingMode::Equilibrium,
              FilmStrainKappaMode::Aggregate, 0.0, 0.5 * Pi_unloaded);
    EXPECT_FALSE(s_unloaded.loaded_branch);
    EXPECT_DOUBLE_EQ(s_unloaded.w_eff, st.n_l);
    EXPECT_DOUBLE_EQ(s_unloaded.dw_eff_dnl, 1.0);
}

// Anchor: derived identity — the equilibrium inversion with the exponential
// augmentation active (exercises the Newton off the pure cubic seed).
TEST(RichardsMechanicsStrainedFilm, EquilibriumInversionWithAugmentation)
{
    StrainedFilmSampleState st;
    // Augmentation amplitude/decay from the prior approved dd1600 PRJ family
    // (potential_augmentation_prefactor 103879 J/kg, exponent 7.5e-7 m) —
    // citation source: Tests/Data/.../ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj.
    st.K_aug = 103879.0;
    st.lambda_aug = 7.5e-7;

    double const Pi_unloaded = barePi(st, st.n_l);
    ASSERT_GT(Pi_unloaded, 0.0);
    double const p_loaded = 3.0 * Pi_unloaded;
    auto const s =
        state(st, FilmStrainCouplingMode::Equilibrium,
              FilmStrainKappaMode::Aggregate, 0.0, p_loaded);
    ASSERT_TRUE(s.loaded_branch);
    EXPECT_NEAR(barePi(st, s.w_eff), p_loaded, 1e-9 * p_loaded);
}

// Anchor: analytical limit — at zero strain and zero confining pressure both
// strained modes reduce EXACTLY to the frozen-geometry evaluation point.
TEST(RichardsMechanicsStrainedFilm, ZeroStrainZeroLoadReducesToBareState)
{
    StrainedFilmSampleState st;
    for (auto const mode : {FilmStrainCouplingMode::Kinematic,
                            FilmStrainCouplingMode::Equilibrium})
    {
        auto const s =
            state(st, mode, FilmStrainKappaMode::Aggregate, 0.0, 0.0);
        EXPECT_DOUBLE_EQ(s.w_eff, st.n_l);
        EXPECT_DOUBLE_EQ(s.dw_eff_dnl, 1.0);
        EXPECT_DOUBLE_EQ(s.dw_eff_deps_v, mode ==
                             FilmStrainCouplingMode::Kinematic
                             ? st.n_l * st.active_nS
                             : 0.0);
        EXPECT_FALSE(s.loaded_branch);
    }
}

// Anchor: physical limit (sign only) — through the fold point, raising the
// confining pressure at fixed water content must RAISE mu_lR (the Derjaguin
// load term: squeezing confined liquid raises its chemical potential; the
// expulsion channel). No magnitude asserted.
TEST(RichardsMechanicsStrainedFilm, LoadRaisesPotentialAtFixedWaterContent)
{
    StrainedFilmSampleState st;
    PotentialExchangeParameters params;
    params.enabled = true;
    params.hamaker_constant = st.hamaker;
    params.specific_surface = st.Sa;
    params.micro_solid_density_reference = st.rho_SR;
    params.micro_solid_volume_fraction_reference = st.active_nS;
    params.micro_potential_convention =
        MicroPotentialConvention::NegativeAttractive;
    params.film_pressure_coupling = true;
    params.film_strain_coupling = FilmStrainCouplingMode::Kinematic;
    params.film_strain_kappa = FilmStrainKappaMode::Aggregate;

    auto const mu_at = [&](double const p_conf)
    {
        auto out = computeVanDerWaalsMicroPotential(
            st.n_l, st.rho_lR, st.active_nS, st.rho_SR, st.hamaker, st.Sa,
            st.sign, 0.0, 0.0, 0.0, 0.0);
        PotentialExchangeLocalSolveContext ctx;
        ctx.phi = 0.4;
        ctx.volumetric_strain = -0.01;  // compression, fixed
        ctx.volumetric_strain_prev = 0.0;
        ctx.confining_pressure_p_conf = p_conf;
        ctx.biot_coefficient = 1.0;
        applyFilmPressureMicroPotential(out, st.n_l, st.rho_lR, st.active_nS,
                                        ctx, params);
        return out.mu_lR;
    };

    double const mu_low = mu_at(1.0e5);
    double const mu_high = mu_at(1.0e6);
    EXPECT_GT(mu_high, mu_low);
}

// Anchor: approved baseline — with the strain coupling Off, the fold point
// must follow the existing (shipped) path: the integrable partner is active
// and the result differs from the bare law only by that partner; with the
// coupling ON the shipped partner must NOT also be applied (no double
// counting; replacement is exclusive). Verified structurally: Off and
// Kinematic disagree at finite strain (different mechanisms), while at zero
// strain and zero load Kinematic equals the bare law exactly but Off equals
// bare law + (zero) partner = bare law as well.
TEST(RichardsMechanicsStrainedFilm, ReplacementIsExclusiveAtZeroStrain)
{
    StrainedFilmSampleState st;
    PotentialExchangeParameters params;
    params.enabled = true;
    params.hamaker_constant = st.hamaker;
    params.specific_surface = st.Sa;
    params.micro_solid_density_reference = st.rho_SR;
    params.micro_solid_volume_fraction_reference = st.active_nS;
    params.micro_potential_convention =
        MicroPotentialConvention::NegativeAttractive;
    params.film_pressure_coupling = true;

    auto const mu_for = [&](FilmStrainCouplingMode const mode,
                            double const eps_v, double const p_conf)
    {
        auto out = computeVanDerWaalsMicroPotential(
            st.n_l, st.rho_lR, st.active_nS, st.rho_SR, st.hamaker, st.Sa,
            st.sign, 0.0, 0.0, 0.0, 0.0);
        PotentialExchangeLocalSolveContext ctx;
        ctx.phi = 0.4;
        ctx.volumetric_strain = eps_v;
        ctx.volumetric_strain_prev = 0.0;
        ctx.confining_pressure_p_conf = p_conf;
        ctx.biot_coefficient = 1.0;
        ctx.drained_bulk_modulus = 1.0e9;  // structural sample stiffness
        params.film_strain_coupling = mode;
        applyFilmPressureMicroPotential(out, st.n_l, st.rho_lR, st.active_nS,
                                        ctx, params);
        return out.mu_lR;
    };

    double const mu_bare = bareMu(st, st.n_l);

    // Zero strain, zero load: both reduce to the bare law (the Off path's
    // partner vanishes at eps_v = 0 and p_conf = 0; the strained path's
    // evaluation point and load term are identities there).
    EXPECT_DOUBLE_EQ(mu_for(FilmStrainCouplingMode::Off, 0.0, 0.0), mu_bare);
    EXPECT_DOUBLE_EQ(mu_for(FilmStrainCouplingMode::Kinematic, 0.0, 0.0),
                     mu_bare);

    // Finite strain + load: the two mechanisms must DIFFER (if the shipped
    // partner were still added on top of the strained law, the strained value
    // would carry both and this distinction would collapse).
    double const eps_v = -0.02;
    double const p_conf = 1.0e6;
    EXPECT_NE(mu_for(FilmStrainCouplingMode::Off, eps_v, p_conf),
              mu_for(FilmStrainCouplingMode::Kinematic, eps_v, p_conf));
}

// ── Live K(rho_d) helper (K_OF_RHO_D_LIVE.md; Vinay 2026-06-10) ────────────
// Physics anchor (CLAUDE.md §3a): analytical limit / derived identity — the
// helper must reproduce the piecewise-linear table exactly at its knots and
// hold the endpoint values outside the range; off-mode must return the
// parse-time scalar bit-for-bit. The knot values below are STRUCTURAL
// in-test constants (not physical material parameters); expected values are
// derived in-file from the linear-interpolation identity.
namespace
{
PotentialExchangeParameters liveKSampleParams()
{
    PotentialExchangeParameters params;
    params.micro_solid_density_reference = 2650.0;  // rho_SR, mirrors the
        // prior approved sample state above (kg/m^3).
    params.potential_augmentation_prefactor = 7.0;  // structural scalar K
    // Structural knots: K(1000) = 10, K(2000) = 30 (J/kg vs kg/m^3).
    params.potential_augmentation_prefactor_vs_dry_density =
        std::make_shared<MathLib::PiecewiseLinearInterpolation const>(
            std::vector<double>{1000.0, 2000.0},
            std::vector<double>{10.0, 30.0});
    return params;
}

// phi such that rho_d = rho_SR*(1-phi) equals the requested dry density.
double phiForDryDensity(PotentialExchangeParameters const& params,
                        double const rho_d)
{
    return 1.0 - rho_d / params.micro_solid_density_reference;
}
}  // namespace

TEST(RichardsMechanicsLiveKOfRhoD, OffModeReturnsScalar)
{
    auto params = liveKSampleParams();
    params.potential_augmentation_prefactor_live_dry_density = false;
    // Off mode: scalar, regardless of table presence and finite phi.
    EXPECT_DOUBLE_EQ(7.0, effectiveAugmentationPrefactor(
                              params, phiForDryDensity(params, 1500.0)));

    // No table: scalar even when the mode flag is on.
    auto params_no_table = liveKSampleParams();
    params_no_table.potential_augmentation_prefactor_live_dry_density = true;
    params_no_table.potential_augmentation_prefactor_vs_dry_density = nullptr;
    EXPECT_DOUBLE_EQ(7.0,
                     effectiveAugmentationPrefactor(
                         params_no_table, phiForDryDensity(params, 1500.0)));
}

TEST(RichardsMechanicsLiveKOfRhoD, LiveModeEvaluatesTable)
{
    auto params = liveKSampleParams();
    params.potential_augmentation_prefactor_live_dry_density = true;

    // At the knots: exact knot values.
    EXPECT_DOUBLE_EQ(10.0, effectiveAugmentationPrefactor(
                               params, phiForDryDensity(params, 1000.0)));
    EXPECT_DOUBLE_EQ(30.0, effectiveAugmentationPrefactor(
                               params, phiForDryDensity(params, 2000.0)));

    // Interior points (two distinct phis): linear-interpolation identity
    // K(rho_d) = 10 + 20*(rho_d - 1000)/1000, derived in-file.
    double const rho_d_1 = 1500.0;
    double const rho_d_2 = 1750.0;
    double const expected_1 = 10.0 + 20.0 * (rho_d_1 - 1000.0) / 1000.0;
    double const expected_2 = 10.0 + 20.0 * (rho_d_2 - 1000.0) / 1000.0;
    EXPECT_DOUBLE_EQ(expected_1,
                     effectiveAugmentationPrefactor(
                         params, phiForDryDensity(params, rho_d_1)));
    EXPECT_DOUBLE_EQ(expected_2,
                     effectiveAugmentationPrefactor(
                         params, phiForDryDensity(params, rho_d_2)));

    // Non-finite phi (the context sentinel): fall back to the scalar.
    EXPECT_DOUBLE_EQ(7.0,
                     effectiveAugmentationPrefactor(
                         params, std::numeric_limits<double>::infinity()));
    EXPECT_DOUBLE_EQ(7.0,
                     effectiveAugmentationPrefactor(
                         params, std::numeric_limits<double>::quiet_NaN()));
}

TEST(RichardsMechanicsLiveKOfRhoD, ClampsAtTableRangeEnds)
{
    auto params = liveKSampleParams();
    params.potential_augmentation_prefactor_live_dry_density = true;

    // Below rho_d_min and above rho_d_max the endpoint values are held
    // (PiecewiseLinearInterpolation::getValue endpoint hold).
    EXPECT_DOUBLE_EQ(10.0, effectiveAugmentationPrefactor(
                               params, phiForDryDensity(params, 500.0)));
    EXPECT_DOUBLE_EQ(30.0, effectiveAugmentationPrefactor(
                               params, phiForDryDensity(params, 2600.0)));
}
