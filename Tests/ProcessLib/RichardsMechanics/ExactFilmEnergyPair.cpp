// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause
//
// EXACT one-Psi strained-film energy pair (film_energy_route = exact) —
// unit tests. Design: ProcessLib/RichardsMechanics/DSM/
// PI_OF_NL_EV_IMPLEMENTATION.md (tests T-2..T-5, T-7; T-1 is run-level,
// T-6/T-8 blocked on decisions Q3/Q4 -> GTEST_SKIP).
//
// Physics anchors (CLAUDE.md §3): analytical limits (zero strain, kappa->0
// reduction to the shipped integrable partner), derived identities
// (FD-vs-analytic chains; the Maxwell cross identity), conservation law
// (loop integral of a gradient field vanishes). No fitted expected values;
// tolerances derive from the FD step / quadrature order in this file.
//
// Sample-state parameter values mirror the prior approved unit tests in
// StrainedFilmPotential.cpp / DSMMicroMacroSingleIntegrationPoint.cpp
// (hamaker 6.0e-20 J, Sa 1000 m^2/kg, rho_SR 2650 kg/m^3, nS 0.6,
// rho_lR 1100 kg/m^3, biot 1.0, K_drained 1.5e8 Pa) — citation source:
// prior user-approved test code in this repository.

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

#include "ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h"

using namespace ProcessLib::RichardsMechanics;

namespace
{
struct ExactPairSampleState
{
    double n_l = 0.3;
    double rho_lR = 1100.0;  // confined micro-liquid density scale
    double active_nS = 0.6;
    double rho_SR = 2650.0;
    double hamaker = 6.0e-20;
    double Sa = 1000.0;
    // NegativeAttractive convention (the MS33 PRJ family): mu_lR < 0,
    // Pi = -rho*mu_lR > 0.
    double sign = -1.0;
    double K_aug = 0.0;     // augmentation off by default; tests switch it on
    double lambda_aug = 0.0;
    double floor = 0.0;
    double biot = 1.0;          // prior approved test value
    double K_drained = 1.5e8;   // Pa, prior approved test value
};

StrainedFilmEnergyPairData pairAt(ExactPairSampleState const& st,
                                  double const n_l, double const eps_v,
                                  double const kappa,
                                  bool const include_S = true)
{
    return computeStrainedFilmEnergyPair(
        n_l, eps_v, kappa, st.biot, st.K_drained, include_S, st.rho_lR,
        st.active_nS, st.rho_SR, st.hamaker, st.Sa, st.sign, st.K_aug,
        st.lambda_aug, 0.0, st.floor);
}

VanDerWaalsMicroPotentialData bare(ExactPairSampleState const& st,
                                   double const w)
{
    return computeVanDerWaalsMicroPotential(w, st.rho_lR, st.active_nS,
                                            st.rho_SR, st.hamaker, st.Sa,
                                            st.sign, st.K_aug, st.lambda_aug,
                                            0.0, st.floor);
}

// Augmentation parameters for the with-aug variants: K_aug magnitude of the
// Model-I fitted family order (cited: prior approved calibration record,
// K = 7654.9 J/kg at dd1400, memory ogs_rm_dsm_potential_physics); lambda
// from h-scale: h0 = n_l/(nS*rho_SR*Sa) ~ 1.9e-7 m at the sample state, take
// lambda = h0 so xi0 ~ 1 (structural probe placement, not a material claim).
ExactPairSampleState withAug()
{
    ExactPairSampleState st;
    st.K_aug = 7654.9;  // J/kg
    st.lambda_aug = st.n_l / (st.active_nS * st.rho_SR * st.Sa);  // m
    return st;
}
}  // namespace

// Anchor: approved baseline — default route must stay Operational so every
// existing PRJ is bit-for-bit unaffected (T-1 run-level counterpart).
TEST(RichardsMechanicsExactFilmPair, DefaultRouteIsOperational)
{
    PotentialExchangeParameters params;
    EXPECT_EQ(params.film_energy_route, FilmEnergyRoute::Operational);
}

// T-7 (helper level): the §3 mode-matrix predicate.
TEST(RichardsMechanicsExactFilmPair, FilmEnergyRouteCombinationMatrix)
{
    using M = FilmStrainCouplingMode;
    using R = FilmEnergyRoute;
    EXPECT_TRUE(isValidFilmEnergyRouteCombination(M::Off, R::Operational));
    EXPECT_TRUE(isValidFilmEnergyRouteCombination(M::Kinematic, R::Operational));
    EXPECT_TRUE(isValidFilmEnergyRouteCombination(M::Equilibrium, R::Operational));
    EXPECT_TRUE(isValidFilmEnergyRouteCombination(M::Kinematic, R::Exact));
    EXPECT_FALSE(isValidFilmEnergyRouteCombination(M::Off, R::Exact));
    EXPECT_FALSE(isValidFilmEnergyRouteCombination(M::Equilibrium, R::Exact));
}

// T-2. Anchor: analytical limit — at eps_v = 0 the strain integrals vanish
// identically (algebraic identity, no FD): mu_mech = 0, Psi_film = 0, and the
// eigenstress half reduces to the unstrained -nS*n_l*Pi(n_l).
TEST(RichardsMechanicsExactFilmPair, ZeroStrainReduction)
{
    for (auto const& st : {ExactPairSampleState{}, withAug()})
    {
        for (double const kappa : {st.active_nS, 1.0})
        {
            auto const p = pairAt(st, st.n_l, 0.0, kappa);
            EXPECT_DOUBLE_EQ(0.0, p.mu_mech);
            EXPECT_DOUBLE_EQ(0.0, p.Psi_film);
            double const Pi0 = -st.rho_lR * bare(st, st.n_l).mu_lR;  // Pa
            EXPECT_NEAR(-st.active_nS * st.n_l * Pi0, p.sigma_sw_m,
                        1e-12 * std::abs(Pi0));  // tol: algebraic identity at
                                                 // the Pi scale
        }
    }
}

// T-3. Anchor: derived identity — FD-vs-analytic for all derivative blocks
// AND the Maxwell cross identity dsigma_sw/dn_l == nS*rho_lR*dmu_mech/deps_v
// (exact for the one-Psi pair). FD tolerance derived from the central-
// difference step: rel error O(d^2 * f'''/f') -> use 5e-5 relative at
// d = 1e-7 * scale, plus an absolute floor at machine-eps * value scale.
TEST(RichardsMechanicsExactFilmPair, MaxwellIdentityAndFDChains)
{
    for (auto const& st : {ExactPairSampleState{}, withAug()})
    {
        for (double const kappa : {st.active_nS, 1.0})
        {
            for (double const eps_v : {-0.03, -0.005, 0.01})
            {
                auto const p = pairAt(st, st.n_l, eps_v, kappa);

                // Maxwell cross identity (analytic vs analytic — the pair is
                // a gradient by construction; tolerance at rounding scale).
                double const lhs = p.dsigma_sw_dnl;  // Pa per n_l
                double const rhs = st.active_nS * st.rho_lR *
                                   p.dmu_mech_deps_v;  // Pa per n_l
                EXPECT_NEAR(lhs, rhs, 1e-9 * std::max(std::abs(lhs), 1.0))
                    << "kappa=" << kappa << " eps_v=" << eps_v;

                // FD: d(mu_mech)/d(eps_v).
                double const de = 1e-7;
                double const fd_eps =
                    (pairAt(st, st.n_l, eps_v + de, kappa).mu_mech -
                     pairAt(st, st.n_l, eps_v - de, kappa).mu_mech) /
                    (2 * de);
                EXPECT_NEAR(p.dmu_mech_deps_v, fd_eps,
                            5e-5 * std::abs(fd_eps) + 1e-10);

                // FD: d(mu_mech)/d(n_l). Absolute floor = FD-numerator
                // roundoff eps_mach*|f|/(2*dn) with margin (the derivative can
                // nearly cancel while |mu_mech| stays large — the floor must
                // scale with the FUNCTION value, not the derivative).
                double const dn = 1e-7 * st.n_l;
                double const fd_nl =
                    (pairAt(st, st.n_l + dn, eps_v, kappa).mu_mech -
                     pairAt(st, st.n_l - dn, eps_v, kappa).mu_mech) /
                    (2 * dn);
                double const fd_nl_floor =
                    1e-15 * std::abs(p.mu_mech) / dn;  // roundoff floor
                EXPECT_NEAR(p.dmu_mech_dnl, fd_nl,
                            5e-5 * std::abs(fd_nl) + fd_nl_floor + 1e-10);

                // FD: d(sigma_sw)/d(n_l).
                double const fd_sig =
                    (pairAt(st, st.n_l + dn, eps_v, kappa).sigma_sw_m -
                     pairAt(st, st.n_l - dn, eps_v, kappa).sigma_sw_m) /
                    (2 * dn);
                EXPECT_NEAR(p.dsigma_sw_dnl, fd_sig,
                            5e-5 * std::abs(fd_sig) + 1e-6);

                // FD: d(mu_mech)/d(rho_lR).
                double const dr = 1e-7 * st.rho_lR;
                auto stp = st;
                stp.rho_lR = st.rho_lR + dr;
                auto stm = st;
                stm.rho_lR = st.rho_lR - dr;
                double const fd_rho =
                    (pairAt(stp, st.n_l, eps_v, kappa).mu_mech -
                     pairAt(stm, st.n_l, eps_v, kappa).mu_mech) /
                    (2 * dr);
                EXPECT_NEAR(p.dmu_mech_drho_lR, fd_rho,
                            5e-5 * std::abs(fd_rho) + 1e-12);
            }
        }
    }
}

// T-4. Anchor: analytical limit — kappa -> 0 reduces the exact pair EXACTLY
// to the shipped integrable partner (frozen-h limit). Series remainder is
// O(kappa*eps_v) relative -> assert |diff| <= 10*kappa*|value| + abs floor
// (C = 10 bounds the series coefficients 1.5, 2, xi0-products at this state).
// Discriminating: the OPERATIONAL form does NOT have this limit.
TEST(RichardsMechanicsExactFilmPair, FrozenHLimitMatchesShippedPartner)
{
    double const kappa = 1e-9;  // sub-threshold: series branch active
    for (auto const& st : {ExactPairSampleState{}, withAug()})
    {
        for (double const eps_v : {-0.03, 0.01})
        {
            auto const b0 = bare(st, st.n_l);
            double const Pi = -st.rho_lR * b0.mu_lR;             // Pa
            double const dPi = -st.rho_lR * b0.dmu_lR_dnl;       // Pa per n_l
            double const d2Pi = -st.rho_lR * b0.d2mu_lR_dnl2;    // Pa per n_l^2
            auto const shipped = computeIntegrableMechanicalMicroPotential(
                Pi, dPi, d2Pi, st.n_l, eps_v, st.biot, st.K_drained,
                st.rho_lR);
            auto const p = pairAt(st, st.n_l, eps_v, kappa);

            double const tol_rel = 10.0 * kappa;
            EXPECT_NEAR(shipped.mu_lR_mech, p.mu_mech,
                        tol_rel * std::abs(shipped.mu_lR_mech) + 1e-12);
            EXPECT_NEAR(shipped.dmu_lR_mech_dnl, p.dmu_mech_dnl,
                        tol_rel * std::abs(shipped.dmu_lR_mech_dnl) + 1e-9);
            EXPECT_NEAR(shipped.dmu_lR_mech_deps_v, p.dmu_mech_deps_v,
                        tol_rel * std::abs(shipped.dmu_lR_mech_deps_v) + 1e-9);
        }
    }
}

// T-5. Anchor: conservation law — for the EXACT pair the work integral around
// a closed (eps_v, n_l) loop vanishes (gradient field); the OPERATIONAL cut
// has a Maxwell defect O(Pi*eps_v) and must NOT vanish (predicted in
// STRAINED_FILM_IMPLEMENTATION.md §9a — this test measures it; report per
// CLAUDE.md §5.1). Trapezoid quadrature is O(N^-2): the exact loop residual
// must fall below rel_bound = 50/N^2 of the path work scale (50 bounds the
// curvature-to-scale ratio on this path, checked by the N-halving assert);
// the operational defect must sit ABOVE the same bound by x100.
// Loop ranges are structural probes around the sample state (CLAUDE.md §1.2),
// matching the strain magnitudes already used in the approved tests.
TEST(RichardsMechanicsExactFilmPair, ReversibilityLoopClosesExactOnly)
{
    auto const st = withAug();
    double const kappa = st.active_nS;
    double const e0 = -0.03, e1 = 0.0;     // eps_v range (compression leg)
    double const n0 = 0.27, n1 = 0.33;     // n_l range around the sample state

    auto const loop = [&](int const N, bool const exact_route)
    {
        // mu half [J/kg] and sigma half [Pa] of the tested route on the
        // drained line p_conf = -K_drained*eps_v.
        auto const mu_half = [&](double const n_l, double const eps_v)
        {
            if (exact_route)
            {
                return pairAt(st, n_l, eps_v, kappa).mu_mech;  // J/kg
            }
            // operational mech share: bare(w_eff) - bare(n_l) + b*p_conf/rho
            double const w_eff = n_l * std::max(1e-6, 1.0 + kappa * eps_v);
            double const p_conf = -st.K_drained * eps_v;  // Pa, drained line
            return bare(st, w_eff).mu_lR - bare(st, n_l).mu_lR +
                   st.biot * p_conf / st.rho_lR;  // J/kg
        };
        auto const sigma_half = [&](double const n_l, double const eps_v)
        {
            // identical in both routes (Pi at w_eff, transmitted load).
            return pairAt(st, n_l, eps_v, kappa).sigma_sw_m;  // Pa
        };

        // dW = sigma_sw deps_v + nS*rho_lR*mu dn_l, trapezoid on 4 edges.
        double W = 0.0, Wabs = 0.0;
        auto const seg = [&](double const xa, double const ya, double const xb,
                             double const yb)
        {
            // x = eps_v, y = n_l
            for (int i = 0; i < N; ++i)
            {
                double const t0 = static_cast<double>(i) / N;
                double const t1 = static_cast<double>(i + 1) / N;
                double const xA = xa + (xb - xa) * t0, xB = xa + (xb - xa) * t1;
                double const yA = ya + (yb - ya) * t0, yB = ya + (yb - ya) * t1;
                double const dx = xB - xA, dy = yB - yA;
                double const fA = sigma_half(yA, xA) * dx +
                                  st.active_nS * st.rho_lR * mu_half(yA, xA) *
                                      dy;  // J/m^3
                double const fB = sigma_half(yB, xB) * dx +
                                  st.active_nS * st.rho_lR * mu_half(yB, xB) *
                                      dy;  // J/m^3
                W += 0.5 * (fA + fB);
                Wabs += std::abs(0.5 * (fA + fB));
            }
        };
        seg(e0, n0, e1, n0);  // strain up at n0
        seg(e1, n0, e1, n1);  // wet at e1
        seg(e1, n1, e0, n1);  // strain down at n1
        seg(e0, n1, e0, n0);  // dry at e0
        return std::pair{W, Wabs};
    };

    int const N = 400;
    double const rel_bound = 50.0 / (static_cast<double>(N) * N);

    auto const [W_exact, Wabs_exact] = loop(N, true);
    auto const [W_exact2, Wabs_exact2] = loop(2 * N, true);
    auto const [W_op, Wabs_op] = loop(N, false);

    std::cout << "[loop] measured: |W_exact|/scale = "
              << std::abs(W_exact) / Wabs_exact
              << ", |W_operational|/scale = " << std::abs(W_op) / Wabs_op
              << " (N=" << N << ")\n";

    // Exact: below the quadrature bound, and converging (halving check).
    EXPECT_LE(std::abs(W_exact), rel_bound * Wabs_exact);
    EXPECT_LE(std::abs(W_exact2), std::abs(W_exact) + 1e-30);
    // Operational: the §9a defect sits far above the quadrature bound.
    EXPECT_GE(std::abs(W_op), 100.0 * rel_bound * Wabs_op);
}

// T-6 — blocked on decisions Q3 (mass-derivative freeze) and Q4 (K_liq
// value+source); see PI_OF_NL_EV_IMPLEMENTATION.md §8.
TEST(RichardsMechanicsExactFilmPair, LiquidCarrierEnergyPressureConsistency)
{
    GTEST_SKIP() << "TODO(Vinay): blocked on Q3 (fixed-volume vs kinematic "
                    "mass-derivative) and Q4 (K_liq value + source); "
                    "PI_OF_NL_EV_IMPLEMENTATION.md §2.2/§8.";
}

// T-8 — run-level expulsion probe (drained oedometer ramp past the
// crossover); magnitudes TODO(Vinay).
TEST(RichardsMechanicsExactFilmPair, ExpulsionProbeDrainedRamp)
{
    GTEST_SKIP() << "TODO(Vinay): run-level probe (T-8); expected magnitudes "
                    "TODO; PI_OF_NL_EV_IMPLEMENTATION.md §5.";
}
