// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// Unit tests of the v4 switches (branch dsm_mass_conservation_v4_tm_krel_
// 2026-10-02; Vinay's ruling 2026-10-02 "(go with L + drop T_m) x (1a, 1b
// separate)"; design ~/ogs-models/scratch/2026-10-02_kkt_v4_tm_krel/
// DESIGN_V4.md section 3.1). NOT adopted. SUPPLEMENTS only: no existing test
// is edited.
//
// UT-V4-1 Kirchhoff table, UT-V4-2 element mobility with ties, UT-V4-3 the
// fold of the Gauss-point flux against the monotone Kirchhoff flux, UT-V4-4
// drop T_m (books expression, tangent struct, Active identity, interior
// dn_l/dp_L and dn_l/deps_v against central differences of the local KKT
// solve), UT-V4-5 gate truth tables. UT-V4-6 (parser) is run as v4Deck tests
// (REC/tests/parser/), not here.
//
// Every test prints its measured numbers (C7 convention of MicroCeilingKkt).
// Deck values come from the AB decks (DESIGN_V4.md 1.3, READ): SaturationTuller
// (residual 0, max 1, area factor 1, shape factor 0.8584073464102069, pore size
// 1e-5 m, surface tension 0.0715 N/m, pressure tolerance 0, cavitation 1.4e8 Pa)
// and RelativePermeabilityGeneralizedPower (residuals 0, min 1e-2, a = 1,
// exponent 3); the live K(rho_d) table of the AB decks (900 1400 1600 1800 /
// 24241.938995 46000.3269694 104698.192423 265909.813902), exponent 7.5e-7 and
// floor 7.06e-5. Test-only proposals (need Vinay's approval, as in
// MicroCeilingKkt.cpp): kV4FdSafetyFactor = 10, kV4SampleSafety = 1.5 (sampled
// curvature bound), the 1e-10 local-solve floor used for the FD tolerance of
// UT-V4-4d (the solver's own literal, RichardsMechanicsFEM-impl.h), the
// IC-like states of UT-V4-4 (phi 0.40, n_l 0.05, p_L -1 MPa, S_L 0.3).

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <vector>

#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.h"
#include "MaterialLib/MPL/Properties/RelativePermeability/RelPermGeneralizedPower.h"
#include "MaterialLib/MPL/VariableType.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/RichardsMechanics/KirchhoffMobility.h"
#include "ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h"

using namespace ProcessLib::RichardsMechanics;

namespace
{
constexpr double v4_eps_mach = std::numeric_limits<double>::epsilon();
constexpr double kV4FdSafetyFactor = 10.0;
constexpr double kV4SampleSafety = 1.5;

// ── Deck pair (AB decks) ─────────────────────────────────────────────────────
constexpr double kV4Shape = 0.8584073464102069;
constexpr double kV4Sigma = 0.0715;
constexpr double kV4Pore = 1e-5;
constexpr double kV4Cav = 1.4e8;
// c = 4 F sigma^2 / (A L^2) (SaturationTuller.h, C_T).
double const kV4C = 4.0 * kV4Shape * kV4Sigma * kV4Sigma / (1.0 * kV4Pore * kV4Pore);

struct V4DeckPair
{
    MaterialPropertyLib::SaturationTuller sat{
        "saturation", 0.0, 1.0, 1.0, kV4Shape, kV4Pore, kV4Sigma, 0.0, kV4Cav};
    MaterialPropertyLib::RelPermGeneralizedPower kr{
        "relative_permeability", 0.0, 0.0, 1e-2, 1.0, 3.0};
    ParameterLib::SpatialPosition pos;
    double S(double const p_c) const
    {
        MaterialPropertyLib::VariableArray v;
        v.capillary_pressure = p_c;
        return std::get<double>(sat.value(v, pos, 0.0, 0.0));
    }
    double k(double const p_c) const
    {
        MaterialPropertyLib::VariableArray v;
        v.capillary_pressure = p_c;
        v.liquid_saturation = S(p_c);
        return std::get<double>(kr.value(v, pos, 0.0, 0.0));
    }
};

V4DeckPair const& v4Deck()
{
    static V4DeckPair const d;
    return d;
}

KirchhoffMobilityTable const& v4Table2048()
{
    static KirchhoffMobilityTable const t = KirchhoffMobilityTable::build(
        [](double const p) { return v4Deck().k(p); }, 2048);
    return t;
}

// Closed-form Kirchhoff potential of the v4Deck pair (in-file derivation, as
// REC/design_facts/kbar_static.py): int_0^x exp(-m c/s^2) ds = x exp(-m c/x^2)
// - sqrt(pi m c) erfc(sqrt(m c)/x); (1 - E)^3 = 1 - 3E + 3E^2 - E^3, so on
// (0, p_k] Phi = x - 3G1 + 3G2 - G3, slope 1e-2 beyond p_k, Phi = x for x <= 0.
double v4G(double const x, double const m)
{
    if (x <= 0.0)
    {
        return 0.0;
    }
    double const a = m * kV4C;
    return x * std::exp(-a / (x * x)) -
           std::sqrt(M_PI * a) * std::erfc(std::sqrt(a) / x);
}
double v4PKinkClosedForm()
{
    return std::sqrt(kV4C / (-std::log(1.0 - std::cbrt(1e-2))));
}
double v4PhiClosedForm(double const x)
{
    if (x <= 0.0)
    {
        return x;
    }
    double const pk = v4PKinkClosedForm();
    double const xs = std::min(x, pk);
    double v = xs - 3.0 * v4G(xs, 1) + 3.0 * v4G(xs, 2) - v4G(xs, 3);
    if (x > pk)
    {
        v += 1e-2 * (x - pk);
    }
    return v;
}
// d^2 k / d p_c^2 of k = S^3, S = 1 - E, E = exp(-c/p^2) (above the floor).
double v4KSecondDerivative(double const p)
{
    double const E = std::exp(-kV4C / (p * p));
    double const S = 1.0 - E;
    double const dS = -E * 2.0 * kV4C / (p * p * p);
    double const d2S =
        -E * (4.0 * kV4C * kV4C / std::pow(p, 6) - 6.0 * kV4C / std::pow(p, 4));
    return 6.0 * S * dS * dS + 3.0 * S * S * d2S;
}
}  // namespace

// ── UT-V4-1 ─────────────────────────────────────────────────────────────────
TEST(RichardsMechanics, DSMv4KirchhoffTable)
{
    auto const& t = v4Table2048();
    auto const& d = v4Deck();
    // (a) p_k against the closed form, p_sat = the largest double with S == 1.
    double const pk_cf = v4PKinkClosedForm();
    double const ulps_pk = std::abs(t.pK() - pk_cf) /
                           (std::nextafter(pk_cf, 2 * pk_cf) - pk_cf);
    std::cout << "UT-V4-1a p_sat=" << std::setprecision(17) << t.pSat()
              << " p_k=" << t.pK() << " closed form p_k=" << pk_cf
              << " |diff| in ulps=" << ulps_pk << " cells=" << t.numberOfCells()
              << "\n";
    EXPECT_EQ(d.S(t.pSat()), 1.0);
    EXPECT_LT(d.S(std::nextafter(t.pSat(), 1e300)), 1.0);
    EXPECT_EQ(d.k(t.pK()), 1e-2);
    EXPECT_GT(d.k(std::nextafter(t.pK(), 0.0)), 1e-2);
    EXPECT_LE(ulps_pk, 1.0);
    EXPECT_EQ(t.numberOfCells(), 2241u);  // DESIGN_V4.md 1.6 (MEASURED there)

    // (b) k_h and Phi_h against the closed forms; bound = per cell
    // (h_j^2/8) max|k''| (linear interpolation error), max|k''| sampled at 33
    // points per cell times kV4SampleSafety.
    auto const& g = t.nodes();
    std::vector<double> M(t.numberOfCells());
    for (std::size_t j = 0; j < M.size(); ++j)
    {
        double m = 0.0;
        for (int q = 0; q <= 32; ++q)
        {
            m = std::max(m, std::abs(v4KSecondDerivative(
                                g[j] + (g[j + 1] - g[j]) * q / 32.0)));
        }
        M[j] = kV4SampleSafety * m;
    }
    std::mt19937_64 rng(20261002);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    double max_rel_k = 0.0;
    double max_ratio_k = 0.0;
    for (int i = 0; i < 100000; ++i)
    {
        double const x = t.pSat() * std::pow(t.pK() / t.pSat(), u01(rng));
        std::size_t const j = std::min(t.cellOf(x), M.size() - 1);
        double const h = g[j + 1] - g[j];
        double const err = std::abs(t.value(x) - d.k(x));
        double const bound = h * h / 8.0 * M[j] + 4 * v4_eps_mach;
        max_rel_k = std::max(max_rel_k, err / d.k(x));
        max_ratio_k = std::max(max_ratio_k, err / bound);
        EXPECT_LE(err, bound) << "x=" << x;
    }
    double max_rel_Phi = 0.0;
    double cum_bound = 0.0;
    double const Phi0 = v4PhiClosedForm(t.pSat());
    for (std::size_t i = 1; i < g.size(); ++i)
    {
        double const h = g[i] - g[i - 1];
        cum_bound += h * h * h / 12.0 * M[i - 1];
        double const exact = v4PhiClosedForm(g[i]) - Phi0;
        double const err = std::abs(t.potential(g[i]) - exact);
        max_rel_Phi = std::max(max_rel_Phi, err / exact);
        EXPECT_LE(err, cum_bound + 1e-12 * exact) << "node " << i;
    }
    std::cout << "UT-V4-1b max|k_h-k|/k=" << max_rel_k
              << " (design table_resolution.py: 5.06e-06) max err/bound="
              << max_ratio_k << " max|Phi_h-Phi|/Phi at nodes=" << max_rel_Phi
              << " (design: 1.93e-07)\n";

    // (c) exact constant branches.
    for (double const a : {t.pK(), 3e4, 1e6, 1e8, 2e8})
    {
        EXPECT_EQ(t.mean(a, a * 1.5).k_bar, 1e-2);
        EXPECT_EQ(t.mean(a, a).k_bar, 1e-2);
    }
    for (double const b : {t.pSat(), 2000.0, 0.0, -1e5})
    {
        EXPECT_EQ(t.mean(b - 1e3, b).k_bar, 1.0);
        EXPECT_EQ(t.mean(b, b).k_bar, 1.0);
    }

    // (d) convex bound on random [a, b] (k_h monotone: extremes at the ends).
    std::uniform_real_distribution<double> uab(-5e3, 6e4);
    for (int i = 0; i < 20000; ++i)
    {
        double a = uab(rng);
        double b = uab(rng);
        if (a > b)
        {
            std::swap(a, b);
        }
        // Convex combination of the node values up to the round-off of
        // Q/len (a few ulps; MEASURED 1.1e-16 below k_h(b) in run 1).
        double const kb = t.mean(a, b).k_bar;
        EXPECT_LE(kb, t.value(a) + 16 * v4_eps_mach * t.value(a));
        EXPECT_GE(kb, t.value(b) - 16 * v4_eps_mach * t.value(a));
    }

    // (e) d kbar/da, d kbar/db against central differences of the DISCRETE
    // kbar; the step stays below the distance to the nearest grid node or
    // branch point (k_h linear across the stencil); tolerance = kV4FdSafetyFactor
    // * (round-off eps*|kbar|/h + Richardson estimate |D(h) - D(h/2)|).
    auto const dist_to_node = [&](double const x)
    {
        double dmin = std::numeric_limits<double>::infinity();
        auto const it = std::lower_bound(g.begin(), g.end(), x);
        if (it != g.end())
        {
            dmin = std::min(dmin, std::abs(*it - x));
        }
        if (it != g.begin())
        {
            dmin = std::min(dmin, std::abs(*(it - 1) - x));
        }
        return dmin;
    };
    struct Pair
    {
        double a, b;
        char const* name;
    };
    double const c0 = g[1000] + 0.3 * (g[1001] - g[1000]);
    std::vector<Pair> const cases = {
        {c0, g[1000] + 0.7 * (g[1001] - g[1000]), "same cell"},
        {c0, g[1001] + 0.4 * (g[1002] - g[1001]), "adjacent cells"},
        {c0, g[1500] + 0.5 * (g[1501] - g[1500]), "many cells"},
        {c0, 5e4, "straddling p_k"},
        {1000.0, c0, "straddling p_sat"},
        {-2e3, 4e4, "straddling both"},
    };
    double worst = 0.0;
    for (auto const& cs : cases)
    {
        for (int which = 0; which < 2; ++which)
        {
            double const x = which == 0 ? cs.a : cs.b;
            double const h0 = std::min(std::cbrt(v4_eps_mach) *
                                           std::max({std::abs(cs.a),
                                                     std::abs(cs.b), 1.0}),
                                       0.25 * dist_to_node(x));
            auto const D = [&](double const h)
            {
                double const ap = which == 0 ? cs.a + h : cs.a;
                double const am = which == 0 ? cs.a - h : cs.a;
                double const bp = which == 1 ? cs.b + h : cs.b;
                double const bm = which == 1 ? cs.b - h : cs.b;
                return (t.mean(ap, bp).k_bar - t.mean(am, bm).k_bar) /
                       ((which == 0 ? ap - am : bp - bm));
            };
            double const fd = D(h0);
            double const an = which == 0 ? t.mean(cs.a, cs.b).dk_da
                                         : t.mean(cs.a, cs.b).dk_db;
            double const tol =
                kV4FdSafetyFactor *
                (v4_eps_mach * std::abs(t.mean(cs.a, cs.b).k_bar) / h0 +
                 std::abs(D(h0) - D(0.5 * h0)));
            worst = std::max(worst, std::abs(fd - an) / (tol + 1e-300));
            EXPECT_LE(std::abs(fd - an), tol)
                << cs.name << (which == 0 ? " d/da" : " d/db") << " an=" << an
                << " fd=" << fd;
            std::cout << "UT-V4-1e " << cs.name << (which == 0 ? " d/da" : " d/db")
                      << " analytic=" << an << " fd=" << fd << " tol=" << tol
                      << "\n";
        }
    }
    // Degenerate a == b inside a cell: kbar = k_h(a), both partials = s/2.
    {
        double const a = c0;
        auto const r = t.mean(a, a);
        double const s = t.slopeRight(a);
        EXPECT_EQ(r.k_bar, t.value(a));
        EXPECT_EQ(r.dk_da, 0.5 * s);
        EXPECT_EQ(r.dk_db, 0.5 * s);
        // Same-cell continuity: the general formula at a tiny width gives
        // the same partials to O(width).
        auto const r2 = t.mean(a, a + 1e-6 * (g[1001] - g[1000]));
        EXPECT_NEAR(r2.dk_da, 0.5 * s, 1e-9 * std::abs(s));
        EXPECT_NEAR(r2.dk_db, 0.5 * s, 1e-9 * std::abs(s));
    }
    std::cout << "UT-V4-1e worst |fd-an|/tol=" << worst << "\n";
}

// ── UT-V4-2 ─────────────────────────────────────────────────────────────────
TEST(RichardsMechanics, DSMv4KirchhoffElementMobilityTies)
{
    auto const& t = v4Table2048();
    auto const& g = t.nodes();
    double const x0 = g[800] + 0.25 * (g[801] - g[800]);
    double const x1 = g[1200] + 0.5 * (g[1201] - g[1200]);
    double const xm = g[1000] + 0.5 * (g[1001] - g[1000]);
    auto const mob = [&](std::array<double, 4> const& p_c)
    {
        std::array<double, 4> p_L{};
        for (int n = 0; n < 4; ++n)
        {
            p_L[n] = -p_c[n];
        }
        return kirchhoffElementMobility<4>(t, p_L);
    };
    double const h = 1e-3;  // [Pa], far below the cell widths (> 2 Pa) here
    // Distinct values: central FD per node = g_n.
    {
        std::array<double, 4> const pc = {x0, xm, x1, 0.5 * (xm + x1)};
        auto const m = mob(pc);
        for (int n = 0; n < 4; ++n)
        {
            auto pp = pc;
            auto pm = pc;
            pp[n] -= h;  // p_L + h
            pm[n] += h;  // p_L - h
            double const fd = (mob(pp).k_bar - mob(pm).k_bar) / (2 * h);
            EXPECT_NEAR(fd, m.dk_dpL[n],
                        1e-6 * std::abs(m.dk_dpL[n]) + 1e-18)
                << "node " << n;
            std::cout << "UT-V4-2 distinct node " << n << " g=" << m.dk_dpL[n]
                      << " fd=" << fd << "\n";
        }
    }
    // Two-way tie at the max (nodes 2, 3): uniform shift of the tie set gives
    // the exact directional derivative; the sum of g over the tie = -dkbar/db.
    {
        std::array<double, 4> const pc = {x0, xm, x1, x1};
        auto const m = mob(pc);
        EXPECT_EQ(m.argmax_mask, 0b1100u);
        auto pp = pc;
        auto pm = pc;
        pp[2] += h;
        pp[3] += h;
        pm[2] -= h;
        pm[3] -= h;
        double const fd_shift = (mob(pp).k_bar - mob(pm).k_bar) / (2 * h);
        double const g_sum = -(m.dk_dpL[2] + m.dk_dpL[3]);  // d/dp_c
        EXPECT_NEAR(fd_shift, g_sum, 1e-6 * std::abs(g_sum));
        EXPECT_EQ(m.dk_dpL[2], m.dk_dpL[3]);
        // One-sided FD that keeps the tie order (node 2 moves up alone: it
        // becomes the unique max) = full dkbar/db; Clarke split = half of it.
        auto p1 = pc;
        p1[2] += h;
        double const fd_one = (mob(p1).k_bar - m.k_bar) / h;
        EXPECT_NEAR(fd_one, -2.0 * m.dk_dpL[2], 1e-4 * std::abs(fd_one));
        std::cout << "UT-V4-2 tie-max shift fd=" << fd_shift << " sum g=" << g_sum
                  << " one-sided=" << fd_one << "\n";
    }
    // Two-way tie at the min (nodes 0, 1).
    {
        std::array<double, 4> const pc = {x0, x0, xm, x1};
        auto const m = mob(pc);
        EXPECT_EQ(m.argmin_mask, 0b0011u);
        auto pp = pc;
        auto pm = pc;
        pp[0] += h;
        pp[1] += h;
        pm[0] -= h;
        pm[1] -= h;
        double const fd_shift = (mob(pp).k_bar - mob(pm).k_bar) / (2 * h);
        double const g_sum = -(m.dk_dpL[0] + m.dk_dpL[1]);
        EXPECT_NEAR(fd_shift, g_sum, 1e-6 * std::abs(g_sum));
    }
    // All equal: g_n = -k_h'(a)/4 (right cell), uniform shift exact.
    {
        std::array<double, 4> const pc = {xm, xm, xm, xm};
        auto const m = mob(pc);
        double const s = t.slopeRight(xm);
        for (int n = 0; n < 4; ++n)
        {
            EXPECT_EQ(m.dk_dpL[n], -(0.5 * s / 4 + 0.5 * s / 4));
        }
        std::array<double, 4> pp = pc;
        std::array<double, 4> pm = pc;
        for (int n = 0; n < 4; ++n)
        {
            pp[n] += h;
            pm[n] -= h;
        }
        double const fd = (mob(pp).k_bar - mob(pm).k_bar) / (2 * h);
        EXPECT_NEAR(fd, s, 1e-6 * std::abs(s));
        std::cout << "UT-V4-2 all-equal slope=" << s << " shift fd=" << fd
                  << "\n";
    }
}

// ── UT-V4-3 ─────────────────────────────────────────────────────────────────
// Two-node 1D boundary element, p_c,0 = 0 (Dirichlet row), p_c,1 = P. Gauss
// flux F_GP = P (k(g1 P) + k(g2 P))/2 at the 2-point Gauss points g1,2 = (1 -+
// 1/sqrt 3)/2; Kirchhoff flux F_K = P kbar(0, P) = Phi_h(P) - Phi_h(0).
TEST(RichardsMechanics, DSMv4GaussFoldVersusMonotoneKirchhoffFlux)
{
    auto const& t = v4Table2048();
    auto const& d = v4Deck();
    double const g1 = 0.5 * (1.0 - 1.0 / std::sqrt(3.0));
    double const g2 = 0.5 * (1.0 + 1.0 / std::sqrt(3.0));
    auto const FGP = [&](double const P)
    { return P * 0.5 * (d.k(g1 * P) + d.k(g2 * P)); };
    auto const FK = [&](double const P) { return P * t.mean(0.0, P).k_bar; };
    double const pk = t.pK();
    int sign_changes = 0;
    double prev = 0.0;
    double min_dFK = std::numeric_limits<double>::infinity();
    for (int i = 0; i <= 400; ++i)
    {
        double const P = pk * (0.5 + 4.5 * i / 400.0);
        double const h = 1e-6 * P;
        double const dGP = (FGP(P + h) - FGP(P - h)) / (2 * h);
        double const dK = (FK(P + h) - FK(P - h)) / (2 * h);
        if (i > 0 && (dGP > 0) != (prev > 0))
        {
            ++sign_changes;
            std::cout << "UT-V4-3 Gauss-point dF/dP changes sign near P/p_k="
                      << P / pk << "\n";
        }
        prev = dGP;
        min_dFK = std::min(min_dFK, dK);
        // Structural: dF_K/dP = k_h(P) > 0 (Phi_h' = k_h).
        EXPECT_GT(dK, 0.0) << "P=" << P;
        EXPECT_NEAR(dK, t.value(P), 1e-5 * t.value(P)) << "P=" << P;
    }
    std::cout << "UT-V4-3 Gauss-point sign changes=" << sign_changes
              << " min dF_K/dP=" << min_dFK << "\n";
    EXPECT_GE(sign_changes, 1);  // the fold of investigator B (1 - 2.65 < 0)
}

// ── UT-V4-4 ─────────────────────────────────────────────────────────────────
namespace
{
std::shared_ptr<AugmentationPrefactorTable const> v4ABKTable()
{
    return std::make_shared<AugmentationPrefactorTable const>(
        std::vector<double>{900.0, 1400.0, 1600.0, 1800.0},
        std::vector<double>{24241.938995, 46000.3269694, 104698.192423,
                            265909.813902});
}

// Parameters: the committed block of MicroCeilingKkt.cpp (DeckEos regime:
// rho_lR = 100 + rho_LR), KKT, F3 on (s = -1), reference n_S, with the AB live
// K(rho_d) table, exponent and floor.
PotentialExchangeParameters v4Parameters()
{
    PotentialExchangeParameters p;
    p.enabled = true;
    p.pressure_tolerance = 0.0;
    p.hamaker_constant = 6.0e-20;
    p.specific_surface = 1000.0;
    p.micro_solid_density_reference = 2650.0;
    p.micro_solid_volume_fraction_reference = 0.6;
    p.micro_liquid_density_reference = 100.0;
    p.micro_liquid_density_a = 1e-16;
    p.micro_liquid_density_b = 1.0;
    p.micro_potential_convention = MicroPotentialConvention::NegativeAttractive;
    p.local_nonlinear_solve_mode =
        LocalNonlinearSolveMode::ScalarReferenceMassStorage;
    p.macro_porosity_update_mode = MacroPorosityUpdateMode::ReferenceAdditiveRate;
    p.micro_solid_volume_fraction_mode = MicroSolidVolumeFractionMode::Reference;
    p.initial_micro_water_content = 0.05;
    p.potential_augmentation_prefactor = 104698.192423;
    p.potential_augmentation_prefactor_vs_dry_density = v4ABKTable();
    p.potential_augmentation_prefactor_live_dry_density = true;
    p.potential_augmentation_exponent = 7.5e-7;
    p.micro_water_content_floor = 7.06e-5;
    p.micro_ceiling_treatment = MicroCeilingTreatment::Kkt;
    p.micro_mass_strain_term_eulerian = true;
    p.macro_storage_uses_macro_porosity = true;
    p.macro_balance_drops_micro_biot_term = true;
    return p;
}

struct V4Case
{
    PotentialExchangeParameters params = v4Parameters();
    double dt = 100.0;
    double rho_LR = 1000.0;
    double alpha_bar = 1.0e-9;
    double mu = 1.0e-3;
    double p_L = -1.0e6;
    double n_l_prev = 0.05;
    double phi_prev = 0.40;  // rho_d = 2650*0.6 = 1590, inside [1400, 1600]
    double phi_m_prev = 0.0;
    double phi_M_prev = 0.0;
    double rho_lR_prev = 0.0;
    double rho_l_prev = 0.0;
    double eps_prev = 0.205 - 1.0e-3;
    double d_eps = 1.0e-3;
    double E = 52.0e6;
    double nu = 0.3;
    double biot = 1.0;
    double K_d() const { return E / (3.0 * (1.0 - 2.0 * nu)); }
    double eps() const { return eps_prev + d_eps; }
    double phi() const { return (phi_prev + d_eps) / (1.0 + d_eps); }
    PotentialExchangeLocalSolveContext context() const
    {
        PotentialExchangeLocalSolveContext ctx{
            .phi = phi(),
            .phi_M_prev = phi_M_prev,
            .phi_m_prev = phi_m_prev,
            .volumetric_strain = eps(),
            .volumetric_strain_prev = eps_prev,
        };
        ctx.confining_pressure_p_conf = -K_d() * eps();
        ctx.drained_bulk_modulus = K_d();
        ctx.biot_coefficient = biot;
        return ctx;
    }
    void setHistory()
    {
        phi_m_prev = n_l_prev * (1.0 - phi_prev) / (1.0 - n_l_prev);
        phi_M_prev = phi_prev - phi_m_prev;
        auto const prev = computePreviousMicroLiquidDensity(
            n_l_prev, rho_LR, context(), params);
        rho_lR_prev = prev.rho_lR;
        rho_l_prev = phi_m_prev * prev.rho_lR;
    }
    MicroCeilingKktSolveData solve() const
    {
        return solveReferenceMassStorageKktState(
            n_l_prev, rho_l_prev, rho_lR_prev, dt, rho_LR, alpha_bar, mu,
            computeYoungLaplaceMacroPotential(p_L, rho_LR,
                                              params.pressure_tolerance),
            context(), params);
    }
};
}  // namespace

TEST(RichardsMechanics, DSMv4DroppedMicroBiotRateBooksAndTangent)
{
    // (a) bitwise the books' expression (massbook_v3.py:170 order, W = dt = 1).
    std::mt19937_64 rng(4);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (int i = 0; i < 10000; ++i)
    {
        double const sat = u(rng);
        double const rho = 900.0 + 200.0 * u(rng);
        double const phim = 0.5 * u(rng);
        double const phim_prev = 0.5 * u(rng);
        double const de = 1e-3 * (u(rng) - 0.5);
        double const books = sat * rho * ((phim - phim_prev) + phim * de) * 1.0;
        EXPECT_EQ(droppedMicroBiotRate(sat, rho, phim, phim_prev, de, 1.0).T,
                  books);
    }
    // (b) the tangent struct against central differences.
    double const S = 0.3;
    double const rho = 1000.0;
    double const pm = 0.2;
    double const pmp = 0.19;
    double const de = 2e-4;
    double const dt = 50.0;
    auto const T = [&](double s_, double r_, double p_, double d_)
    { return droppedMicroBiotRate(s_, r_, p_, pmp, d_, dt).T; };
    auto const an = droppedMicroBiotRate(S, rho, pm, pmp, de, dt);
    auto const cd = [&](auto f, double x)
    {
        double const h = std::cbrt(v4_eps_mach) * std::max(std::abs(x), 1e-3);
        return (f(x + h) - f(x - h)) / (2 * h);
    };
    double const fS = cd([&](double x) { return T(x, rho, pm, de); }, S);
    double const fr = cd([&](double x) { return T(S, x, pm, de); }, rho);
    double const fp = cd([&](double x) { return T(S, rho, x, de); }, pm);
    double const fd = cd([&](double x) { return T(S, rho, pm, x); }, de);
    EXPECT_NEAR(fS, an.dT_dS_L, 1e-8 * std::abs(an.dT_dS_L));
    EXPECT_NEAR(fr, an.dT_drho_LR, 1e-8 * std::abs(an.dT_drho_LR));
    EXPECT_NEAR(fp, an.dT_dphi_m, 1e-8 * std::abs(an.dT_dphi_m));
    EXPECT_NEAR(fd, an.dT_ddeps, 1e-6 * std::abs(an.dT_ddeps));
    std::cout << "UT-V4-4b dS " << an.dT_dS_L << "/" << fS << " drho "
              << an.dT_drho_LR << "/" << fr << " dphi_m " << an.dT_dphi_m << "/"
              << fp << " ddeps " << an.dT_ddeps << "/" << fd << "\n";
}

TEST(RichardsMechanics, DSMv4DroppedMicroBiotRateActiveIdentity)
{
    // (c) Active branch: T dt == S_L (rho_LR/rho_lR) dt S_s (DESIGN_V4.md 1.1;
    // exact when rho_lR == rho_lR,prev, which the v4Deck EOS gives: rho_lR = 100 +
    // rho_LR). An at-ceiling state (n_l_prev = phi_prev, phi_M_prev = 0) under
    // suction mu_LR = mu_lR(n_max)/2 (the construction of MicroCeilingKkt.cpp,
    // flagged there) with compaction.
    V4Case c;
    c.n_l_prev = c.phi_prev;
    c.d_eps = 0.0;
    c.phi_m_prev = c.phi_prev;
    c.phi_M_prev = 0.0;
    {
        auto const prev = computePreviousMicroLiquidDensity(
            c.n_l_prev, c.rho_LR, c.context(), c.params);
        c.rho_lR_prev = prev.rho_lR;
        c.rho_l_prev = c.phi_m_prev * prev.rho_lR;
    }
    auto const mu_at_wall =
        computeActiveMicroPotential(c.phi(), c.rho_LR, c.context(), c.params)
            .mu_lR;
    c.p_L = 0.5 * c.rho_LR * mu_at_wall;
    int n_active = 0;
    double worst = 0.0;
    for (double const de : {0.0, -1e-5, -1e-4})
    {
        c.d_eps = de;
        auto const r = c.solve();
        if (r.status != MicroCeilingKktStatus::Active)
        {
            continue;
        }
        ++n_active;
        double const S_L = 0.3;
        double const phi_m = r.local.n_l;  // = n_max = phi_s on the Active branch
        auto const T = droppedMicroBiotRate(S_L, c.rho_LR, phi_m, c.phi_m_prev,
                                            de, c.dt);
        double const lhs = T.T * c.dt;
        double const rhs =
            S_L * (c.rho_LR / r.local.rho_lR) * c.dt * r.exchange_received;
        double const scale = S_L * c.rho_LR * (phi_m + c.phi_m_prev);
        worst = std::max(worst, std::abs(lhs - rhs) / scale);
        EXPECT_EQ(r.local.rho_lR, c.rho_lR_prev);
        EXPECT_LE(std::abs(lhs - rhs), 64 * v4_eps_mach * scale)
            << "d_eps " << de;
        std::cout << "UT-V4-4c d_eps=" << de << " T dt=" << lhs
                  << " S(rho_LR/rho_lR) dt S_s=" << rhs << "\n";
    }
    EXPECT_GE(n_active, 1);
    std::cout << "UT-V4-4c Active cases=" << n_active
              << " worst |diff|/scale=" << worst << "\n";
}

TEST(RichardsMechanics, DSMv4InteriorMicroWaterContentTangents)
{
    // (d) Interior: dn_l/dp_L (computeImplicitNlDpL) and dn_l/deps_v
    // (computeImplicitNlDEpsV, with d mu_lR/d eps_v at fixed n_l built exactly
    // as the assembler captures it: integrable partner + live-K chain) against
    // central differences of the local KKT solve re-solved at p_L +- h and
    // eps_v +- h. FD floor from the local solve's residual tolerance 1e-10
    // (RichardsMechanicsFEM-impl.h, test-only use of that literal).
    // An interior state: dt = 1 s and a macro potential slightly above the
    // micro potential of the previous state, mu_LR = 0.99 mu_lR(n_l,prev)
    // (mu_lR < 0: water moves in slowly; test-only construction).
    V4Case c;
    c.dt = 1.0;
    c.setHistory();
    c.p_L = 0.99 * c.rho_LR *
            computeActiveMicroPotential(c.n_l_prev, c.rho_LR, c.context(),
                                        c.params)
                .mu_lR;
    auto const base = c.solve();
    std::cout << "UT-V4-4d p_L=" << c.p_L << " status="
              << static_cast<int>(base.status) << " n_l=" << base.local.n_l
              << " n_l_prev=" << c.n_l_prev << "\n";
    ASSERT_EQ(base.status, MicroCeilingKktStatus::Interior);
    double const n = base.local.n_l;
    double const rho_lR = base.local.rho_lR;
    auto const ctx = c.context();
    auto const& p = c.params;
    auto const macro = computeYoungLaplaceMacroPotential(c.p_L, c.rho_LR,
                                                         p.pressure_tolerance);
    auto const micro_potential =
        computeActiveMicroPotential(n, c.rho_LR, ctx, p);
    double const alpha_M = c.alpha_bar * c.rho_LR / c.mu;
    auto const exchange = computePotentialDrivenMassExchange(
        alpha_M, macro.mu_LR, micro_potential.mu_lR);
    // dn/dp_L (existing helper, beta_LR = 0).
    double const dn_dpL = computeImplicitNlDpL(
        c.n_l_prev, c.p_L, c.dt, c.rho_LR, 0.0, c.alpha_bar, c.mu, macro,
        micro_potential, exchange, ctx, p, n, rho_lR);
    // The assembler capture (film block of assembleWithJacobian).
    double const nS = computeActiveMicroSolidVolumeFraction(n, ctx, p);
    auto const vdw = computeVanDerWaalsMicroPotential(
        n, rho_lR, nS, p.micro_solid_density_reference, p.hamaker_constant,
        p.specific_surface, microPotentialSignFactorFromParameters(p),
        effectiveAugmentationPrefactor(p, ctx.phi),
        p.potential_augmentation_exponent, 0.0, p.micro_water_content_floor);
    auto const mech = computeIntegrableMechanicalMicroPotential(
        -rho_lR * vdw.mu_lR, -rho_lR * vdw.dmu_lR_dnl,
        -rho_lR * vdw.d2mu_lR_dnl2, n, c.eps(), c.biot, c.K_d(), rho_lR);
    double const dphi = porosityDerivativeWrtVolumetricStrain(
        1.0, c.phi(), c.phi_prev, c.d_eps);
    double const dmu_fixed_n =
        mech.dmu_lR_mech_deps_v +
        (vdw.dmu_lR_dK + (vdw.dmu_lR_dK + n * vdw.ddmu_lR_dnl_dK) * c.eps()) *
            effectiveAugmentationPrefactorPhiDerivative(p, ctx.phi) * dphi;
    double const dn_deps = computeImplicitNlDEpsV(
        c.n_l_prev, c.dt, c.rho_LR, c.mu, micro_potential, exchange, ctx, p,
        dphi, dmu_fixed_n, n, rho_lR);
    // dr/dn for the FD floor.
    double const one_m_n = 1.0 - n;
    double const drdn = (1.0 - c.phi()) * rho_lR / (one_m_n * one_m_n);
    // Also an independent check of d mu_lR/d eps_v at fixed n_l: central
    // differences of computeActiveMicroPotential(n fixed, context(eps +- h)).
    double const h_e = std::cbrt(v4_eps_mach) * std::abs(c.eps());
    V4Case cp = c;
    cp.d_eps += h_e;
    V4Case cm = c;
    cm.d_eps -= h_e;
    double const dmu_fd =
        (computeActiveMicroPotential(n, c.rho_LR, cp.context(), p).mu_lR -
         computeActiveMicroPotential(n, c.rho_LR, cm.context(), p).mu_lR) /
        (2 * h_e);
    std::cout << "UT-V4-4d status Interior n=" << n << " dmu/deps|n capture="
              << dmu_fixed_n << " fd=" << dmu_fd
              << " live-K share=" << dmu_fixed_n - mech.dmu_lR_mech_deps_v
              << "\n";
    EXPECT_NEAR(dmu_fixed_n, dmu_fd,
                kV4FdSafetyFactor * std::pow(v4_eps_mach, 2.0 / 3.0) *
                        (std::abs(dmu_fd) +
                         std::abs(micro_potential.mu_lR) / std::abs(c.eps())) +
                    1e-6 * std::abs(dmu_fd));
    for (double const factor : {0.5, 1.0, 2.0})
    {
        double const he = factor * h_e;
        V4Case ep = c;
        ep.d_eps += he;
        V4Case em = c;
        em.d_eps -= he;
        auto const rp = ep.solve();
        auto const rm = em.solve();
        ASSERT_EQ(rp.status, MicroCeilingKktStatus::Interior);
        ASSERT_EQ(rm.status, MicroCeilingKktStatus::Interior);
        double const fd_e = (rp.local.n_l - rm.local.n_l) / (2 * he);
        double const floor_e = 2.0 * 1e-10 / drdn / (2 * he);
        double const tol_e = kV4FdSafetyFactor * (floor_e +
                                                std::pow(v4_eps_mach, 2.0 / 3.0) *
                                                    std::abs(dn_deps));
        EXPECT_NEAR(fd_e, dn_deps, tol_e) << "factor " << factor;
        // phi_m chain (as assembled): dphi_m/deps = (1-phi)/(1-n)^2 dn/deps -
        // n/(1-n) dphi/deps against FD of phi_m = n (1-phi)/(1-n).
        double const pmp = rp.local.n_l * (1.0 - ep.phi()) / (1.0 - rp.local.n_l);
        double const pmm = rm.local.n_l * (1.0 - em.phi()) / (1.0 - rm.local.n_l);
        double const dphim_an = (1.0 - c.phi()) / (one_m_n * one_m_n) * dn_deps -
                                n / one_m_n * dphi;
        double const dphim_fd = (pmp - pmm) / (2 * he);
        std::cout << "UT-V4-4d eps h=" << he << " dn/deps an=" << dn_deps
                  << " fd=" << fd_e << " tol=" << tol_e
                  << " dphi_m/deps an=" << dphim_an << " fd=" << dphim_fd
                  << "\n";
        EXPECT_NEAR(dphim_fd, dphim_an,
                    kV4FdSafetyFactor * (floor_e + 1e-9 * std::abs(dphim_an)) +
                        std::abs(dphim_an) * 1e-6);
    }
    double const h_p = std::cbrt(v4_eps_mach) * std::abs(c.p_L);
    for (double const factor : {0.5, 1.0, 2.0})
    {
        double const hp = factor * h_p;
        V4Case pp = c;
        pp.p_L += hp;
        V4Case pm = c;
        pm.p_L -= hp;
        double const fd_p =
            (pp.solve().local.n_l - pm.solve().local.n_l) / (2 * hp);
        double const floor_p = 2.0 * 1e-10 / drdn / (2 * hp);
        double const tol_p = kV4FdSafetyFactor * (floor_p +
                                                std::pow(v4_eps_mach, 2.0 / 3.0) *
                                                    std::abs(dn_dpL));
        std::cout << "UT-V4-4d p h=" << hp << " dn/dpL an=" << dn_dpL
                  << " fd=" << fd_p << " tol=" << tol_p << "\n";
        EXPECT_NEAR(fd_p, dn_dpL, tol_p) << "factor " << factor;
    }
}

// ── UT-V4-5 ─────────────────────────────────────────────────────────────────
TEST(RichardsMechanics, DSMv4ClosedMacroGateTruthTable)
{
    using L = MicroCeilingClosedMacroGate;
    using F = MicroCeilingSaturationGate;
    auto const chi_deck = [](double const S) { return S >= 1.0 ? 1.0 : 0.0; };
    int checked = 0;
    for (L const level : {L::Off, L::Relperm, L::BishopRelperm})
    {
        for (int bits = 0; bits < 16; ++bits)
        {
            bool const prev_active = bits & 1;
            bool const prev_phiM_zero = bits & 2;
            bool const active = bits & 4;
            bool const acted_prev = bits & 8;
            bool const acts = closedMacroGateActs(level, prev_active,
                                                  prev_phiM_zero, active);
            EXPECT_EQ(acts, level != L::Off && prev_active && prev_phiM_zero &&
                                active);
            // chi gated only at bishop_relperm; chi_prev by the marker rule.
            BishopFactorValues const in{0.25, 0.5, 0.75};
            auto const out =
                closedMacroGatedBishopFactors(level, acts, acted_prev, in, chi_deck);
            if (level == L::BishopRelperm && acts)
            {
                EXPECT_EQ(out.chi, 1.0);
                EXPECT_EQ(out.dchi_dS_L, 0.0);
                EXPECT_EQ(out.chi_prev, acted_prev ? 1.0 : 0.5);
            }
            else
            {
                EXPECT_EQ(out.chi, in.chi);
                EXPECT_EQ(out.chi_prev, in.chi_prev);
                EXPECT_EQ(out.dchi_dS_L, in.dchi_dS_L);
            }
            // Union with Fix B for k (latch_old reuses the acted_prev bit).
            for (F const fixb : {F::Off, F::BishopRelperm, F::Bishop})
            {
                bool const k_gated =
                    saturationGateActsOnRelativePermeability(fixb, acted_prev,
                                                             active) ||
                    acts;
                EXPECT_EQ(k_gated,
                          (fixb == F::BishopRelperm && acted_prev && active) ||
                              acts);
                ++checked;
            }
        }
    }
    // Fix B latch unchanged (regression of nextSaturatedLatch).
    for (int bits = 0; bits < 8; ++bits)
    {
        bool const a = bits & 1;
        bool const l = bits & 2;
        double const chi = (bits & 4) ? 1.0 : 0.0;
        EXPECT_EQ(nextSaturatedLatch(a, l, chi), a && (l || chi == 1.0));
    }
    std::cout << "UT-V4-5 truth-table cases=" << checked << "\n";
}
