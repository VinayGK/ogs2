// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// Unit tests of the v5 probe switch macro_storage_exact_time_levels (branch
// dsm_mass_conservation_v5_P_exact_2026-10-02; design ~/ogs-models/scratch/
// 2026-10-02_kkt_v5_P_exact/DESIGN_V5.md section 5.1). PROBE, NOT adopted.
// SUPPLEMENTS only: no existing test is edited.
//
// UT-V5-1 telescoping (exact dyadic arithmetic + a Higham bound), UT-V5-2
// pass-through, UT-V5-3 secant against the books' product term (a round-off
// IDENTITY, not a test of the assembled residual: review must-fix M1),
// UT-V5-4 lumping sum, UT-V5-5 single-IP tangent (FD and the cancellation of
// the phi_M chains, DESIGN_V5.md 2.3), UT-V5-6 defaults and the runtime
// admissibility predicate. The parser guards are deck tests
// (REC5/tests/parser/), not here.
//
// Every test prints its numbers (C7 convention of MicroCeilingKkt.cpp).
// Round-off bounds use Higham's gamma_n = n u / (1 - n u), u = 2^-53
// (N. J. Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed.,
// SIAM 2002, Lemma 3.1 and eq. (4.4)); each bound is derived next to its
// assertion by counting the roundings of the operation chain.
// Deck values: the SaturationTuller of the AB decks (DESIGN_V4.md 1.3, READ;
// the same constants as MacroBalanceV4.cpp). Test-only proposals (need
// Vinay's approval, as in MicroCeilingKkt.cpp / MacroBalanceV4.cpp):
// kV5FdSafetyFactor = 10 on the evaluation round-off of the FD quotients,
// kV5SampleSafety = 1.5 on the sampled third derivative, the synthetic states
// of UT-V5-3/4/5 (rho_LR = 1000 kg/m3, dt = 3600 s, the phi_M(p_L) closed form
// of UT-V5-5), the seeds.

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.h"
#include "MaterialLib/MPL/VariableType.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/RichardsMechanics/PotentialExchangeParameters.h"

using namespace ProcessLib::RichardsMechanics;

namespace
{
constexpr double kV5U = std::numeric_limits<double>::epsilon() / 2.0;  // u
double gammaN(double const n)
{
    return n * kV5U / (1.0 - n * kV5U);
}
constexpr double kV5FdSafetyFactor = 10.0;  // test-only proposal
constexpr double kV5SampleSafety = 1.5;     // test-only proposal
constexpr double kV5Rho = 1000.0;           // [kg/m3] test-only state
constexpr double kV5Dt = 3600.0;            // [s] test-only state

// AB deck SaturationTuller (MacroBalanceV4.cpp, DESIGN_V4.md 1.3).
constexpr double kV5Shape = 0.8584073464102069;
constexpr double kV5Sigma = 0.0715;
constexpr double kV5Pore = 1e-5;
constexpr double kV5Cav = 1.4e8;
double const kV5C =
    4.0 * kV5Shape * kV5Sigma * kV5Sigma / (1.0 * kV5Pore * kV5Pore);

struct V5Tuller
{
    MaterialPropertyLib::SaturationTuller sat{
        "saturation", 0.0, 1.0, 1.0, kV5Shape, kV5Pore, kV5Sigma, 0.0, kV5Cav};
    ParameterLib::SpatialPosition pos;
    double S(double const p_c) const
    {
        MaterialPropertyLib::VariableArray v;
        v.capillary_pressure = p_c;
        return std::get<double>(sat.value(v, pos, 0.0, 0.0));
    }
    double dS(double const p_c) const
    {
        MaterialPropertyLib::VariableArray v;
        v.capillary_pressure = p_c;
        v.liquid_saturation = S(p_c);
        return std::get<double>(
            sat.dValue(v, MaterialPropertyLib::Variable::capillary_pressure,
                       pos, 0.0, 0.0));
    }
};
V5Tuller const& v5Tuller()
{
    static V5Tuller const t;
    return t;
}
// Closed forms on 0 < p < cavitation (in-file derivation): E = exp(-c/p^2),
// S = 1 - E, S' = -2c p^-3 E, S'' = -E (4c^2 p^-6 - 6c p^-4),
// S''' = -E (8c^3 p^-9 - 36c^2 p^-7 + 24c p^-5).
double v5S1(double const p)
{
    double const E = std::exp(-kV5C / (p * p));
    return -2.0 * kV5C * E / (p * p * p);
}
double v5S2(double const p)
{
    double const E = std::exp(-kV5C / (p * p));
    return -E * (4.0 * kV5C * kV5C / std::pow(p, 6) - 6.0 * kV5C / std::pow(p, 4));
}
double v5S3(double const p)
{
    double const E = std::exp(-kV5C / (p * p));
    return -E * (8.0 * std::pow(kV5C, 3) / std::pow(p, 9) -
                 36.0 * kV5C * kV5C / std::pow(p, 7) +
                 24.0 * kV5C / std::pow(p, 5));
}

std::uint64_t bits(double const x)
{
    std::uint64_t b;
    std::memcpy(&b, &x, sizeof b);
    return b;
}
}  // namespace

// ── UT-V5-1 telescoping ─────────────────────────────────────────────────────
TEST(RichardsMechanics, DSMv5MacroStorageTelescoping)
{
    // (i) Exact arithmetic: S_n = i/2^20 in [0, 1], phi_n = j/2^20 in
    // [0, 0.5], rho = 1. Every difference is an exact multiple of 2^-20;
    // every product an exact multiple of 2^-40 below 2^1 (41 bits); the off
    // variant's partial sums are below 2^11 (K = 2048 steps of |dS dphi| <= 1)
    // with granularity 2^-40, i.e. at most 51 bits < 53: no operation
    // rounds, so the assertions are exact equalities.
    std::mt19937_64 rng(20261002);
    constexpr int K = 2048;
    constexpr double q = 1.0 / 1048576.0;  // 2^-20
    std::uniform_int_distribution<int> iS(0, 1048576), iP(0, 524288);
    std::vector<double> S(K + 1), P(K + 1);
    for (int n = 0; n <= K; ++n)
    {
        S[n] = iS(rng) * q;
        P[n] = iP(rng) * q;
    }
    double sum_on = 0.0, sum_off = 0.0, sum_prod = 0.0;
    for (int n = 1; n <= K; ++n)
    {
        double const dS = S[n] - S[n - 1];
        double const dP = P[n] - P[n - 1];
        double const c_on = macroStorageCoefficient(true, P[n], P[n - 1]);
        double const c_off = macroStorageCoefficient(false, P[n], P[n - 1]);
        sum_on += c_on * dS + S[n] * dP;
        sum_off += c_off * dS + S[n] * dP;
        sum_prod += 1.0 * (S[n] - S[n - 1]) * (P[n] - P[n - 1]);  // facts_post order
    }
    double const telescoped = S[K] * P[K] - S[0] * P[0];
    std::cout << std::setprecision(17) << "UT-V5-1i K=" << K
              << " sum_on=" << sum_on << " telescoped=" << telescoped
              << " sum_off-sum_on=" << sum_off - sum_on
              << " sum_dS_dphi=" << sum_prod << "\n";
    EXPECT_EQ(sum_on, telescoped);
    EXPECT_EQ(sum_off - sum_on, sum_prod);

    // (ii) Random reals, rho = kV5Rho, K2 = 10000 steps. Computed term
    // fl(rho * fl(fl(c dS^) + fl(S dphi^))) with dS^, dphi^ rounded: each of
    // the two products carries <= 4 roundings (difference, product, sum,
    // rho), so |t^ - t| <= gamma_4 a_n, a_n = rho(|c dS| + |S dphi|).
    // Recursive summation adds gamma_{K2-1} sum|t^| (Higham (4.4)). The
    // reference fl(rho fl(fl(S_K phi_K) - fl(S_0 phi_0))) carries <= 3
    // roundings: gamma_3 B, B = rho(|S_K phi_K| + |S_0 phi_0|). With
    // gamma_a + gamma_b (1 + gamma_a) <= gamma_{a+b}:
    // |sum^ - ref^| <= gamma_{K2+3} (A + B), A = sum a_n.
    constexpr int K2 = 10000;
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    std::vector<double> S2(K2 + 1), P2(K2 + 1);
    for (int n = 0; n <= K2; ++n)
    {
        S2[n] = u01(rng);
        P2[n] = 0.5 * u01(rng);
    }
    double s2 = 0.0, A = 0.0;
    for (int n = 1; n <= K2; ++n)
    {
        double const c = macroStorageCoefficient(true, P2[n], P2[n - 1]);
        double const dS = S2[n] - S2[n - 1];
        double const dP = P2[n] - P2[n - 1];
        s2 += kV5Rho * (c * dS + S2[n] * dP);
        A += kV5Rho * (std::abs(c * dS) + std::abs(S2[n] * dP));
    }
    double const ref2 = kV5Rho * (S2[K2] * P2[K2] - S2[0] * P2[0]);
    double const B = kV5Rho * (std::abs(S2[K2] * P2[K2]) + std::abs(S2[0] * P2[0]));
    double const bound2 = gammaN(K2 + 3) * (A + B);
    std::cout << "UT-V5-1ii K=" << K2 << " |sum - telescoped|="
              << std::abs(s2 - ref2) << " bound gamma_{K+3}(A+B)=" << bound2
              << " ratio=" << std::abs(s2 - ref2) / bound2 << "\n";
    EXPECT_LE(std::abs(s2 - ref2), bound2);
}

// ── UT-V5-2 pass-through ────────────────────────────────────────────────────
TEST(RichardsMechanics, DSMv5MacroStorageCoefficientPassThrough)
{
    std::vector<double> v = {0.0,
                             -0.0,
                             std::numeric_limits<double>::denorm_min(),
                             -std::numeric_limits<double>::denorm_min(),
                             std::numeric_limits<double>::min(),
                             1.0 - 1e-12,
                             0.4,
                             std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()};
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    for (int i = 0; i < 1000; ++i)
    {
        v.push_back(u(rng));
    }
    int n = 0;
    for (double const a : v)
    {
        for (double const b : v)
        {
            EXPECT_EQ(bits(macroStorageCoefficient(false, a, b)), bits(a));
            EXPECT_EQ(bits(macroStorageCoefficient(true, a, b)), bits(b));
            ++n;
        }
    }
    std::cout << "UT-V5-2 pairs checked bitwise: " << n << "\n";
}

// ── UT-V5-3 secant against the books' product term (an identity) ───────────
TEST(RichardsMechanics, DSMv5SecantEqualsBooksProductTerm)
{
    // C = fl(fl(fl(rho dphi^) s^) dpc^) with s^ = fl(dS^/dpc^): s^ dpc^ =
    // dS^ (1 + d) exactly in the reals, so C = X (1 + theta_4) with
    // X = rho dS^ dphi^; P = fl(fl(rho dS^) dphi^) = X (1 + theta_2). Hence
    // |C - P| <= gamma_6 |X| <= gamma_6 / (1 - gamma_2) |P|. This is
    // round-off of an algebraic identity, NOT a measurement of the assembled
    // residual (review must-fix M1); the time levels the assembly actually
    // uses are tested by the per-step identity on V5 runs (TESTS_V5.md).
    auto const& T = v5Tuller();
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> lp(std::log(2.0e3), std::log(1.0e6));
    std::uniform_real_distribution<double> uphi(0.0, 0.3);
    double max_rel = 0.0;
    int n_nonzero = 0;
    double const rel_bound = gammaN(6) / (1.0 - gammaN(2));
    for (int i = 0; i < 200000; ++i)
    {
        double const pc = std::exp(lp(rng));
        double const pcp = std::exp(lp(rng));
        double const phiM = uphi(rng), phiMp = uphi(rng);
        double const S = T.S(pc), Sp = T.S(pcp);
        if (pc == pcp)
        {
            continue;
        }
        double const s = (S - Sp) / (pc - pcp);
        double const C = kV5Rho * (phiM - phiMp) * s * (pc - pcp);
        double const Pb = kV5Rho * (S - Sp) * (phiM - phiMp);
        if (Pb != 0.0)
        {
            ++n_nonzero;
            max_rel = std::max(max_rel, std::abs(C - Pb) / std::abs(Pb));
        }
        EXPECT_LE(std::abs(C - Pb), rel_bound * std::abs(Pb));
    }
    // dp_c == 0: S is a function of p_c alone (SaturationTuller), so the
    // saturation change is exactly 0 and the secant term (code: tangent
    // branch, times dp_c = 0) is 0.
    for (double const pc : {2.0e3, 2.2e4, 6.0e4, 1.0e6, 1.0e8})
    {
        EXPECT_EQ(T.S(pc) - T.S(pc), 0.0);
    }
    std::cout << std::setprecision(6) << "UT-V5-3 samples with P != 0: "
              << n_nonzero << " max |C - P|/|P| = " << max_rel
              << " bound gamma_6/(1-gamma_2) = " << rel_bound << "\n";
}

// ── UT-V5-4 lumping sum ─────────────────────────────────────────────────────
TEST(RichardsMechanics, DSMv5LumpedSecantResidualSum)
{
    // Bilinear quad, 2x2 Gauss, own shape functions N_a = (1 +- xi)(1 +- eta)/4.
    // M = sum_q N_q^T c_q w_q N_q; lumped L = diag(colsum M); the sum of the
    // lumped residual entries sum_j colsum_j dp_j against the IP form
    // sum_q c_q w_q (N_q . dp). Exact identity when sum_a N_qa = 1.
    // Bound: lumped chain <= 13 roundings per elementary term (3 in the
    // 4-factor product, 3 over q, 3 in the column sum, 1 times dp_j, 3 over
    // j): gamma_13 T, T = sum_q sum_i sum_j |N_qi c_q w_q N_qj dp_j|; IP chain
    // <= 9 (dot product of 4: gamma_4; times c, w; 3 over q): gamma_9 T'
    // (T' = sum_q |c_q w_q| sum_j |N_qj dp_j| <= T / min_q sum_i N_qi, here
    // T' <= T(1 + gamma_3)); partition-of-unity defect of the computed N:
    // sum_q |c_q w_q (N_q . dp)| (|d_q| + gamma_3 sum_i N_qi), d_q = fl(sum N) - 1.
    double const g = 1.0 / std::sqrt(3.0);
    std::array<std::array<double, 2>, 4> const gp = {
        {{-g, -g}, {-g, g}, {g, -g}, {g, g}}};
    std::mt19937_64 rng(11);
    std::uniform_real_distribution<double> u(-1.0, 1.0), up(0.1, 2.0);
    double max_ratio = 0.0;
    for (int trial = 0; trial < 10000; ++trial)
    {
        Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
        Eigen::Vector4d dp;
        for (int a = 0; a < 4; ++a)
        {
            dp[a] = 1e5 * u(rng);
        }
        double ip_form = 0.0, T = 0.0, defect = 0.0;
        for (auto const& [xi, eta] : gp)
        {
            Eigen::RowVector4d N;
            N << 0.25 * (1 - xi) * (1 - eta), 0.25 * (1 + xi) * (1 - eta),
                0.25 * (1 + xi) * (1 + eta), 0.25 * (1 - xi) * (1 + eta);
            double const c = 1e3 * u(rng) * 1e-6;  // rho a_S s [kg/m3/Pa]
            double const w = up(rng);
            M.noalias() += N.transpose() * c * N * w;
            double const Ndp = N.dot(dp);
            ip_form += c * w * Ndp;
            double const sN = N.sum();
            for (int i = 0; i < 4; ++i)
            {
                for (int j = 0; j < 4; ++j)
                {
                    T += std::abs(N[i] * c * w * N[j] * dp[j]);
                }
            }
            defect += std::abs(c * w * Ndp) *
                      (std::abs(sN - 1.0) + gammaN(3) * std::abs(sN));
        }
        Eigen::Vector4d const L = M.colwise().sum().transpose();
        double const lumped_sum = (L.asDiagonal() * dp).sum();
        double const bound =
            gammaN(13) * T + gammaN(9) * T * (1.0 + gammaN(3)) + defect;
        max_ratio = std::max(max_ratio, std::abs(lumped_sum - ip_form) / bound);
        EXPECT_LE(std::abs(lumped_sum - ip_form), bound);
    }
    std::cout << "UT-V5-4 10000 random elements: max |lumped sum - IP form| / "
                 "bound = "
              << max_ratio << "\n";
}

// ── UT-V5-5 tangent of a single IP ──────────────────────────────────────────
namespace
{
// Synthetic Interior IP (test-only closed form): phi_M(p_L) = phi0 +
// a tanh((p_L - x0)/b); phi_M,prev fixed; S_L,prev = S(p_c,prev).
struct V5Ip
{
    double phi0 = 0.08, a = 0.01, x0 = -6.0e4, b = 2.0e4;
    double phiMp = 0.10;
    double pcp = 7.0e4;
    double phi(double x) const { return phi0 + a * std::tanh((x - x0) / b); }
    double dphi(double x) const
    {
        double const sc = 1.0 / std::cosh((x - x0) / b);
        return a / b * sc * sc;
    }
    double d2phi(double x) const
    {
        double const sc = 1.0 / std::cosh((x - x0) / b);
        return -2.0 * a / (b * b) * sc * sc * std::tanh((x - x0) / b);
    }
    double d3phi(double x) const
    {
        double const sc = 1.0 / std::cosh((x - x0) / b);
        double const th = std::tanh((x - x0) / b);
        return -2.0 * a / (b * b * b) * (sc * sc * sc * sc - 2.0 * sc * sc * th * th);
    }
    double Sp() const { return v5Tuller().S(pcp); }
    // residual parts in p_L (x), p_c = -x
    double r5(double x) const
    {
        return kV5Rho * phiMp * (v5Tuller().S(-x) - Sp()) / kV5Dt;
    }
    double r4(double x) const
    {
        return kV5Rho * phi(x) * (v5Tuller().S(-x) - Sp()) / kV5Dt;
    }
    // third derivatives (D(x) = S(-x) - S_p: D' = -S', D'' = S'', D''' = -S''')
    double r5_3(double x) const { return -kV5Rho * phiMp * v5S3(-x) / kV5Dt; }
    double r4_3(double x) const
    {
        double const D = v5Tuller().S(-x) - Sp();
        return kV5Rho / kV5Dt *
               (d3phi(x) * D - 3.0 * d2phi(x) * v5S1(-x) +
                3.0 * dphi(x) * v5S2(-x) - phi(x) * v5S3(-x));
    }
};

struct FdResult
{
    double fd, bound, h;
};
// Central difference with the step derived from the bound
// |FD - f'| <= h^2/6 M3 + eta/h, eta = kV5FdSafetyFactor u |f|_scale
// (evaluation round-off of each f, test-only factor): h = cbrt(3 eta / M3);
// M3 = kV5SampleSafety * max |f'''| sampled at 65 points of [x-h, x+h].
template <typename F, typename F3>
FdResult centralFd(F&& f, F3&& f3, double const x, double const f_scale)
{
    double const eta = kV5FdSafetyFactor * kV5U * f_scale;
    double M3 = std::abs(f3(x));
    double h = std::cbrt(3.0 * eta / M3);
    double m = 0.0;
    for (int k = 0; k <= 64; ++k)
    {
        m = std::max(m, std::abs(f3(x - h + 2.0 * h * k / 64.0)));
    }
    M3 = kV5SampleSafety * m;
    double const fd = (f(x + h) - f(x - h)) / (2.0 * h);
    return {fd, h * h / 6.0 * M3 + eta / h + 4.0 * kV5U * std::abs(fd), h};
}
}  // namespace

TEST(RichardsMechanics, DSMv5SingleIpTangentAndCancellation)
{
    V5Ip const ip;
    auto const& T = v5Tuller();
    std::cout << std::setprecision(10);
    for (double const pc : {4.0e4, 6.0e4, 9.0e4, 1.5e5})
    {
        double const x = -pc;
        double const D = T.S(pc) - ip.Sp();
        double const dSdpc = T.dS(pc);           // MPL dValue, as the code
        double const dDdx = -dSdpc;
        // closed form S' against MPL dValue: printed, not asserted
        std::cout << "UT-V5-5 p_c=" << pc << " S=" << T.S(pc)
                  << " dS/dp_c MPL=" << dSdpc << " closed form=" << v5S1(pc)
                  << "\n";
        double const scale5 = kV5Rho * ip.phiMp * (std::abs(T.S(pc)) + std::abs(ip.Sp())) / kV5Dt;
        double const scale4 = kV5Rho * (ip.phi0 + std::abs(ip.a)) * (std::abs(T.S(pc)) + std::abs(ip.Sp())) / kV5Dt;

        // (a) v5: the S-chain alone is the exact derivative of r5.
        double const an5 = kV5Rho * ip.phiMp * dDdx / kV5Dt;
        auto const fd5 = centralFd([&](double y) { return ip.r5(y); },
                                   [&](double y) { return ip.r5_3(y); }, x, scale5);
        std::cout << "  (a) v5 S-chain=" << an5 << " FD=" << fd5.fd
                  << " |diff|=" << std::abs(fd5.fd - an5)
                  << " bound=" << fd5.bound << " h=" << fd5.h << "\n";
        EXPECT_LE(std::abs(fd5.fd - an5), fd5.bound + 4.0 * kV5U * std::abs(an5));

        // (b) v4: its analytic S-chain misses rho dS dphi_M/dp_L / dt.
        double const an4 = kV5Rho * ip.phi(x) * dDdx / kV5Dt;
        double const miss_expected = kV5Rho * D * ip.dphi(x) / kV5Dt;
        auto const fd4 = centralFd([&](double y) { return ip.r4(y); },
                                   [&](double y) { return ip.r4_3(y); }, x, scale4);
        double const miss = fd4.fd - an4;
        std::cout << "  (b) v4 S-chain=" << an4 << " FD=" << fd4.fd
                  << " miss=" << miss << " rho dS dphi/dp_L/dt=" << miss_expected
                  << " |diff|=" << std::abs(miss - miss_expected)
                  << " bound=" << fd4.bound << "\n";
        EXPECT_LE(std::abs(miss - miss_expected),
                  fd4.bound + 8.0 * kV5U * (std::abs(an4) + std::abs(miss_expected)));

        // (c) [correction derivative incl. phi_M chain] + [v4 omitted chain]
        //     == [v5 S-chain] - [v4 S-chain]  (algebra; bound: each side has
        //     <= 6 roundings per product term and one subtraction / addition:
        //     gamma_8 times the sum of the absolute terms).
        double const corr_deriv = kV5Rho / kV5Dt *
                                  ((ip.phiMp - ip.phi(x)) * dDdx - ip.dphi(x) * D);
        double const omitted = kV5Rho * ip.dphi(x) * D / kV5Dt;
        double const lhs = corr_deriv + omitted;
        double const rhs = an5 - an4;
        double const abs_terms =
            kV5Rho / kV5Dt *
            (std::abs(ip.phiMp * dDdx) + std::abs(ip.phi(x) * dDdx) +
             2.0 * std::abs(ip.dphi(x) * D));
        std::cout << "  (c) lhs=" << lhs << " rhs=" << rhs
                  << " |diff|=" << std::abs(lhs - rhs)
                  << " bound=" << gammaN(8) * abs_terms << "\n";
        EXPECT_LE(std::abs(lhs - rhs), gammaN(8) * abs_terms);
    }
}

// ── UT-V5-6 defaults and the runtime admissibility predicate ────────────────
TEST(RichardsMechanics, DSMv5DefaultsAndAdmissibility)
{
    PotentialExchangeParameters const p;
    EXPECT_FALSE(p.macro_storage_exact_time_levels);
    EXPECT_FALSE(isMacroStorageExactTimeLevels(&p));
    EXPECT_FALSE(isMacroStorageExactTimeLevels(nullptr));
    PotentialExchangeParameters q;
    q.macro_storage_exact_time_levels = true;
    EXPECT_TRUE(isMacroStorageExactTimeLevels(&q));
    // The v4 guard function is unchanged by the new member.
    EXPECT_FALSE(anyV4SwitchOn(q));

    EXPECT_TRUE(macroStorageExactTimeLevelsAdmissible(0.0, 0.0));
    EXPECT_TRUE(macroStorageExactTimeLevelsAdmissible(-0.0, 0.0));  // -0 == 0
    EXPECT_FALSE(macroStorageExactTimeLevelsAdmissible(
        std::numeric_limits<double>::denorm_min(), 0.0));
    EXPECT_FALSE(macroStorageExactTimeLevelsAdmissible(
        0.0, std::numeric_limits<double>::denorm_min()));
    EXPECT_FALSE(macroStorageExactTimeLevelsAdmissible(4.5e-10, 0.0));
    EXPECT_FALSE(macroStorageExactTimeLevelsAdmissible(
        std::numeric_limits<double>::quiet_NaN(), 0.0));
    std::cout << "UT-V5-6 defaults false; admissibility predicate: 6 cases\n";
}
