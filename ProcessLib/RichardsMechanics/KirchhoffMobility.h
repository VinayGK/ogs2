// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

// Variant 1a of the v4 round (branch dsm_mass_conservation_v4_tm_krel_
// 2026-10-02; Vinay's ruling 2026-10-02 "(go with L + drop T_m) x (1a, 1b
// separate)"; design ~/ogs-models/scratch/2026-10-02_kkt_v4_tm_krel/
// DESIGN_V4.md 2.3). NOT adopted.
//
// Kirchhoff element-mean relative permeability. One mobility per element:
//
//   kbar_e = (1/(b - a)) * int_a^b k_h(s) ds   (b > a),   kbar_e = k_h(a) (b == a),
//
// a = min and b = max of the element's nodal capillary pressures, k_h the
// piecewise-linear interpolant of the deck law k_rel(S_L(p_c)) on a log-uniform
// p_c grid with constant extensions. Phi_h = int k_h is a Kirchhoff potential
// with Phi_h' = k_h > 0, so the flux of an element whose nodal p_c varies along
// one direction is monotone in that direction (DESIGN_V4.md 2.3.5; the range
// mean over [min, max] of the nodal values is the DESIGN's reading for 2D quads,
// DESIGN_V4.md 2.3.1). Literature (locators not re-checked here): Bastian 1999
// (habilitation, Kiel), Forsyth, Wu & Pruess 1995, Adv. Water Resour. 18:25-38.
//
// No tolerance literal is used: the constant branches are exact comparisons,
// p_sat and p_k are found by bisection on the bit patterns of non-negative
// doubles (exactly the last double with k == k_sat, the first with k ==
// k_flat), and the derivatives are written without differences of nearly equal
// numbers (DESIGN_V4.md 2.3.3).

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace ProcessLib::RichardsMechanics
{
struct KirchhoffMeanResult
{
    double k_bar;  // [-] element mobility
    double dk_da;  // [1/Pa] d kbar / d a (a = min nodal p_c)
    double dk_db;  // [1/Pa] d kbar / d b (b = max nodal p_c)
};

class KirchhoffMobilityTable
{
public:
    // Largest x in [lo, hi) with pred(x) true, given pred(lo) true, pred(hi)
    // false, 0 <= lo < hi (bisection on the bit patterns of non-negative
    // doubles, whose integer order is the value order). Exact: the result r has
    // pred(r) true and pred(nextafter(r, +inf)) false when pred is monotone.
    template <typename Pred>
    static double largestTrue(Pred&& pred, double const lo, double const hi)
    {
        auto lb = std::bit_cast<std::uint64_t>(lo);
        auto hb = std::bit_cast<std::uint64_t>(hi);
        while (hb - lb > 1)
        {
            std::uint64_t const mb = lb + (hb - lb) / 2;
            if (pred(std::bit_cast<double>(mb)))
            {
                lb = mb;
            }
            else
            {
                hb = mb;
            }
        }
        return std::bit_cast<double>(lb);
    }

    // Build the table of the deck composite k_of_pc(p_c) = k_rel(S_L(p_c)),
    // which must be non-increasing in p_c with constant branches k_sat on
    // (-inf, p_sat] and k_flat on [p_k, +inf) (checked by the caller through
    // the property types, DESIGN_V4.md 2.3.2 item 1).
    //   p_sat: the largest double with k == k_sat = k_of_pc(0);
    //   p_k:   the smallest double with k == k_flat = k_of_pc(DBL_MAX) (the
    //          floor of the k law, or the frozen post-cavitation value).
    //   grid:  ceil(cells_per_decade * log10(p_k/p_sat)) log-uniform cells,
    //          node values = the exact deck law at the nodes.
    template <typename KOfPc>
    static KirchhoffMobilityTable build(KOfPc&& k_of_pc,
                                        int const cells_per_decade)
    {
        KirchhoffMobilityTable t;
        t.k_sat_ = k_of_pc(0.0);
        double const big = std::numeric_limits<double>::max();
        t.k_flat_ = k_of_pc(big);
        if (t.k_flat_ == t.k_sat_)
        {
            // Constant law: no table, kbar = k_sat everywhere.
            t.constant_ = true;
            t.g_ = {0.0};
            t.k_ = {t.k_sat_};
            t.s_ = {};
            return t;
        }
        auto const is_sat = [&](double const p) { return k_of_pc(p) == t.k_sat_; };
        // Upper bracket of p_sat: doubling from 1 Pa (k(big) != k_sat).
        double hi = 1.0;
        while (is_sat(hi) && hi < big / 2)
        {
            hi *= 2.0;
        }
        if (is_sat(hi))
        {
            hi = big;
        }
        double const p_sat = largestTrue(is_sat, 0.0, hi);
        auto const is_flat = [&](double const p)
        { return k_of_pc(p) == t.k_flat_; };
        // p_k = the smallest double with is_flat: the successor of the largest
        // double in [p_sat, big) with !is_flat.
        double const last_not_flat =
            largestTrue([&](double const p) { return !is_flat(p); }, p_sat, big);
        double const p_k = std::nextafter(last_not_flat, big);

        int const n_cells = std::max(
            1, static_cast<int>(std::ceil(cells_per_decade *
                                          std::log10(p_k / p_sat))));
        t.g_.resize(n_cells + 1);
        t.k_.resize(n_cells + 1);
        t.s_.resize(n_cells);
        double const ratio = p_k / p_sat;
        for (int i = 0; i <= n_cells; ++i)
        {
            t.g_[i] = (i == 0)         ? p_sat
                      : (i == n_cells) ? p_k
                                       : p_sat * std::pow(ratio,
                                                          static_cast<double>(
                                                              i) /
                                                              n_cells);
            t.k_[i] = k_of_pc(t.g_[i]);
        }
        for (int i = 0; i < n_cells; ++i)
        {
            t.s_[i] = (t.k_[i + 1] - t.k_[i]) / (t.g_[i + 1] - t.g_[i]);
        }
        return t;
    }

    bool isConstant() const { return constant_; }
    double pSat() const { return g_.front(); }
    double pK() const { return g_.back(); }
    double kSat() const { return k_sat_; }
    double kFlat() const { return k_flat_; }
    std::size_t numberOfCells() const { return s_.size(); }
    std::vector<double> const& nodes() const { return g_; }
    std::vector<double> const& nodeValues() const { return k_; }

    // Index of the cell [g_j, g_{j+1}) that contains x, for g_0 <= x < g_N
    // (at a grid node: the cell to the right).
    std::size_t cellOf(double const x) const
    {
        auto const it = std::upper_bound(g_.begin(), g_.end(), x);
        return static_cast<std::size_t>(std::distance(g_.begin(), it)) - 1;
    }

    // k_h(x): k_sat for x <= g_0, k_flat for x >= g_N, linear in between.
    double value(double const x) const
    {
        if (constant_ || x <= g_.front())
        {
            return k_sat_;
        }
        if (x >= g_.back())
        {
            return k_flat_;
        }
        std::size_t const j = cellOf(x);
        return k_[j] + s_[j] * (x - g_[j]);
    }

    // Slope of k_h in the cell to the right of x (0 on the extensions).
    double slopeRight(double const x) const
    {
        if (constant_ || x < g_.front() || x >= g_.back())
        {
            return 0.0;
        }
        return s_[cellOf(x)];
    }

    // int_{g_0}^{x} k_h(s) ds (negative for x < g_0), summed cell by cell.
    // For tests and the creation log; not used in the assembly.
    double potential(double const x) const
    {
        if (constant_)
        {
            return k_sat_ * x;
        }
        if (x <= g_.front())
        {
            return k_sat_ * (x - g_.front());
        }
        double q = 0.0;
        std::size_t const n = s_.size();
        for (std::size_t j = 0; j < n && g_[j] < x; ++j)
        {
            double const u = std::min(x, g_[j + 1]);
            q += 0.5 * (k_[j] + value(u)) * (u - g_[j]);
        }
        if (x > g_.back())
        {
            q += k_flat_ * (x - g_.back());
        }
        return q;
    }

    // kbar(a, b) and its exact partial derivatives (DESIGN_V4.md 2.3.3), a <= b.
    //   dkbar/db = (1/D^2) sum_j s_j (u_j - l_j)((u_j - a) + (l_j - a))/2,
    //   dkbar/da = (1/D^2) sum_j s_j (u_j - l_j)((b - u_j) + (b - l_j))/2,
    // D = the sum of the piece lengths, Q = the sum of the piece integrals
    // (exact trapezoids of the linear k_h), kbar = Q/D (a convex combination of
    // the node values: min k_h <= kbar <= max k_h, no 0/0).
    KirchhoffMeanResult mean(double const a, double const b) const
    {
        if (constant_)
        {
            return {k_sat_, 0.0, 0.0};
        }
        if (a >= g_.back())
        {
            return {k_flat_, 0.0, 0.0};
        }
        if (b <= g_.front())
        {
            return {k_sat_, 0.0, 0.0};
        }
        if (a == b)
        {
            double const half_slope = 0.5 * slopeRight(a);
            return {value(a), half_slope, half_slope};
        }
        double q = 0.0;       // [Pa] integral of k_h over [a, b]
        double len = 0.0;     // [Pa] sum of the piece lengths
        double sum_b = 0.0;   // [Pa] sum s (u-l)((u-a)+(l-a))/2
        double sum_a = 0.0;   // [Pa] sum s (u-l)((b-u)+(b-l))/2
        if (a < g_.front())
        {
            double const u = std::min(b, g_.front());
            q += k_sat_ * (u - a);
            len += u - a;
        }
        std::size_t const n = s_.size();
        std::size_t const j_begin = a <= g_.front() ? 0 : cellOf(a);
        for (std::size_t j = j_begin; j < n && g_[j] < b; ++j)
        {
            double const l = std::max(a, g_[j]);
            double const u = std::min(b, g_[j + 1]);
            if (!(u > l))
            {
                continue;
            }
            double const k_l = k_[j] + s_[j] * (l - g_[j]);
            double const k_u = k_[j] + s_[j] * (u - g_[j]);
            double const w = u - l;
            q += 0.5 * (k_l + k_u) * w;
            len += w;
            sum_b += s_[j] * w * ((u - a) + (l - a)) * 0.5;
            sum_a += s_[j] * w * ((b - u) + (b - l)) * 0.5;
        }
        if (b > g_.back())
        {
            double const l = std::max(a, g_.back());
            q += k_flat_ * (b - l);
            len += b - l;
        }
        double const len2 = len * len;
        return {q / len, sum_a / len2, sum_b / len2};
    }

private:
    bool constant_ = false;
    double k_sat_ = 1.0;
    double k_flat_ = 1.0;
    std::vector<double> g_;  // [Pa] nodes g_0 = p_sat < ... < g_N = p_k
    std::vector<double> k_;  // [-] deck law at the nodes
    std::vector<double> s_;  // [1/Pa] cell slopes
};

// Element mobility from the nodal liquid pressures p_L (p_c = -p_L) and its
// gradient d kbar / d p_L,n (DESIGN_V4.md 2.3.4): with the tie sets T_a, T_b of
// the minimum and maximum (exact equality),
//   g_n = -[dkbar/da * 1{n in T_a}/|T_a| + dkbar/db * 1{n in T_b}/|T_b|],
// an element of Clarke's generalised gradient that splits the derivative
// equally among tied nodes; exact along a uniform shift of the tied nodes.
template <int NPOINTS>
struct KirchhoffElementMobility
{
    double k_bar = 1.0;
    std::array<double, NPOINTS> dk_dpL{};
    double a = 0.0;  // [Pa] min nodal p_c
    double b = 0.0;  // [Pa] max nodal p_c
    unsigned argmin_mask = 0;  // bit n set: node n in T_a
    unsigned argmax_mask = 0;  // bit n set: node n in T_b
    KirchhoffMeanResult mean{1.0, 0.0, 0.0};
};

template <int NPOINTS, typename PVector>
KirchhoffElementMobility<NPOINTS> kirchhoffElementMobility(
    KirchhoffMobilityTable const& table, PVector const& p_L)
{
    KirchhoffElementMobility<NPOINTS> out;
    std::array<double, NPOINTS> p_c{};
    for (int n = 0; n < NPOINTS; ++n)
    {
        p_c[n] = -p_L[n];
    }
    out.a = *std::min_element(p_c.begin(), p_c.end());
    out.b = *std::max_element(p_c.begin(), p_c.end());
    unsigned n_a = 0;
    unsigned n_b = 0;
    for (int n = 0; n < NPOINTS; ++n)
    {
        if (p_c[n] == out.a)
        {
            out.argmin_mask |= 1u << n;
            ++n_a;
        }
        if (p_c[n] == out.b)
        {
            out.argmax_mask |= 1u << n;
            ++n_b;
        }
    }
    out.mean = table.mean(out.a, out.b);
    out.k_bar = out.mean.k_bar;
    for (int n = 0; n < NPOINTS; ++n)
    {
        double const share_a =
            (out.argmin_mask >> n) & 1u ? out.mean.dk_da / n_a : 0.0;
        double const share_b =
            (out.argmax_mask >> n) & 1u ? out.mean.dk_db / n_b : 0.0;
        out.dk_dpL[n] = -(share_a + share_b);  // dp_c/dp_L = -1
    }
    return out;
}
}  // namespace ProcessLib::RichardsMechanics
