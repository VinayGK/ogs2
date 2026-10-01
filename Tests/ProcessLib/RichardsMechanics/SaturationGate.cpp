// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// Unit tests of the latched saturation gate (potential_exchange
// micro_ceiling_saturation_gate; branch
// dsm_mass_conservation_v3_kkt_vii_gate_2026-10-01; design part B.4 of
// ~/ogs-models/scratch/2026-10-01_kkt_iv_vii_fixes/DESIGN_FIXES.md; ruled by
// Vinay 2026-10-01 ~15:15 CEST; NOT adopted). Every test is a SUPPLEMENT: no
// existing test is edited.
//
// Expected values: the latch rule and the gate predicates are the definition of
// DESIGN_FIXES.md B.4 (in-file derivation: truth tables below). The air-entry
// edge of the deck's Tuller retention is derived in the file, not tuned:
//   S_L(p_c) = 1 - exp(-C_T/p_c^2) rounds to exactly 1 in double precision when
//   exp(-C_T/p_c^2) <= 2^-54 (half the spacing 2^-53 of the doubles just below
//   1), i.e. for p_c <= p_edge = sqrt(C_T / (54 ln 2)),
//   C_T = 4 F_gamma sigma^2 / (A_n L^2).
// The parameters are those of the Model VII deck (area_factor_tuller 1.0,
// pore_area_shapefactor_tuller 0.8584073464102069, characteristic_pore_size
// 1e-5 m, surface_tension 0.0715 N/m, no pressure tolerance; relative
// permeability: floor 1e-2, enhancement 1.0, exponent 3; Bishop cutoff 1);
// their provenance is the deck's own header comment (TODO(Vinay) of that
// comment: the source of the Tuller geometry is outstanding; the test uses the
// values as the deck consumes them, it does not endorse them). p_edge = 2165.57
// Pa is the value MEASURED in the stall states of DESIGN_FIXES.md B.1 (2.166
// kPa).
//
// Nothing here is a fit-and-verify: no parameter is calibrated; every assertion
// concerns the logic of the gate and the FD-consistency of the gated function.

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <variant>
#include <vector>

#include "MaterialLib/MPL/Properties/BishopsSaturationCutoff.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.h"
#include "MaterialLib/MPL/Properties/RelativePermeability/RelPermGeneralizedPower.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/RichardsMechanics/PotentialExchangeParameters.h"

using namespace ProcessLib::RichardsMechanics;
namespace MPL = MaterialPropertyLib;

namespace
{
constexpr auto kOff = MicroCeilingSaturationGate::Off;
constexpr auto kBishop = MicroCeilingSaturationGate::Bishop;
constexpr auto kBishopRelperm = MicroCeilingSaturationGate::BishopRelperm;

// The deck's laws, as the Model VII deck consumes them (see the header).
struct DeckLaws
{
    MPL::SaturationTuller saturation{"saturation", 0.0, 1.0, 1.0,
                                     0.8584073464102069, 1e-5, 0.0715, 0.0,
                                     1.4e8};
    MPL::RelPermGeneralizedPower relperm{"relative_permeability", 0.0, 0.0,
                                         1e-2, 1.0, 3.0};
    MPL::BishopsSaturationCutoff bishop{"bishops_effective_stress", 1.0};
    ParameterLib::SpatialPosition pos;

    double S(double const p_c) const
    {
        MPL::VariableArray v;
        v.capillary_pressure = p_c;
        return std::get<double>(saturation.value(v, pos, 0.0, 0.0));
    }
    double chi(double const S_L) const
    {
        MPL::VariableArray v;
        v.liquid_saturation = S_L;
        return std::get<double>(bishop.value(v, pos, 0.0, 0.0));
    }
    double krel(double const S_L) const
    {
        MPL::VariableArray v;
        v.liquid_saturation = S_L;
        return std::get<double>(relperm.value(v, pos, 0.0, 0.0));
    }
};

// p_edge of the header derivation.
double airEntryEdge()
{
    constexpr double F_gamma = 0.8584073464102069;
    constexpr double sigma = 0.0715;
    constexpr double A_n = 1.0;
    constexpr double L = 1e-5;
    double const C_T = 4.0 * F_gamma * sigma * sigma / (A_n * L * L);  // Pa^2
    return std::sqrt(C_T / (54.0 * std::log(2.0)));                    // Pa
}
}  // namespace

// UT-B1a. The latch rule L_new = Active AND (L_old OR chi_deck(S_L) == 1),
// complete truth table (2 x 2 x 2 inputs).
TEST(SaturationGateLatch, TruthTable)
{
    // columns: status_active, latch_old, chi_deck_at_S_L -> expected L_new
    struct Row
    {
        bool active;
        bool latch_old;
        double chi;
        bool expected;
        char const* what;
    };
    std::array<Row, 8> const rows{{
        {true, false, 1.0, true, "set: Active, saturated, not yet latched"},
        {true, false, 0.0, false, "stay off: Active but hydration front"},
        {true, true, 1.0, true, "keep: latched, still saturated"},
        {true, true, 0.0, true,
         "keep: latched, S_L < 1 again (load-induced suction)"},
        {false, false, 1.0, false, "never set off Active (chi 1)"},
        {false, false, 0.0, false, "never set off Active (chi 0)"},
        {false, true, 1.0, false, "release: latched point left Active"},
        {false, true, 0.0, false, "release: latched point left Active"},
    }};
    for (auto const& r : rows)
    {
        EXPECT_EQ(nextSaturatedLatch(r.active, r.latch_old, r.chi), r.expected)
            << r.what;
    }
}

// UT-B1b. A sequence of converged steps: set at the first saturated Active
// step, kept through a drop of S_L below 1 (the gated case), released when the
// point leaves Active, and NOT re-set by a later Active step with chi 0.
TEST(SaturationGateLatch, SequenceSetKeepReleaseResetRule)
{
    struct Step
    {
        bool active;
        double chi_deck;
        bool expected_latch;
    };
    std::vector<Step> const steps{
        {false, 0.0, false},  // hydration, interior
        {true, 0.0, false},   // front reaches the ceiling, S_L < 1: no latch
        {true, 0.0, false},   // still the front (unlatched points keep chi_deck)
        {true, 1.0, true},    // plateau reached: set
        {true, 1.0, true},    // keep
        {true, 0.0, true},    // keep through chi_deck = 0 (load-induced suction)
        {false, 0.0, false},  // release
        {true, 0.0, false},   // Active again with chi_deck 0: not re-latched
        {true, 1.0, true},    // saturated again: set
    };
    bool latch = false;
    for (std::size_t i = 0; i < steps.size(); ++i)
    {
        latch = nextSaturatedLatch(steps[i].active, latch, steps[i].chi_deck);
        EXPECT_EQ(latch, steps[i].expected_latch) << "step " << i;
    }
}

// UT-B1c. The gate predicates: the gate acts only at level != off AND latched
// at the end of the previous step AND Active in THIS iterate; k_rel is gated at
// bishop_relperm only (bishop is the labelled probe).
TEST(SaturationGateLatch, PredicatesByLevel)
{
    for (bool const latch_old : {false, true})
    {
        for (bool const active : {false, true})
        {
            bool const both = latch_old && active;
            EXPECT_FALSE(saturationGateActs(kOff, latch_old, active));
            EXPECT_EQ(saturationGateActs(kBishop, latch_old, active), both);
            EXPECT_EQ(saturationGateActs(kBishopRelperm, latch_old, active),
                      both);
            EXPECT_FALSE(
                saturationGateActsOnRelativePermeability(kOff, latch_old, active));
            EXPECT_FALSE(saturationGateActsOnRelativePermeability(
                kBishop, latch_old, active));
            EXPECT_EQ(saturationGateActsOnRelativePermeability(
                          kBishopRelperm, latch_old, active),
                      both);
        }
    }
}

// UT-B1d. Switch off or an ungated point: the deck values pass through
// bitwise; at a gated point the output does not depend on the deck values at
// all (S is a constant there: dchi/dS = 0, no derivative of the latch).
TEST(SaturationGateLatch, GatedFactorsPassThroughOrAreConstant)
{
    DeckLaws const laws;
    auto const chi_deck = [&](double S) { return laws.chi(S); };
    std::vector<BishopFactorValues> const deck_inputs{
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.3, 0.7, 12.5},  // arbitrary values: the gate must ignore them
    };
    for (auto const level : {kOff, kBishop, kBishopRelperm})
    {
        for (bool const latch_old : {false, true})
        {
            for (bool const active : {false, true})
            {
                bool const acts = saturationGateActs(level, latch_old, active);
                for (auto const& in : deck_inputs)
                {
                    auto const out = saturationGatedBishopFactors(
                        level, latch_old, active, in, chi_deck);
                    if (!acts)
                    {
                        EXPECT_EQ(out.chi, in.chi);
                        EXPECT_EQ(out.chi_prev, in.chi_prev);
                        EXPECT_EQ(out.dchi_dS_L, in.dchi_dS_L);
                    }
                    else
                    {
                        EXPECT_EQ(out.chi, laws.chi(1.0));
                        EXPECT_EQ(out.chi_prev, laws.chi(1.0));
                        EXPECT_EQ(out.dchi_dS_L, 0.0);
                    }
                }
            }
        }
    }
}

// UT-B2. One latched, Active integration point whose suction crosses the
// air-entry edge of the deck's Tuller retention (p_edge = 2165.57 Pa). The old
// rule has the jump chi: 1 -> 0 there (the cliff of the VII stall); the gated
// p_FR = -chi p_c is continuous, chi and k_rel have no jump, and the tangent
// of the gated function agrees with central differences to round-off.
TEST(SaturationGateEdge, ChiAndKrelContinuousAcrossTheAirEntryEdge)
{
    DeckLaws const laws;
    auto const chi_deck = [&](double S) { return laws.chi(S); };
    double const p_edge_analytic = airEntryEdge();

    // Locate the edge numerically on the deck's own retention (S_L == 1 for
    // p_c <= edge): bisection on the predicate S_L(p_c) == 1.
    double lo = 1.0e3;  // Pa, S_L == 1
    double hi = 1.0e4;  // Pa, S_L < 1
    ASSERT_EQ(laws.S(lo), 1.0);
    ASSERT_LT(laws.S(hi), 1.0);
    for (int i = 0; i < 200; ++i)
    {
        double const mid = 0.5 * (lo + hi);
        (laws.S(mid) == 1.0 ? lo : hi) = mid;
    }
    double const p_edge = lo;
    std::cout << "air-entry edge (numerical) = " << p_edge
              << " Pa, analytic = " << p_edge_analytic
              << " Pa, DESIGN_FIXES B.1 MEASURED 2165.57 Pa\n";
    EXPECT_NEAR(p_edge, p_edge_analytic, 1e-3);

    // The old rule: chi(S_L(p_c)) jumps 1 -> 0 and p_FR = -chi p_c jumps by
    // p_edge across the edge.
    double const below = std::nextafter(p_edge, 0.0);
    double const above = hi;  // first bisection point with S_L < 1
    EXPECT_EQ(laws.chi(laws.S(below)), 1.0);
    EXPECT_EQ(laws.chi(laws.S(above)), 0.0);
    double const old_jump = std::abs(-laws.chi(laws.S(above)) * above -
                                     -laws.chi(laws.S(below)) * below);
    EXPECT_NEAR(old_jump, p_edge, 1e-3 * p_edge);

    // The gated point: latch_old true, Active.
    auto const gated_at = [&](double const p_c)
    {
        double const S_L = laws.S(p_c);
        return saturationGatedBishopFactors(
            kBishopRelperm, true, true, {laws.chi(S_L), 1.0, 0.0}, chi_deck);
    };
    auto const gated_p_FR = [&](double const p_c)
    { return -gated_at(p_c).chi * p_c; };  // Pa

    // Sweep from 1 kPa to 1 MPa with a log grid that includes points one ulp
    // either side of the edge; nodes straddling the edge are the point.
    std::vector<double> p_grid;
    for (int i = 0; i <= 400; ++i)
    {
        p_grid.push_back(1.0e3 * std::pow(10.0, 3.0 * i / 400.0));  // Pa
    }
    p_grid.push_back(below);
    p_grid.push_back(p_edge);
    p_grid.push_back(above);
    std::sort(p_grid.begin(), p_grid.end());

    double const k_unit = laws.krel(1.0);
    EXPECT_EQ(k_unit, 1.0);  // the deck's law at S = 1 (enhancement 1.0)
    for (auto const p_c : p_grid)
    {
        auto const g = gated_at(p_c);
        EXPECT_EQ(g.chi, 1.0) << "p_c = " << p_c;
        EXPECT_EQ(g.chi_prev, 1.0);
        EXPECT_EQ(g.dchi_dS_L, 0.0);
    }
    // continuity: -p_FR is p_c itself, no jump anywhere on the grid (the
    // increments of p_FR equal the increments of -p_c to the rounding of the
    // subtraction).
    for (std::size_t i = 1; i < p_grid.size(); ++i)
    {
        double const d_p_FR = gated_p_FR(p_grid[i]) - gated_p_FR(p_grid[i - 1]);
        double const d_p_c = -(p_grid[i] - p_grid[i - 1]);
        EXPECT_EQ(d_p_FR, d_p_c) << "i = " << i;
    }

    // FD-consistency of the gated residual function p_FR(p_c): analytic
    // derivative d p_FR / d p_c = -(chi + dchi/dS * p_c * dS/dp_c) = -chi at a
    // gated point (dchi/dS = 0); central difference with h = 1e-3 p_c; the error
    // bound is the round-off of the difference quotient, 16 ulp * p_c / h (the
    // 16 is the unit of the existing KKT tests, kRoundoffUlps).
    constexpr double kRoundoffUlps = 16.0;
    for (auto const p_c : p_grid)
    {
        double const h = 1e-3 * p_c;
        double const fd =
            (gated_p_FR(p_c + h) - gated_p_FR(p_c - h)) / (2.0 * h);
        double const analytic = -gated_at(p_c).chi;  // [-]
        double const tol =
            kRoundoffUlps * std::numeric_limits<double>::epsilon() * p_c / h;
        EXPECT_NEAR(fd, analytic, tol) << "p_c = " << p_c;
    }
    // ... and the SAME check with the old rule fails at the edge: the old
    // function has no derivative there, the central difference across the cliff
    // is far from the analytic -chi (documents what the gate removes).
    {
        double const h = 1e-3 * p_edge;
        auto const old_p_FR = [&](double const p_c)
        { return -laws.chi(laws.S(p_c)) * p_c; };
        double const fd_old =
            (old_p_FR(p_edge + h) - old_p_FR(p_edge - h)) / (2.0 * h);
        double const analytic_old = -laws.chi(laws.S(p_edge));
        EXPECT_GT(std::abs(fd_old - analytic_old), 100.0);
    }

    // The old k_rel above the edge (documentation of the second carrier): it
    // sits on the deck's floor 1e-2 at 1 MPa, the gated value stays k_rel(1).
    EXPECT_EQ(laws.krel(laws.S(1.0e6)), 1e-2);
}

// UT-B3. The level bishop is the probe: chi is gated, k_rel is not.
TEST(SaturationGateEdge, ProbeLevelGatesChiOnly)
{
    EXPECT_TRUE(saturationGateActs(kBishop, true, true));
    EXPECT_FALSE(saturationGateActsOnRelativePermeability(kBishop, true, true));
    EXPECT_TRUE(saturationGateActsOnRelativePermeability(kBishopRelperm, true,
                                                         true));
}
