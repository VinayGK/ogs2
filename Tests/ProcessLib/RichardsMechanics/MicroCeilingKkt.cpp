// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// Unit tests of the KKT treatment of the micro water-content ceiling
// (potential_exchange micro_ceiling_treatment = kkt; branch
// dsm_mass_conservation_v3_kkt_ceiling_2026-09-30; NOT adopted).
//
// Batch 1 of DESIGN.md section 4.1 (record folder
// ~/ogs-models/scratch/2026-09-30_kkt_ceiling_impl/): UT-1, UT-2, UT-3, UT-4N,
// UT-5. Every test is a SUPPLEMENT: no existing test is edited. Expected
// values come from DERIVATION.md (D-n) and WEAK_FORMS.md (W-n) of that folder
// (in-file derivation) and from independent evaluations written in this file;
// none is tuned. The independent evaluations (oracle below) use the vdW and
// Maxwell helper functions as the DEFINITION of mu_lR (these are tested
// elsewhere); they do not call the KKT solver.
//
// Literals that are test-only proposals and need Vinay's approval (DESIGN.md
// D3, guardrail 1.2): kFdSafetyFactor, kRoundoffUlps, kDenseScanNodesPerDecade,
// the factor 1/2 in suctionForActiveState, the 0.1 fraction of the bracket
// first-order root in activeStrainSteps, eps_v = 0.55 of UT-4N block 3, the
// synthetic compressibility of UT-5, the hump constructions of UT-4N block 1.

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h"

using namespace ProcessLib::RichardsMechanics;

namespace
{
constexpr double eps_mach = std::numeric_limits<double>::epsilon();

// Named constants (flagged, see the header comment).
constexpr double kFdSafetyFactor = 10.0;  // C of DESIGN.md 4.2 (SCOPE.md l.154)
constexpr double kRoundoffUlps = 16.0;    // "to round-off", in units of eps_mach
constexpr int kDenseScanNodesPerDecade = 200;  // independent of N_dec = 8

enum class Regime
{
    CommittedBlock,  // (i)  a = 1.3, rho_l0 = 1300: rho_lR depends on n_l
    DeckEos          // (ii) a = 1e-16, rho_l0 = 100: rho_lR = 100 + rho_LR
};

enum class Context
{
    Production,  // (P) finite p_conf, finite K_d, nonzero eps_v
    Sentinel     // (S) NaN p_conf and K_d (film term self-disables)
};

PotentialExchangeParameters makeParameters(Regime const regime, double const s)
{
    // The parameter block of DSMMicroMacroMassStorageCoupledSolveResiduals
    // (DSMMicroMacroSingleIntegrationPoint.cpp), re-used as a block.
    PotentialExchangeParameters p;
    p.enabled = true;
    p.pressure_tolerance = 0.0;
    p.hamaker_constant = 6.0e-20;
    p.specific_surface = 1000.0;
    p.micro_solid_density_reference = 2650.0;
    p.micro_solid_volume_fraction_reference = 0.6;
    p.micro_liquid_density_reference = 1300.0;
    p.micro_liquid_density_a = 1.3;
    p.micro_liquid_density_b = 1.0;
    p.micro_potential_convention = MicroPotentialConvention::NegativeAttractive;
    p.local_nonlinear_solve_mode =
        LocalNonlinearSolveMode::ScalarReferenceMassStorage;
    p.macro_porosity_update_mode = MacroPorosityUpdateMode::ReferenceAdditiveRate;
    p.micro_solid_volume_fraction_mode =
        MicroSolidVolumeFractionMode::CurrentPorositySplit;
    p.initial_micro_water_content = 0.05;
    if (regime == Regime::DeckEos)
    {
        // The EOS block of every suite deck (VII deck l.239-241): rho_lR =
        // 100 + rho_LR = 1100 for rho_LR = 1000.
        p.micro_liquid_density_a = 1e-16;
        p.micro_liquid_density_reference = 100.0;
        p.micro_liquid_density_b = 1.0;
    }
    p.micro_ceiling_treatment = MicroCeilingTreatment::Kkt;
    p.micro_mass_strain_term_eulerian = (s < 0.0);
    return p;
}

// One local-solve case: everything the solver needs plus the history.
struct Case
{
    PotentialExchangeParameters params;
    Context context = Context::Production;
    double s = 1.0;  // strain-term sign of the parameters

    double dt = 100.0;       // block of the committed test
    double rho_LR = 1000.0;  // block
    double alpha_bar = 1.0e-9;
    double mu = 1.0e-3;
    double p_L = -1.0e6;

    // History (kept fixed when a perturbation is applied).
    double n_l_prev = 0.05;
    double phi_M_prev = 0.18;
    double phi_m_prev = 0.08;
    double phi_prev = 0.26;
    double rho_lR_prev = 0.0;
    double rho_l_prev = 0.0;

    // Kinematics.
    double eps_prev = 1.0e-3 - 1.0e-3;  // set by makers
    double d_eps = 1.0e-3;
    double phi = 0.26;  // total porosity at the end of the step

    // Elastic constants of the production-like context: VII deck l.475-476
    // (E flagged there as an uncited working value, Vinay-approved
    // 2026-06-15), nu = 0.3.
    double E = 52.0e6;
    double nu = 0.3;
    double biot = 1.0;  // MEASURED: deck biot_coefficient 1.0

    double K_d() const { return E / (3.0 * (1.0 - 2.0 * nu)); }
    double eps() const { return eps_prev + d_eps; }
    double alpha_M() const { return alpha_bar * rho_LR / mu; }
};

// Porosity law of PorosityFromMassBalance with alpha = 1 and beta_SR = 0:
// phi = (phi_prev + d_eps)/(1 + d_eps).
double porosityLaw(double const phi_prev, double const d_eps)
{
    return (phi_prev + d_eps) / (1.0 + d_eps);
}

PotentialExchangeLocalSolveContext makeContext(Case const& c)
{
    PotentialExchangeLocalSolveContext ctx{
        .phi = c.phi,
        .phi_M_prev = c.phi_M_prev,
        .phi_m_prev = c.phi_m_prev,
        .volumetric_strain = c.eps(),
        .volumetric_strain_prev = c.eps_prev,
    };
    if (c.context == Context::Production)
    {
        ctx.confining_pressure_p_conf = -c.K_d() * c.eps();  // drained line
        ctx.drained_bulk_modulus = c.K_d();
        ctx.biot_coefficient = c.biot;
    }
    return ctx;
}

// History of the committed fixture: rho_l_prev = phi_m_prev * rho_lR_prev.
void setHistory(Case& c)
{
    PotentialExchangeLocalSolveContext const ctx = makeContext(c);
    auto const prev = computePreviousMicroLiquidDensity(c.n_l_prev, c.rho_LR,
                                                        ctx, c.params);
    c.rho_lR_prev = prev.rho_lR;
    c.rho_l_prev = c.phi_m_prev * prev.rho_lR;
}

// The state of the committed block: micro far below the ceiling.
Case makeBelowCeilingCase(Regime const regime, Context const context,
                          double const s)
{
    Case c;
    c.params = makeParameters(regime, s);
    c.context = context;
    c.s = s;
    c.p_L = -1.0e6;  // committed DSMMicroMacroScalarMassStorageLocalSolve
                     // ReferencePath: "1 MPa suction"
    c.n_l_prev = 0.05;
    c.phi_M_prev = 0.18;
    c.phi_m_prev = 0.08;
    c.phi_prev = 0.26;
    c.phi = 0.26;  // block value (explicit, as the committed test)
    c.d_eps = 1.0e-3;
    // (S): the committed test's strains (0 -> 1e-3). (P): MEASURED mean total
    // strain of the free-swelling sample, +0.205 (diag_V1_VII_offset_
    // 2026-09-30/README.md), reached with the same increment.
    c.eps_prev = (context == Context::Production) ? 0.205 - c.d_eps : 0.0;
    setHistory(c);
    return c;
}

// The micro at its ceiling: n_l_prev = phi_prev, phi_M_prev = 0, total
// strain eps0 (P), increment d_eps (porosity law).
Case makeAtCeilingCase(Regime const regime, double const s,
                       double const eps0 = 0.205, double const d_eps = 0.0)
{
    Case c;
    c.params = makeParameters(regime, s);
    c.context = Context::Production;
    c.s = s;
    c.phi_prev = 0.26;  // block value
    c.n_l_prev = c.phi_prev;
    c.phi_M_prev = 0.0;
    c.phi_m_prev = c.phi_prev;
    c.eps_prev = eps0;
    c.d_eps = 0.0;
    c.phi = c.phi_prev;
    setHistory(c);  // the history belongs to the state at d_eps = 0
    c.d_eps = d_eps;
    c.phi = porosityLaw(c.phi_prev, d_eps);
    return c;
}

struct SolveResult
{
    MicroCeilingKktSolveData kkt;
    MicroMacroMassStorageCoupledSolveData base;
};

YoungLaplaceMacroPotentialData macroPotential(Case const& c)
{
    return computeYoungLaplaceMacroPotential(c.p_L, c.rho_LR,
                                             c.params.pressure_tolerance);
}

MicroCeilingKktSolveData solveKkt(Case const& c)
{
    return solveReferenceMassStorageKktState(
        c.n_l_prev, c.rho_l_prev, c.rho_lR_prev, c.dt, c.rho_LR, c.alpha_bar,
        c.mu, macroPotential(c), makeContext(c), c.params);
}

MicroMacroMassStorageCoupledSolveData solveBase(Case const& c)
{
    return solveReferenceMassStorageCoupledState(
        c.n_l_prev, c.rho_l_prev, c.rho_lR_prev, c.dt, c.rho_LR, c.alpha_bar,
        c.mu, macroPotential(c), makeContext(c), c.params);
}

// ── Independent evaluation of the local problem (the oracle) ────────────────
// f(n) = R_m(n, rho_lR_EOS(n), lambda = 0) as in DERIVATION.md D-2.3, written
// out here from the definitions (not through the solver's lambda).
struct Oracle
{
    Case const& c;
    PotentialExchangeLocalSolveContext ctx;
    double mu_LR;

    explicit Oracle(Case const& c_)
        : c(c_), ctx(makeContext(c_)), mu_LR(macroPotential(c_).mu_LR)
    {
    }

    double nS(double const n) const
    {
        return computeActiveMicroSolidVolumeFraction(n, ctx, c.params);
    }
    double rho_lR_eos(double const n) const
    {
        return computeReducedMicroLiquidDensity(n, c.rho_LR, nS(n), c.params)
            .rho_lR;
    }
    // mu_lR = mu_vdW + mu_mech (production) or mu_vdW (sentinel); [J/kg].
    double mu_lR(double const n, double const rho_lR) const
    {
        auto const vdw = computeVanDerWaalsMicroPotential(
            n, rho_lR, nS(n), c.params.micro_solid_density_reference,
            c.params.hamaker_constant, c.params.specific_surface,
            microPotentialSignFactorFromParameters(c.params),
            effectiveAugmentationPrefactor(c.params, ctx.phi),
            c.params.potential_augmentation_exponent, 0.0,
            c.params.micro_water_content_floor);
        if (c.context == Context::Sentinel)
        {
            return vdw.mu_lR;
        }
        double const Pi = -rho_lR * vdw.mu_lR;
        double const dPi = -rho_lR * vdw.dmu_lR_dnl;
        double const d2Pi = -rho_lR * vdw.d2mu_lR_dnl2;
        auto const mech = computeIntegrableMechanicalMicroPotential(
            Pi, dPi, d2Pi, n, c.eps(), c.biot, c.K_d(), rho_lR);
        return vdw.mu_lR + mech.mu_lR_mech;
    }
    double rho_l(double const n, double const rho_lR) const
    {
        double const phi_cs = std::clamp(ctx.phi, 0.0, 1.0 - 1e-12);
        return (1.0 - phi_cs) / std::max(1e-12, 1.0 - n) * n * rho_lR;
    }
    double rhohat_pot(double const n) const
    {
        return c.alpha_M() * (mu_LR - mu_lR(n, rho_lR_eos(n)));
    }
    double eps_dot() const { return c.d_eps / c.dt; }
    // Booked storage rate S_s(n) = (c_s rho_l - rho_l_prev)/dt.
    double S_s(double const n) const
    {
        double const cs = 1.0 - c.s * c.dt * eps_dot();
        return (cs * rho_l(n, rho_lR_eos(n)) - c.rho_l_prev) / c.dt;
    }
    double f(double const n) const
    {
        double const rl = rho_lR_eos(n);
        double const rho = rho_l(n, rl);
        return rho - c.rho_l_prev - c.dt * c.alpha_M() * (mu_LR - mu_lR(n, rl)) -
               c.s * c.dt * rho * eps_dot();
    }
};

// True iff f < 0 at every node of a log-uniform grid on [n_floor, n_max]
// (kDenseScanNodesPerDecade nodes per decade, independent of N_dec).
bool denseScanAllNegative(Oracle const& o, double const n_floor,
                          double const n_max, double* max_f = nullptr)
{
    int const cells = static_cast<int>(std::ceil(
        kDenseScanNodesPerDecade * std::log10(n_max / n_floor)));
    bool all_negative = true;
    double fmax = -std::numeric_limits<double>::infinity();
    for (int j = 0; j <= cells; ++j)
    {
        double const n =
            (j == cells) ? n_max
                         : n_floor * std::pow(n_max / n_floor,
                                              static_cast<double>(j) / cells);
        double const fv = o.f(n);
        fmax = std::max(fmax, fv);
        if (!(fv < 0.0))
        {
            all_negative = false;
        }
    }
    if (max_f)
    {
        *max_f = fmax;
    }
    return all_negative;
}

// Suction of the active-state cases: mu_LR := mu_lR(n_max)/2, i.e. half the
// micro potential at the ceiling, so that mu_LR - mu_lR(n_max) =
// -mu_lR(n_max)/2 > 0 whenever mu_lR(n_max) < 0 (water wants IN while the
// micro is full). Test-only construction (the factor 1/2 needs approval).
void setSuctionForActiveState(Case& c)
{
    // The oracle needs the mu_LR of the case only through c.p_L; mu_lR does
    // not depend on it, so a first oracle with any p_L gives mu_lR(n_max).
    Oracle const o(c);
    double const n_max = c.phi;
    double const mu_at_ceiling = o.mu_lR(n_max, o.rho_lR_eos(n_max));
    c.p_L = 0.5 * c.rho_LR * mu_at_ceiling;  // unsaturated: p_L < 0
}

double dtFirstOrderRoot(Case const& c)
{
    // d_eps at which rho_lR d_eps/dt = rhohat_pot (first-order estimate of
    // the strain at which S_-1 = rho_lR d_eps/dt crosses rhohat_pot, D-2.4).
    Oracle const o(c);
    double const n_max = c.phi;
    return c.dt * o.rhohat_pot(n_max) / o.rho_lR_eos(n_max);
}

// Relative round-off tolerance for bookkeeping identities.
double roundoffTolerance(double const scale)
{
    return kRoundoffUlps * eps_mach * scale;
}

// ── Central differences (DESIGN.md 4.2) ─────────────────────────────────────
double fdStep(double const x, double const x_scale, double const factor = 1.0)
{
    return factor * std::cbrt(eps_mach) * std::max(std::abs(x), x_scale);
}

// Tolerance C * eps^(2/3) * (|analytic| + scale_term).
double fdTolerance(double const analytic, double const scale_term)
{
    return kFdSafetyFactor * std::pow(eps_mach, 2.0 / 3.0) *
           (std::abs(analytic) + scale_term);
}

std::vector<Regime> const regimes = {Regime::CommittedBlock, Regime::DeckEos};
std::array<double, 2> const signs = {+1.0, -1.0};
}  // namespace

// ── UT-1 ────────────────────────────────────────────────────────────────────
// Below the ceiling (f(n_max) >= 0) the KKT solve is the base solve, bit for
// bit, in both EOS regimes, both strain signs and both context variants.
TEST(RichardsMechanics, DSMMicroCeilingKktBelowCeilingIsBitwiseClamp)
{
    for (auto const regime : regimes)
    {
        for (auto const context : {Context::Production, Context::Sentinel})
        {
            for (double const s : signs)
            {
                SCOPED_TRACE(testing::Message()
                             << "regime " << static_cast<int>(regime)
                             << " context " << static_cast<int>(context)
                             << " s " << s);
                Case const c = makeBelowCeilingCase(regime, context, s);
                Oracle const o(c);
                double const n_max = c.phi;
                // Premise, from the independent evaluation (no skip).
                ASSERT_GE(o.f(n_max), 0.0)
                    << "fixture is not below the ceiling";

                auto const kkt = solveKkt(c);
                auto const base = solveBase(c);

                EXPECT_EQ(kkt.status, MicroCeilingKktStatus::Interior);
                EXPECT_EQ(kkt.local.n_l, base.n_l);
                EXPECT_EQ(kkt.local.rho_lR, base.rho_lR);
                EXPECT_EQ(kkt.local.exchange.rho_l_hat, base.exchange.rho_l_hat);
                EXPECT_EQ(kkt.local.converged, base.converged);
                EXPECT_EQ(kkt.multiplier, 0.0);
                EXPECT_EQ(kkt.rejected_exchange, 0.0);
                EXPECT_EQ(kkt.exchange_received, kkt.local.exchange.rho_l_hat);
            }
        }
    }
}

// ── UT-2 ────────────────────────────────────────────────────────────────────
// Micro at the ceiling, suction draws water in: the active branch. Complemen-
// tarity, multiplier definition and mass bookkeeping from the returned state.
TEST(RichardsMechanics,
     DSMMicroCeilingKktActiveBranchComplementarityAndBookkeeping)
{
    for (auto const regime : regimes)
    {
        for (double const s : signs)
        {
            SCOPED_TRACE(testing::Message() << "regime "
                                            << static_cast<int>(regime)
                                            << " s " << s);
            Case c = makeAtCeilingCase(regime, s);
            setSuctionForActiveState(c);
            Oracle const o(c);
            constexpr double n_floor = 1e-16;  // as the solver (existing literal)
            double const n_max = c.phi;

            // Premises, from independent evaluations (no skip).
            double const rhohat_pot_wall = o.rhohat_pot(n_max);
            ASSERT_GT(rhohat_pot_wall, 0.0);
            double max_f = 0.0;
            ASSERT_TRUE(denseScanAllNegative(o, n_floor, n_max, &max_f))
                << "f >= 0 somewhere on [n_floor, n_max], max f = " << max_f;

            auto const kkt = solveKkt(c);
            ASSERT_EQ(kkt.status, MicroCeilingKktStatus::Active);
            // Memoryless: the solver has no state argument, so two calls with
            // identical inputs return identical output (no branch flip under
            // a fixed residual; the sweep over p_L of UT-4 is held, DESIGN D2).
            {
                auto const again = solveKkt(c);
                EXPECT_EQ(again.status, kkt.status);
                EXPECT_EQ(again.local.n_l, kkt.local.n_l);
                EXPECT_EQ(again.local.rho_lR, kkt.local.rho_lR);
                EXPECT_EQ(again.multiplier, kkt.multiplier);
                EXPECT_EQ(again.exchange_received, kkt.exchange_received);
                EXPECT_EQ(again.rejected_exchange, kkt.rejected_exchange);
            }
            EXPECT_EQ(kkt.local.n_l, n_max);
            EXPECT_EQ(kkt.n_max, n_max);
            EXPECT_GE(kkt.multiplier, 0.0);

            // lambda = rho_lR (mu_LR - mu_lR(n_max)) (S_s = 0 here, D-2.4).
            double const rho_lR_wall = o.rho_lR_eos(n_max);
            double const lambda_def =
                rho_lR_wall * (o.mu_LR - o.mu_lR(n_max, rho_lR_wall));
            EXPECT_NEAR(kkt.multiplier, lambda_def,
                        roundoffTolerance(std::abs(lambda_def)) +
                            roundoffTolerance(rho_lR_wall *
                                              std::abs(o.mu_LR)));

            // The micro takes what the storage balance gives: nothing here.
            double const scale_rate = rhohat_pot_wall;
            EXPECT_NEAR(kkt.exchange_received, 0.0,
                        roundoffTolerance(scale_rate));
            EXPECT_NEAR(kkt.rejected_exchange, rhohat_pot_wall,
                        roundoffTolerance(scale_rate));
            EXPECT_NEAR(kkt.local.exchange.rho_l_hat, rhohat_pot_wall,
                        roundoffTolerance(scale_rate));

            // Complementarity (exact: n_l == n_max).
            EXPECT_EQ(kkt.multiplier * (n_max - kkt.local.n_l), 0.0);

            // Mass bookkeeping from the returned state alone (independent of
            // the solver's S_s): c_s rho_l - rho_l_prev == dt * received and
            // received + rejected == rhohat_pot.
            double const cs = 1.0 - c.s * c.dt * o.eps_dot();
            double const rho_l_back = o.rho_l(kkt.local.n_l, kkt.local.rho_lR);
            double const lhs = cs * rho_l_back - c.rho_l_prev;
            double const rhs = c.dt * kkt.exchange_received;
            EXPECT_NEAR(lhs, rhs, roundoffTolerance(rho_l_back) +
                                      roundoffTolerance(c.rho_l_prev));
            EXPECT_NEAR(kkt.exchange_received + kkt.rejected_exchange,
                        kkt.local.exchange.rho_l_hat,
                        roundoffTolerance(scale_rate));
        }
    }
}

// ── UT-3 ────────────────────────────────────────────────────────────────────
// Bisection on the strain increment between an active and an inactive point:
// the received exchange is continuous at the kink, for both strain signs.
//
// Scenario (test-only construction, flagged): the kink is placed at the strain
// step of the Model VII top probe, d* = 0.230 - 0.205 (MEASURED live-K top
// probe and mean free-swelling strain, DESIGN.md D-2.4 / UT-8 line,
// diag_V1_VII_offset_2026-09-30/README.md), by choosing the suction so that
// f(n_max; d*) = 0 exactly:  mu_LR = mu_lR(n_max(d*), eps0 + d*) + S_s(d*)/alpha_M.
// The direction (sigma = +1 or -1) is the one in which the independent f turns
// from negative (active) at d = 0 to positive at 2 d*.
namespace
{
// Sets c.p_L so that f(n_max; d*) = 0 (DERIVED, see above). Returns false if
// the resulting p_L is not unsaturated (p_L >= 0).
bool setSuctionForKinkAt(Case& c, double const d_star)
{
    Case c_star = c;
    c_star.d_eps = d_star;
    c_star.phi = porosityLaw(c_star.phi_prev, d_star);
    Oracle const o(c_star);
    double const n_max = c_star.phi;
    double const mu_lR_star = o.mu_lR(n_max, o.rho_lR_eos(n_max));
    double const mu_LR = mu_lR_star + o.S_s(n_max) / c.alpha_M();
    c.p_L = c.rho_LR * mu_LR;
    return c.p_L < 0.0;
}
}  // namespace

TEST(RichardsMechanics,
     DSMMicroCeilingKktLeavingTheCeilingAndContinuityAtTheKink)
{
    double const d_star_magnitude = 0.230 - 0.205;  // see above
    int const grid_cells = 20;  // coarse scan of d in [0, 2 d*] for a kink

    for (auto const regime : regimes)
    {
        for (double const s : signs)
        {
            SCOPED_TRACE(testing::Message() << "regime "
                                            << static_cast<int>(regime)
                                            << " s " << s);
            Case base;
            double lo = 0.0;
            double hi = 0.0;
            double grid_d_max = 0.0;
            bool found = false;
            for (double const sigma : {+1.0, -1.0})
            {
                base = makeAtCeilingCase(regime, s);
                double const d_star = sigma * d_star_magnitude;
                grid_d_max = 2.0 * d_star;
                if (!setSuctionForKinkAt(base, d_star))
                {
                    continue;
                }
                // First adjacent pair (Active, Interior) on the grid.
                double prev_d = 0.0;
                auto prev_status = MicroCeilingKktStatus::Interior;
                for (int j = 0; j <= grid_cells && !found; ++j)
                {
                    double const d = 2.0 * d_star * j / grid_cells;
                    Case c = base;
                    c.d_eps = d;
                    c.phi = porosityLaw(c.phi_prev, d);
                    auto const status = solveKkt(c).status;
                    if (j == 0 && status != MicroCeilingKktStatus::Active)
                    {
                        break;  // not active at d = 0: the other direction
                    }
                    if (j > 0 &&
                        prev_status == MicroCeilingKktStatus::Active &&
                        status == MicroCeilingKktStatus::Interior)
                    {
                        lo = prev_d;
                        hi = d;
                        found = true;
                    }
                    prev_d = d;
                    prev_status = status;
                }
                if (found)
                {
                    break;
                }
            }
            ASSERT_TRUE(found) << "no (Active, Interior) pair in either "
                                  "direction";
            auto const at_d = [&](double const d_eps)
            {
                Case c = base;
                c.d_eps = d_eps;
                c.phi = porosityLaw(c.phi_prev, d_eps);
                return c;
            };
            // The Interior grid node with the largest independent f(n_max)
            // (the root then sits well below the wall).
            double d_inside = hi;
            {
                double best = -std::numeric_limits<double>::infinity();
                for (int j = 0; j <= grid_cells; ++j)
                {
                    double const d = grid_d_max * j / grid_cells;
                    Case const c = [&]
                    {
                        Case cc = base;
                        cc.d_eps = d;
                        cc.phi = porosityLaw(cc.phi_prev, d);
                        return cc;
                    }();
                    if (solveKkt(c).status != MicroCeilingKktStatus::Interior)
                    {
                        continue;
                    }
                    double const margin = Oracle(c).f(c.phi);
                    if (margin > best)
                    {
                        best = margin;
                        d_inside = d;
                    }
                }
            }

            // Bisection between lo (active) and hi (inactive) down to
            // adjacent doubles (count-limited).
            for (int it = 0; it < 200; ++it)
            {
                double const mid = 0.5 * (lo + hi);
                if (mid == lo || mid == hi)
                {
                    break;
                }
                if (solveKkt(at_d(mid)).status ==
                    MicroCeilingKktStatus::Active)
                {
                    lo = mid;
                }
                else
                {
                    hi = mid;
                }
            }
            double const h = std::abs(hi - lo);

            auto const plus = solveKkt(at_d(hi));   // inactive side of the kink
            auto const minus = solveKkt(at_d(lo));  // active side
            EXPECT_EQ(plus.status, MicroCeilingKktStatus::Interior);
            EXPECT_EQ(minus.status, MicroCeilingKktStatus::Active);
            EXPECT_GE(minus.multiplier, 0.0);
            EXPECT_EQ(plus.multiplier, 0.0);

            // The status switch coincides with the sign change of the
            // INDEPENDENT f(n_max) (to round-off: the two evaluations of f
            // may differ in the last bits at the bracket).
            Case const c_lo = at_d(lo);
            Case const c_hi = at_d(hi);
            Oracle const o_lo(c_lo);
            Oracle const o_hi(c_hi);
            double const f_scale =
                kRoundoffUlps * eps_mach *
                (std::abs(o_lo.rho_l(c_lo.phi, o_lo.rho_lR_eos(c_lo.phi))) +
                 c_lo.rho_l_prev +
                 c_lo.dt * std::abs(o_lo.rhohat_pot(c_lo.phi)));
            EXPECT_LE(o_lo.f(c_lo.phi), f_scale);
            EXPECT_GE(o_hi.f(c_hi.phi), -f_scale);

            // Continuity at the bracketed kink: |rhohat(+) - rhohat(-)| <=
            // L h + round-off, L = max(rho_lR/dt, alpha_M |d mu_lR/d eps along
            // n = phi(eps)|) (D-5.1, D-5.2; the inactive side carries the
            // Maxwell slope). For either sign s (D-4.3).
            double const rho_lR_lo = o_lo.rho_lR_eos(c_lo.phi);
            auto const mu_path = [&](double const d_eps)
            {
                Case const c = at_d(d_eps);
                Oracle const o(c);
                return o.mu_lR(c.phi, o.rho_lR_eos(c.phi));
            };
            double const h_fd = fdStep(lo, 1.0);
            double const dmu_path =
                (mu_path(lo + h_fd) - mu_path(lo - h_fd)) / (2.0 * h_fd);
            double const L = std::max(rho_lR_lo / c_lo.dt,
                                      c_lo.alpha_M() * std::abs(dmu_path));
            double const jump =
                std::abs(plus.exchange_received - minus.exchange_received);
            double const scale_rate = std::abs(minus.exchange_received) +
                                      std::abs(plus.exchange_received) +
                                      rho_lR_lo * c_lo.phi / c_lo.dt;
            EXPECT_LE(jump, L * h + roundoffTolerance(scale_rate))
                << "jump " << jump << " L*h " << L * h << " h " << h;
            GTEST_LOG_(INFO) << "MEASURED UT-3 regime " << static_cast<int>(regime)
                             << " s " << s << ": kink at d_eps = " << lo
                             << ", bracket h = " << h << ", |jump| = " << jump
                             << ", L*h = " << L * h << ", round-off allowance "
                             << roundoffTolerance(scale_rate);

            // At the Interior grid node: interior state, no multiplier,
            // independent residual of the 2x2 base solve small (the committed
            // test's measure).
            Case const c_in = at_d(d_inside);
            auto const in = solveKkt(c_in);
            Oracle const o_in(c_in);
            ASSERT_EQ(in.status, MicroCeilingKktStatus::Interior);
            EXPECT_EQ(in.multiplier, 0.0);
            EXPECT_LT(in.local.n_l, c_in.phi);
            double const rho_l_in = o_in.rho_l(in.local.n_l, in.local.rho_lR);
            double const R_m =
                rho_l_in - c_in.rho_l_prev -
                c_in.dt * in.local.exchange.rho_l_hat -
                c_in.s * c_in.dt * rho_l_in * o_in.eps_dot();
            double const R_rho =
                in.local.rho_lR - o_in.rho_lR_eos(in.local.n_l);
            double const residual_norm =
                std::abs(R_m) / std::max(1.0, std::abs(c_in.rho_l_prev)) +
                std::abs(R_rho) / std::max(1.0, std::abs(in.local.rho_lR));
            EXPECT_EQ(in.exchange_received, in.local.exchange.rho_l_hat);
            if (regime == Regime::DeckEos)
            {
                // The production EOS (rho_lR constant): the base solve
                // converges; its residual is small (precedent: committed test
                // DSMMicroMacroMassStorageCoupledSolveResiduals).
                EXPECT_TRUE(in.local.converged);
                EXPECT_LE(residual_norm, 1e-8);
            }
            else if (!in.local.converged)
            {
                // MEASURED limitation of the BASE solve (not of the KKT rule):
                // with rho_lR(n_l) (regime i) the 1D predictor of
                // solveReferenceMassStorageCoupledState does not converge for
                // a root close below the ceiling (its tangent omits
                // dmu_lR/drho_lR * drho_lR/dn) and the base returns it with
                // converged = false. The KKT Interior branch returns the base
                // result bit for bit, so the residual bound is not asserted
                // here; the state is recorded.
                RecordProperty("base_predictor_not_converged_regime_i_s" +
                                   std::to_string(static_cast<int>(s)),
                               std::to_string(residual_norm));
                GTEST_LOG_(WARNING)
                    << "UT-3 regime (i) s " << s
                    << ": base solve not converged at the Interior node, "
                       "normalised residual "
                    << residual_norm << " (pre-existing base behaviour)";
            }
        }
    }
}

// ── UT-4N ───────────────────────────────────────────────────────────────────
// The bracketed scan: non-monotone f is detected (also a thin hump the nodes
// miss), a monotone f raises no alarm, the verdict does not depend on the
// amplitude of f; the production fixture is Active with an independent dense
// scan in agreement; above the cubic-core threshold the premise check fires
// and the result is the base solve.
TEST(RichardsMechanics, DSMMicroCeilingKktScanDetectsNonMonotoneAndFallsBack)
{
    using Status = MicroCeilingKktStatus;
    int const nodes_per_decade = 8;  // the default (proposal N_dec)
    double const n_floor = 1e-6;

    // Block 1: synthetic functions (theory_fix_check.py Part B construction):
    // f = k (n-a)(n-b)(n-c) with roots a < b < c, n_max in (b, c): f(n_floor)
    // < 0, hump f > 0 on (a, b), f(n_max) < 0.
    struct Hump
    {
        double a, b, c, n_max;
        char const* name;
    };
    std::array<Hump, 3> const humps = {
        Hump{1e-4, 3e-4, 1e-2, 5e-3, "wide hump"},
        Hump{1e-3, 1e-3 * (1.0 + 1e-3), 1e-2, 5e-3, "thin hump"},
        Hump{1e-5, 1e-5 * 1.5, 1e-1, 5e-2, "low hump"}};
    for (auto const& h : humps)
    {
        for (double const k : {1e-12, 1.0, 1e12})
        {
            SCOPED_TRACE(testing::Message() << h.name << " k " << k);
            auto const f = [&](double const n)
            { return k * (n - h.a) * (n - h.b) * (n - h.c); };
            auto const fn = [&](double const n)
            {
                return k * ((n - h.b) * (n - h.c) + (n - h.a) * (n - h.c) +
                            (n - h.a) * (n - h.b));
            };
            ASSERT_LT(f(n_floor), 0.0);
            ASSERT_LT(f(h.n_max), 0.0);  // the bare rule would say "active"
            EXPECT_EQ(scanForInteriorCeilingRoot(f, fn, n_floor, h.n_max,
                                                 nodes_per_decade),
                      Status::NonMonotone);
        }
    }
    // Monotone increasing f with f(n_max) < 0: no alarm.
    for (double const k : {1e-12, 1.0, 1e12})
    {
        double const a = 1e-3;
        auto const f = [&](double const n) { return k * (n - 3.0 * a); };
        auto const fn = [&](double const) { return k; };
        EXPECT_EQ(scanForInteriorCeilingRoot(f, fn, n_floor, 0.5 * a,
                                             nodes_per_decade),
                  Status::Active);
    }
    // Premise violations: f(n_floor) >= 0 and a non-finite value.
    {
        auto const f = [](double const n) { return 1e-3 - n; };  // f(n_floor) > 0
        auto const fn = [](double const) { return -1.0; };
        EXPECT_EQ(scanForInteriorCeilingRoot(f, fn, n_floor, 1e-5,
                                             nodes_per_decade),
                  Status::PremiseViolated);
        auto const f_nan = [](double const) {
            return std::numeric_limits<double>::quiet_NaN();
        };
        EXPECT_EQ(scanForInteriorCeilingRoot(f_nan, fn, n_floor, 1e-5,
                                             nodes_per_decade),
                  Status::PremiseViolated);
    }

    for (auto const regime : regimes)
    {
        // Block 2: production fixture eps_v = +0.205 (MEASURED mean total
        // strain): Active, and the independent dense scan agrees.
        {
            SCOPED_TRACE(testing::Message() << "block 2, regime "
                                            << static_cast<int>(regime));
            Case c = makeAtCeilingCase(regime, +1.0, 0.205);
            setSuctionForActiveState(c);
            Oracle const o(c);
            double max_f = 0.0;
            EXPECT_TRUE(denseScanAllNegative(o, 1e-16, c.phi, &max_f))
                << "max f " << max_f;
            EXPECT_EQ(solveKkt(c).status, Status::Active);
        }
        // Block 3: eps_v above the cubic-core threshold 1/2 (test-only, no
        // physical claim): f(n_floor) > 0, the premise check fires, and the
        // returned state is the base solve, bit for bit.
        {
            SCOPED_TRACE(testing::Message() << "block 3, regime "
                                            << static_cast<int>(regime));
            Case c = makeAtCeilingCase(regime, +1.0, 0.55);
            // p_L as in block 2 (half the micro potential at the ceiling at
            // THIS strain): the candidate must be active (f(n_max) < 0).
            setSuctionForActiveState(c);
            Oracle const o(c);
            ASSERT_LT(o.f(c.phi), 0.0);
            ASSERT_GT(o.f(1e-16), 0.0);
            auto const kkt = solveKkt(c);
            auto const base = solveBase(c);
            EXPECT_EQ(kkt.status, Status::PremiseViolated);
            EXPECT_EQ(kkt.local.n_l, base.n_l);
            EXPECT_EQ(kkt.local.rho_lR, base.rho_lR);
            EXPECT_EQ(kkt.local.exchange.rho_l_hat, base.exchange.rho_l_hat);
            EXPECT_EQ(kkt.multiplier, 0.0);
            EXPECT_EQ(kkt.exchange_received, base.exchange.rho_l_hat);
        }
    }
}

// ── UT-5 ────────────────────────────────────────────────────────────────────
// The active-branch tangents of the received exchange against central
// differences of the SOLVER (not of the closed form), with the step rule of
// DESIGN.md 4.2 and a step-halving scan.
TEST(RichardsMechanics, DSMMicroCeilingKktActiveTangentsVersusCentralDifference)
{
    double ut5_worst_ratio = 0.0;  // max over all cases of (FD error)/(tolerance)
    int ut5_cases = 0;
    for (auto const regime : regimes)
    {
        for (double const s : signs)
        {
            for (bool const compressible : {false, true})
            {
                for (double const frac : {0.0, 0.1})
                {
                    SCOPED_TRACE(testing::Message()
                                 << "regime " << static_cast<int>(regime)
                                 << " s " << s << " compressible "
                                 << compressible << " frac " << frac);
                    ++ut5_cases;
                    Case c0 = makeAtCeilingCase(regime, s);
                    setSuctionForActiveState(c0);
                    // d_eps: 0 or a fixed fraction of the first-order root
                    // (still active, asserted below).
                    double const d_eps_base = frac * dtFirstOrderRoot(c0);
                    Case c = c0;
                    c.d_eps = d_eps_base;
                    c.phi = porosityLaw(c.phi_prev, d_eps_base);
                    // Synthetic compressibility beta_LR = 1/(1e3 |p_L|)
                    // (test-only, needs approval): rho_LR(p) = rho_LR0 (1 +
                    // beta (p - p_L0)).
                    double const beta_LR =
                        compressible ? 1.0 / (1.0e3 * std::abs(c.p_L)) : 0.0;
                    double const rho_LR0 = c.rho_LR;
                    double const p_L0 = c.p_L;

                    auto const solve_at = [&](double const dp, double const de)
                    {
                        Case cp = c;
                        cp.p_L = p_L0 + dp;
                        cp.rho_LR = rho_LR0 * (1.0 + beta_LR * dp);
                        cp.d_eps = d_eps_base + de;
                        cp.phi = porosityLaw(cp.phi_prev, cp.d_eps);
                        auto const r = solveKkt(cp);
                        EXPECT_EQ(r.status, MicroCeilingKktStatus::Active);
                        return r;
                    };

                    auto const base = solve_at(0.0, 0.0);
                    ASSERT_EQ(base.status, MicroCeilingKktStatus::Active);

                    // Analytic tangents.
                    Oracle const o(c);
                    double const n_max = c.phi;
                    auto const eos = computeReducedMicroLiquidDensity(
                        n_max, c.rho_LR, o.nS(n_max), c.params);
                    double const dphi = porosityDerivativeWrtVolumetricStrain(
                        1.0, c.phi, c.phi_prev, c.d_eps);
                    auto const tan = computeCeilingKktActiveExchangeTangents(
                        c.s, c.d_eps, c.phi, base.local.rho_lR, eos.drho_lR_dnl,
                        eos.drho_lR_drho_LR, rho_LR0 * beta_LR, dphi, c.dt);

                    // Scale of the difference S_s: |c_s rho_l| + rho_l_prev.
                    double const cs = 1.0 - c.s * c.d_eps;
                    double const rho_l_wall = c.phi * base.local.rho_lR;
                    double const scale_S =
                        (std::abs(cs * rho_l_wall) + c.rho_l_prev) / c.dt;

                    // d/dp_L, with a step-halving scan.
                    double const x_p = std::abs(c.p_L);
                    double const tol_p =
                        fdTolerance(tan.drhohat_dpL, scale_S / x_p);
                    std::vector<double> errors_p;
                    for (double const factor : {0.5, 1.0, 2.0})
                    {
                        double const hh = fdStep(c.p_L, x_p, factor);
                        double const fd =
                            (solve_at(+hh, 0.0).exchange_received -
                             solve_at(-hh, 0.0).exchange_received) /
                            (2.0 * hh);
                        errors_p.push_back(std::abs(fd - tan.drhohat_dpL));
                    }
                    for (double const e : errors_p)
                    {
                        EXPECT_LE(e, tol_p) << "dp_L";
                        ut5_worst_ratio = std::max(ut5_worst_ratio, e / tol_p);
                    }

                    // d/deps_v.
                    double const x_e = std::abs(c.eps());
                    double const tol_e =
                        fdTolerance(tan.drhohat_deps_v, scale_S / x_e);
                    std::vector<double> errors_e;
                    for (double const factor : {0.5, 1.0, 2.0})
                    {
                        double const hh = fdStep(c.eps(), x_e, factor);
                        double const fd =
                            (solve_at(0.0, +hh).exchange_received -
                             solve_at(0.0, -hh).exchange_received) /
                            (2.0 * hh);
                        errors_e.push_back(std::abs(fd - tan.drhohat_deps_v));
                    }
                    for (double const e : errors_e)
                    {
                        EXPECT_LE(e, tol_e) << "deps_v";
                        ut5_worst_ratio = std::max(ut5_worst_ratio, e / tol_e);
                    }

                    // dn/dp_L = 0 and dn/deps_v = phi' by differencing n_l.
                    {
                        double const hp = fdStep(c.p_L, x_p);
                        double const dn_dp = (solve_at(+hp, 0.0).local.n_l -
                                              solve_at(-hp, 0.0).local.n_l) /
                                             (2.0 * hp);
                        EXPECT_EQ(dn_dp, 0.0);
                        double const he = fdStep(c.eps(), x_e);
                        double const dn_de = (solve_at(0.0, +he).local.n_l -
                                              solve_at(0.0, -he).local.n_l) /
                                             (2.0 * he);
                        EXPECT_NEAR(dn_de, dphi,
                                    fdTolerance(dphi, c.phi / x_e));
                    }

                    // s = -1, alpha = 1, constant rho_lR: drhohat/deps_v =
                    // rho_lR/dt to round-off (DERIVED, D-2.4 identity).
                    if (s < 0.0 && regime == Regime::DeckEos)
                    {
                        EXPECT_NEAR(tan.drhohat_deps_v,
                                    base.local.rho_lR / c.dt,
                                    roundoffTolerance(base.local.rho_lR / c.dt));
                    }
                }
            }
        }
    }
    GTEST_LOG_(INFO) << "MEASURED UT-5: " << ut5_cases
                     << " cases; largest (FD error)/(tolerance C eps^(2/3) "
                        "(|analytic| + scale)) over d/dp_L and d/deps_v and "
                        "the three step factors: "
                     << ut5_worst_ratio;
}
