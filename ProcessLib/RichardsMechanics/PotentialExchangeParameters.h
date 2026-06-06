// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>

namespace ProcessLib::RichardsMechanics
{
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

    // Optional consistency switches for the hierarchical DSM branch.
    // Default micro-pressure density is the confined micro-liquid density.
    bool use_micro_liquid_density_for_micro_pressure = true;

    // ── Film-pressure coupling (maxwell beamer sec.5) ──────────────────────
    // Master gate for the equipresent stress-coupled micro potential
    // mu_lR(p_film = p_disj + sigma'). Default OFF preserves the existing path
    // bit-for-bit. When true: mu_lR carries the effective-stress (film) term in
    // ALL local solves and the macro exchange, the swelling stress becomes the
    // eigenstrain form (S1 < 0 -> compression drains), and the sharp gate is
    // replaced by a C1 activation of width film_pressure_gate_width.
    bool film_pressure_coupling = false;
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
};
}  // namespace ProcessLib::RichardsMechanics
