// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateRichardsMechanicsProcess.h"

#include <algorithm>
#include <cassert>
#include <memory>
#include <vector>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "NumLib/CreateNewtonRaphsonSolverParameters.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Common/HydroMechanics/CreateInitialStress.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "RichardsMechanicsProcess.h"
#include "RichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace RichardsMechanics
{
namespace
{
MicroPotentialConvention parseMicroPotentialConvention(
    std::string const& convention)
{
    if (convention == "positive_reduced")
    {
        return MicroPotentialConvention::PositiveReduced;
    }
    if (convention == "negative_attractive")
    {
        return MicroPotentialConvention::NegativeAttractive;
    }

    OGS_FATAL(
        "RichardsMechanics: unsupported potential_exchange "
        "micro_potential_convention '{}'. Currently supported: "
        "'positive_reduced', 'negative_attractive'.",
        convention);
}

LocalNonlinearSolveMode parseLocalNonlinearSolveMode(
    std::string const& mode)
{
    if (mode == "scalar_exchange")
    {
        return LocalNonlinearSolveMode::ScalarExchange;
    }
    if (mode == "scalar_microstate_storage_mode")
    {
        return LocalNonlinearSolveMode::ScalarReferenceStorage;
    }
    if (mode == "scalar_micro_macro_mass_storage_mode")
    {
        return LocalNonlinearSolveMode::ScalarReferenceMassStorage;
    }

    OGS_FATAL(
        "RichardsMechanics: unsupported potential_exchange "
        "local_nonlinear_solve_mode '{}'. Currently supported: "
        "'scalar_exchange', 'scalar_microstate_storage_mode', "
        "'scalar_micro_macro_mass_storage_mode'.",
        mode);
}

FilmStrainCouplingMode parseFilmStrainCouplingMode(std::string const& mode)
{
    if (mode == "off")
    {
        return FilmStrainCouplingMode::Off;
    }
    if (mode == "kinematic")
    {
        return FilmStrainCouplingMode::Kinematic;
    }
    if (mode == "equilibrium")
    {
        return FilmStrainCouplingMode::Equilibrium;
    }

    OGS_FATAL(
        "RichardsMechanics: unsupported potential_exchange "
        "film_strain_coupling '{}'. Currently supported: 'off', 'kinematic', "
        "'equilibrium'. (DSM/STRAINED_FILM_IMPLEMENTATION.md)",
        mode);
}

FilmStrainKappaMode parseFilmStrainKappaMode(std::string const& mode)
{
    if (mode == "aggregate")
    {
        return FilmStrainKappaMode::Aggregate;
    }
    if (mode == "unity")
    {
        return FilmStrainKappaMode::Unity;
    }

    OGS_FATAL(
        "RichardsMechanics: unsupported potential_exchange "
        "film_strain_kappa '{}'. Currently supported: 'aggregate' "
        "(kappa = 1 - phi_M, the integrable completion of the eigenstress "
        "scale), 'unity' (kappa = 1, naive geometric reading).",
        mode);
}

MacroPorosityUpdateMode parseMacroPorosityUpdateMode(
    std::string const& mode)
{
    if (mode == "algebraic_split")
    {
        return MacroPorosityUpdateMode::AlgebraicSplit;
    }
    if (mode == "additive_macro_porosity_rate_mode")
    {
        return MacroPorosityUpdateMode::ReferenceAdditiveRate;
    }

    OGS_FATAL(
        "RichardsMechanics: unsupported potential_exchange "
        "macro_porosity_update_mode '{}'. Currently supported: "
        "'algebraic_split', 'additive_macro_porosity_rate_mode'.",
        mode);
}

MicroSolidVolumeFractionMode parseMicroSolidVolumeFractionMode(
    std::string const& mode)
{
    if (mode == "reference")
    {
        return MicroSolidVolumeFractionMode::Reference;
    }
    if (mode == "current_porosity_split")
    {
        return MicroSolidVolumeFractionMode::CurrentPorositySplit;
    }

    OGS_FATAL(
        "RichardsMechanics: unsupported potential_exchange "
        "micro_solid_volume_fraction_mode '{}'. Currently supported: "
        "'reference', 'current_porosity_split'.",
        mode);
}

}  // namespace

void checkMPLProperties(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    std::array const required_medium_properties = {
        MaterialPropertyLib::reference_temperature,
        MaterialPropertyLib::bishops_effective_stress,
        MaterialPropertyLib::relative_permeability,
        MaterialPropertyLib::saturation,
        MaterialPropertyLib::porosity,
        MaterialPropertyLib::biot_coefficient};
    std::array const required_liquid_properties = {
        MaterialPropertyLib::viscosity, MaterialPropertyLib::density};
    std::array const required_solid_properties = {MaterialPropertyLib::density};

    for (auto const& m : media)
    {
        checkRequiredProperties(*m.second, required_medium_properties);
        checkRequiredProperties(m.second->phase(MaterialPropertyLib::PhaseName::AqueousLiquid),
                                required_liquid_properties);
        checkRequiredProperties(m.second->phase(MaterialPropertyLib::PhaseName::Solid),
                                required_solid_properties);
    }
}

void validateMicroPorosityAndPotentialExchangeConfiguration(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    std::optional<MicroPorosityParameters> const& micro_porosity_parameters,
    std::optional<PotentialExchangeParameters> const&
        potential_exchange_parameters,
    std::map<int, PotentialExchangeParameters> const&
        potential_exchange_parameters_by_material)
{
    namespace MPL = MaterialPropertyLib;

    bool const micro_porosity_enabled = micro_porosity_parameters.has_value();
    bool any_saturation_micro = false;
    bool const any_dsm_exchange_enabled =
        (potential_exchange_parameters &&
         potential_exchange_parameters->enabled) ||
        std::any_of(potential_exchange_parameters_by_material.begin(),
                    potential_exchange_parameters_by_material.end(),
                    [](auto const& item) { return item.second.enabled; });

    for (auto const& [material_id, medium] : media)
    {
        bool const has_saturation_micro =
            medium->hasProperty(MPL::PropertyType::saturation_micro);
        any_saturation_micro = any_saturation_micro || has_saturation_micro;

        if (has_saturation_micro && !micro_porosity_enabled)
        {
            OGS_FATAL(
                "RichardsMechanics: medium {} defines 'saturation_micro' but "
                "the process has no <micro_porosity> block. Define "
                "<micro_porosity> or remove 'saturation_micro'.",
                material_id);
        }
    }

    if (micro_porosity_enabled && !any_saturation_micro && !any_dsm_exchange_enabled)
    {
        OGS_FATAL(
            "RichardsMechanics: <micro_porosity> is configured, but no medium "
            "defines 'saturation_micro'. Define 'saturation_micro' in at least "
            "one medium or remove <micro_porosity>.");
    }

    if (any_dsm_exchange_enabled && !micro_porosity_enabled)
    {
        OGS_FATAL(
            "RichardsMechanics: potential_exchange.enabled=true requires "
            "a <micro_porosity> process block.");
    }

    for (auto const& [material_id, potential_exchange_params] :
         potential_exchange_parameters_by_material)
    {
        if (media.find(material_id) == media.end())
        {
            OGS_FATAL(
                "RichardsMechanics: potential_exchange medium override "
                "references unknown material id {}.",
                material_id);
        }

    }
}

PotentialExchangeParameters parsePotentialExchangeParameters(
    BaseLib::ConfigTree const& config,
    std::optional<PotentialExchangeParameters> const& defaults,
    std::string const& context)
{
    auto const enabled =
        config.getConfigParameter<bool>("enabled",
                                        defaults ? defaults->enabled : false);

    auto const pressure_tolerance = config.getConfigParameter<double>(
        "pressure_tolerance",
        defaults ? defaults->pressure_tolerance : 0.0);
    if (pressure_tolerance < 0.0)
    {
        OGS_FATAL(
            "RichardsMechanics: {} pressure_tolerance must be >= 0, got {:g}.",
            context, pressure_tolerance);
    }

    auto const micro_potential_convention = parseMicroPotentialConvention(
        config.getConfigParameter<std::string>(
            "micro_potential_convention",
            defaults ? toString(defaults->micro_potential_convention)
                     : "positive_reduced"));
    auto const local_nonlinear_solve_mode = parseLocalNonlinearSolveMode(
        config.getConfigParameter<std::string>(
            "local_nonlinear_solve_mode",
            defaults ? toString(defaults->local_nonlinear_solve_mode)
                     : "scalar_exchange"));
    auto const macro_porosity_update_mode = parseMacroPorosityUpdateMode(
        config.getConfigParameter<std::string>(
            "macro_porosity_update_mode",
            defaults ? toString(defaults->macro_porosity_update_mode)
                     : "algebraic_split"));
    auto const micro_solid_volume_fraction_mode =
        parseMicroSolidVolumeFractionMode(
            config.getConfigParameter<std::string>(
                "micro_solid_volume_fraction_mode",
                defaults
                    ? toString(defaults->micro_solid_volume_fraction_mode)
                    : "reference"));
    auto get_positive_required_or_default =
        [&](char const* const key, double const fallback)
    {
        auto const value = config.getConfigParameterOptional<double>(key);
        double const selected = value ? *value : fallback;
        if (!(selected > 0.0))
        {
            OGS_FATAL(
                "RichardsMechanics: {} {} must be > 0, got {:g}.", context,
                key, selected);
        }
        return selected;
    };

    auto get_positive_optional_or_default =
        [&](char const* const key, std::optional<double> const fallback)
            -> std::optional<double>
    {
        auto const value = config.getConfigParameterOptional<double>(key);
        std::optional<double> selected = value ? std::optional<double>{*value}
                                               : fallback;
        if (selected && !(*selected > 0.0))
        {
            OGS_FATAL(
                "RichardsMechanics: {} {} must be > 0 if provided, got {:g}.",
                context, key, *selected);
        }
        return selected;
    };

    double const default_hamaker =
        defaults ? defaults->hamaker_constant : 0.0;
    double const default_surface =
        defaults ? defaults->specific_surface : 0.0;
    double const default_rho_sr =
        defaults ? defaults->micro_solid_density_reference : 0.0;
    double const default_ns =
        defaults ? defaults->micro_solid_volume_fraction_reference : 0.0;
    double const default_rho_l0 =
        defaults ? defaults->micro_liquid_density_reference : 0.0;
    double const default_a_rho =
        defaults ? defaults->micro_liquid_density_a : 0.0;
    double const default_b_rho =
        defaults ? defaults->micro_liquid_density_b : 0.0;

    double hamaker_constant = 0.0;
    double specific_surface = 0.0;
    double micro_solid_density_reference = 0.0;
    double micro_solid_volume_fraction_reference = 0.0;
    double micro_liquid_density_reference = 0.0;
    double micro_liquid_density_a = 0.0;
    double micro_liquid_density_b = 0.0;
    bool const uses_micro_liquid_density_eos =
        local_nonlinear_solve_mode ==
        LocalNonlinearSolveMode::ScalarReferenceMassStorage;

    if (enabled)
    {
        hamaker_constant =
            get_positive_required_or_default("hamaker_constant",
                                             default_hamaker);
        specific_surface =
            get_positive_required_or_default("specific_surface",
                                             default_surface);
        micro_solid_density_reference = get_positive_required_or_default(
            "micro_solid_density_reference", default_rho_sr);
        micro_solid_volume_fraction_reference =
            get_positive_required_or_default(
                "micro_solid_volume_fraction_reference", default_ns);
        if (uses_micro_liquid_density_eos)
        {
            micro_liquid_density_reference =
                get_positive_required_or_default(
                    "micro_liquid_density_reference", default_rho_l0);
            micro_liquid_density_a = get_positive_required_or_default(
                "micro_liquid_density_a", default_a_rho);
            micro_liquid_density_b = get_positive_required_or_default(
                "micro_liquid_density_b", default_b_rho);
        }
        else
        {
            micro_liquid_density_reference =
                get_positive_optional_or_default(
                    "micro_liquid_density_reference",
                    defaults ? std::optional<double>{
                                   defaults->micro_liquid_density_reference}
                             : std::nullopt)
                    .value_or(0.0);
            micro_liquid_density_a = get_positive_optional_or_default(
                                         "micro_liquid_density_a",
                                         defaults
                                             ? std::optional<double>{
                                                   defaults->micro_liquid_density_a}
                                             : std::nullopt)
                                         .value_or(0.0);
            micro_liquid_density_b = get_positive_optional_or_default(
                                         "micro_liquid_density_b",
                                         defaults
                                             ? std::optional<double>{
                                                   defaults->micro_liquid_density_b}
                                             : std::nullopt)
                                         .value_or(0.0);
        }
    }
    else
    {
        hamaker_constant = get_positive_optional_or_default(
                               "hamaker_constant",
                               defaults ? std::optional<double>{
                                              defaults->hamaker_constant}
                                        : std::nullopt)
                               .value_or(0.0);
        specific_surface = get_positive_optional_or_default(
                               "specific_surface",
                               defaults ? std::optional<double>{
                                              defaults->specific_surface}
                                        : std::nullopt)
                               .value_or(0.0);
        micro_solid_density_reference =
            get_positive_optional_or_default(
                "micro_solid_density_reference",
                defaults ? std::optional<double>{
                               defaults->micro_solid_density_reference}
                         : std::nullopt)
                .value_or(0.0);
        micro_solid_volume_fraction_reference =
            get_positive_optional_or_default(
                "micro_solid_volume_fraction_reference",
                defaults ? std::optional<double>{
                               defaults->micro_solid_volume_fraction_reference}
                         : std::nullopt)
                .value_or(0.0);
        micro_liquid_density_reference =
            get_positive_optional_or_default(
                "micro_liquid_density_reference",
                defaults ? std::optional<double>{
                               defaults->micro_liquid_density_reference}
                         : std::nullopt)
                .value_or(0.0);
        micro_liquid_density_a =
            get_positive_optional_or_default(
                "micro_liquid_density_a",
                defaults
                    ? std::optional<double>{defaults->micro_liquid_density_a}
                    : std::nullopt)
                .value_or(0.0);
        micro_liquid_density_b =
            get_positive_optional_or_default(
                "micro_liquid_density_b",
                defaults
                    ? std::optional<double>{defaults->micro_liquid_density_b}
                    : std::nullopt)
                .value_or(0.0);
    }

    auto const initial_micro_water_content = get_positive_optional_or_default(
        "initial_micro_water_content",
        defaults ? defaults->initial_micro_water_content : std::nullopt);

    auto const use_fd_jacobian_for_exchange = config.getConfigParameter<bool>(
        "fd_jacobian_for_exchange",
        defaults ? defaults->use_fd_jacobian_for_exchange : false);

    auto const fd_jacobian_perturbation = config.getConfigParameter<double>(
        "fd_jacobian_perturbation",
        defaults ? defaults->fd_jacobian_perturbation : 1e-8);
    if (!(fd_jacobian_perturbation > 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} fd_jacobian_perturbation must be > 0, got {:g}.",
            context, fd_jacobian_perturbation);
    }

    auto const local_jacobian_perturbation = config.getConfigParameter<double>(
        "local_jacobian_perturbation",
        defaults ? defaults->local_jacobian_perturbation : 1e-8);
    if (!(local_jacobian_perturbation > 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} local_jacobian_perturbation must be > 0, got {:g}.",
            context, local_jacobian_perturbation);
    }

    // ── Augmentation prefactor K, optionally as a function of dry density ──
    // Two mutually exclusive ways to set K (J/kg):
    //   (1) scalar  <potential_augmentation_prefactor>
    //   (2) table   <potential_augmentation_prefactor_vs_dry_density> with
    //       child lists <dry_densities>/<prefactors>, evaluated at the
    //       material's <dry_density> (initial/target rho_d, kg/m^3).
    // Path (2) resolves to a scalar K = K(rho_d) at parse time; K is constant
    // in time (initial/target rho_d, Vinay 2026-06-08) so no Jacobian term is
    // introduced. The shared table inherits into per-<medium id> overrides via
    // `defaults`, while each medium supplies its own <dry_density>.
    std::shared_ptr<MathLib::PiecewiseLinearInterpolation const>
        potential_augmentation_prefactor_vs_dry_density =
            defaults ? defaults->potential_augmentation_prefactor_vs_dry_density
                     : nullptr;
    if (auto k_curve_config = config.getConfigSubtreeOptional(
            "potential_augmentation_prefactor_vs_dry_density"))
    {
        auto dry_densities =
            k_curve_config->getConfigParameter<std::vector<double>>(
                "dry_densities");
        auto prefactors =
            k_curve_config->getConfigParameter<std::vector<double>>(
                "prefactors");
        if (dry_densities.size() < 2 ||
            dry_densities.size() != prefactors.size())
        {
            OGS_FATAL(
                "RichardsMechanics: {} "
                "potential_augmentation_prefactor_vs_dry_density requires "
                "<dry_densities> and <prefactors> of equal length >= 2, got "
                "{} and {}.",
                context, dry_densities.size(), prefactors.size());
        }
        potential_augmentation_prefactor_vs_dry_density =
            std::make_shared<MathLib::PiecewiseLinearInterpolation const>(
                std::move(dry_densities), std::move(prefactors));
    }

    std::optional<double> dry_density =
        config.getConfigParameterOptional<double>("dry_density");
    if (!dry_density && defaults)
    {
        dry_density = defaults->dry_density;
    }

    // ── LIVE K(rho_d) (K_OF_RHO_D_LIVE.md; Vinay 2026-06-10) ───────────────
    // When true, the parse-time freeze below is SKIPPED for the live
    // evaluation path: the table stays live and K is re-evaluated at the
    // evolving rho_d = rho_SR*(1-phi) at run time (see
    // effectiveAugmentationPrefactor). The scalar stored into
    // potential_augmentation_prefactor then only serves as the FALLBACK for
    // evaluation sites without a porosity in scope.
    auto const potential_augmentation_prefactor_live_dry_density =
        config.getConfigParameter<bool>(
            "potential_augmentation_prefactor_live_dry_density",
            defaults
                ? defaults->potential_augmentation_prefactor_live_dry_density
                : false);
    if (potential_augmentation_prefactor_live_dry_density &&
        !potential_augmentation_prefactor_vs_dry_density)
    {
        OGS_FATAL(
            "RichardsMechanics: {} "
            "potential_augmentation_prefactor_live_dry_density=true requires "
            "a <potential_augmentation_prefactor_vs_dry_density> table.",
            context);
    }

    auto const potential_augmentation_prefactor_scalar =
        config.getConfigParameterOptional<double>(
            "potential_augmentation_prefactor");

    double potential_augmentation_prefactor;
    if (potential_augmentation_prefactor_vs_dry_density)
    {
        if (potential_augmentation_prefactor_scalar)
        {
            OGS_FATAL(
                "RichardsMechanics: {} both a scalar "
                "potential_augmentation_prefactor and a "
                "potential_augmentation_prefactor_vs_dry_density table were "
                "given; they are mutually exclusive.",
                context);
        }
        if (!dry_density && !potential_augmentation_prefactor_live_dry_density)
        {
            OGS_FATAL(
                "RichardsMechanics: {} "
                "potential_augmentation_prefactor_vs_dry_density requires a "
                "<dry_density> (rho_d, kg/m^3) to evaluate K(rho_d).",
                context);
        }
        // Live mode: this is NOT a freeze — the table stays live; the value
        // stored here is only the fallback K for phi-less evaluation sites
        // (initial/target rho_d if given, else inherited scalar / 0).
        potential_augmentation_prefactor =
            dry_density
                ? potential_augmentation_prefactor_vs_dry_density->getValue(
                      *dry_density)
                : (defaults ? defaults->potential_augmentation_prefactor
                            : 0.0);
    }
    else
    {
        potential_augmentation_prefactor =
            potential_augmentation_prefactor_scalar
                ? *potential_augmentation_prefactor_scalar
                : (defaults ? defaults->potential_augmentation_prefactor : 0.0);
    }
    if (!(potential_augmentation_prefactor >= 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} potential_augmentation_prefactor must be >= 0, got {:g}.",
            context, potential_augmentation_prefactor);
    }

    auto const potential_augmentation_exponent = config.getConfigParameter<double>(
        "potential_augmentation_exponent",
        defaults ? defaults->potential_augmentation_exponent : 0.0);
    if (potential_augmentation_prefactor > 0.0 &&
        !(potential_augmentation_exponent > 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} potential_augmentation_exponent must be > 0 when "
            "potential_augmentation_prefactor > 0, got {:g}.",
            context, potential_augmentation_exponent);
    }

    auto const use_micro_liquid_density_for_micro_pressure =
        config.getConfigParameter<bool>(
            "use_micro_liquid_density_for_micro_pressure",
            defaults ? defaults->use_micro_liquid_density_for_micro_pressure
                     : true);
    if (!use_micro_liquid_density_for_micro_pressure)
    {
        WARN(
            "RichardsMechanics: {} use_micro_liquid_density_for_micro_pressure=false selected; micro pressure will use bulk rho_LR instead of confined rho_lR.",
            context);
    }

    // ── Film-pressure coupling (maxwell beamer sec.5); default ON ──────────
    // RETIRED OFF path 2026-06-08 (Vinay): the model is consolidated on the film
    // coupling (biot=alpha). The bare-Pi OFF formulation is no longer selectable;
    // if a PRJ requests false it is overridden to true with a warning.
    bool film_pressure_coupling = config.getConfigParameter<bool>(
        "film_pressure_coupling",
        defaults ? defaults->film_pressure_coupling : true);
    if (!film_pressure_coupling)
    {
        WARN(
            "RichardsMechanics: {} film_pressure_coupling=false requested, but "
            "the bare-Pi OFF DSM path is RETIRED (consolidated on the film "
            "coupling, 2026-06-08). Overriding to film_pressure_coupling=true.",
            context);
        film_pressure_coupling = true;
    }
    // NOTE: the eigenstrain Biot b is unified with the poroelastic
    // biot_coefficient MPL property (no separate film_pressure_biot_b param).
    auto const film_pressure_gate_width = config.getConfigParameter<double>(
        "film_pressure_gate_width",
        defaults ? defaults->film_pressure_gate_width : 0.0);
    if (!(film_pressure_gate_width >= 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} film_pressure_gate_width must be >= 0, got {:g}.",
            context, film_pressure_gate_width);
    }
    // DEPRECATED 2026-06-06: swelling stress is now (1-phi_M)*p_film; this modulus is unused.
    auto const film_pressure_swelling_modulus =
        config.getConfigParameter<double>(
            "film_pressure_swelling_modulus",
            defaults ? defaults->film_pressure_swelling_modulus : 0.0);
    if (!(film_pressure_swelling_modulus >= 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} film_pressure_swelling_modulus must be >= 0, got {:g}.",
            context, film_pressure_swelling_modulus);
    }

    // ── Strained-film disjoining law h(w_m, eps_v) ──────────────────────────
    // (DSM/STRAINED_FILM_IMPLEMENTATION.md; Vinay 2026-06-09.) PRJ-selectable
    // variants: 'off' (default, frozen geometry, bit-for-bit), 'kinematic'
    // (variant A, spacing follows the volumetric strain), 'equilibrium'
    // (variant B, spacing tracks the film force balance Pi = p_conf). When ON,
    // the strained law REPLACES the shipped integrable mechanical partner (its
    // frozen-h truncation) — never both (no double counting).
    auto const film_strain_coupling = parseFilmStrainCouplingMode(
        config.getConfigParameter<std::string>(
            "film_strain_coupling",
            defaults ? toString(defaults->film_strain_coupling) : "off"));
    auto const film_strain_kappa = parseFilmStrainKappaMode(
        config.getConfigParameter<std::string>(
            "film_strain_kappa",
            defaults ? toString(defaults->film_strain_kappa) : "aggregate"));

    // Macro-porosity floor phi_M,min: keeps the macro pore from collapsing into
    // the interlayer (n_l capped at (phi-floor)/(1-floor)); 0 -> no floor.
    auto const macro_porosity_floor = config.getConfigParameter<double>(
        "macro_porosity_floor",
        defaults ? defaults->macro_porosity_floor : 0.0);
    if (!(macro_porosity_floor >= 0.0 && macro_porosity_floor < 1.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} macro_porosity_floor must be in [0, 1), got {:g}.",
            context, macro_porosity_floor);
    }
    auto const macro_floor_cutoff_width = config.getConfigParameter<double>(
        "macro_floor_cutoff_width",
        defaults ? defaults->macro_floor_cutoff_width : 0.0);
    if (!(macro_floor_cutoff_width >= 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} macro_floor_cutoff_width must be >= 0, got {:g}.",
            context, macro_floor_cutoff_width);
    }

    // Disjoining-pressure floor via a micro-water-content lower bound: clamps the
    // water content used in the vdW disjoining law so Pi = Pi(max(n_l, floor)),
    // capping Pi instead of diverging as n_l -> 0. 0 (default) -> no floor -> the
    // disjoining evaluation is byte-identical to before.
    auto const micro_water_content_floor = config.getConfigParameter<double>(
        "micro_water_content_floor",
        defaults ? defaults->micro_water_content_floor : 0.0);
    if (!(micro_water_content_floor >= 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} micro_water_content_floor must be >= 0, got "
            "{:g}.",
            context, micro_water_content_floor);
    }
    // (The former "experimental -- verify before trusting" WARN was dropped:
    // the film coupling is now the consolidated standard path, not opt-in.)
    return PotentialExchangeParameters{
        enabled,
        pressure_tolerance,
        hamaker_constant,
        specific_surface,
        micro_solid_density_reference,
        micro_solid_volume_fraction_reference,
        micro_liquid_density_reference,
        micro_liquid_density_a,
        micro_liquid_density_b,
        micro_potential_convention,
        local_nonlinear_solve_mode,
        macro_porosity_update_mode,
        micro_solid_volume_fraction_mode,
        initial_micro_water_content,
        use_fd_jacobian_for_exchange,
        fd_jacobian_perturbation,
        local_jacobian_perturbation,
        potential_augmentation_prefactor,
        potential_augmentation_exponent,
        micro_water_content_floor,
        use_micro_liquid_density_for_micro_pressure,
        film_pressure_coupling,
        film_pressure_gate_width,
        film_pressure_swelling_modulus,
        macro_porosity_floor,
        macro_floor_cutoff_width,
        film_strain_coupling,
        film_strain_kappa,
        potential_augmentation_prefactor_vs_dry_density,
        dry_density,
        potential_augmentation_prefactor_live_dry_density};
}

template <int DisplacementDim>
std::unique_ptr<Process> createRichardsMechanicsProcess(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "RICHARDS_MECHANICS");
    DBUG("Create RichardsMechanicsProcess.");

    auto const coupling_scheme =
        //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(coupling_scheme && (*coupling_scheme == "staggered"));

    /// \section processvariablesrm Process Variables

    //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_p;
    ProcessVariable* variable_u;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        /// Primary process variables as they appear in the global component
        /// vector:
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__RICHARDS_MECHANICS__process_variables__pressure}
             "pressure",
             //! \ogs_file_param_special{prj__processes__process__RICHARDS_MECHANICS__process_variables__displacement}
             "displacement"});
        variable_p = &per_process_variables[0].get();
        variable_u = &per_process_variables[1].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (auto const& variable_name : {"pressure"s, "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_p = &process_variables[0][0].get();
        variable_u = &process_variables[1][0].get();
    }

    DBUG("Associate displacement with process variable '{:s}'.",
         variable_u->getName());

    if (variable_u->getNumberOfGlobalComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '{:s}' is different "
            "from the displacement dimension: got {:d}, expected {:d}",
            variable_u->getName(),
            variable_u->getNumberOfGlobalComponents(),
            DisplacementDim);
    }

    DBUG("Associate pressure with process variable '{:s}'.",
         variable_p->getName());
    if (variable_p->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '{:s}' is not a scalar variable but has "
            "{:d} components.",
            variable_p->getName(),
            variable_p->getNumberOfGlobalComponents());
    }

    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, materialIDs(mesh), config);

    /// \section parametersrm Process Parameters
    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != DisplacementDim)
        {
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is {:d}, displacement "
                "dimension is {:d}",
                b.size(), DisplacementDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);
    DBUG("Check the media properties of RichardsMechanics process ...");
    checkMPLProperties(media);
    DBUG("Media properties verified.");

    // Initial stress conditions
    auto const initial_stress =
        ProcessLib::createInitialStress<DisplacementDim>(config, parameters,
                                                         mesh);

    std::optional<MicroPorosityParameters> micro_porosity_parameters;
    if (auto const micro_porosity_config =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__micro_porosity}
        config.getConfigSubtreeOptional("micro_porosity"))
    {
        micro_porosity_parameters = MicroPorosityParameters{
            NumLib::createNewtonRaphsonSolverParameters(
                //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__micro_porosity__nonlinear_solver}
                micro_porosity_config->getConfigSubtree("nonlinear_solver")),
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__micro_porosity__mass_exchange_coefficient}
            micro_porosity_config->getConfigParameter<double>(
                "mass_exchange_coefficient")};
    }

    std::optional<PotentialExchangeParameters> potential_exchange_parameters;
    std::map<int, PotentialExchangeParameters>
        potential_exchange_parameters_by_material;
    if (auto const potential_exchange_config =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__potential_exchange}
            config.getConfigSubtreeOptional("potential_exchange"))
    {
        potential_exchange_parameters = parsePotentialExchangeParameters(
            *potential_exchange_config, std::nullopt,
            "potential_exchange");

        for (auto medium_config :
             potential_exchange_config->getConfigSubtreeList("medium"))
        {
            int const material_id = medium_config.getConfigAttribute<int>("id");
            if (!potential_exchange_parameters_by_material
                     .emplace(material_id,
                              parsePotentialExchangeParameters(
                                  medium_config,
                                  potential_exchange_parameters,
                                  fmt::format(
                                      "potential_exchange medium id {}",
                                      material_id)))
                     .second)
            {
                OGS_FATAL(
                    "RichardsMechanics: duplicate potential_exchange medium override for material id {}.",
                    material_id);
            }
        }
    }

    validateMicroPorosityAndPotentialExchangeConfiguration(
        media, micro_porosity_parameters, potential_exchange_parameters,
        potential_exchange_parameters_by_material);

    auto const mass_lumping =
        //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__mass_lumping}
        config.getConfigParameter<bool>("mass_lumping", false);

    auto const explicit_hm_coupling_in_unsaturated_zone =
        //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__explicit_hm_coupling_in_unsaturated_zone}
        config.getConfigParameter<bool>(
            "explicit_hm_coupling_in_unsaturated_zone", false);

    auto const is_linear =
        //! \ogs_file_param{prj__processes__process__linear}
        config.getConfigParameter("linear", false);

    bool const use_numerical_jacobian =
        jacobian_assembler->isPerturbationEnabled();

    RichardsMechanicsProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(media_map),
        std::move(solid_constitutive_relations),
        initial_stress,
        specific_body_force,
        micro_porosity_parameters,
        potential_exchange_parameters,
        potential_exchange_parameters_by_material,
        mass_lumping,
        explicit_hm_coupling_in_unsaturated_zone,
        use_numerical_jacobian};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<RichardsMechanicsProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme, is_linear);
}

template std::unique_ptr<Process> createRichardsMechanicsProcess<2>(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

template std::unique_ptr<Process> createRichardsMechanicsProcess<3>(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

}  // namespace RichardsMechanics
}  // namespace ProcessLib
