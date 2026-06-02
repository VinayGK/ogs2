// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BaseLib/ConfigTree.h"
#include "SaturationTuller.h"

#include <cmath>
#include <limits>

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationTuller> createSaturationTuller(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    auto const property_type = config.getConfigParameter<std::string>("type");
    if (property_type != "SaturationTuller" && property_type != "TullerRetention")
    {
        OGS_FATAL(
            "createSaturationTuller expected property type "
            "'SaturationTuller' or alias 'TullerRetention' but got '{:s}'.",
            property_type);
    }

    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create SaturationTuller medium property {:s}.", property_name);

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__SaturationTuller__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");

    auto const maximum_liquid_saturation_opt =
        //! \ogs_file_param{properties__property__SaturationTuller__maximum_liquid_saturation}
        config.getConfigParameterOptional<double>("maximum_liquid_saturation");
    auto const residual_gas_saturation_opt =
        //! \ogs_file_param{properties__property__SaturationTuller__residual_gas_saturation}
        config.getConfigParameterOptional<double>("residual_gas_saturation");

    if (!maximum_liquid_saturation_opt && !residual_gas_saturation_opt)
    {
        OGS_FATAL(
            "SaturationTuller requires either 'maximum_liquid_saturation' or "
            "'residual_gas_saturation' in the property configuration.");
    }

    double maximum_liquid_saturation = maximum_liquid_saturation_opt
                                           ? *maximum_liquid_saturation_opt
                                           : 1.0 - *residual_gas_saturation_opt;

    if (maximum_liquid_saturation_opt && residual_gas_saturation_opt)
    {
        double const maximum_from_residual_gas = 1.0 - *residual_gas_saturation_opt;
        if (std::abs(maximum_liquid_saturation - maximum_from_residual_gas) > 1e-12)
        {
            OGS_FATAL(
                "Inconsistent SaturationTuller bounds: "
                "maximum_liquid_saturation = {:g} but 1 - "
                "residual_gas_saturation = {:g}.",
                maximum_liquid_saturation, maximum_from_residual_gas);
        }
    }

    auto const area_factor_tuller =
        //! \ogs_file_param{properties__property__SaturationTuller__area_factor_tuller}
        config.getConfigParameter<double>("area_factor_tuller");
    auto const pore_area_shapefactor_tuller =
        //! \ogs_file_param{properties__property__SaturationTuller__pore_area_shapefactor_tuller}
        config.getConfigParameter<double>("pore_area_shapefactor_tuller");
    auto const characteristic_pore_size =
        //! \ogs_file_param{properties__property__SaturationTuller__characteristic_pore_size}
        config.getConfigParameter<double>("characteristic_pore_size");
    auto const surface_tension =
        //! \ogs_file_param{properties__property__SaturationTuller__surface_tension}
        config.getConfigParameter<double>("surface_tension");

    auto const pressure_tolerance =
        //! \ogs_file_param{properties__property__SaturationTuller__pressure_tolerance}
        config.getConfigParameterOptional<double>("pressure_tolerance")
            .value_or(0.0);

    // Optional Frydman--Baker cavitation cutoff. When omitted, the cutoff is
    // disabled (infinity) and the original asymptotic curve is recovered.
    auto const cavitation_pressure =
        //! \ogs_file_param{properties__property__SaturationTuller__cavitation_pressure}
        config.getConfigParameterOptional<double>("cavitation_pressure")
            .value_or(std::numeric_limits<double>::infinity());

    return std::make_unique<SaturationTuller>(
        std::move(property_name), residual_liquid_saturation,
        maximum_liquid_saturation, area_factor_tuller,
        pore_area_shapefactor_tuller, characteristic_pore_size, surface_tension,
        pressure_tolerance, cavitation_pressure);
}
}  // namespace MaterialPropertyLib
