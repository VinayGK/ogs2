// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SaturationTuller.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

namespace MaterialPropertyLib
{
SaturationTuller::SaturationTuller(
    std::string name, double const residual_liquid_saturation,
    double const maximum_liquid_saturation, double const area_factor_tuller,
    double const pore_area_shapefactor_tuller,
    double const characteristic_pore_size, double const surface_tension,
    double const pressure_tolerance, double const cavitation_pressure)
    : S_L_res_(residual_liquid_saturation),
      S_L_max_(maximum_liquid_saturation),
      coefficient_(4.0 * pore_area_shapefactor_tuller * surface_tension *
                   surface_tension /
                   (area_factor_tuller * characteristic_pore_size *
                    characteristic_pore_size)),
      pressure_tolerance_(pressure_tolerance),
      cavitation_pressure_(cavitation_pressure),
      saturation_at_cavitation_(saturationUncut(cavitation_pressure))
{
    name_ = std::move(name);

    if (!(0.0 <= S_L_res_ && S_L_res_ <= S_L_max_ && S_L_max_ <= 1.0))
    {
        OGS_FATAL(
            "SaturationTuller bounds must satisfy 0 <= S_L_res <= S_L_max <= "
            "1, but got S_L_res = {:g}, S_L_max = {:g}.",
            S_L_res_, S_L_max_);
    }

    if (!(area_factor_tuller > 0.0))
    {
        OGS_FATAL("SaturationTuller requires area_factor_tuller > 0, got {:g}.",
                  area_factor_tuller);
    }
    if (!(pore_area_shapefactor_tuller > 0.0))
    {
        OGS_FATAL(
            "SaturationTuller requires pore_area_shapefactor_tuller > 0, got "
            "{:g}.",
            pore_area_shapefactor_tuller);
    }
    if (!(characteristic_pore_size > 0.0))
    {
        OGS_FATAL(
            "SaturationTuller requires characteristic_pore_size > 0, got "
            "{:g}.",
            characteristic_pore_size);
    }
    if (!(surface_tension > 0.0))
    {
        OGS_FATAL("SaturationTuller requires surface_tension > 0, got {:g}.",
                  surface_tension);
    }
    if (!(pressure_tolerance_ >= 0.0))
    {
        OGS_FATAL(
            "SaturationTuller requires pressure_tolerance >= 0, got {:g}.",
            pressure_tolerance_);
    }
    if (!(coefficient_ > 0.0))
    {
        OGS_FATAL("SaturationTuller internal coefficient must be > 0, got {:g}.",
                  coefficient_);
    }
    if (!(cavitation_pressure_ > pressure_tolerance_))
    {
        OGS_FATAL(
            "SaturationTuller requires cavitation_pressure > pressure_tolerance "
            "({:g}), got {:g}. The cavitation cutoff must lie strictly inside "
            "the drained branch of the curve.",
            pressure_tolerance_, cavitation_pressure_);
    }
}

double SaturationTuller::saturationUncut(double const p_cap) const
{
    if (p_cap <= pressure_tolerance_)
    {
        return S_L_max_;
    }

    double const e = std::exp(-coefficient_ / (p_cap * p_cap));
    double const S_eff = 1.0 - e;
    double const S = S_L_res_ + (S_L_max_ - S_L_res_) * S_eff;

    return std::clamp(S, S_L_res_, S_L_max_);
}

PropertyDataType SaturationTuller::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p_cap = variable_array.capillary_pressure;

    // Beyond cavitation the macro meniscus has drained completely into the
    // micro structure: freeze the macro saturation (zero storage response).
    if (p_cap >= cavitation_pressure_)
    {
        return saturation_at_cavitation_;
    }

    return saturationUncut(p_cap);
}

PropertyDataType SaturationTuller::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (variable != Variable::capillary_pressure)
    {
        OGS_FATAL(
            "SaturationTuller::dValue is implemented for derivatives with "
            "respect to capillary pressure only.");
    }

    double const p_cap = variable_array.capillary_pressure;
    // Plateau (and the dry corner): zero storage response.
    if (p_cap <= pressure_tolerance_ || p_cap >= cavitation_pressure_)
    {
        return 0.0;
    }

    double const e = std::exp(-coefficient_ / (p_cap * p_cap));
    double const dS_eff_dp_cap = -(2.0 * coefficient_ / (p_cap * p_cap * p_cap)) * e;
    return dS_eff_dp_cap * (S_L_max_ - S_L_res_);
}

PropertyDataType SaturationTuller::d2Value(
    VariableArray const& variable_array, Variable const variable1,
    Variable const variable2, ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/, double const /*dt*/) const
{
    (void)variable1;
    (void)variable2;
    assert((variable1 == Variable::capillary_pressure) &&
           (variable2 == Variable::capillary_pressure) &&
           "SaturationTuller::d2Value is implemented for derivatives with "
           "respect to capillary pressure only.");

    double const p_cap = variable_array.capillary_pressure;
    if (p_cap <= pressure_tolerance_ || p_cap >= cavitation_pressure_)
    {
        return 0.0;
    }

    double const p2 = p_cap * p_cap;
    double const p4 = p2 * p2;
    double const p6 = p4 * p2;
    double const e = std::exp(-coefficient_ / p2);
    double const d2S_eff_dp_cap2 =
        2.0 * coefficient_ * e * (3.0 * p2 - 2.0 * coefficient_) / p6;
    return d2S_eff_dp_cap2 * (S_L_max_ - S_L_res_);
}
}  // namespace MaterialPropertyLib

