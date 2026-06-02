// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

/**
 * \brief Tuller/Young-Laplace-inspired saturation model used in the DSM
 *        micro-macro reference branch
 * macro retention law.
 *
 * The implemented saturation relation is (for capillary pressure p_c):
 * \f[
 * S_L(p_c)=
 * \begin{cases}
 * S_{L,\max}, & p_c \le p_\mathrm{tol},\\
 * S_{L,\mathrm{res}} +
 * (S_{L,\max}-S_{L,\mathrm{res}})
 * \left(1-\exp\left(-\frac{C_T}{p_c^2}\right)\right),
 *     & p_\mathrm{tol} < p_c < p_\mathrm{cav},\\
 * S_L(p_\mathrm{cav}), & p_c \ge p_\mathrm{cav},
 * \end{cases}
 * \f]
 * with
 * \f[
 * C_T = \frac{4 F_\gamma \sigma^2}{A_n L^2}.
 * \f]
 *
 * The parameters \f$A_n\f$, \f$F_\gamma\f$, \f$L\f$, \f$\sigma\f$ correspond
 * to the DSM parameter names `AreaFactorTuller`, `PoreAreaShapefactorTuller`,
 * `CharacteristicPoreSize`, and `SurfaceTension`.
 *
 * \f$p_\mathrm{cav}\f$ is an optional cavitation pressure (Frydman--Baker
 * dual-porosity reading of the DSM closure): once the capillary pressure
 * reaches it, the macro capillary meniscus has cavitated and the macro pore
 * water has drained completely into the micro structure. Beyond
 * \f$p_\mathrm{cav}\f$ the macro saturation is frozen at its value there (zero
 * storage response), so the unphysical post-cavitation branch of the capillary
 * closure is never evaluated. \f$p_\mathrm{cav}\f$ should be chosen in the
 * drained tail of the curve where the slope is already negligible, so the
 * \f$C^1\f$ kink introduced by the freeze is mild. If \f$p_\mathrm{cav}\f$ is
 * not given, the original asymptotic (uncut) curve is used.
 */
class SaturationTuller final : public Property
{
public:
    SaturationTuller(std::string name, double residual_liquid_saturation,
                     double maximum_liquid_saturation,
                     double area_factor_tuller,
                     double pore_area_shapefactor_tuller,
                     double characteristic_pore_size, double surface_tension,
                     double pressure_tolerance, double cavitation_pressure);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'SaturationTuller' is implemented on the "
                "'media' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& pos,
                             double const t, double const dt) const override;

    //! Saturation evaluated on the uncut capillary branch at \f$p_c\f$.
    double saturationUncut(double p_cap) const;

private:
    double const S_L_res_;
    double const S_L_max_;
    double const coefficient_;
    double const pressure_tolerance_;
    //! Cavitation pressure; std::numeric_limits<double>::infinity() disables
    //! the cutoff (original asymptotic curve).
    double const cavitation_pressure_;
    //! Frozen saturation value on the post-cavitation plateau, i.e.
    //! saturationUncut(cavitation_pressure_). When the cutoff is disabled
    //! (cavitation_pressure_ == +inf) this evaluates to S_L_res_, but it is
    //! never returned because the plateau branch p_c >= +inf is unreachable.
    double const saturation_at_cavitation_;
};
}  // namespace MaterialPropertyLib
