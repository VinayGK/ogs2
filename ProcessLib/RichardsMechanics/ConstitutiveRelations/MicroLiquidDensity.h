// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

#include "BaseLib/StrongType.h"

namespace ProcessLib::RichardsMechanics
{
using MicroLiquidDensity =
    BaseLib::StrongType<double, struct MicroLiquidDensityTag>;

constexpr std::string_view ioName(struct MicroLiquidDensityTag*)
{
    return "micro_liquid_density";
}
}  // namespace ProcessLib::RichardsMechanics
