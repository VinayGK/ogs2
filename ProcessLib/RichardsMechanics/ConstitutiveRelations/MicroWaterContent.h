// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

#include "BaseLib/StrongType.h"

namespace ProcessLib::RichardsMechanics
{
using MicroWaterContent =
    BaseLib::StrongType<double, struct MicroWaterContentTag>;

constexpr std::string_view ioName(struct MicroWaterContentTag*)
{
    return "micro_water_content";
}

using MicroExchangeSource =
    BaseLib::StrongType<double, struct MicroExchangeSourceTag>;

constexpr std::string_view ioName(struct MicroExchangeSourceTag*)
{
    return "micro_exchange_source";
}
}  // namespace ProcessLib::RichardsMechanics
