// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

#include "BaseLib/StrongType.h"

namespace ProcessLib::RichardsMechanics
{
using MicroPorosity = BaseLib::StrongType<double, struct MicroPorosityTag>;

constexpr std::string_view ioName(struct MicroPorosityTag*)
{
    return "micro_porosity";
}
}  // namespace ProcessLib::RichardsMechanics
