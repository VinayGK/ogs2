// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <string_view>

#include "BaseLib/StrongType.h"

// Per-integration-point state and output fields of the KKT micro-water ceiling
// (branch dsm_mass_conservation_v3_kkt_ceiling_2026-09-30; DESIGN.md 3.3 of the
// record folder ~/ogs-models/scratch/2026-09-30_kkt_ceiling_impl/). NOT adopted.
//
// All twelve are written ONLY when potential_exchange micro_ceiling_treatment =
// kkt (otherwise they stay 0 and nothing reads them). The local KKT solve never
// reads any of them (memoryless); they are outputs and diagnostics, plus the
// memory of the group-3 iteration diagnostics. They become integration-point
// writers and secondary variables through the StatefulData reflection.
namespace ProcessLib::RichardsMechanics
{
// 0 Interior, 1 Active, 2 NonMonotone, 3 PremiseViolated (DESIGN.md 3.2). The
// assembler reads status == 1 as the active flag (no detector, no tolerance).
using MicroCeilingStatus =
    BaseLib::StrongType<double, struct MicroCeilingStatusTag>;
constexpr std::string_view ioName(struct MicroCeilingStatusTag*)
{
    return "micro_ceiling_status";
}

// lambda >= 0 of this evaluation [Pa]; 0 unless Active. (The code identifier is
// never "lambda": that symbol is also the vdW film spacing.)
using MicroCeilingMultiplier =
    BaseLib::StrongType<double, struct MicroCeilingMultiplierTag>;
constexpr std::string_view ioName(struct MicroCeilingMultiplierTag*)
{
    return "micro_ceiling_multiplier";
}

// rhohat the micro gained and the macro sink uses [kg/(m3 s)]: S_s when Active,
// rhohat_pot otherwise. The field the books integrate in KKT runs
// (DERIVATION.md 0.2).
using MicroExchangeReceived =
    BaseLib::StrongType<double, struct MicroExchangeReceivedTag>;
constexpr std::string_view ioName(struct MicroExchangeReceivedTag*)
{
    return "micro_exchange_received";
}

// r = rhohat_pot - S_s >= 0 on the active branch, else 0 [kg/(m3 s)]. Identity:
// micro_exchange_source - micro_exchange_received = r.
using MicroCeilingRejectedExchange =
    BaseLib::StrongType<double, struct MicroCeilingRejectedExchangeTag>;
constexpr std::string_view ioName(struct MicroCeilingRejectedExchangeTag*)
{
    return "micro_ceiling_rejected_exchange";
}

// Cumulative number of changes of (status == 1) between two successive
// evaluations at this IP [count]; includes step boundaries and rejected
// attempts.
using MicroCeilingFlips =
    BaseLib::StrongType<double, struct MicroCeilingFlipsTag>;
constexpr std::string_view ioName(struct MicroCeilingFlipsTag*)
{
    return "micro_ceiling_flips";
}

// Cumulative number of status-2 evaluations [count].
using MicroCeilingNonMonotone =
    BaseLib::StrongType<double, struct MicroCeilingNonMonotoneTag>;
constexpr std::string_view ioName(struct MicroCeilingNonMonotoneTag*)
{
    return "micro_ceiling_nonmonotone";
}

// Cumulative number of status-3 evaluations [count].
using MicroCeilingPremise =
    BaseLib::StrongType<double, struct MicroCeilingPremiseTag>;
constexpr std::string_view ioName(struct MicroCeilingPremiseTag*)
{
    return "micro_ceiling_premise";
}

// Group 3, iteration diagnostics (memory): total eps_v at the previous distinct
// evaluation [-].
using MicroCeilingEpsSeen =
    BaseLib::StrongType<double, struct MicroCeilingEpsSeenTag>;
constexpr std::string_view ioName(struct MicroCeilingEpsSeenTag*)
{
    return "micro_ceiling_eps_seen";
}

// Last nonzero increment eps_v^k - eps_v^(k-1) of the current attempt (0 at
// attempt start) [-].
using MicroCeilingIncLast =
    BaseLib::StrongType<double, struct MicroCeilingIncLastTag>;
constexpr std::string_view ioName(struct MicroCeilingIncLastTag*)
{
    return "micro_ceiling_inc_last";
}

// The time t of the attempt the memory belongs to [s]. A retried attempt has a
// different t, because the retry has a smaller dt.
using MicroCeilingAttemptT =
    BaseLib::StrongType<double, struct MicroCeilingAttemptTTag>;
constexpr std::string_view ioName(struct MicroCeilingAttemptTTag*)
{
    return "micro_ceiling_attempt_t";
}

// Cumulative number of pairs of successive nonzero increments of opposite sign
// with status 1 at both evaluations [count] (the parity signature,
// DERIVATION.md 4.7).
using MicroCeilingIncAlt =
    BaseLib::StrongType<double, struct MicroCeilingIncAltTag>;
constexpr std::string_view ioName(struct MicroCeilingIncAltTag*)
{
    return "micro_ceiling_inc_alt";
}

// The same with equal sign [count].
// Tag of one trace evaluation (DESIGN.md 3.8): which element and integration
// point, at which time, and whether the evaluation is the output re-evaluation.
// Built by the assembler only for the elements of micro_ceiling_trace_elements.
struct MicroCeilingTraceTag
{
    std::size_t element_id = 0;
    std::size_t integration_point = 0;
    double t = 0.0;
    double dt = 0.0;
    bool output_reevaluation = false;
};

using MicroCeilingIncSame =
    BaseLib::StrongType<double, struct MicroCeilingIncSameTag>;
constexpr std::string_view ioName(struct MicroCeilingIncSameTag*)
{
    return "micro_ceiling_inc_same";
}

// Saturated latch L of the gate micro_ceiling_saturation_gate (branch
// dsm_mass_conservation_v3_kkt_vii_gate_2026-10-01; design part B.4 of
// ~/ogs-models/scratch/2026-10-01_kkt_iv_vii_fixes/DESIGN_FIXES.md): 1 when the
// point is KKT-active and was, at this or an earlier accepted step, active with
// the deck's Bishop factor equal to 1; else 0. Written only with
// micro_ceiling_saturation_gate != off (otherwise it stays 0 and nothing reads
// it). A stateful field: its previous-step copy is the L_old of the next step.
using MicroSaturatedLatch =
    BaseLib::StrongType<double, struct MicroSaturatedLatchTag>;
constexpr std::string_view ioName(struct MicroSaturatedLatchTag*)
{
    return "micro_saturated_latch";
}

// 1b closed-macro gate marker (branch dsm_mass_conservation_v4_tm_krel_
// 2026-10-02; DESIGN_V4.md 2.1, 2.4): 1 when the gate
// micro_ceiling_closed_macro_gate acted at this evaluation (previous step
// Active with phi_M == 0, and Active now), else 0. Written at every evaluation
// with potential exchange on, at every level (identically 0 with the switch
// off). Read only by the bishop_relperm level (chi_prev rule, DESIGN_V4.md
// 2.4.4) through its previous-step copy; otherwise an output map of where 1b
// acted (micro_closed_macro_gate_ip).
using MicroClosedMacroGateActed =
    BaseLib::StrongType<double, struct MicroClosedMacroGateActedTag>;
constexpr std::string_view ioName(struct MicroClosedMacroGateActedTag*)
{
    return "micro_closed_macro_gate";
}

// DIAGNOSTIC, NOT FOR PRODUCTION (swelling-stress fix (b), 2026-10-04): the
// level L = -n_S n_l [Pi + b sigma'_mean] [Pa] that the LAST evaluation of the
// swelling stress used (level form only; identically 0 with the default step
// form). Stateful: its previous-step copy is the L_prev of the next step, so
// that sigma_sw = L - L_ref holds exactly, with L_ref the level of the initial
// state, and the lag of sigma'_mean within a step does not accumulate over the
// steps. 0 = not yet set (the level of the previous state is then recomputed
// from the previous-step state variables).
using SwellingLevelUsed =
    BaseLib::StrongType<double, struct SwellingLevelUsedTag>;
constexpr std::string_view ioName(struct SwellingLevelUsedTag*)
{
    return "swelling_level_used";
}
}  // namespace ProcessLib::RichardsMechanics
