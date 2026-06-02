# DSM Native Hierarchical Porosity Swap Plan

## Goal
Replace the additive/naive micro-macro porosity split in the RM DSM-native path with the hierarchical split used in the PINION derivation.

## Target model

- State variable kept as-is:
  - `n_l` = micro water content (`MicroWaterContent`), aggregate-referenced.
- Porosity split law (REV-referenced):
  - `phi = phi_M + (1 - phi_M) * n_l`
  - `phi_M = (phi - n_l) / (1 - n_l)`
  - `phi_m = (1 - phi_M) * n_l`

## Implementation scope

File:
- `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`

Changes:
1. Replace `computeTransportPorosityUpdate(...)` with hierarchical split evaluation.
2. Keep `macro_porosity_update_mode="additive_macro_porosity_rate_mode"` as a config alias, but evaluate the hierarchical law.
3. Use previous micro-porosity (`phi_m_prev`) consistently in local solve context (instead of previous `n_l`).
4. Initialize `n_l` from `(phi, phi_M)` via hierarchical inversion, not via `phi - phi_M`.
5. Initialize `phi_m` from hierarchical relation `phi_m = (1 - phi_M) * n_l`.
6. Update split state consistently for all DSM-local modes through `computeTransportPorosityUpdate(...)`.

## Expected behavior changes

- `micro_porosity` now represents REV micro porosity from the hierarchical split, not directly `n_l`.
- Macro/micro split remains bounded and consistent with total porosity.
- Legacy project files using `additive_macro_porosity_rate_mode` run on this branch with hierarchical split behavior.

## Swelling slope correction (2026-05-19) — SUPERSEDED (2026-06-01)

> **HISTORICAL / SUPERSEDED.** This section documents the removed parameter
> `micro_water_content_swelling_slope` (and the companion
> `accumulate_swelling_contributions` switch). Both were deleted from
> `PotentialExchangeParameters` and `CreateRichardsMechanicsProcess` when
> `computeReferenceMicroPorositySwellingStressIncrement` was rewritten to derive
> the swelling-stress increment from the micro-liquid-density / micro-porosity
> state directly (see its current signature), rather than from an externally
> supplied slope. Retained for provenance only; the slope parameter no longer
> exists in the code or in project files.

**Issue identified**: `computeReferenceMicroPorositySwellingStressIncrement` was using
`delta_phi_m = phi_m - phi_m_prev` as the swelling driver. In the hierarchical split
`phi_m = (1 - phi_M) * n_l`, so for a fixed-boundary confined test (phi = const):
  `phi_m ≈ nS_ref * n_l`  (nS_ref ≈ 0.5 for typical density)
This made the effective swelling slope ≈ 0.5× the calibrated value, reducing the
achievable swelling pressure compared to the pre-hierarchical model for the same slope.

**Fix**: Changed the swelling driver to `delta_n_l = n_l - n_l_prev` (the LOCAL micro
water content, which is the `MicroWaterContent` state variable). Rationale:
- The parameter is named `micro_water_content_swelling_slope` — it should couple to n_l.
- Before the hierarchical split, `n_l ≡ phi_m`, so old calibrations are preserved exactly.
- After the split, the slope now couples to the LOCAL change (film-thickness-proportional),
  not the REV-scale porosity change.

**Pressure cap note**: Even with this fix the LinearElastic model (E=52 MPa, ν=0.3)
caps swelling at `P_sw_max = K_bulk × slope × phi0 ≈ 2 MPa`. Matching high-density
Villar targets (up to 22.6 MPa) requires a pressure-dependent constitutive model
(e.g. ModCamClay via MFront).

## Validation checklist

1. Compile RM process.
2. Run DSM micromacro Richards tests, especially:
   - `beacon_1a01_dsm_micromacro_inflowprobe`
   - `beacon_1a01_dsm_micromacro_stressprobe`
3. Verify:
   - no negative porosity states,
   - `phi_m <= phi`,
   - stable local implicit `n_l` solve,
   - expected trend shifts vs additive baseline.
