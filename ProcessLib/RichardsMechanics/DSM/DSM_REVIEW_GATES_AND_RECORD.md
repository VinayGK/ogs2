# DSM Review — Gates and Record (+ unit-test recommendations)

Scope: deep review of the DSM native (p_disj augmentation + Tuller) port on
branch `claude/piofnlev-review-tests-slj9pd`, focused on (1) the **gating /
clamping logic** ("gates") and (2) the **recorded-baseline test
infrastructure** ("record"), with concrete unit-test recommendations.

Files reviewed:
- `ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h`
- `ProcessLib/RichardsMechanics/PotentialExchangeParameters.h`
- `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h` (gate clamps + active-fraction helpers)
- `MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.{h,cpp}`
- `Tests/ProcessLib/RichardsMechanics/DSMMicroMacroSingleIntegrationPoint.cpp`

---

## 1. Gates

### 1.1 Porosity-split clamp ("Gate 1/2 fix") — duplicated and lossy

The clamp appears **twice, verbatim**, at `RichardsMechanicsFEM-impl.h:2967-2987`
(constitutive-state path) and `:4033-4056` (output-field path):

```cpp
phi_m = std::min(phi_m, phi_total);
phi_M = phi_total - phi_m;   // >= 0
```

Findings:
- **DRY**: the two copies must stay in lockstep; a future edit to one is a
  silent-divergence hazard. Extract a single helper
  (`clampHierarchicalSplit(phi_total, phi_m) -> {phi_M, phi_m}`) and call it
  from both sites.
- **Mass is silently dropped when the clamp bites.** Capping `phi_m` at
  `phi_total` discards the excess micro porosity (and the water it carries)
  with no rebalancing into `phi_M` storage and no diagnostic. The code comment
  correctly flags this as an output clamp over a real constitutive gap (missing
  micro-swelling saturation), but operationally there is no way to know *how
  often* or *how hard* the clamp fires in a run. Add at minimum a one-time
  warning / counter when `phi_m > phi_total` so a silently-clamped run is
  distinguishable from a healthy one.
- The clamp guarantees `phi_M >= 0` and `phi_m + phi_M == phi_total` by
  construction — good — but those invariants are currently only checked
  indirectly via the porosity-recomposition assertions in the unit tests, and
  never on the clamp branch itself (the clamp branch needs `phi_m > phi_total`
  to trigger, which no committed test forces).

### 1.2 Young–Laplace branch gate (`PotentialExchange.h:41`)

`saturated_branch = (p_LR > -pressure_tolerance)` returns `mu_LR = 0` on the
saturated side and `p_LR / rho_LR` otherwise. This is **continuous only for
`pressure_tolerance == 0`**: at the boundary `p_LR = -ptol` the unsaturated
value is `-ptol / rho_LR ≠ 0`, so for `ptol > 0` there is a finite jump of
`ptol / rho_LR` in `mu_LR`. With the values used in the suite (`ptol = 1e-12`)
the jump is ~1e-15 and harmless, but it is a latent Newton-convergence trap if
anyone raises `ptol`. Worth a continuity assertion across the branch.

### 1.3 vdW guards (`PotentialExchange.h:105-156`)

Comprehensive `OGS_FATAL` guards on every denominator/positive input
(`n_l, rho_lR, nS, rho_SR, A, Sa > 0`; `K >= 0`; `lambda > 0` when `K > 0`).
This is solid defensive code. None of these guard paths is currently exercised
by a test (see §3).

### 1.4 SaturationTuller branch gates (`SaturationTuller.cpp`)

- Constructor validates bounds and `cavitation_pressure > pressure_tolerance`
  with clear `OGS_FATAL`s — good.
- `value()`/`dValue()`/`d2Value()` analytic derivatives were checked by hand
  and are **correct** (`dS_eff/dp = -(2c/p³)e`,
  `d²S_eff/dp² = 2c·e·(3p² − 2c)/p⁶`).
- **Plateau-edge discontinuity for `ptol > 0`**: `saturationUncut` returns
  `S_L_max` for `p_cap <= ptol`, but for `p_cap` just above `ptol` it returns
  `S_L_res + (S_L_max − S_L_res)(1 − e^{−c/p²})`, which is strictly below
  `S_L_max`. The two only meet as `ptol → 0`. Continuous at the cavitation
  freeze in value (by construction, `saturation_at_cavitation_`) but the
  derivative is intentionally kinked there.

### 1.5 active micro-solid-fraction clamps (`impl.h:396-427`)

`computeActiveMicroSolidVolumeFraction` / `computePreviousMicroSolidVolumeFraction`
return `1 − n_l` (CurrentPorositySplit) with `clamp(n_l, 0, 1−1e-12)` and a
`1e-16` floor. Consistent with the documented "active_nS = 1 − n_l" fix. See
§2.3 for a stale *test-side* copy of the old definition.

---

## 2. Record (recorded-baseline test infrastructure)

### 2.1 The baseline CSVs are missing from the branch — two tests cannot pass

`DSMMicroMacroOverlapTransferBaselineHistory` and
`DSMMicroMacroStrainCoupledOverlapBaselineHistory` load:

```
<data_path>/RichardsMechanics/DSMMicroMacroOverlapTransferBaseline.csv
<data_path>/RichardsMechanics/DSMMicroMacroStrainCoupledOverlapBaseline.csv
```

Neither file exists anywhere in the branch (or its history), and nothing
generates them at build/test time. Both tests then hit
`ASSERT_EQ(baseline_rows.size(), 5)` on an empty vector and **fail on a fresh
checkout / in CI**. This is almost certainly a port-omission: the data lived on
the original author machine and was not carried onto the fresh-upstream port.

Action: commit the two CSVs to `Tests/Data/RichardsMechanics/`, or make the
loader `GTEST_SKIP()` with a clear message when the file is absent (skipping is
strictly worse than committing the data — prefer committing).

### 2.2 The baseline-history tests barely use the record they load

Even with the CSVs present, these tests use only `row.pressure` (as input) and
`row.saturation` (checked `== 1`). Every other recorded column — `n_l`,
`mu_lR`, `phi_m`, `phi_M`, `rho_lR`, `rho_l_hat`, `epsilon_sw`, `stress_xx` —
is loaded into the struct and then **never compared**. The rest of each
assertion block is `EXPECT_TRUE(std::isfinite(...))` plus a
`phi_m + phi_M == phi` recomposition check. As written these are smoke tests,
not regression guards: a sign flip or coefficient drift in the exchange would
sail through as long as the result stays finite.

Action: compare the computed `n_l`, `mu_lR`, `rho_l_hat`, `phi_m`, `phi_M`
(and `epsilon_sw`/`stress_xx` for the strain-coupled case) against the recorded
columns with `comparisonTolerance(...)`, the same way
`DSMMicroMacroSingleIntegrationPointReferencePath` already does against its
independent oracle.

### 2.3 Stale test-side oracle for `active_nS` (latent)

`referenceMicroSolidVolumeFraction` (test file, ~line 157) still returns
`1 − phi_M − phi_m` (= `1 − phi`) for `CurrentPorositySplit` mode, i.e. the
*old* definition. Production (`impl.h:408`) now returns `1 − n_l`. This helper
feeds `solveDsmMicromacroReferenceSinglePoint`, but every committed
reference-oracle test runs in the default `Reference` mode, so the stale branch
is currently dead. It is a trap: the first reference-oracle test written in
CurrentPorositySplit mode will see the oracle disagree with production. Fix the
helper to `1 − n_l` to match the production definition.

### 2.4 What is good in the record machinery

`DSMMicroMacroSingleIntegrationPointReferencePath` compares the production
helper `solveImplicitMicroWaterContent` against a genuinely **independent**
reimplementation (its own Newton loop + FD Jacobian in
`solveDsmMicromacroReferenceSinglePoint`), including an analytic-vs-FD
`dn_l/dp_L` tangent check. This is the right pattern and should be the template
for the strengthened history tests.

---

## 3. Recommended unit tests

Priority order. All are GP-local, no mesh/solver needed, and belong in
`DSMMicroMacroSingleIntegrationPoint.cpp` or a new
`SaturationTuller`-dedicated test file.

**P0 — restore the record**
1. Commit `DSMMicroMacroOverlapTransferBaseline.csv` and
   `DSMMicroMacroStrainCoupledOverlapBaseline.csv`; confirm the two history
   tests pass.
2. Upgrade both history tests to assert computed-vs-recorded for every
   physically meaningful column (§2.2).

**P1 — gate invariants (currently untested branches)**
3. Porosity-split clamp: drive `n_l → phi` under confinement so
   `phi_m > phi_total` and assert post-clamp `phi_M >= 0`,
   `phi_m <= phi_total`, `phi_m + phi_M == phi_total`. Add a sibling test that
   the clamp fires (so the diagnostic from §1.1, once added, is covered).
4. `computeVanDerWaalsMicroPotential` guard deaths: one
   `EXPECT_DEATH`/`OGS_FATAL`-equivalent per guard (`n_l<=0`, `rho_lR<=0`,
   `nS<=0`, `rho_SR<=0`, `A<=0`, `Sa<=0`, `K<0`, `lambda<=0 while K>0`).
5. `SaturationTuller` constructor guard deaths: bad bounds
   (`S_L_res > S_L_max`, `> 1`), non-positive `area_factor`/`shapefactor`/
   `pore_size`/`surface_tension`, `cavitation <= pressure_tolerance`.

**P2 — analytic-derivative correctness (FD cross-checks)**
6. `computeVanDerWaalsMicroPotential`: `dmu_lR_dnl`, `dmu_lR_drho_lR`,
   `dmu_lR_dnS`, `dmu_lR_drho_SR` (and the `omega_l` partials) vs central
   finite differences, with and without the exponential augmentation.
7. `SaturationTuller::dValue`/`d2Value` vs central FD across a pressure sweep
   spanning plateau, drained branch, and cavitation freeze.
8. `computeImplicitNlDpL` tangent vs FD for **all three**
   `LocalNonlinearSolveMode` values and both `MacroPorosityUpdateMode` values
   (ReferencePath currently covers only ScalarExchange + AlgebraicSplit; the
   mass-storage modes are only touched by the broken history tests).

**P3 — closed-form / invariant properties**
9. vdW augmentation additivity: with `K = 0` the result is byte-identical to
   the pure-vdW form; with `K > 0` the increment equals
   `sign · K · exp(−h/λ)` at a known state.
10. Sign convention: `NegativeAttractive` flips the sign of `mu_lR` relative to
    `PositiveReduced` at an identical state.
11. `computePotentialDrivenMassExchange`: `rho_L_hat == −rho_l_hat` and the
    three exact derivatives at representative `(alpha_M, mu_LR, mu_lR)`.
12. Young–Laplace branch continuity at `p_LR = -ptol` for `ptol = 0`, and the
    documented finite jump for `ptol > 0` (lock the magnitude `ptol/rho_LR`).
13. After fixing §2.3, an end-to-end reference-oracle test in
    `CurrentPorositySplit` mode that locks `active_nS = 1 − n_l` through the
    exchange (guards against regressing the documented active-fraction fix).

---

## 4. Summary

- The **gates** are individually well-guarded (vdW, Tuller, branch tolerances)
  and the analytic derivatives that were checkable are correct. The two real
  issues are the **duplicated, silently-lossy porosity-split clamp** (§1.1) and
  minor `ptol > 0` discontinuities (§1.2, §1.4).
- The **record** has one blocking defect — **two committed tests depend on CSV
  baselines that are not in the branch** (§2.1) — plus a soft defect: the
  history tests **don't actually assert against the recorded values** (§2.2),
  and a dormant stale oracle (§2.3).
- Highest-value next steps: commit the missing CSVs and make the history tests
  real regression checks, then backfill the untested gate/guard branches and
  the multi-mode tangent checks.
