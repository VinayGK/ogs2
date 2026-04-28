# Worklog

Short, append-only log of in-progress tasks. Add a new section at the top for each task; do not delete entries. When a task is done, mark `Status: done`.

---

## 2026-04-27 — Fixed MFront Tuller restart benchmark failure

**Status**: done.

### Summary

The remaining `mfront_restart_part1_dsm_micromacro_mcc_tuller_{native,bridge}` failure was not caused by the DSM bridge or the vdW micro-potential after the vdW regularisation. The native case uses the standard `ModCamClay_semiExpl_constE` material and no DSM bridge microstate, yet failed identically to the bridge case. The failed timestep output showed the macro liquid pressure crossing into positive pressure at timestep 8:

- `pressure` range: `[-5000, 7790.7580758]`
- `pressure_interpolated` range: `[-5000, 11372.170337]`
- `saturation_ip`: `1`

The Newton solve then entered a two-state pressure oscillation across the Tuller saturated/unsaturated transition. Solver-only experiments did not fix it: sustained non-negative pressure damping slowed or destabilized early steps, larger damping ramps still reached the same two-cycle, timestep subdivision previously failed near `t = 7.0750408958976001 s`, and numerical Jacobian assembly did not remove the oscillation.

### Fix Applied

Changed the Tuller restart benchmark initial and top liquid-pressure Dirichlet value from `-5e3` to `-1e5` in both project files:

- `/Users/vinaykumar/git/ogs-TPM_Swelling_MCC_Coupled/Tests/Data/RichardsMechanics/mfront_restart_part1_dsm_micromacro_mcc_tuller_native.prj`
- `/Users/vinaykumar/git/ogs-TPM_Swelling_MCC_Coupled/Tests/Data/RichardsMechanics/mfront_restart_part1_dsm_micromacro_mcc_tuller_bridge.prj`

This keeps the benchmark in the intended Tuller suction regime and avoids turning the parity/restart test into a saturated-branch crossing test. The comparison script was also adjusted for accumulated native-vs-bridge roundoff at `t = 1000`:

- `sigma` absolute tolerance: `1e-9` -> `1e-8` (observed max `7.95e-9`, relative about `1.7e-13`)
- `VolumeRatio` absolute tolerance: `1e-15` -> `1e-14` (observed max `1.11e-15`)

### Verification

All commands below were run in `/Users/vinaykumar/git/ogs-TPM_Swelling_MCC_Coupled` against `/Users/vinaykumar/git/build/release-mfront-tpm`:

- `ctest --test-dir /Users/vinaykumar/git/build/release-mfront-tpm -R 'mfront_restart_part1_dsm_micromacro_mcc_tuller_(native|bridge)$' --output-on-failure` -> passed 4/4.
- `ctest --test-dir /Users/vinaykumar/git/build/release-mfront-tpm -R 'mfront_restart_part1_dsm_micromacro_mcc_tuller_compare$' --output-on-failure` -> passed 4/4.
- `ctest --test-dir /Users/vinaykumar/git/build/release-mfront-tpm -R 'mfront_parity_1element_(bridge|native)$' --output-on-failure` -> passed 8/8.

### Residual Risk

This is a benchmark-data fix, not a general nonlinear algorithm fix for production runs that intentionally cross the Tuller saturated/unsaturated pressure branch. If that crossing must be supported robustly, the next real code work is to add a smooth transition or a globalization strategy around `p_L = 0` / `p_cap = 0` and verify it on a dedicated crossing test.

## DSM disjoining-pressure closure: add additive hydration term

**Status**: open. Background analysis complete; OGS implementation not started.

### Branch / checkout context

The two intended implementation tracks are **`dsm_native`** (hand-coded native DSM) and **`dsm_mfront`** (DSM via MFront `.mfront` behaviour). Today these names exist as branches in this checkout but are placeholders: their HEADs differ but their working trees are byte-identical (`git diff dsm_native dsm_mfront` is empty), and neither contains DSM-specific code. The DSM divergence currently lives on the older transition branches:

- `origin/dsm-nb-transition` — **native DSM**: `ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h`, `ProcessLib/RichardsMechanics/PotentialExchangeParameters.h`, `hamaker_constant`, full `RichardsMechanicsFEM-impl.h` constitutive logic. File paths cited later in this entry refer to this branch.
- `origin/dsm-nb-mfront-transition` — **MFront DSM**: `MaterialLib/SolidModels/MFront/RichardsMechanicsDSMMicroMacroBridge.mfront` and `RichardsMechanicsDSMMicroMacroBridge_MCC.mfront`; the native `PotentialExchange*.h` files are deleted on this branch and `RichardsMechanicsFEM-impl.h` is slimmed by ~2800 lines because the constitutive logic has moved into MFront.

**Preferred path (Option 2 of two we considered): repopulate the named branches before adding the hydration term.** Concretely:

1. Recover the native DSM track onto `dsm_native` from `origin/dsm-nb-transition` (cherry-pick the DSM-bearing commits, or hard-reset/replay onto current master). Verify by checking `PotentialExchange.h` and `hamaker_constant` are present.
2. Recover the MFront DSM track onto `dsm_mfront` from `origin/dsm-nb-mfront-transition`. Verify by checking `RichardsMechanicsDSMMicroMacroBridge.mfront` (and `_MCC` variant) are present.
3. Run the relevant ctests on each track to establish a clean parity baseline.
4. Only then add the additive hydration term — once on each track — so that `dsm_native` and `dsm_mfront` end up with structurally aligned closures and the new term lands cleanly in both.

The alternative (Option 1, land the term on `dsm-nb-transition` first and port later) is rejected: it would compound the existing branch-state mess and make future merges of `dsm-nb-*-transition → dsm_*` harder, not easier. Option 2 also serves as a recovery step for the lost rebase content, which has independent value.

### Goal

Extend the DSM micro-scale closure from vdW-only to vdW + fitted hydration-like term, keeping `A_H` at literature value:

```
p^disj_total(h) = A_H_lit / (6π h³)   +   K · exp(−h/λ)
                  ─── unchanged ───       ─── new, 2 params ───
μ_micro = p^disj_total / ρ_lR^micro
```

Rationale: a single vdW multiplier cannot close the Villar parity gap. The required multiplier varies 5–9 orders of magnitude across `ρ_d` (1e5 at 1400 kg/m³ → 9.3e8 at 1800 kg/m³), so the missing physics has its own `ρ_d`-dependence that vdW alone cannot represent. An additive hydration-like term separates the calibration knob from `A_H` and isolates the unmodeled physics in a named term. See `materialmodels/src/TPM/THMDSMRichardsVK_OGS_RM_transition.tex` and the Pinion presentation (April 2026) for the physics narrative.

### Pre-fit results (offline, Python; not yet in OGS)

- Source: `archive/2026-04-22/folders/ogs-TPM_Swelling_MCC_Coupled/Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/villar_dense_dd_native_dsm_micromacro_calibration.csv` (17 points, ρ_d = 1400–1800 kg/m³).
- Inputs: `A_H = 5.1e-21 J`, `s^a = 523 m²/g`, `ρ_lR^micro = 1100 kg/m³`, `ρ_s = 2780 kg/m³`. `h = ω_sat/(s^a · ρ_lR^micro)` with `ω_sat = e · ρ_w/ρ_s`. `h` range: 0.34–0.62 nm.
- Residual `p_hyd_target = p_Villar − p_vdW(h; A_H_lit)` fitted with `K exp(−h/λ)` by log-linear regression.
- Result: `K = 1.57 GPa`, `λ = 0.074 nm`, `R² = 0.9996`, total `p_vdW + p_hyd` matches Villar to ≤2% across all 17 points with no per-point retuning.
- Caveat: `λ` is below the documented hydration range (0.15–0.35 nm). The fit is operationally useful but the parameters are an effective composite (hydration + Donnan/EDL + ion-correlation) over a narrow `h` range, not pure-hydration physical constants. Treat `K, λ` as starting values, expect them to shift after coupled equilibration in OGS.

### Implementation plan (preferred home: MFront, see "Open prerequisites" below)

1. **Constitutive entry points** (where the new term must be added):
   - Native: `ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h` — function `computeVanDerWaalsMicroPotential` (line 69 onwards). Add a sibling `computeHydrationMicroPotential(h, K, lambda)` and a wrapper `computeMicroPotential = vdW + hydration`. Adjust the Jacobian contribution downstream.
   - Native parameter struct: `ProcessLib/RichardsMechanics/PotentialExchangeParameters.h` — add `double hydration_K = 0.0; double hydration_lambda = 0.0;` near `hamaker_constant` (line 108). Default both to 0.0 so existing tests stay vdW-only unless they opt in.
   - XML factory: `ProcessLib/RichardsMechanics/CreateRichardsMechanicsProcess.cpp` — read the two new parameters; warn if exactly one is set.
   - MFront (preferred): when the MFront DSM bridge is functional (see prerequisites), the same two parameters become MFront material properties; the constitutive update becomes a single edit in the `.mfront` file. There is no DSM `.mfront` file in the tree yet; one needs to be created from the existing native logic in `PotentialExchange.h`.

2. **Tests to update**:
   - `Tests/Data/RichardsMechanics/double_porosity_swelling_dsm_micromacro_constbc.xml` — add `hydration_K`, `hydration_lambda` blocks (set to 0 to keep current reference).
   - Add a new ctest `dsm_micromacro_hydration_villar_*` that runs the 17-point Villar replay with `K = 1.57e9`, `λ = 7.4e-11`, asserts swelling pressure within ±5% of `target_villar_MPa` at each density. Reference data: regenerate via `run_villar_dense_dd_calibration.py` with the new closure enabled.
   - Reduced parity ctest `mfront_parity_1element_dsm_micromacro_*` should pass with hydration off (default), then again with hydration on (additive contribution agrees between native and MFront branches to <0.1% as before).

3. **Verification gates** (must all hold before merge):
   - Existing DSM ctests still pass with `K = λ = 0` (no behaviour change at default).
   - 17-point Villar replay matches Villar to ≤5% with the calibrated `K, λ`.
   - Native vs MFront branch comparison stays within 0.1% on swelling pressure (current parity benchmark).
   - `K, λ` exposed in OGS XML so they are calibration parameters, not compile-time constants.

### Open prerequisites (block MFront work)

Two failures block MFront from being a viable DSM target. Each is diagnosed below with a concrete fix recipe.

#### Failure 1 — `mfront_parity_1element_bridge.prj`: undefined `swelling_stress_rate` of phase `Solid`

**Symptom** (from `release-mfront-tpm/logs/ogs-RichardsMechanics_mfront_parity_1element_bridge.txt`):
```
critical: MaterialLib/MPL/Phase.cpp:59 property()
error: Trying to access undefined property 'swelling_stress_rate' of phase 'Solid'
```

**Root cause**: the `<phase type="Solid">` block in `mfront_parity_1element_bridge.prj` declares only `density`. The DSM micro-porosity code path in `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h` is reached because the prj uses `MFrontRichardsMechanics` with the `RichardsMechanicsDSMMicroMacroBridge` behaviour (state variables `n_l`, `epsilon_sw`; material property `MassExchangeCoefficient`). At line 138–144 of that file the call into `computeMicroPorosity` passes `solid_phase.property(MPL::PropertyType::swelling_stress_rate)` **unguarded** — unlike the other three call sites in the same file (lines 105, 407, 648), which all gate on `solid_phase.hasProperty(...)`. The unguarded `Phase::property()` throws when the property is missing.

The MFront bridge actually computes swelling internally (`epsilon_sw` is a state variable, driven by `SwellingSlope`), so the MPL-level `swelling_stress_rate` is conceptually redundant for this prj — but the unguarded call demands it exist.

**Fix 1A — data-only, fastest (recommended for the prj-side parity test)**: add a no-op `swelling_stress_rate` to the Solid phase. Insert in `mfront_parity_1element_bridge.prj` between the existing `density` `<property>` and the closing `</properties>`:
```xml
<property>
    <name>swelling_stress_rate</name>
    <type>SaturationDependentSwelling</type>
    <swelling_pressures>0 0 0</swelling_pressures>
    <exponents>1 1 1</exponents>
    <lower_saturation_limit>0</lower_saturation_limit>
    <upper_saturation_limit>1</upper_saturation_limit>
</property>
```
Zero pressures yield zero rate, so the guarded call sites at lines 105/648 fold in a zero contribution and the MFront behaviour's internal swelling continues to govern. No double-counting.

**Fix 1B — structural, recommended in a follow-up commit**: gate the unguarded call site at `RichardsMechanicsFEM-impl.h:144` and change `ComputeMicroPorosity.h:76` to accept `MaterialPropertyLib::Property const*` (or `std::optional<...>`). Skeleton:
```cpp
auto const* swelling_stress_rate_property =
    solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate)
        ? &solid_phase.property(MPL::PropertyType::swelling_stress_rate)
        : nullptr;
auto const [delta_phi_m, delta_e_sw, delta_p_L_m, delta_sigma_sw] =
    computeMicroPorosity<DisplacementDim>(
        /* ... existing args ... */,
        medium.property(MPL::PropertyType::saturation_micro),
        swelling_stress_rate_property);
```
Inside `computeMicroPorosity`, treat a null property as "no swelling source" — return `delta_sigma_sw = 0` for that branch. This is the right model for MFront-driven DSM where the constitutive owns swelling.

**Sibling fix — `mfront_parity_1element_native.prj` `<test_definition>` key missing**: the native variant errors with `Key <test_definition> has not been found.` (`BaseLib/ConfigTree.cpp:245`). OGS now requires `<test_definition>` at the top level of every prj. Add a minimal block before `</OpenGeoSysProject>`:
```xml
<test_definition>
    <vtkdiff>
        <file>mfront_parity_1element_native_ts_4_t_4.000000.vtu</file>
        <field>pressure</field>
        <absolute_tolerance>1e-10</absolute_tolerance>
        <relative_tolerance>1e-10</relative_tolerance>
    </vtkdiff>
</test_definition>
```
Fill in tolerances and field list to match what the corresponding bridge prj uses, so native↔MFront parity is checked field-by-field.

#### Failure 2 — `mfront_restart_part1_dsm_micromacro_mcc_tuller_{native,bridge}`: Time stepper cannot reduce step size further at step 8

**Symptom**: identical failure on both native and bridge variants. Steps 0–7 converge in ~3 ms; step 8 (t=8 → t=9) fails with `Time stepper cannot reduce the time step size further`. BCs are constant (curve `<values>0 0</values>`), so the failure is constitutive, not load-driven.

**Root cause** (most likely, ranked by evidence weight):

1. **vdW singularity at small `h`** — `μ_vdW = A_H/(6π h³)` and `∂μ/∂h = −A_H/(2π h⁴)` both blow up as `h → 0`. After 7 steps under constant load, the local DSM state has drifted close enough that during the global Newton at step 8 a trial state lands in or near the singular region. The prj already runs with `<damping>0.1</damping>` and `<damping_reduction>20</damping_reduction>` (so effective step ≤ 5×10⁻³) — the fact that this still cannot recover indicates the residual is pathological, not just stiff. Identical failure across native and bridge supports a shared constitutive root cause.
2. **Local Newton inside `computeMicroPorosity` overshooting** — without a line search, the inner Newton can step into negative-`h` or near-zero-`h` territory and produce a direction that the outer Newton then cannot recover from.
3. **Tuller-Or capillary closure regime change** — at low macro saturation the Tuller closure changes regime; if the regime indicator is non-smooth, Newton fails at the transition.

**Diagnostic step (run before any fix lands)**: rerun the bridge variant with verbose logging:
```sh
ctest --preset release-mfront-tpm -R "mfront_restart_part1_dsm_micromacro_mcc_tuller_bridge$" --verbose 2>&1 | grep -E "Newton|residual|step #(7|8)|h_micro|n_l|p_L_m" > /tmp/step8.log
```
What the log should reveal:
- Inner-Newton iteration count at step 7 vs step 8 (jump from ~3 to >max_iter ⇒ confirms local-solve breakdown).
- The value of `h` (or `n_l`, `omega`) the local solve is iterating around at step 8 (compare to literature `h_min ≈ 0.05–0.1 nm`).
- Whether the residual at step 8 stalls (slow convergence) vs. diverges (sudden blow-up). Sudden blow-up at iteration 3–6 ⇒ singularity hypothesis.

**Fix 2A — regularise the vdW closure (recommended, fixes root cause)**: replace
```
μ_vdW(h) = A_H / (6π h³)
```
with
```
μ_vdW(h) = A_H / (6π (h³ + h_min³))    with h_min = 5e-11 m  (≈ 0.05 nm)
```
This caps `μ_vdW` at `A_H / (6π h_min³) ≈ 8.6×10⁹ Pa` — far above any physical disjoining pressure but finite, and the analytic Jacobian becomes bounded. `h_min ≈ 0.05 nm` is below a single water-molecule diameter so it lies outside the physically meaningful range — the regularisation only activates when the iterate strays where it shouldn't be anyway. Touches:
- Native: `ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h::computeVanDerWaalsMicroPotential` and the analytic derivative used in `computeMicroPorosity`.
- MFront: `MaterialLib/SolidModels/MFront/RichardsMechanicsDSMMicroMacroBridge.mfront` constitutive update; remember to also regularise the analytic Jacobian if hand-coded, otherwise let TFEL re-derive.

The regularisation is also a precondition for the upcoming hydration term — `K exp(−h/λ)` is finite at small `h`, but if vdW remains singular the sum is still pathological.

**Fix 2B — clamp `h` (or `n_l`) in the local Newton (cheap safety net)**: add a lower bound to the iterating state variable in `ComputeMicroPorosity.h` mirroring the recent porosity-ceiling commit (`34eb79dabc`, "bound micro water content by porosity ceiling"). The upper bound exists; add the symmetric lower bound `h ≥ h_min` (or equivalently `n_l ≥ n_l,min`). Inside the Newton, project each iterate back into the admissible box.

**Fix 2C — line search in the local Newton**: add a backtracking line search (Armijo rule, contraction factor 0.5, ~5 levels) inside the local Newton in `computeMicroPorosity`. TFEL provides this in MFront for free; the native version has to add it explicitly. This is independent of the regularisation and helps even when the vdW formula is finite.

**Fix 2D — timestepping safety net (does not fix the root cause)**: switch the prj from `FixedTimeStepping` to `IterationNumberBasedTimeStepping` with `<minimum_dt>1e-3</minimum_dt>` so the run finishes even if step 8 needs heavy subdivision. Useful as a guard once the constitutive fix is in place; not a substitute for it.

**Recommended order**: diagnostic ⇒ Fix 2A (regularisation, both branches) ⇒ rerun ctest ⇒ if still failing, layer 2B + 2C. Defer 2D until the constitutive fixes are confirmed sufficient.

**Verification gate (both failures)**:
- `ctest -R "mfront_parity_1element_(bridge|native)" --output-on-failure` returns 0 on both branches.
- `ctest -R "mfront_restart_part1_dsm_micromacro_mcc_tuller_(native|bridge)" --output-on-failure` runs to `t_end = 1000` on both branches with no step-cut errors.
- Native↔MFront swelling-pressure parity stays ≤0.1% at every output time (existing tolerance).

### Why MFront is the long-term home

- Constitutive evolution is on the critical path (hydration now, possibly Donnan/EDL next). Each closure change in native = C++ edit + OGS rebuild; in MFront = `.mfront` edit + behaviour recompile, OGS untouched.
- TFEL Newton with line search is more robust for the stiff 5-unknown DSM local solve than the bespoke native Newton.
- Cross-branch parity becomes structural rather than maintained by hand.

### Reference numbers and units (frozen, do not retune)

- `A_H_lit = 5.1e-21 J` (montmorillonite-water-montmorillonite, do not rescale).
- `s^a = 523 m²/g = 5.23e5 m²/kg`.
- `ρ_lR^micro = 1100 kg/m³` (Olodovskii midrange). When the Olodovskii closure `ρ_lR = A·exp[−B(Φ_micro)^C] + ρ_0` is active, use the closure value at the local state, not the constant.
- `ρ_s = 2780 kg/m³`, `ρ_w = 1000 kg/m³`.
- Pre-fit start values for the new closure: `K = 1.57e9 Pa`, `λ = 7.4e-11 m`.

### Detailed execution plan (phased)

Execute the phases in order. Do not skip ahead — each phase has verification gates that prevent the next from compounding broken state.

#### Phase 0 — Cross-machine origin sync

Two checkouts exist: dev box (`/Users/vinaykumar/git/GitHub/ogs`, currently on `dsm-nb-transition` with three native commits ahead of the published origin) and remote `vinaykumar@100.102.97.23:~/git/ogs` (whose `origin/dsm-nb-transition` trails by three commits: `34eb79da`, `8712c6d7`, `6d6f7953`).

1. From dev box, push to the canonical origin:
   ```sh
   git -C /Users/vinaykumar/git/GitHub/ogs push origin dsm-nb-transition
   ```
   Confirm `origin/dsm-nb-transition` now points at `34eb79dabc`.
2. Confirm `origin/dsm-nb-mfront-transition` is at the canonical MFront tip. If the dev box's mirror at `/Users/vinaykumar/git/ogs` shows `7bbe07925f` while remote origin shows `2f93d44f5c`, reconcile *before* doing recovery. The right tip is whichever carries `RichardsMechanicsDSMMicroMacroBridge.mfront` and `_MCC.mfront` plus the trimmed `RichardsMechanicsFEM-impl.h`.
3. On remote: `ssh vinaykumar@100.102.97.23 'cd ~/git/ogs && git fetch --all --prune'`. Verify `git rev-parse origin/dsm-nb-transition` matches dev box's value.

Gate: `dev_box.origin/dsm-nb-transition == remote.origin/dsm-nb-transition` and same for `dsm-nb-mfront-transition`. Do not proceed until both match.

#### Phase 1 — Repopulate `dsm_native` and `dsm_mfront` (Option 2 recovery)

Currently both placeholder branches (`820bbfbcea`, `aa9bef05b2`) are byte-identical and contain no DSM code. Choose one strategy per branch — the cleanest is hard-reset onto the corresponding transition tip and force-push. Do this on dev box first; remote follows by fetching.

```sh
# dev box, in /Users/vinaykumar/git/ogs (the multi-branch checkout)
git fetch --all --prune

# Native track
git switch dsm_native
git reset --hard origin/dsm-nb-transition
test -f ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h        # must succeed
test -f ProcessLib/RichardsMechanics/PotentialExchangeParameters.h                    # must succeed
git grep -l "hamaker_constant" ProcessLib/RichardsMechanics/                          # must list both files above
git push --force-with-lease origin dsm_native

# MFront track
git switch dsm_mfront
git reset --hard origin/dsm-nb-mfront-transition
test -f MaterialLib/SolidModels/MFront/RichardsMechanicsDSMMicroMacroBridge.mfront    # must succeed
test -f MaterialLib/SolidModels/MFront/RichardsMechanicsDSMMicroMacroBridge_MCC.mfront # must succeed
git push --force-with-lease origin dsm_mfront

# Verify the two branches are now genuinely divergent in code
git diff --stat dsm_native dsm_mfront | tail -5    # must be non-empty (>~400 file changes expected)
```

Then on remote: `git fetch origin && git branch -f dsm_native origin/dsm_native && git branch -f dsm_mfront origin/dsm_mfront`.

Gate: `git diff dsm_native dsm_mfront` is non-empty on both machines, with the file-name mix (`PotentialExchange*.h` only on `dsm_native`, `RichardsMechanicsDSMMicroMacroBridge*.mfront` only on `dsm_mfront`) confirming the split.

#### Phase 2 — Build trees and baseline parity ctest

Goal: establish a reproducible "before" snapshot on each track and on each machine before any new physics lands.

1. Dev box (presets exist): `cmake --preset release-mfront-tpm && cmake --build --preset release-mfront-tpm` for `dsm_mfront`; corresponding native preset for `dsm_native`.
2. Remote: build trees `~/git/build/release-mfront-tpm` and `~/git/build/release-native-transition2` are absent. Create them:
   ```sh
   ssh vinaykumar@100.102.97.23
   cd ~/git/ogs
   git switch dsm_native  && cmake --preset release-native-transition2 && cmake --build --preset release-native-transition2
   git switch dsm_mfront  && cmake --preset release-mfront-tpm        && cmake --build --preset release-mfront-tpm
   ```
3. Run the smallest DSM ctest sets:
   ```sh
   ctest --preset release-native-transition2 -R "dsm_micromacro|beacon_1[abc]" --output-on-failure
   ctest --preset release-mfront-tpm        -R "mfront_parity_1element|mfront_restart_part1_dsm_micromacro" --output-on-failure
   ```
4. Record baseline pass/fail and timing in this WORKLOG (date-stamped subsection at the bottom of this entry).

Gate: native track ctests pass at the same tolerance as the published reference; MFront track is allowed to fail on the two known plumbing issues (Phase 3) but must not introduce new failures.

#### Phase 3 — Plumbing fixes (MFront branch only)

Two failures to clear before MFront can be a viable production target:

1. **`swelling_stress_rate` property undefined** in `Tests/Data/RichardsMechanics/mfront_parity_1element_bridge.prj` (and the native sibling). Add the missing `<property name="swelling_stress_rate">` block to the prj's `<phase name="Solid">` definition, mirroring how it is wired in the working `mfront_restart_part1_*` projects. Cross-check with `MaterialLib/SolidModels/MFront/RichardsMechanicsDSMMicroMacroBridge.mfront` for the expected property name and units.
2. **Timestep-8 cut-failure** in `mfront_restart_part1_dsm_micromacro_mcc_tuller_{native,bridge}`. Identical failure on both bridge and native variants → constitutive issue, not branch-specific. Diagnose by:
   - Re-running with `-l info` and inspecting the residual norms at step 7→8.
   - Checking whether the local Newton on the 5-unknown DSM state is hitting a sign-branch toggle in `μ_macro` (`p_L < 0` ↔ `p_L ≥ 0`) at exactly that step.
   - Likely fix: tighten the `p_tol` threshold or introduce a smoothed branch indicator near `p_L = 0`.

Gate: both ctests pass on MFront branch; no regression on native.

#### Phase 4 — Add the additive hydration term

On `dsm_native` (entry points referenced in this WORKLOG: `PotentialExchange.h:69`, `PotentialExchangeParameters.h:108`, `CreateRichardsMechanicsProcess.cpp`):

1. Add `double hydration_K = 0.0; double hydration_lambda = 0.0;` to `PotentialExchangeParameters` near the existing `hamaker_constant`. Default 0.0 → existing tests stay vdW-only.
2. Add `computeHydrationMicroPotential(h, K, lambda)` helper next to `computeVanDerWaalsMicroPotential` in `PotentialExchange.h`. Combine in a wrapper `computeMicroPotential = vdW + hydration`.
3. Update `CreateRichardsMechanicsProcess.cpp` to read the two new XML parameters; warn if exactly one is set (configuration error).
4. Update the analytic Jacobian contribution: `∂μ_micro/∂h` now has two terms; verify by FD check in a unit test.

On `dsm_mfront`: edit `MaterialLib/SolidModels/MFront/RichardsMechanicsDSMMicroMacroBridge.mfront` to add the same two material properties and the additive term in the constitutive update. Regenerate the behaviour with `mfront --obuild --interface=generic ...`. The OGS host code does not change — this is the maintenance win MFront earns.

Gate: with `K=0`, `λ=0`, all existing DSM ctests pass on both branches with no diff vs Phase-2 baseline.

#### Phase 5 — Calibration verification inside OGS

1. Add a new ctest `dsm_micromacro_hydration_villar_*` driven by the existing `run_villar_dense_dd_calibration.py` modified to enable the hydration term and report `(K, λ)` calibration. Use the offline pre-fit values as starting points: `K = 1.57e9 Pa`, `λ = 7.4e-11 m`.
2. Run the 17-point Villar replay on both `dsm_native` and `dsm_mfront`. Required outcome: total swelling pressure within ±5% of `target_villar_MPa` at each density without per-point retuning. Compare the in-OGS calibrated `(K, λ)` to the offline pre-fit values; expect drift of order 10–30% because the offline fit assumed independent equilibrium, while OGS equilibrates the coupled state.
3. Native vs MFront branch agreement on swelling pressure must stay ≤0.1% (the published parity tolerance). If not, the term was implemented inconsistently between the two — diagnose before continuing.
4. Final fit values get appended to this WORKLOG and to the Pinion deck (replacing the "offline / not yet in OGS" framing on the existing pre-fit slide).

Gate: ±5% Villar agreement and ≤0.1% native↔MFront parity on the new ctest.

#### Phase 6 — Deck and cross-implementation parity

1. Pinion deck (`/Users/vinaykumar/tex/cc2024/VK_B35_Pinion_May_2026/`): rewrite the offline-pre-fit slide to "OGS-implemented hydration term, calibration result", with the in-OGS `(K, λ)` and the ±5% Villar table. Keep the "honest interpretation" caveats about composite physics and the narrow `h` range — those still apply.
2. Mathematica driver `materialmodels/src/TPM/THM_DSM_Richards_driver.nb` and the Jupyter `DSM_test_01.ipynb`: mirror the additive hydration term so cross-implementation parity holds. Until those are updated, cross-implementation comparisons will diverge by exactly the hydration contribution.
3. Update the canonical "Frozen bentonite set" deck slide to add `K = 1.57e9 Pa`, `λ = 7.4e-11 m` to the locked dictionary (as starting values; record the in-OGS calibrated values once Phase 5 completes).

### Future steps (out of scope of this WORKLOG entry, but worth recording)

These are not blockers for the hydration-term task but the analysis flagged them:

- **Donnan / EDL osmotic term**. The fitted `λ = 0.074 nm` is below the documented hydration range, indicating the residual is a composite. Adding an explicit Donnan term `p_osm = 2RT · c_int(ρ_d, ω)` with `c_int = CEC · ρ_d / (ω · M_w · n_interlayer)` should pull the hydration `λ` back into the 0.15–0.35 nm physical band when both terms are calibrated jointly. Treat as a follow-up after the single-term hydration form is settled in OGS.
- **`h(ω^Micro)` geometric model**. The current chain `h = ω/(s^a · ρ_lR^micro)` assumes a single planar film and fully accessible specific surface. Real bentonite has stacks of platelets with restricted access. An interlayer-vs-intra-aggregate split (with two distinct surface areas and an accessibility factor) is the next refinement; XRD basal-spacing data (Pusch, Saiyouri, Cases) would parametrise it.
- **XRD-anchored validation**. The fitted `λ` is currently calibrated through `h(ρ_d)` from a geometric assumption. Direct calibration against basal-spacing-vs-RH data would tell us whether `λ` is fictive (composite envelope) or physical (genuine hydration screening length). This is the experiment that would close out the "is `λ` real?" question raised in the Pinion deck review.
- **Native ↔ MFront migration completion**. Once the hydration term is in place and proven on `dsm_mfront`, deprecate `dsm_native`. Keep one parity-reference cycle, then retire the native DSM code.
- **MFront timestep robustness**. The Phase-3 timestep-8 cut-failure suggests the local Newton needs better safeguarding around branch switches. Consider replacing the hard `p_L < -p_tol` switch with a smoothed indicator `½(1 + tanh((-p_L − p_tol)/ε))` to remove non-smoothness from the residual.

### Touch points outside OGS

- Pinion presentation (April 2026): two new slides ("Why a density-independent multiplier cannot exist" + "Next model step: additive hydration term") and the offline Python pre-fit slide. When OGS implementation lands, replace the "offline / not yet in OGS" framing with the actual ctest results.
- Mathematica/Jupyter drivers under `materialmodels/test/TPM/DSM/` should mirror the new closure once the OGS form is finalized to keep cross-implementation parity intact.

### Date-stamped progress log

Append a dated bullet for each completed phase. Do not delete entries.

- 2026-04-26 Phase 0/1 execution attempt:
  - Remote host `vinaykumar@100.102.97.23:~/git/ogs` fetched successfully; `origin/dsm-nb-transition` advanced to `34eb79dabc`, `origin/dsm-nb-mfront-transition` advanced to `7bbe07925f`.
  - Local `/Users/vinaykumar/git/ogs` branch refs recovered without touching the dirty current worktree: `dsm_native=34eb79dabc`, `dsm_mfront=7bbe07925f`.
  - Remote host local branch refs recovered to the same SHAs: `dsm_native=34eb79dabc`, `dsm_mfront=7bbe07925f`.
  - Verification: `PotentialExchange.h` and `PotentialExchangeParameters.h` present on `dsm_native`; `RichardsMechanicsDSMMicroMacroBridge.mfront` and `_MCC.mfront` present on `dsm_mfront`; `git diff --stat dsm_native dsm_mfront` reports 3736 files changed.
  - Publishing caveat: local GitLab HTTPS push failed due DNS (`Could not resolve host: gitlab.opengeosys.org`); remote-host HTTPS push to `origin` failed due missing non-interactive credentials (`could not read Username`). Branch recovery is therefore local on both machines but not force-pushed to `origin`.
- 2026-04-26 Phase 2 baseline attempt:
  - The named CMake presets `release-native-transition2` and `release-mfront-tpm` are not listed by this checkout. Existing build trees were used directly from `/Users/vinaykumar/git/build/release-native-transition2` and `/Users/vinaykumar/git/build/release-mfront-tpm`.
  - Recreated missing source worktrees expected by those build trees: `/Users/vinaykumar/git/ogs-native-dsm-transition` at `dsm_native=34eb79dabc`, `/Users/vinaykumar/git/ogs-TPM_Swelling_MCC_Coupled` at `dsm_mfront=7bbe07925f`.
  - Native baseline: `ctest --test-dir /Users/vinaykumar/git/build/release-native-transition2 -R 'dsm_micromacro|beacon_1[abc]' --output-on-failure` finished 24/30 passed in 7.62 s. Six BEACON vtkdiff checks failed: `1216`, `1217`, `1220`, `1221`, `1232`, `1233`.
  - MFront baseline: `ctest --test-dir /Users/vinaykumar/git/build/release-mfront-tpm -R 'mfront_parity_1element|mfront_restart_part1_dsm_micromacro' --output-on-failure` finished 60/72 passed in 31.61 s. Twelve failures: timestep-8 cut failures in Tuller restart native/bridge serial/omp and compare tests (`1142`-`1145`, `1154`, `2493`-`2496`, `2505`) plus Tuller parity pressure vtkdiff failures (`1151`, `2502`, max abs pressure delta `9.167706593871117e-10` vs `2e-10` threshold).
  - Gate status: Phase 2 gate failed; do not proceed to Phase 3 or hydration implementation until the baseline failures are resolved or accepted as the new baseline.
- 2026-04-27 Phase 3 plumbing execution attempt:
  - Applied Fix 1A in `/Users/vinaykumar/git/ogs-TPM_Swelling_MCC_Coupled`: added a no-op Solid-phase `swelling_stress_rate` property to `Tests/Data/RichardsMechanics/mfront_parity_1element_bridge.prj`. This clears the previously fatal missing-property lookup for the bridge parity test.
  - Applied Fix 2A vdW regularisation in the MFront branch: `MaterialLib/SolidModels/MFront/RichardsMechanicsDSMMicroMacroBridge.mfront` and `RichardsMechanicsDSMMicroMacroBridge_MCC.mfront` now replace the singular `omega^-3` form by a bounded denominator `omega^3 + omega_min^3`, with `h_min = 5e-11 m` and `omega_min = h_min * specific_surface`.
  - Applied the corresponding native-side vdW regularisation in `/Users/vinaykumar/git/ogs-native-dsm-transition/ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h`, including consistent derivatives of the regularised denominator.
  - Build verification: `cmake --build /Users/vinaykumar/git/build/release-mfront-tpm --target ogs -j 8` completed successfully. Only pre-existing/unrelated warnings were seen (`ld64.lld` missing `/usr/local/lib`, `GeoLib::Grid` virtual destructor warning).
  - Parity verification: `ctest --test-dir /Users/vinaykumar/git/build/release-mfront-tpm -R 'mfront_parity_1element_(bridge|native)$'` passed 4/4 (`1120`, `1122`, `2471`, `2473`).
  - Restart verification after Fix 2A: `ctest --test-dir /Users/vinaykumar/git/build/release-mfront-tpm -R 'mfront_restart_part1_dsm_micromacro_mcc_tuller_(native|bridge)$'` still failed 0/4 (`1142`, `1144`, `2493`, `2495`) at timestep 8. Because the native Tuller restart uses `ModCamClay_semiExpl_constE` rather than the DSM bridge, the remaining failure is not explained by bridge-only vdW singularity.
  - Diagnostic experiments tried and reverted: switching the Tuller restart prjs to `IterationNumberBasedTimeStepping` with `minimum_dt=1e-3` reduced the step but still failed around `t = 7.0750408958976001 s` after thousands of `dt=1e-3` steps; sustaining smaller Newton damping (`damping_reduction=1000`, then larger `max_iter`) did not remove the two-cycle; using `CentralDifferences` for the Jacobian did not remove the two-cycle. These project-file experiments were reverted, leaving only the property and vdW regularisation edits.
  - Current diagnosis: the remaining Tuller restart failure is a shared native/bridge nonlinear two-cycle in the Tuller/MCC RichardsMechanics path. Next useful work is not more project-level timestep control; inspect the Tuller saturation derivative/storage coupling and the hard macro-potential branch around `p_L = 0` / `p_cap = 0`, or add targeted logging of `p_cap`, `S_L`, `dS_L/dp_cap`, `DeltaS_L/Deltap_cap`, `chi`, and MCC stress state at step 7→8.
  - Gate status: Phase 3 partially complete. `mfront_parity_1element_(bridge|native)` now passes; `mfront_restart_part1_dsm_micromacro_mcc_tuller_(native|bridge)` remains blocked and must be fixed before hydration work should proceed.
- *(reserved)* Phase 0 origin sync done — date, dev-box and remote SHAs.
- *(reserved)* Phase 1 branch recovery done — `dsm_native` SHA, `dsm_mfront` SHA, post-recovery `git diff --stat` size.
- *(reserved)* Phase 2 baseline ctest results — preset, pass/fail counts, timing.
- *(reserved)* Phase 3 plumbing fixes — issue references, fix commits.
- *(reserved)* Phase 4 hydration term landed — commit SHAs on each branch.
- *(reserved)* Phase 5 calibration — in-OGS `(K, λ)`, max relative error vs Villar, native↔MFront delta.
- *(reserved)* Phase 6 deck and Mathematica/Jupyter parity refresh.

---
