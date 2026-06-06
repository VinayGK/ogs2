# DSM Hierarchical Porosity — Patch Recipe

**Purpose:** This file contains everything needed to re-implement the
`dsm_native_hierarchical` functionality on a fresh `master` checkout. It is
the reconstruction recipe. `AGENTS.md` is the status tracker; `HIERARCHICAL_POROSITY.md`
is the physics narrative. This file is the **how-to-apply** reference.

**Maintenance rule:** Every time code, tests, or benchmark PRJ files change on
`dsm_native_hierarchical`, update the relevant section of this file. Do not
leave this file stale. See the maintenance checklist at the bottom.

---

## Files changed (complete list)

| File | Nature of change |
|------|-----------------|
| `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h` | Core: hierarchical storage, Jacobian, swelling driver, viscosity guards, dead-code removal |
| `ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h` | vdW dimensional fix (`/rho_lR`), `dmu_lR_drho_lR` non-zero, literature A comment |
| `Tests/ProcessLib/RichardsMechanics/DSMMicroMacroSingleIntegrationPoint.cpp` | Unit test refactor — 13/13 passing after hierarchical changes |
| `Tests/Data/RichardsMechanics/ANCHORS_MS33_Model*/` | PRJ files: `hamaker_constant` 5.1e-21 → 2.2e-20 J, K recalibrated; VTU reference outputs updated |

---

## Build instructions (fresh master)

```bash
# Assumes OGS source at /Users/vinaykumar/git/ogs
# Build dir convention used on this machine:
BUILD=/Users/vinaykumar/git/build/release-omp-mfront

cmake -B $BUILD -S /Users/vinaykumar/git/ogs \
  -DCMAKE_BUILD_TYPE=Release \
  -DOGSEnabledElements="HM;Richards;TRM" \
  -DOGS_USE_MFront=ON \
  -DOGS_USE_PETSC=OFF

cmake --build $BUILD --target ogs testrunner -j$(nproc)
```

---

## Verification commands

```bash
BUILD=/Users/vinaykumar/git/build/release-omp-mfront

# Unit tests — must all pass before any commit
$BUILD/bin/testrunner --gtest_filter='*DSMMicroMacro*'

# Benchmark smoke — zero rejected steps on canonical LE set
$BUILD/bin/ogs Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj
$BUILD/bin/ogs Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj
$BUILD/bin/ogs Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1800.prj
$BUILD/bin/ogs Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj
$BUILD/bin/ogs Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj
$BUILD/bin/ogs Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj

# ctest for registered DSM tests
ctest --test-dir $BUILD --output-on-failure -R 'DSMMicroMacro'
```

**Passing criteria:**
- `testrunner` DSMMicroMacro: 13/13 pass.
- All canonical LE runs: 0 rejected steps.
- Model I Villar MAE ≤ 0.02% (currently 0.016%).

---

## Step-by-step patch diffs

Apply in order. Each section is one logical commit. The commit hash on
`dsm_native_hierarchical` is given for reference — use it for `git show` to see
the full diff including PRJ/VTU changes.

---

### Step 1 — REV-scale storage + hierarchical split consistency
**Commit:** `0d7a9edd64`  
**File:** `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`

**What changes and why:**

1. `computeActiveMicroSolidVolumeFraction`: the active solid volume fraction is
   `1 - phi_M` (macro skeleton), not `1 - phi_M - phi_m`. The micro pore space
   belongs to the pore fluid, not the solid.

```diff
-    return std::max(1e-16, 1.0 - split.phi_M - split.phi_m);
+    return std::max(1e-16, 1.0 - split.phi_M);
```

2. Previous active solid (fallback path from prev-step state): use only `phi_M_prev`,
   not `phi_M_prev + phi_m_prev`.

```diff
-    double const total_prev_porosity = std::clamp(
-        std::max(0.0, local_context.phi_M_prev) +
-            std::max(0.0, local_context.phi_m_prev),
-        0.0, 1.0 - 1e-12);
-    return std::max(1e-16, 1.0 - total_prev_porosity);
+    double const phi_M_prev_safe =
+        std::clamp(std::max(0.0, local_context.phi_M_prev), 0.0, 1.0 - 1e-12);
+    return std::max(1e-16, 1.0 - phi_M_prev_safe);
```

3. REV-scale liquid apparent density in local Newton residual. Previously `n_l *
   rho_lR` (aggregate-referenced). Now `phi_m * rho_lR` where
   `phi_m = (1-phi)/(1-n_l) * n_l` (hierarchical, REV-referenced).

```diff
-        double const rho_l = n_l * micro_liquid_density.rho_lR;
+        double const phi_h = std::isfinite(local_context.phi)
+            ? std::clamp(local_context.phi, 0.0, 1.0 - 1e-12)
+            : std::clamp(local_context.phi_M_prev + local_context.phi_m_prev,
+                         0.0, 1.0 - 1e-12);
+        double const one_minus_n_l_h = std::max(1e-12, 1.0 - n_l);
+        double const rho_l =
+            (1.0 - phi_h) / one_minus_n_l_h * n_l * micro_liquid_density.rho_lR;
```

4. Jacobian updated consistently: `d(rho_l_REV)/dn_l` now includes the
   `(1-phi)/(1-n_l)^2` factor.

```diff
-        double const jacobian = micro_density.drho_l_dn_l -
+        double const drho_l_REV_dn_l =
+            one_minus_phi_M_jac / one_minus_n_l_jac * micro_density.rho_lR +
+            one_minus_phi_M_jac * micro_density.drho_lR_dnl;
+        double const jacobian = drho_l_REV_dn_l -
                                 dt_safe * drho_l_hat_dn_l -
-                                dt_safe * micro_density.drho_l_dn_l *
+                                dt_safe * drho_l_REV_dn_l *
                                     volumetric_strain_rate;
```

5. `rho_l_prev` now uses `phi_m_prev * rho_lR_prev` (REV-scale, hierarchical):

```diff
     double const rho_l_prev =
         prev_micro_liquid_density
-            ? n_l_prev * prev_micro_liquid_density->rho_lR
+            ? local_context.phi_m_prev * prev_micro_liquid_density->rho_lR
             : 0.0;
```

---

### Step 2 — Thermodynamic swelling stress `sigma_sw = -phi_m * Pi`
**Commit:** `88d42c98fd`  
**File:** `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`

**What changes and why:**  
The Pi-path swelling stress is thermodynamically `sigma_sw = -phi_m * Pi`. In
the naive era, `n_l` was absorbed into `K`, hiding the explicit `n_l` factor.
Now `n_l` is explicit and the sign convention is made clear.

```diff
-            delta_sigma_sw.noalias() += n_S * (Pi_curr - Pi_prev) * identity2;
+            // Thermodynamic form: sigma_sw = -phi_m * Pi = -n_S * n_l * Pi
+            // Sign: n_l*Pi increases during hydration, so
+            // n_l_prev*Pi_prev - n_l*Pi_curr < 0 -> compressive increment.
+            delta_sigma_sw.noalias() +=
+                n_S * (n_l_prev * Pi_prev - n_l * Pi_curr) * identity2;
```

Both the `if/else` branches (vdW+aug path and legacy slope path) receive the
same change. After this step, K must be recalibrated — see PRJ parameter notes.

---

### Step 3 — Pi-path uses bulk `rho_LR` (Gibbs–Duhem consistency)
**Commit:** `c4888b6db4`  
**File:** `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`

**What changes and why:**  
`Pi = rho_LR * K * exp(-xi)`. The Gibbs–Duhem relation ties Pi to the bulk
liquid chemical potential (through `rho_LR`), not the micro liquid density
`rho_lR`. Switching the EOS density also removes the first-step spike: at `t=0`,
`rho_lR_prev ≈ 1e-6` while `rho_LR ≈ 1000`, so using `rho_lR` caused a `1e6×`
tensile `sigma_sw` spike.

Add `rho_LR` parameter to `computeReferenceMicroPorositySwellingStressIncrement`
and propagate it through `updateSwellingState` and both call sites:

```diff
 computeReferenceMicroPorositySwellingStressIncrement(
     double const n_l_prev, double const n_l,
-    double const n_S, double const rho_lR, double const rho_lR_prev,
+    double const n_S, double const rho_lR, double const rho_lR_prev,
+    double const rho_LR,
     ...
```

Inside the function replace micro density with bulk density for Pi:

```diff
-            double const rho_curr = rho_lR;
-            double const rho_prev = ... rho_lR_prev ...;
-            double const Pi_curr = rho_curr * K * std::exp(-xi_curr);
-            double const Pi_prev = rho_prev * K * std::exp(-xi_prev);
+            double const Pi_curr = rho_LR * K * std::exp(-xi_curr);
+            double const Pi_prev = rho_LR * K * std::exp(-xi_prev);
```

Also add `(void)rho_lR; (void)rho_lR_prev;` since those args are now unused in
the Pi branch (kept in signature for the legacy slope path).

---

### Step 4 — Remove `use_micro_liquid_density_for_pi` dead flag
**Commit:** `ce9178fa96`  
**File:** `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`

The flag `use_micro_liquid_density_for_pi` was removed after Step 3 made
`rho_LR` the fixed Pi density. Update the comment at the first-step spike
explanation to remove the flag reference:

```diff
-                // while rho_lR = 1000.  With use_micro_liquid_density_for_pi=true,
+                // while rho_lR = 1000.  When micro density enters the Pi-path
+                // swelling stress (Pi = rho_lR * mu_lR), the density mismatch
```

---

### Step 5 — vdW dimensional fix and literature Hamaker constant
**Commit:** `0d579e8aeb`  
**File:** `ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h`

**What changes and why:**  
The vdW base potential `mu_lR_vdW` must have units J/kg. The original formula
was missing a `/rho_lR` factor, giving Pa instead of J/kg. Fix:

```diff
-            (rho_SR * rho_SR * rho_SR) / (n_l * n_l * n_l);
+            (rho_SR * rho_SR * rho_SR) / (n_l * n_l * n_l * rho_lR);
```

The derivative `dmu_lR_drho_lR` is now non-zero (was erroneously 0.0):

```diff
-    out.dmu_lR_drho_lR = 0.0;
+    out.dmu_lR_drho_lR = -out.mu_lR / rho_lR;
```

**PRJ parameter change:** All PRJ files update `hamaker_constant`:
```
5.1e-21 J  →  2.2e-20 J
```
`2.2e-20 J` is the Israelachvili & Adams (1978) SFA mica-water-mica value,
standard smectite proxy. **Do NOT calibrate `hamaker_constant`** — it is a
material constant. Only `potential_augmentation_prefactor` (K) is calibrated.

**K recalibration after this step** (Villar bisection):
```
dd1400: K = 7654.9 J/kg
dd1600: K = 29984.9 J/kg
dd1800: K = 118582.6 J/kg
```

---

### Step 6 — Viscosity guards and micro-pressure density default
**Commit:** `66b782afa1`  
**File:** `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`

Add `requirePositiveViscosity` helper and call it at every local-solve entry
point to give a clear error instead of a silent NaN:

```cpp
inline void requirePositiveViscosity(char const* caller, double const mu)
{
    if (!(std::isfinite(mu) && mu > 0.0))
        OGS_FATAL("{} requires finite mu > 0, got {:g}.", caller, mu);
}
```

Called at: `computePotentialExchangeUpdate`, `solveReferenceMassStoragePredictorState`,
`solveReferenceMassStorageCoupledState`, `solveImplicitMicroWaterContent`,
`computeImplicitNlDpL`, and both local Jacobian assembly sites.

---

### Step 7 — Dead-code removal
**Commit:** `4d47efff55`  
**File:** `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`

Remove the 3-argument overload of `computeCompatibilityMicroHydraulicOutput`
(used bulk `rho_LR` in the vdW denominator — ~10% error — and was dead code
since the 4-argument version with `local_context` was the only active path).
Replace with a documentation comment explaining which overload is active and why.

```diff
-inline CompatibilityMicroHydraulicOutputData
-computeCompatibilityMicroHydraulicOutput(
-    double const n_l, double const rho_LR, ...) { ... }
+// NOTE: the only active overload is the 4-argument version below.
+// The 3-argument overload was removed 2026-05-22 (dead code, ~10% density error).
```

---

### Step 8 — Unit test refactor
**Commit:** `3ac6b7de1f`  
**File:** `Tests/ProcessLib/RichardsMechanics/DSMMicroMacroSingleIntegrationPoint.cpp`

The unit tests must be updated to use the hierarchical `phi_m_from_n_l` lambda
for the REV-scale reference density, and the coupled-state solve must pass `phi`
explicitly. Key changes:

1. Add `phi_m_from_n_l` lambda inside test functions:
```cpp
auto const phi_m_from_n_l = [&](double const n_l_eval)
{
    double const one_minus_n_l = std::max(1e-12, 1.0 - n_l_eval);
    return (1.0 - phi_safe) / one_minus_n_l * n_l_eval;
};
```

2. Replace `n_l_prev * rho_lR` with `phi_m_from_n_l(n_l_prev) * rho_lR` for
   `rho_l_prev`.

3. The reference scalar solve must pass the `local_context` struct with `phi`:
```cpp
return solveImplicitMicroWaterContent(
    n_l_prev, dt, rho_LR, alpha_bar, mu, macro_potential_ref,
    {.phi = phi,
     .volumetric_strain = volumetric_strain,
     .volumetric_strain_prev = volumetric_strain_prev},
    potential_exchange_params).n_l;
```

4. Bound direction fix — `EXPECT_LE(scalar_update.n_l, phi + ...)` replaces
   `EXPECT_GT(scalar_update.n_l, phi)` (hierarchical split makes `n_l` bounded
   differently from the naive case).

---

## PRJ parameter summary

Every `ANCHORS_MS33_Model*/` PRJ file must have:

```xml
<hamaker_constant>2.2e-20</hamaker_constant>  <!-- Israelachvili & Adams 1978, DO NOT calibrate -->
```

K values per dry density (calibrated to Villar targets, Model I):
```
dd1400: <potential_augmentation_prefactor>7654.9</potential_augmentation_prefactor>
dd1600: <potential_augmentation_prefactor>29984.9</potential_augmentation_prefactor>
dd1800: <potential_augmentation_prefactor>118582.6</potential_augmentation_prefactor>
```

MCC pre-consolidation cap `pc` (MCC models only: `*_mcc_native`, `*_mcc`) — parameter
`pc_char_mcc` (MCC_NATIVE) / `InitialPreConsolidationPressure` (MCC_SUITE). It sets both
the initial yield-cap size and the residual-normalization `pc_char`. Values are the
CIMNE-UPC saturated pre-consolidation p0* (EURAD-2 MS report, Table tab:CIMNE_BBM;
calibrated to Eriksson 2017 compaction pressures):
```
dd1400: 6.0e6    (6 MPa)
dd1600: 12.0e6   (12 MPa)  <- also ModelIII / ModelIV(bentonite) / ModelVII (dd1600 reference)
dd1800: 24.0e6   (24 MPa)
```
**K consistency (2026-05-29):** MCC_NATIVE K was aligned to the LE submission's agreed
EMDD≡rho_d calibration (`potential_augmentation_prefactor` = 26950/71900/214400 for dd1400/1600/1800;
III/IV/VII use the dd1600 value 71900). It previously carried the stale EMDD=0.8*rho_d K
(5500/13050/31280 -> 1.12/2.61/6.09 MPa). With the consistent EMDD≡rho_d K the developed
swelling stress vs the CIMNE caps gives:
- dd1400: p=4.91 MPa < 6 -> elastic (coincides with LE = Dixon).
- dd1600: p=15.2 MPa > 12 -> **yields and hardens** (pc 12->15.8); LE elastic 14.1 (+7.8%).
- dd1800: ~40 MPa >> 24 -> **fails** (plastic excursion +67% exceeds the MFront-MCC integrator;
  not fixable by dt reduction).
- III gap / IV pellets: elastic (gap/pellet stress relief keeps p ~7 / ~4.5 MPa < 12 cap).
- VII free: fails on the MCC supercritical dry side at the wetting front -> use the LE variant.

ModelIV's pellet zone (rho_d=900) shares the single `pc_char_mcc`; CIMNE has no p0* at 0.9 g/cm3
(open: needs a separate pellet-zone pc parameter).

**Per-medium DSM swelling override (already in the code, no change needed):** the
`<potential_exchange>` block accepts per-material `<medium id="N">` sub-blocks that override the
base DSM params (parsed in `CreateRichardsMechanicsProcess.cpp` via `parsePotentialExchangeParameters`
with the base as defaults; the local assembler picks per material id in `LocalAssemblerInterface.h`).
ModelIV now uses this to differentiate the pellet zone (material id 1): `n_s=0.324` (=1-phi0_pellet),
`initial_micro_water_content=6.59e-4`, lower `potential_augmentation_prefactor`. Result: heterogeneous
swelling (clay ~2.45 vs pellet ~0.1 MPa) + partial density homogenisation. NOTE: ModelIV currently
runs at a **demo** base K=13050 (EMDD=0.8, runnable); the EMDD≡rho_d K=71900 makes the multi-element
gap/pellet models (III, IV) impractically slow (III crawled to step 9371). Pellet K_pellet=1200 is a
demo value, to be calibrated to Dixon EMDD≡rho_d sigma_sat(0.9)~0.35 MPa for production.

**Model III (gap) stays elastic — by design:** the gap is a soft LE spring (E_gap=1e7 Pa, phi0=0.985)
so the bentonite equilibrates at ~7 MPa < 12 cap (no yielding, no evolving yield surface, no plastic
gap closure). A stiffer/contact-like gap would build stress to the cap but suppress closure; true
contact mechanics is absent in OGS RM. Documented limitation, not a calibration error.

All PRJ files must also have:
```xml
<use_micro_liquid_density_for_micro_pressure>true</use_micro_liquid_density_for_micro_pressure>
```

Bishop effective-stress coefficient (group decision 2026-05-29): **all DSM models** use
```xml
<property><name>bishops_effective_stress</name><type>BishopsSaturationCutoff</type><cutoff_value>1</cutoff_value></property>
```
i.e. chi=0 while S_L<1, chi=1 only at full saturation. Molecular suction is carried solely
by the Pi-path disjoining-pressure swelling source (no double counting). Rationale: with
`BishopsPowerLaw` (chi=S_L) the high initial suction (~100 MPa, S_L~0.33) gives chi*s~33 MPa
effective stress that busts the MCC cap at t=0 (instant yield -> divergence). With the cutoff
this vanishes and III/IV/VII integrate (previously failed at step 1/82).

Initial stress: confined models use `sigma0 = -1.5e5` Pa (seating). Free-boundary models
(ModelVII free swelling) **also** use -1.5e5 now — the former -33 MPa balanced the PowerLaw
Bishop force (sigma0_eff = alpha*chi*s); with chi=0 that balance gives 0, so -33 MPa is obsolete.

ModelVII (free swelling) + MCC: documented limitation. Even with cutoff+sigma0 fixed it diverges
near the wetting front (~step 74, t~19 d): unconfined + sharp front -> deviatoric stress on the
MCC supercritical/dry side (p<<pc/2) -> dilatant softening -> local return-map instability.
Decision: use the **LE** variant (`ANCHORS_MS33_LE_PER_DD/ModelVII`) for the free-swelling case.

Optional macro-retention `SaturationTuller` (implemented, not yet wired into MS33 prjs which use
`SaturationVanGenuchten`, p_b=27 MPa) supports `cavitation_pressure` (Frydman-Baker cutoff: freezes
S_L beyond p_cav, zero macro storage). Tuller/Or-faithful value: `cavitation_pressure=1.4e8` Pa
(~140 MPa, homogeneous-nucleation tensile strength of water; Or & Tuller 2002 WRR 38(5) 1061, via
Zheng 1991 / Speedy 1982) — NOT the tens-of-MPa Lu figure.

---

## Maintenance checklist

**When any of the following change, update this file before committing:**

| Trigger | Section to update |
|---------|-------------------|
| Any hunk in `RichardsMechanicsFEM-impl.h` changes | Relevant step diff above |
| Any hunk in `PotentialExchange.h` changes | Step 5 diff |
| DSMMicroMacro unit test changes | Step 8 description + passing criteria |
| PRJ `hamaker_constant`, K, or MCC `pc` (`pc_char_mcc`/`InitialPreConsolidationPressure`) values change | PRJ parameter summary table |
| New benchmark model added to canonical LE set | Verification commands + passing criteria |
| Build system flags change | Build instructions |
| New step added beyond Step 8 | New section following the same template |

**Rule:** This file must always be sufficient for an agent starting from clean
`master` to reproduce `dsm_native_hierarchical` without access to the branch
history. If a diff here is stale, the file has failed its purpose.
