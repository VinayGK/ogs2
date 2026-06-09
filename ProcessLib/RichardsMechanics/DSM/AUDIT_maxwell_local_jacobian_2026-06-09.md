# Audit — Maxwell DSM local Jacobian + parallelization (2026-06-09)

Branch: `dsm_maxwell_jac_parallel` (off `dsm_native_maxwell_conjugate` @ d98f5f8324).
Method: three independent read-only audits (Jacobian-block inventory; per-GP
micro-solve cost; assembly parallelism / thread-safety). Line numbers refer to
`ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h` unless noted.
Scope: the local element Jacobian on the DSM `potential_exchange` /
`film_pressure_coupling` (maxwell-conjugate) path. Read-only; no edits.

---

## Headline finding

**The element-level parallel Jacobian assembly already exists and RichardsMechanics
already uses it.** `RichardsMechanicsProcess` inherits `AssemblyMixin`, which drives
`Assembly::ParallelVectorMatrixAssembler`; the element loop is OpenMP work-shared
(`ParallelVectorMatrixAssembler.cpp:172`, `#pragma omp for nowait`), with
thread-local scratch + a per-thread Jacobian-assembler copy
(`cpp:356-360`) and a `#pragma omp critical` scatter into the global matrix
(`MatrixElementCache.h:234`). The Maxwell local assembler is already
per-element-reentrant.

**It just defaults to serial:** thread count = `OGS_ASM_THREADS`, default **1**
when unset (`BaseLib/OgsAsmThreads.cpp:7`). So the parallelism is present but
dormant. "Parallelize the local Jacobian" is therefore mostly **enable + verify +
two cheap hardening fixes**, not a new parallel loop. Details in §3.

---

## Part A — Jacobian block inventory (assembleWithJacobian, body 3697–4678)

Flag defaults (`PotentialExchangeParameters.h`): `use_fd_jacobian_for_exchange =
false` (analytic by default); `film_pressure_coupling = true` (ON by default).

| # | Block | Lines | Term | Analytic/FD | Gate + default | Notes |
|---|-------|-------|------|-------------|----------------|-------|
| 1 | K[p,u] | 4200–4206 | exchange↔ε_v integrable partner (strain-view) | analytic | `!film_pressure_coupling` | fires only when film OFF (mutually exclusive with #2) |
| 2 | K[p,u] | 4287–4293 | integrable Maxwell p–u tangent (film ON) | analytic | `film_pressure_coupling` (ON) | **live**, complete, smooth in ε_v |
| 3 | K[p,p] | 4619–4623 | direct macro exchange-source tangent | analytic by default; FD if flag | always assembled; value FD iff `use_fd_jacobian_for_exchange` (default false) | **live**; analytic `dmu_lR_vdw_dpL` (4354) omits the direct p_conf(p_L) Bishop channel by design |
| 4 | K[u,p] | 4497–4502 | swelling-eigenstress u–p (σ_sw/n_l·dn_l/dpL) | analytic | `enable_dsm_swelling_up_jacobian && film` — **constexpr false (4394)** | **never assembled**; comment: "tiny, does not improve convergence" on dd1400 |
| 5 | K[u,u] | 4517–4522 | dσ_sw/dε_v = +(1−φ_M)·n_l·b·K_d | analytic | same constexpr-OFF gate | **never assembled** |
| 6 | K[u,p] | 4594–4599 | legacy swelling u–p (film-OFF form) | analytic | `else if enable_dsm_swelling_up_jacobian` — constexpr false | **never assembled** |
| 7 | K[p,p] | 4637–4642 | lagged secant micro-pressure term | secant | `!use_vdw_micro_potential...` | **dead on the Maxwell path** (vdW path sets the flag true) |

FD path: in `computePotentialExchangeUpdate` (179–230), central-difference of the
exchange residual; feeds **only block #3**; triggered by `use_fd_jacobian_for_exchange`.

**Audit conclusion (completeness):** on the default path the **pressure row**
(K[p,p] #3, K[p,u] #2) is analytic and live. The **displacement-side swelling
coupling** (K[u,p] #4/#6, K[u,u] #5) is **compile-time disabled** — the consistent
u-side tangent is never assembled. Benign for linear-elastic / weakly-coupled
cases (the existing note calls it convergence-neutral on dd1400), but it is a
genuine gap and is exactly the block that was missing for route-R.

---

## Part B — per-GP cost (the real expense)

ONE fused integration-point loop at **3763**; residual + Jacobian in the same pass.
Cross-GP **data-parallel: yes** — all state indexed by `ip`
(`current_states_[ip]`, `prev_states_[ip]`, `material_states_[ip]`), micro solve
reads frozen `prev_states_[ip]`, writes only `current_states_[ip]`; no neighbour
coupling. Caveat: the loop `+=`-accumulates into shared element matrices
(`local_Jac`, `Kup`, `Kpu`, …) — independent compute, but the scatter needs
per-thread partials, not a naive shared write.

Per-GP heavy work (parallelization / optimization targets):
- **Micro local solves** (dominant): scalar n_l Newton (1280–1495, cap 25) and a
  **2×2 (n_l, ρ_lR) Newton with a finite-difference Jacobian** (928–1149, cap 60):
  5 `evaluate()` per iteration (1 + 4 FD), each running a **nested EOS Newton**
  (`computeReducedMicroLiquidDensity`, 495–602, cap 30). So one outer 2×2 iteration
  ≈ 5 EOS solves; dozens of vdW evals per GP funnel through here (call at 3512).
- **Elastic tangent reconstructions** `computeElasticTangentStiffness`: 2/GP (film
  OFF) or 3/GP (film ON) — 3354, 4099, 4171, 4273.
- `C_el.inverse()` once/GP (3650); `updateConstitutiveRelation` once/GP (3673).
- `createConstitutiveModels` **inside** the ip-loop (3768) — redundant per-GP
  reconstruction, hoistable to once/element (pure win).

Note (separate from parallelism): the 2×2 micro Jacobian is **finite-difference**,
not analytic. An analytic micro-Jacobian would cut per-GP cost (no 4 FD evals/iter)
and likely the iteration count — a bigger lever than threading for single-element cost.

---

## Part C — parallelism state & thread-safety

1. Driver: `RichardsMechanicsProcess.cpp:307` → `AssemblyMixin.h:435`
   (`pvma_.assembleWithJacobian`) → `ParallelVectorMatrixAssembler.cpp:345`
   (`#pragma omp parallel`), loop at `:172`.
2. Partition: per-element work-sharing (no coloring); thread-local `local_Jac_data`,
   a `jacobian_assembler_.copy()` per thread, thread-local element cache; deferred
   scatter under `#pragma omp critical`. Composes with PETSc/MPI domain decomposition.
3. Thread-safety: **already per-element-reentrant.** One assembler instance per
   element; all assembly writes are per-ip/per-element. Shared objects
   (`solid_material_`, media map) are read-only / `const`. The only real (benign)
   race: `MaterialLib/MPL/Property.h:300` `mutable bool property_used` written
   `=true` on the shared Property — idempotent, harmless, but ThreadSanitizer-flagged.
   No statics, no element-loop-shared scratch, no nested omp.
4. Recommendation: the element loop is the right granularity and already exists.
   **Per-GP OpenMP is NOT recommended** (would nest inside the element-level
   `omp parallel`; GP counts too small to amortize fork/join; collides with the
   element scheme).

---

## Recommended actions

1. **Enable + verify** element-level assembly threading for the DSM Maxwell path:
   run a Maxwell benchmark with `OGS_ASM_THREADS>1` and confirm bit-stability vs the
   serial baseline (and a speedup). This is the actual "parallelize" deliverable.
2. **Hardening (cheap, pure wins):**
   - make `Property::property_used` an `std::atomic<bool>` (or set-once) for a clean
     TSan run under threaded assembly;
   - hoist `createConstitutiveModels` out of the ip-loop (3768) — once/element.
3. **NOT** per-GP nested OpenMP (see §C.4).
4. **Optional, separate from parallelism (bigger single-element lever):** replace the
   2×2 FD micro-Jacobian (928–1149) with an analytic tangent. Owner decision.
5. **Audit gap to flag (not parallelism):** the u-side swelling Jacobian
   (#4/#5/#6) is compile-time OFF — the default Maxwell Jacobian is p-row-complete
   but u-side-incomplete. Owner decision whether to enable/complete it.

---

## Gap-closing implementation + verification (2026-06-09, branch dsm_maxwell_jac_parallel)

Implemented per Vinay's AceGen derivation (`THM_DSM_Richards_maxwell_web.wl`).
Tangent-only: NO residual was changed. Granted autonomy on this branch.

### What changed (files + lines)

- **`RichardsMechanicsFEM-impl.h` — analytic micro 2×2 Jacobian** (item a). In
  `solveReferenceMassStorageCoupledState` (the `ScalarReferenceMassStorage` /
  `scalar_micro_macro_mass_storage_mode` 2×2 (n_l, ρ_lR) Newton):
  - Added `evaluate_analytic_jacobian` lambda (~after the `evaluate` lambda)
    assembling `J = d(mass_residual, density_residual)/d(n_l, ρ_lR)` in closed
    form: `d ρ_l/d n_l = (1-φ)·ρ_lR/(1-n_l)²`, `d ρ_l/d ρ_lR = (1-φ)/(1-n_l)·n_l`;
    `J11 = d ρ_l/d n_l·(1-dt·ε̇_v) + dt·α_M·d μ_lR/d n_l`;
    `J12 = d ρ_l/d ρ_lR·(1-dt·ε̇_v) + dt·α_M·d μ_lR/d ρ_lR`;
    `J21 = -drho_lR_dnl` (from the EOS, carries the live-nS chain in
    CurrentPorositySplit); `J22 = 1` exact (ρ_lR_EOS independent of the outer
    2×2 unknown ρ_lR). The `μ_lR` derivatives reuse the shared helper chain
    (macro-floor cutoff + film-pressure coupling + augmentation); the live-nS
    chain `d μ_lR/d n_l` is recovered by re-running the vdW helper with the
    correct `dnS_dnl` (-1 split, 0 reference) then re-applying the film coupling,
    reproducing the TOTAL n_l-derivative the FD of `evaluate` sees.
  - Branched the Newton-loop Jacobian on `use_fd_jacobian_for_exchange` (PRJ key
    `fd_jacobian_for_exchange`): **analytic by default**, the prior 4-FD-eval path
    kept reachable as the opt-in fallback.

- **`RichardsMechanicsFEM-impl.h` — u-side swelling Jacobian** (item b). Flipped
  `constexpr bool enable_dsm_swelling_up_jacobian` (~4394) from `false` to `true`.
  The film-ON block (K[u,p] ~4497–4502, K[u,u] ~4517–4522) implements the Maxwell
  identity: `d σ_sw/d ε_v = +(1-φ_M)·n_l·b·K_drained` and
  `d σ_sw/d n_l = -(1-φ_M)·(p_film + n_l·Π')` routed through `dn_l/dpL` — the exact
  transpose of the eigenstress block (`.wl` L54, L135). Block formulas already
  matched the identity; no formula fix was needed.

SPLICE B (algebraic `R_φM = φ_M - (φ - n_S·n_l)` storage form) was **NOT** needed:
full quadratic convergence was achieved without the φ_M formulation change.

### Verification — Model I dd1400 (Maxwell film path), maxjac vs floor (pre-edit)

Binaries: pre-edit `build/floor_20260609/bin/ogs`, new `build/maxjac_20260609/bin/ogs`.

(i) **SOLUTION-UNCHANGED — PASS.** Final VTU (`ts_99`) field-by-field max
relative diff = **6.748e-15** (round-off). displacement/pressure/saturation/
transport_porosity bit-identical (0.0); sigma, swelling_stress, micro_pressure,
dry_density_solid ≤ 1.7e-15 rel. Tangent-only confirmed.

(ii) **QUADRATIC CONVERGENCE — YES.** Per-iteration local 2×2 residual norm of a
representative converged GP solve (measured via a temporary env-gated trace,
since removed):

| iter | \|R\| (analytic) | \|R_{k+1}\|/\|R_k\|² |
|-----:|-----------------:|---------------------:|
| 0 | 4.838e-02 | — |
| 1 | 5.566e-04 | 0.238 |
| 2 | 7.227e-08 | 0.233 |
| 3 | 1.30e-15  | 0.250 |

`\|R_{k+1}\|/\|R_k\|²` is constant (~0.24) → quadratic (converged digits double
each iteration). The analytic cascade equals the FD cascade to the 1e-15 floor
(FD: 0.238/0.233/0.234), so the analytic Jacobian = the FD Jacobian to machine
precision while removing 4 `evaluate()` calls per iteration.

(iii) **iters/step (measured).** Global Newton iters/step IDENTICAL pre vs post:
mean **1.338**, max 2 (308 steps, 412 global solves) for both analytic and FD.
This benchmark's GLOBAL problem is near-linear per step (the p-row tangent was
already complete and `micro_liquid_density_a=1e-16` decouples the EOS, CLAUDE.md
§2 EOS-bypass regime), so the global iteration count cannot move; the consistency
gain is visible only in the LOCAL 2×2 cascade above. Mean LOCAL converged-solve
length = 3.137 iters, identical analytic vs FD.

Note (pre-existing, NOT introduced): ~30 GP local solves hit the 60-iter cap and
fall back to the predictor (residual plateau at \|R\|=1.0908e-02, mass=-5.904).
The plateau value and count are **byte-identical between analytic and FD**, so the
analytic Jacobian neither caused nor changed them; under the degenerate `a=1e-16`
2×2 (`density_residual≡0`) this is a property of the existing backtracking/`det`
guard, and it does not affect the converged solution (verified by (i)).

Predicted (not yet verified): on a benchmark where the EOS is active (`a` not
bypassed) and the global problem is genuinely two-way coupled, enabling the u-side
blocks + analytic micro tangent is expected to reduce global iters/step; this was
not exercised here.
