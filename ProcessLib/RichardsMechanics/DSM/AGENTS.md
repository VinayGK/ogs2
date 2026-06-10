# AGENTS.md — ANCHORS MS33 DSM (streamlined)

## Scope
- Models: I / III / IV / VII (plus V_LE/V_MCC where present).
- Core implementation: `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`.
- Branch: `dsm_native_hierarchical`.
- Binary: `/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs`.

## Roadmap (one-line commit refs)
- Step 1 REV-scale storage + split consistency: `0d7a9edd64`.
- Step 2 thermodynamic swelling stress + K recalibration: `88d42c98fd`.
- Step 3 Pi-path Gibbs–Duhem consistency + flag cleanup: `c4888b6db4`, `ce9178fa96`.
- Step 5 vdW dimensional fix (`/rho_lR`) + literature A lock: `0d579e8aeb`.
- Step 6 DSM hardening (viscosity guards, micro-pressure density default true): `66b782afa1`.
- Step 7 dead-code removal (compatibility overload/unused flag): `4d47efff55`, `ce9178fa96`.
- Step 8 DSM micro-macro test refactor (13/13 passing): `3ac6b7de1f`.

## Key physics/implementation invariants
- Porosity split: `phi = phi_M + phi_m`, with micro state carried by `n_l`.
- Storage is REV-scale: `phi_m * rho_lR`.
- Swelling stress uses thermodynamic Pi-path tied to `rho_LR` for Gibbs–Duhem consistency.
- vdW potential terms are additive; never replace additive update with assignment.
- `hamaker_constant` is literature/material-fixed; calibration target is K (not A).

## Execution instructions
- Keep committed runs benchmark-spec compliant.
- After physics changes, require:
  1. Model-I Villar target check within tolerance,
  2. canonical LE reruns with zero rejected steps,
  3. `./bin/testrunner --gtest_filter='*DSMMicroMacro*'` passing.

## Current summary
- Production path stable under latest DSM fixes.
- Canonical LE outcomes unchanged in accepted/rejected-step sense.
- Open benchmark-quality work is primarily calibration/interpretation side (not immediate solver-break state).

## DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md maintenance rule

`ProcessLib/RichardsMechanics/DSM/DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md` is the reconstruction
recipe for this branch from a fresh `master`. It must stay current.

**Update DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md before committing whenever:**
- Any hunk in `RichardsMechanicsFEM-impl.h` or `PotentialExchange.h` changes.
- The DSMMicroMacro unit tests change (step 8 section + passing count).
- Any PRJ `hamaker_constant`, `potential_augmentation_prefactor` (K), or pre-consolidation
  pressure (`pc_char_mcc` / `InitialPreConsolidationPressure`, the MCC cap `pc`) value changes.
- A new benchmark model is added to the canonical LE set.
- Build flags or the verification `ctest` invocation changes.
- A new step beyond Step 8 is added (add a new numbered section).

Do not mark a step done in AGENTS.md unless DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md already reflects it.

## Known limitations (logged 2026-05-27)

### Forgotten Maxwell pair — mean stress absent from μ_lR  (NEXT IMPLEMENTATION GOAL)

**Location:** `PotentialExchange.h` — `μ_lR = μ_lR(n_l, ρ_lR)` (aggregate fraction
`1−n_l`, no stress/strain); swelling enters σ one-way as the Π-eigenstress
`σ_sw = −φ_m·Π(n_l)` (`RichardsMechanicsFEM-impl.h`, commit `72f4f3a192`).

**Issue:** from one `Ψ(ε,n_l)`, `∂σ/∂n_l = ∂μ_lR/∂ε` (Maxwell). The Π-eigenstress
gives `∂σ/∂n_l ≠ 0` (mean-stress, isotropic); `μ_lR` has no stress dependence so
`∂μ_lR/∂ε = 0` → **broken pair** → not derivable from a single `Ψ` →
non-conservative (`∮≠0`) past the gate. (= the "unlicensed equipresence deletion"
of the T&N lecture's Maxwell-pairs unit.)

**Gate:** exact **iff `σ_n < Π`** — below the disjoining pressure the missing
mechanical-expulsion term is physically zero; swelling-driven monotonic loading is
in-domain. Breaks past `Π` (over-compaction): `Π(n_l)=σ_n` is never closed (the
exchange does only the chemical balance `ψ_M=ψ_m`) → `φ_M→0` crash. Current claim:
restricted-domain admissibility (A).

**Right fix = the next goal:** give `μ_lR` a mean-effective-stress dependence — the
Maxwell partner of the swelling eigenstress, derived from one `Ψ` (coefficient
fixed by the existing swelling closure, **no new constant**; re-verify `K`). Full
spec — diagnosis, code plan, tiers (B1 sharp / B2 smear `Π`), verification — in
**`MAXWELL_PAIR_RESTORATION.md`**.

### Hydraulic-side double-counting of suction in Darcy flux

**Location:** `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h:2627` (and
the parallel branches at 3368 and 3991) — assembly of the macro Darcy flux:

```cpp
laplace_p += dNdx_p^T · (k_intr·k_rel·ρ_LR/μ) · dNdx_p · w
// drives:  q_L = -(k_intr·k_rel/μ) · ∇p_L
```

**Issue:** `p_L` is the primary `pressure` process variable, which the boundary
conditions impose as *total* suction (lab-measured, capillary + molecular).
The full ∇p_L therefore drives the macro Darcy flux, including the molecular
(disjoining-pressure Π) component. Physically the molecular component should
drive the micro→macro mass exchange via the DSM source term ρ̇_micro→macro,
NOT bulk advection in the macro pore.

**Manifestation:** boundary suction ramps of order 100 MPa create huge ∇p_L
gradients that propagate the wetting front much faster than the material
physically would. Affects all transient III/IV/VII results at finite t. Does
NOT affect Model I (no spatial gradient) and does NOT affect t→∞ equilibria
of III/IV/VII (saturated swelling pressure, asymptotic gap closure, equilibrium
void ratio).

**Mechanical-side companion (already fixed):** the same total p_L was being
fed into Bishop's effective stress via the `BishopsPowerLaw(exponent=1)` path
(χ = S_L → χ·p_L = molecular component leaks into σ_eff during unsaturated
phase). Fixed by switching all 4 MS33 PRJs to `BishopsSaturationCutoff(cutoff=1)`
so χ=0 below S_L=1; DSM Π then carries all swelling-source work in the
unsaturated regime.

**Right fix (open):** split p_L into p_L_macro (capillary, bounded by
Young–Laplace ~3 MPa for compacted bentonite) and p_L_micro (disjoining, Π).
Use ∇p_L_macro only in Darcy; route the (p_L − p_L_macro) residual through
the existing DSM micro↔macro exchange. Requires a new constitutive law for
the split and a process-variable refactor in `RichardsMechanicsProcess`.
Estimate: 1–2 weeks careful work + tests.

### MCC tension-apex non-convergence (ModelVII AND ModelIV) — use LE

**Location:** `MaterialLib/SolidModels/MFront/ModCamClay_semiExpl_constE.mfront`
(yield `f = q² + M²·p·(p−pc)`, `p = −trace(σ)/3 + pamb`).

**Issue:** in the free-swelling ModelVII, differential swelling at the wetting
front drives integration points to the tensile APEX of the cam-clay ellipse —
the failing Gauss point sits at `p ≈ 0` with a residual deviatoric `q ≈ 0.16 MPa`
and `pc ≈ 12 MPa` (healthy). At p=0 the ellipse pinches to the single point
(0,0); a state with q≠0 has no admissible return → MFront `status -1`.

**ModelIV joins this class (2026-05-29).** Once the clay–pellet ModelIV uses the
*physical* soft pellet modulus (`E_mcc_pellet = 9.2549 MPa = C·ρ_d³` at ρ_d=900,
parity with the LE variant), the compliant pellet lets differential swelling drive
the pellet/clay interface to the apex (`p → 0.03 MPa`, `q → 4.8 MPa`, `eqpl=0`,
`status -1` at ~22 d). The earlier MCC ModelIV that "completed" (1977 steps) used an
unphysically stiff pellet (`E = 52 MPa = E_clay`) that suppressed the interface
deformation. Per the user decision (physical params preferred), **ModelIV is now
submitted with the LE variant** (`ANCHORS_MS33_LE_PER_DD/ModelIV`, soft pellet,
converges: clay 13.9 / pellet ~0 MPa). ModelIII does NOT fail because its soft gap
is a per-material LinearElastic zone (id=1) — no MCC integration there, no apex.

**Tested (2026-05-29) and rejected as the fix:** added a gated `TensionCutoff`
(`pt_cutoff`) @Parameter (default −1 = OFF; converging models I/III/IV stay
byte-identical, verified dd1400 elastic 4.908 MPa / eqpl=0). `pt_cutoff=0`
(strict no-tension) routes p<0 elastic predictions to a volumetric-only return
that caps the mean stress at 0. It works (nodal min mean −63 → +52 kPa) but VII
still fails at the same step: the cutoff caps p but leaves q, so `f = q² > 0`
persists at the apex. Completing the Jacobian did not change the outcome →
structural, not a numerical bug. The earlier "negative mean stress" reading was
a nodal-extrapolation artifact; the true obstruction is the p→0 / q≠0 apex
coincidence. `pc_min` and `pamb` were separately ruled out (see
`project_dsm_mcc_bishop_cutoff.md` memory).

**Verdict:** the gated cutoff is left disabled in the .mfront as the recorded
experiment, NOT a production lever. ModelVII free swelling and ModelIV clay–pellet
both use the LE variant (`ANCHORS_MS33_LE_PER_DD/ModelVII` and `.../ModelIV`).
A genuine cure would relax BOTH p→0 and q→0 (full apex/fissuring collapse), which
zeroes wetting-front stiffness and is not globally solvable in OGS RM.

**Pragmatic interim:** submit current results with explicit caveat in the
deliverable (deck frame "Known Limitation — Hydraulic Side of OGS RM-DSM").

**2026-06-01 -- full-Pi closure; beta_sw retired; EMDD=rho_d calibration.**
The disjoining-pressure eigenstress (Delta sigma_sw = n_S(n_l_prev Pi_prev - n_l Pi) I,
Pi = -rho_lR psi_Micro) is the implemented micro-swelling closure; the legacy
micro_water_content_swelling_slope (beta_sw) branch and its
<accumulate_swelling_contributions> PRJ tag are removed. K re-fit (dt-converged P*
basis) to the Dixon (2023) MX-80 anchor under the EMDD=rho_d ANCHORS-groups
agreement: targets 4.92/14.16/40.86 MPa, K=35625.4/85312.6/224610 J/kg at
dd1400/1600/1800. CLAUDE.md 12.1 updated; 12.2 provenance synced across the LE
suite; 12.2 blocks added to the ModelVII experimental BC variants. Verified:
ogs+testrunner build clean; 14/14 DSM unit tests pass (incl. corrected
active_nS=1-n_l, section-2 incident); full MS33 suite (10/11 to t_end; dd1600
documented corner crash; endpoints ~2% high, first-order dt error).

**2026-06-08 -- K(rho_d): augmentation prefactor as a function of dry density.**
Branch `dsm_native_pdisj_maxwell_kofdd` (off `dsm_native_pdisj_maxwell`); build
`/Users/vinaykumar/git/build/kofdd_20260608/bin/ogs`. Implemented Option A
(parse-time table): `potential_augmentation_prefactor` may be set by a
`<potential_augmentation_prefactor_vs_dry_density>` block (child lists
`<dry_densities>`/`<prefactors>`) evaluated at each material's `<dry_density>`
(rho_d). K is resolved to a scalar at parse time -> *initial/target* rho_d,
constant in time, so NO Jacobian/tangent change downstream (Vinay's call
2026-06-08). The shared table inherits into per-`<medium id>` overrides via the
existing defaults mechanism; the scalar key and the table are mutually exclusive
(OGS_FATAL if both); a table requires a `<dry_density>`; `getValue` clamps
outside the node range. Files: `PotentialExchangeParameters.h` (+2 fields,
forward-decl MathLib::PiecewiseLinearInterpolation), `CreateRichardsMechanics
Process.cpp` (parse + resolve). DONE.
- Test (pellet block, Model IV): curve-vs-scalar equivalence. New PRJs
  `ms33_modelIV_pellets_kref100x.prj` (scalar K) and `..._kofdd.prj` (table K),
  both at k0 x100 spec (test acceleration, Vinay 2026-06-08; rate-only, endpoint
  unchanged), differing ONLY in how K is specified. Table nodes {(900, 20600),
  (1600, 103879.0)} J/kg carried verbatim from `ms33_modelIV_pellets.prj` so the
  table reproduces each material's existing per-material K. VERIFIED 2026-06-08:
  both run to t=200 d (ts_689); vtkdiff abs-max diff = 0 on all 14 output fields
  -> bit-for-bit identical. Registered both run-only in `Tests.cmake`.
- Back-compat VERIFIED: unmodified `ms33_modelI_dd1600.prj` runs clean through the
  refactored parser (scalar path).
- OPEN (deferred, Vinay): (i) closed-form vs table interpolation shape for
  intermediate rho_d (linear/log) is a modelling choice, not yet decided;
  (ii) the *current/evolving* rho_d variant (K riding porosity, needs a tangent
  term, double-count risk vs the exp(-xi) porosity dependence) is NOT
  implemented -- this delivery is initial/target rho_d only.

---

## Full consistent tangent for the Maxwell local Jacobian (2026-06-09, branch dsm_maxwell_jac_parallel)

Tangent-only gap-closing (Vinay's AceGen derivation, `THM_DSM_Richards_maxwell_web.wl`); residual UNCHANGED.

- (a) Analytic micro 2×2 Jacobian in `solveReferenceMassStorageCoupledState`
  (RichardsMechanicsFEM-impl.h): replaced the 4 FD `evaluate()` calls with a
  closed-form `J = d(mass_res, dens_res)/d(n_l, ρ_lR)` (`evaluate_analytic_jacobian`
  lambda). Reuses the helper-chain μ_lR derivatives; recovers the live-nS chain by
  re-running the vdW helper with the right `dnS_dnl`. Branch on
  `fd_jacobian_for_exchange` (default analytic, FD path kept as fallback). J22=1 exact.
- (b) Enabled the u-side swelling Jacobian: `enable_dsm_swelling_up_jacobian`
  false→true (~L4394). The film-ON K[u,p]/K[u,u] block already matched the Maxwell
  identity (`d σ_sw/d ε_v = +(1-φ_M)·n_l·b·K_drained`; `d σ_sw/d n_l = -(1-φ_M)·(p_film+n_l·Π')`),
  no formula fix needed.
- VERIFIED 2026-06-09 (Model I dd1400, Maxwell film path; maxjac vs pre-edit floor binary):
  (i) solution-unchanged max rel diff = 6.748e-15 (PASS); (ii) local 2×2 cascade
  quadratic, `|R_{k+1}|/|R_k|²` ~0.24 const (4.84e-2 → 5.57e-4 → 7.23e-8 → 1.3e-15),
  analytic = FD to round-off; (iii) global iters/step identical 1.338 (max 2) pre vs
  post (this benchmark's global problem is near-linear per step, `a=1e-16` EOS-bypass —
  global count cannot move; gain is in the local cascade). SPLICE B not needed.
  Full numbers + table in `AUDIT_maxwell_local_jacobian_2026-06-09.md`.
- OPEN: a two-way-coupled benchmark (EOS active, `a` not bypassed) to measure the
  predicted global iters/step reduction — not yet exercised.

---

## Phase A — park analytic micro + u-side OFF; restore FD micro; fix dd1800 FMA fragility (2026-06-09, branch dsm_maxwell_jac_parallel)

The 2026-06-09 "full consistent tangent" delivery above did NOT survive an
at-scale 6-model MS33 check. Root cause established:

1. **dd1800 broke from FMA fragility, not math.** The `if (use_fd_jacobian) {…}
   else {analytic}` boundary changed clang's FMA fusion choices under the build
   default `-ffp-contract=fast`; on the dd1800 near-singular tangent that tipped
   the global Newton path.
2. **Analytic micro 2×2 has a real J11/J12 error on the dense / EOS-active case.**
   The §VERIFIED dd1400 result above is **solution-invariant ONLY under the
   `a=1e-16` EOS-bypass** (`density_residual≡0` degenerates the 2×2). The audit's
   (i)/(ii) claim is RELABELED accordingly (see CORRECTION note in
   `AUDIT_maxwell_local_jacobian_2026-06-09.md`). Not solution-invariant on dd1800.
3. **u-side blocks singularize** dd1800 and ModelIII gap2mm (SparseLU failure).

Phase A (this delivery) — ship-safe non-regression, analytic + u-side RETAINED
but parked OFF by default:

- `solveReferenceMassStorageCoupledState` (RichardsMechanicsFEM-impl.h): added
  `constexpr bool use_analytic_micro_jacobian = false` (parked off); gate now takes
  the FD micro 2×2 path when `use_fd_jacobian_for_exchange || !use_analytic_micro_jacobian`
  → **FD micro = parent**. `evaluate_analytic_jacobian` retained, opt-in (Phase B).
  Decoupled from `use_fd_jacobian_for_exchange` (default false) so block #3 stays
  analytic exactly as parent.
- Localized FP-contraction guard around the function: file-scope
  `#pragma STDC FP_CONTRACT OFF` + body `#pragma clang fp contract(off)` (clang) —
  removes the dd1800 FD-reassociation fragility.
- `enable_dsm_swelling_up_jacobian` back to `false` (~L4488).
- `ParallelVectorMatrixAssembler.cpp` copy()-guard kept (math-neutral).

VERIFIED 2026-06-09 (maxjac_omp NEW vs mxconj_omp parent OTHER, OGS_ASM_THREADS=4,
fresh runs on identical inputs): all 6 MS33 (I dd1400/1600/1800, III gap2mm,
IV pellets_kref20x, VII freeswelling) **complete on both**; identical accepted-step
counts (308/311/308/438/636/507); final-VTU fields bit-identical to parent to
round-off (every field ≤ ~1e-12 rel-to-scale; mostly 1e-14–1e-16). **dd1800 now
completes** (308 steps, 0 rejects). Full table in the audit doc. Compare workspace
`~/ogs-models/maxjac_compare_2026-06-09/{*/phaseA_new,*/phaseA_other}`.

- OPEN (Phase B): correct the analytic micro 2×2 J11/J12 on the EOS-active/dense
  case; re-derive/condition the u-side blocks so they don't singularize stiff cases;
  then re-verify on a two-way-coupled (EOS-active) benchmark before flipping either
  constexpr ON.

---

## Phase B — diagnose analytic micro 2×2 on dd1800; NO Jacobian error found, root cause is global-solver fragility (2026-06-09, branch dsm_maxwell_jac_parallel)

Goal was to find and fix the "real J11/J12 error" on the dense / EOS-active dd1800
case that Phase A point 2 (above) predicted. **That prediction is NOT supported by
the measurements below and is relabeled accordingly (CLAUDE.md §5).**

Method (temporary diagnostic, since removed; tree clean): set
`use_analytic_micro_jacobian=true` and added an env-gated trace
(`DSM_MICRO_JAC_TRACE`) that, at every micro Newton iterate, computed BOTH the
analytic `evaluate_analytic_jacobian` J and an independent central-difference FD J
of `evaluate` (the exact numerical derivative of the unchanged residual = ground
truth), logging per-entry value + relative difference. Ran
`Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1800.prj`
(staged in /tmp), single-thread, with analytic ON and (separately) FD ON.

MEASURED (dd1800, maxjac_omp, OMP_NUM_THREADS=1, 2026-06-09):
- **Analytic J == FD J to round-off at every iterate.** Over all 58 676 traced
  analytic-run iterates: max rel diff J11 = 3.5e-9, J12 = 1.2e-8 (= the FD
  central-difference truncation floor at h=1e-8), J22 = 0 (exact). J21: analytic
  ≈6.1e-15 vs FD 0.0 — both numerically zero under `a=1e-16` (EOS inert ⇒
  `drho_lR_dnl≈0`); the rel=1.0 is a 0/0 artifact, not a Jacobian error.
- **No det sign flips, no singular dets** (analytic vs FD) over the whole trace.
- **Converged micro states identical FD-vs-analytic to ≤3.2e-12 in n_l** over the
  entire common range (1392 micro-solves, up to abort); ρ_lR diff 0.0. The
  analytic path reproduces the FD-path solution entry-for-entry.
- Many micro-solves (2448/2884 in the FD run; 956/1392 in the analytic run) run to
  `max_iterations=60` pinned at `n_l_ceiling` — but this is **tolerated/normal**
  (the FD parent run does it too and completes with **0 rejected steps**).
- **Why analytic-ON fails dd1800:** with analytic ON the *global* run diverges at
  time step #110 ("Newton: the linear solver failed in the compute() step"), step
  size driven to 0.1 s. The FD parent sails through step #110 (Δt≈5720 s, 308
  accepted / 0 rejected). The micro 2×2 J is verified-correct, so the only thing the
  flag changes is FP accumulation order (analytic skips the 4 `evaluate()` calls);
  at dd1800's near-singular *global* tangent that ~1e-12 perturbation tips the
  brittle global linear solve into non-factorizability. This is the SAME dd1800
  fragility Phase A point 1 documented for the FMA boundary — it is a global-solver
  conditioning issue, NOT a micro-tangent error.
- **Premise check:** every MS33 PRJ (all Models I/III/IV/V/VII) sets
  `micro_liquid_density_a=1e-16`, so the micro EOS is bypassed everywhere; J21≈0
  throughout the registered suite. dd1800 differs from dd1400 only in
  `micro_solid_volume_fraction_reference` (0.6475 vs 0.5036), i.e. it is denser
  (stiffer global problem), NOT "EOS-active." There is no EOS-active MS33 case in
  the suite where a missing micro chain term could surface.

CONCLUSION: there is **no missing chain term to add** — the analytic micro 2×2
J11/J12 already matches the FD ground truth to round-off, including the live-nS
chain. Per the task STOP condition ("discrepancy deeper than a missing chain term"),
NO Jacobian change was committed. Diagnostic reverted; `use_analytic_micro_jacobian`
left at its Phase A default (`false`, opt-in). maxjac_omp rebuilt clean; dd1800
re-verified complete (308 accepted, 0 rejected) on the reverted FD-default binary.

Relabeling (CLAUDE.md §5): Phase A point 2's "analytic micro 2×2 has a real J11/J12
error on the dense/EOS-active case" was a *plausible-but-unverified* consequence
claim; the Phase B trace measures analytic==FD to round-off and equal converged
states, so the dd1800 break is reattributed to global-solver fragility (Phase A
point 1's mechanism), not a micro-tangent error.

- OPEN (Phase B, remaining): the analytic micro path is correct but NOT robust on
  dd1800 because it perturbs FP order on a brittle global solve. Candidate
  directions (none implemented, none verified): (a) keep analytic OFF by default
  (current state) — the local cascade gain is real but the global solve is
  near-singular at dd1800 regardless; (b) condition the *global* tangent / time
  stepper at dd1800 so it is no longer 1e-12-fragile, then analytic-ON is safe;
  (c) exercise a genuine EOS-active (`a`≠1e-16) two-way-coupled benchmark — none
  exists in the registered MS33 suite — to measure any global iters/step gain the
  consistent micro tangent could buy. u-side singularization (Phase A point 3) is
  untouched by this diagnosis and remains OPEN.

---

## Phase C — dd1800 conditioning DIAGNOSED + fix found (2026-06-10) — DONE

Resolves Phase B OPEN direction (b) ("condition the global tangent so dd1800 is no
longer 1e-12-fragile, then analytic-ON is safe").

DIAGNOSIS (measured, env-gated SVD probe, since reverted): the single-element MS33
tangent is 12×12; the **pressure block is intrinsically near-singular on every
step** (4 pressure diagonals ~4e-17 vs displacement ~1e6; **cond ≈ 5.77e22**;
null-space = a pure pressure DOF). Root of the #110 break: Eigen `IterScaling`
(`<scaling>true>`) overflows the ~0 pressure row to **NaN**; the bare un-scaled
matrix factorizes fine. Analytic-ON only nudges a pressure off-diagonal across the
IterScaling overflow boundary; FD parent stays just under. So: **global
conditioning (the scaling step), not a tangent error.**

FIX (verified, no literal, no recompile): set the `<linear_solver>` to
**`<scaling>false</scaling>`** (keep SparseLU). With analytic-ON + scaling=false,
ALL 6 MS33 complete with parent-identical step counts (308/311/308/438/636/507) and
parent-identical fields (≤6e-12 rel-to-scale). dd1800 #110 now passes (Δt=5720,
0 rejects). iters/step byte-identical to parent (no global iteration gain — EOS
bypass; the analytic tangent's value is correctness + per-GP cost, not convergence).
scaling=false is a no-op on the current FD default (verified parent-identical on the
clean binary), so the PRJ change is safe but only meaningful once analytic-ON ships.

u-side blocks (`enable_dsm_swelling_up_jacobian`) under the fix: dd1400/1600/1800
become parent-identical, BUT Model III gap2mm still singularizes and **Model IV /
Model VII are NOT solution-invariant** (dry_density 12% / sigma 0.25% shifts) — the
u-side blocks remain OPEN/unsafe, separate from this conditioning fix.

Status: NO change committed pending owner decision (the fix is meaningful only
bundled with the analytic-ON enablement, a numerical-method call → present to Vinay
per §9). Tree clean; analytic flag + diagnostic reverted; maxjac_omp rebuilt clean
(3b64bf9e). Full numbers in DSM/AUDIT_maxwell_local_jacobian_2026-06-09.md Phase C.

DONE 2026-06-10: LANDED. Vinay chose "land it" (GUARDRAIL EXEMPTION §9/§12.3,
user-approved). `use_analytic_micro_jacobian` flipped false->true (analytic micro
2x2 = default); all 9 registered DSM ctests carrying `<potential_exchange>` set
`<scaling>false</scaling>` (Eigen SparseLU; per-PRJ inline block; solver-only, no
§12.2 material change). Rebuilt maxjac_omp. VERIFIED (measured): all 9 ctests
complete to identical final ts and parent-identical to round-off vs FD baseline
(mxconj_omp, scaling=true) — max rel diff <= 6e-12 (table in AUDIT Phase D). MS LE
standard (ModelI/III/IV/VII) passes. Run-only ctests (no reference VTU) => no
reference-VTU refresh, no §3/§12.5 flag. u-side blocks STILL parked OFF (unsafe;
mIII singularizes, mIV/mVII solution-shift — separate work). See AUDIT Phase D.
