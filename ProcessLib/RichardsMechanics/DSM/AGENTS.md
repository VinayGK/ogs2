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

## Strained-film disjoining law h(w_m, eps_v) (2026-06-09, branch dsm_native_h_of_eps)

Goal (Vinay): reversible de-swelling/expulsion under load — film thickness
varies with compressive stress so the potential reverses; the dissipative
residual stays a future flow rule. Design + decision record:
**STRAINED_FILM_IMPLEMENTATION.md** (this directory).

- DONE 2026-06-09: enums + params (film_strain_coupling off|kinematic|
  equilibrium, film_strain_kappa aggregate|unity), computeStrainedFilmState +
  invertDisjoiningPressure (PotentialExchange.h), fold-point rewiring
  (applyFilmPressureMicroPotential REPLACES the shipped integrable partner when
  ON — D3 provisional), eigenstress threading (eps_v sentinel args), PRJ
  parsing, unit tests (StrainedFilmPotential.cpp). Off = bit-for-bit.
- STRUCTURAL FINDING [D]: a pure geometric squeeze of any repulsive Pi(h) can
  never reverse the potential (Pi'(h)<0 ⇒ imbibition); the reversal lives in
  the Derjaguin load term +b*p_conf/rho_lR, made h-live here. Emergent gate
  b*K_drained > 3*kappa*Pi(h) (kinematic) / min()-branch at p_conf = Pi(w_m)
  (equilibrium) — no bolted-on Macaulay gate.
- HONESTY NOTE: implemented cut = operational Derjaguin form, NOT yet
  Maxwell-exact; exact one-Psi closed forms derived in the design doc §9a,
  AWAITING Vinay's review before coding. Do not cite the branch as
  "Maxwell-exact".
- TODO: build + unit tests + dd1400 off-mode regression (in progress
  2026-06-09); §9a exact forms; confined expulsion probe; K re-calibration
  [PRED: saturated swelling-pressure equilibrium shifts in both modes].

## 2026-06-11 — LIVE K(rho_d) variant (Vinay: "K(rho_d) try it")

- DONE (2026-06-11): live (evolving dry-density) K(rho_d) implemented per
  Vinay's order; see K_OF_RHO_D_LIVE.md. New PRJ bool
  `potential_augmentation_prefactor_live_dry_density` (default false =
  parse-time freeze, bit-for-bit); helper effectiveAugmentationPrefactor
  (PotentialExchangeParameters.h) evaluates K_table(rho_SR*(1-phi)) at all
  FEM sites with porosity in scope (context phi / new defaulted
  total_porosity arg on the swelling increment / assembly phi); scalar
  fallback where no phi exists. Endpoint-clamped (getValue endpoint hold).
- Verified: 31/31 RichardsMechanics unit tests (3 new
  RichardsMechanicsLiveKOfRhoD, structural knots); dd1400 off-mode
  regression sigma_zz = -4.9218 MPa = recorded baseline
  (runs/2026-06-10_0841_dsm_native_h_of_eps_successful).
- PROVISIONAL: linear knot interpolation (shape undecided); dK tangent
  OMITTED first cut [PRED: extra Newton iterations, not benchmarked]; no
  live-mode production run yet (behavior under live K = predicted only).

## Equipresent Pi(n_l, eps_v) + compressible-liquid carrier D2 (2026-06-11, branch dsm_native_Pi_fofnlev)

Goal (Vinay, 2026-06-11): make the load coupling energetically compliant —
"Pi does not carry p_conf; the p_conf is still a bolt-on to the micro, not an
energetically compliant bolt-on; equipresence says Pi(n_l, <mech. state>)."
Equipresent argument is eps_v (configuration), not p_conf (force). Two
deliverables, PRJ-selectable, default off (bit-for-bit):
(E) exact one-Psi film pair (closes STRAINED_FILM_IMPLEMENTATION.md §9a) and
(L) compressible-liquid carrier (D2 proper — over-pressure from Psi_liq with
confined K_liq instead of the bolted +b*p_conf/rho_lR).

Design + decision record: **PI_OF_NL_EV_IMPLEMENTATION.md** (this directory) —
complete implementation/test/docs/beamer plan, written for an implementing
agent; decision queue Q1–Q5 (Vinay) inside.

- DONE 2026-06-11: branch + worktree created off 7ff8861847; design doc
  written; memory file project_dsm_pi_fofnlev + MEMORY.md pointer; beamer
  maxwell_from_psi.tex Step 19/23 "in design" markers.
- DONE 2026-06-11 (was: blocked on Q1; resolved by Vinay's "implement that
  now"): (E) implemented — FilmEnergyRoute enum/param/predicate,
  computeStrainedFilmEnergyPair (x_over_kappa-stable closed forms), exact
  fold branch (bare at TRUE n_l + one-Psi partner, g-cutoff product rule,
  eigenstress site unchanged), film_energy_route parsing + mode-matrix
  OGS_FATAL, ExactFilmEnergyPair.cpp (8 tests: 6 active + 2 Q3/Q4 skips).
  VERIFIED: 36/36 unit tests; T-5 loop measured |∮|/scale 8.4e-9 (exact) vs
  0.93 (operational — §9a "small" prediction corrected by measurement);
  T-1 dd1400 off-mode bitwise-identical (12/12 VTUs, parent-head binary vs
  new, one input). Build ~/git/build/Pi_fofnlev_20260611. UNCOMMITTED.
- PARKED 2026-06-11 (Vinay: "park it, this is getting very intricate"):
  (L) computeMicroLiquidCompression + tests T-6/T-8. Self-contained decision
  brief (routes (a)/(b), K_liq candidate + magnitudes, Q3<->Q4 coupling,
  unpark conditions) is in PI_OF_NL_EV_IMPLEMENTATION.md §8 — read THAT
  before any (L) work; do not unpark without Vinay's Q3+Q4 answers.
- DONE 2026-06-11 (partial): STRAINED_FILM_IMPLEMENTATION.md §9a annotated
  (measured correction); beamer Step 21/23/7b updated to implemented+measured
  status. STILL TODO: "Step 24 first numbers" frame (gated on T-8 / MS33 VII
  runs); Doxygen tag doc for film_energy_route (joint TODO with the
  undocumented film_strain tags).
