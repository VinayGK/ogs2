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

## Form (a) vs Form (b) paired comparison (2026-06-12) — READ BEFORE TOUCHING THE FILM ROUTES

Definitions (the two FORMS of the micro potential; beamer maxwell_from_psi.tex Step 21b):
- FORM (a): mu(Pi(n_l), eps_v) — Pi frozen at n_l, strain a separate argument
  (the ch.1 maxwell partner; default/off mode; maxwell_conjugate lineage).
- FORM (b): mu(Pi(n_l, eps_v)) — strain enters THROUGH the film state
  (film_strain_coupling=kinematic + film_energy_route=exact on THIS branch).
- [D] (a) = (b) Taylor-truncated at eps_v->0. Mutually exclusive at runtime
  (create-time guard); coexistent in one build as limit test + baseline.

Paired runs (common pre-recalibration K base, two binaries:
mc_20260608 d98f5f8324 for (a), pi_fofnlev_20260611 @4c7a5d03b2 for (b)),
record: ogs/formAB_2026-06-12/FORMAB_RESULTS.md + eurad-anchors snapshot. MEASURED:
1. dd1600 control (eps_v~0): BITWISE identical across forms AND binaries —
   the truncation identity holds exactly in running code.
2. VII discriminator: e_end 1.4995 (a) vs 1.3530 (b) — self-relaxing drive
   removes 0.146 of the over-swell; exact within 3e-4 of operational here.
3. OPPOSITE PULL: Task-13 DD 1.4592 (a) vs 1.4839 (b) vs expt 1.4139 —
   (b) helps VII but hurts Task-13. No reversible form satisfies both =>
   Task-13 residual is the IRREVERSIBLE channel (see TASK13_MCC_BLOCKAGE.md
   + MCC_INTERNAL_SWELLING_IMPLEMENTATION.md design).
4. (b) also better-conditioned on the Task-13 wetting front (1256 steps/138 s
   vs 2673/713 s) and ends force-balanced (Pi +0.55 MPa ~ load vs -1.85).
Variant PRJs committed next to the base files (suffix _formB_piexact_2026-06-12).

## Live-K(rho_d) analytic tangent completion (2026-06-12)

- DONE 2026-06-12: analytic dK/dphi = -rho_SR*(table segment slope) tangent
  (Vinay-approved completion of the live-K first cut; Jacobian-only,
  residual untouched) wired into the live p-u augmentation Jacobian block.
  New `AugmentationPrefactorTable::getSegmentSlope` (exact clamped
  piecewise-linear slope; zero outside/at edges, left slope at interior
  knots), `effectiveAugmentationPrefactorPhiDerivative`, mu-level exact
  K-partials `dmu_lR_dK`/`ddmu_lR_dnl_dK` (mu_aug linear in K), and the
  `PorosityFromMassBalance` dphi/deps_v = (alpha-phi)/(1+w) chain. Details +
  measured verification: K_OF_RHO_D_LIVE.md "Analytic tangent completion".
- MEASURED: 41 RM unit tests (39 pass + 2 designed skips; 2 new FD-vs-
  analytic tangent tests); dd1400 off-mode bitwise-identical vs the
  pre-tangent h_of_eps_20260609 binary (12/12 VTUs); truncated 1a_robin_A_Kl
  live-K sanity converges with iterations equal to before (17/2/2).
- CURE TEST (task42 1b *_Kl step-1 singularity): NOT CURED — both 1b_A_Kl
  and 1b_B_Kl still die in step #1, but the Newton trajectory measurably
  changed (contraction to |dx|_uz=1.09 over 7 its before the it.8 blow-up,
  vs pre-tangent monotonic divergence). Second mechanism suspected
  (hypothesis, not verified); evidence in task42_case1_2026-06-12/
  out_1b_{A,B}_Kl/run.log + _diagnostics_1bKl/README_DIAG.md addendum.
  OPEN: Vinay's call on the next probe; no further patching done.

## Ultracode-review fixes (2026-06-14) — branch dsm_native_Pi_fofnlev_review_fixes_2026-06-14

Implements the fixes in DSM/ULTRACODE_REVIEW_2026-06-14.md, one commit per
fix. Branch off dsm_native_Pi_fofnlev tip 9795f252e1. Build:
~/git/build/pi_fofnlev_fixes_20260614. Baseline: 39 RM unit tests pass + 2
designed skips. Every Jacobian-only fix re-verified the off-mode dd1400
final-VTU SHA256 = 91f404a5...577 (the STEP-0 reference) BYTE-FOR-BYTE.

- DONE 2026-06-14 (H2, JAC-only, commit 7e2eee0190): exact-route fold passed
  the parse-time scalar K into computeStrainedFilmEnergyPair while the bare
  out.mu_lR used the live effectiveAugmentationPrefactor(phi); g_cut corrupted
  under live K. Now uses effectiveAugmentationPrefactor(params, local_context
  .phi) at L730 (matches :768/:1364/:2084). Scalar-mode bit-for-bit; off-mode
  bitwise verified.
- DONE 2026-06-14 (M1, JAC-only, commit 6708ef1d98): the p-u Maxwell exchange
  tangent always used the integrable-partner dmu_lR_mech_deps_v; now dispatches
  dmu_lR/deps_v on the route (Off -> integrable; operational-strained ->
  d(bare(w_eff))/deps_v + b*(dp_conf/deps_v)/rho with dp_conf/deps_v=-K_drained;
  exact -> g_cut*pair.dmu_mech_deps_v). Off-mode bitwise verified.
- DONE 2026-06-14 (M2+L2, JAC-only, commit cbe3e11ed6): wired the displacement-
  side live-K swelling-eigenstress tangent d(delta_sigma_sw)/dK*dK/dphi*
  (dphi/deps_v, dphi/dp) into K[u,u]/K[u,p] (the 1b compliant-top cure
  candidate). SCOPE DECISION (announced): wired ONLY the live-K chain, NOT the
  pre-existing swelling u-p/u-u term Vinay set OFF 2026-06-01 (enable_dsm_
  swelling_up_jacobian left at its default; independent term). Gated on
  film_pressure_coupling && dK/dphi != 0 && PorosityFromMassBalance. Off-mode
  bitwise verified.
- DONE 2026-06-14 (L1, JAC-only, commit 78a71ae6df): exact-fold dg_dnl mixed
  live-nS out.dmu_lR_dnl with frozen-nS pair.dmu_bare_dnl_pre under
  CurrentPorositySplit; recompute the bare-pre derivative with the caller's
  dnS_dnl. Reference mode bit-for-bit. Off-mode bitwise verified.
- DONE 2026-06-14 (N1, JAC-only, commit 6a03a58213): zero the live-K
  dphi/deps_v=(alpha-phi)/(1+w) chain when PorosityFromMassBalance clamps phi
  (detect via unclamped-vs-stored phi; bounds are private). Applied at the
  live-K p-u block and the new M2 block. Off-mode bitwise verified.
- SUPERSEDED 2026-06-14 (L3, DOC-only, commit f23f69c5b4): documented the
  live-K dn_l/dK local-solve strain channel in ScalarReferenceMassStorage mode
  as a DELIBERATE PARTIAL TANGENT (not wired). Judgment call (prompt-authorized):
  wiring it would risk the converged forward solve for a LOW-severity gap; the
  dominant cure is M2. Off-mode bitwise verified. NOTE: the "risks the converged
  forward solve" concern was the basis for not wiring; it is resolved below by
  keeping the wiring strictly Jacobian-only (the local solve / residual / the
  computeImplicitNlDpL return are all untouched). See the WIRED entry next.
- DONE 2026-06-14 (L3, JAC-only, WIRES the above): the L3 implicit-n_l(K) strain
  channel is now wired into the global displacement Jacobian, RESIDUAL-SAFE.
  - New sibling helper computeImplicitNlDK (RichardsMechanicsFEM-impl.h, right
    after computeImplicitNlDpL): returns dn_l/dK = -(dr/dK)/(dr/dn_l) for
    ScalarReferenceMassStorage (0 in every other mode and at dt<=0). It rebuilds
    dr_dn_l by the SAME REV-mass reduction as computeImplicitNlDpL and uses
    dr/dK = -dt*exchange.drho_l_hat_dmu_lR*micro_potential.dmu_lR_dK (K enters r
    only through the exchange/mu_lR; the mass term carries no K). Reads ONLY the
    already-converged (n_l, rho_lR, micro_potential, exchange) the caller
    threads in -- it never re-solves and never mutates forward state.
  - Wired at the M2 swelling-eigenstress site as the implicit half of M2's
    explicit-K chain: d(delta_sigma_sw)/dn_l * dn_l/dK * dK/dphi * dphi/d(.),
    folded into the SAME dsig_sw_deps_v_scalar / dsig_sw_dp_scalar that ride the
    existing C*C_el^-1*identity2 map into K[u,u]/K[u,p]. d(delta_sigma_sw)/dn_l
    = -n_S*(Pi_curr - n_l*rho_curr*dmu_lR/dw_eval*dw_eval/dn_l) (dw_eval/dn_l=1
    on the OFF branch, film_state.dw_eff_dnl on the strained branch), matching
    the residual increment's argument chain and its dnS_dnl=0.
  - SCOPE (residual-consistency): the implicit chain fires only where the
    telescoped eigenstress IS the assembled residual (OFF + film-ON-operational);
    on the EXACT route (H1, residual = one-Psi pair) it is gated OFF -> the
    cured 1b_B tangent is left bit-for-bit untouched. Off / frozen K / clamped
    table edge / non-mass-storage -> all factors 0 -> Jacobian bit-for-bit.
  - JACOBIAN-ONLY: the local forward n_l solve, the residual, and the
    computeImplicitNlDpL return value are all unchanged. The prior agent's note
    inside computeImplicitNlDpL is replaced with a WIRED note.
  - PREDICTED (not yet verified, §5; build is a later phase): this completes the
    live-K mass-storage displacement tangent and is the candidate cure for the
    1b_A form-(a) step-1 divergence the 1b cure verdict left open; the converged
    root is unaffected by construction. Off-mode bitwise: NOT re-run here (build
    deferred) -> the off/frozen/clamped/non-mass-storage gates make the change a
    no-op there by construction, to be confirmed at build.
- DONE 2026-06-14 (H1, RESIDUAL-CHANGING, Vinay-authorized, commit 6391e357a2):
  under film_energy_route=Exact, source the eigenstress increment from the
  one-Psi pair.sigma_sw_m (drained-line, telescoped curr-prev) instead of the
  operational Pi(w_eff)-b*p_conf. VALIDATED by the new assembled loop-closure
  probe (AssembledExactPairClosesOperationalSigmaDoesNot): |W|/scale = 8.41e-09
  with exact sigma (post-H1, CLOSES) vs 0.0675 with operational sigma (pre-H1,
  does NOT close). H1 is CORRECT (not reverted). PHYSICS TRADEOFF (predicted):
  drained-line p_conf vs actual GP p_conf; they agree on the drained line,
  differ off it (e.g. fully confined). Off-mode unaffected (film off -> branch
  not entered): bitwise verified.
- DONE 2026-06-14 (NEW TEST §8, commit 9732b46498): added
  RichardsMechanicsLiveKOfRhoD.AssembledDisplacementTangentExactKinematicLiveK
  (StrainedFilmPotential.cpp) — helper-level FD-vs-analytic identity on the
  exact-route mu_lR eps_v tangent (H2/M1) and the live-K eigenstress eps_v
  tangent (M2), anchor (d), scale-derived tolerances, no Vinay expected value.
  Scope-noted: no run-level assembleWithJacobian harness exists in the unit-
  test dir, so this checks the tangent FORMULAE the assembly reconstructs.
- DONE 2026-06-14 (N2/N3, DOC-only, commit 9c8423766c): N2 labelled the loop-
  test 100x separation as the MEASURED conservative floor per §5.1 (defect/
  bound ~3.0e3 at N=400, grows like N^2; not a derived 100). N3 re-confirmed
  computeMaxwellConjugateMicroPotential fully dead (zero live callers) and
  strengthened the RETIRED banner; kept on disk per §6.3.
- L5 — NOT FIXED (deliberate, guardrail §1.1/§12.2). The two formB_piexact
  PRJs on disk (ms33_modelI_dd1600_formB_piexact_2026-06-12.prj,
  ms33_modelVII_freeswelling_formB_piexact_2026-06-12.prj) carry TODO(Vinay)
  §12.2 provenance locators (E/nu, micro-EOS, Tuller geom, specific_surface,
  lambda) inherited verbatim from the base PRJs. They are NOT in Tests.cmake
  and STAY OUT pending Vinay's cited source locators (cannot be invented).
  Their calibration K=103879 J/kg is §12.1-clean (Dixon 2023 Fig.1).

OPEN (for Vinay): (1) the 1b cure verdict (M2) — run task42 1b_*_Kl with the
new binary; (2) whether to also flip enable_dsm_swelling_up_jacobian (the
pre-existing OFF swelling tangent), kept untouched here; (3) §12.2 locators to
register the piexact PRJs.

## 1b cure verdict (2026-06-14, the M2 headline)

Ran task42 1b_A_Kl.prj and 1b_B_Kl.prj (live K(rho_d) table, film_pressure_
coupling OFF, ScalarReferenceMassStorage; 1b_A = form (a) Off coupling, 1b_B =
form (b) kinematic + exact route) with the review-fixed binary (commit
936488482c), scratch ~/git/build/pi_fofnlev_fixes_20260614/cure_1b_test.
Pre-fix baseline (out_1b_*_Kl/run.log in task42): BOTH died in time step #1
("Newton: the linear solver failed"). MEASURED with the new binary:
 - 1b_B_Kl (exact route, H1+H2+M1+M2 all active): CURED — passes step #1 and
   keeps stepping (reached Time step #200+ with ZERO step-1 failures and zero
   nonlinear failures; observed VTUs out to t=2.59e6 s in an earlier run).
 - 1b_A_Kl (form (a), Off film_strain_coupling): NOT cured — still fails in
   time step #1 (terminated with error). The M2 live-K eigenstress tangent now
   fires for it (gate = live-K flag), so this confirms the review §3 "second,
   still unidentified mechanism" is NOT (only) the M2 eigenstress tangent for
   the form-(a) ScalarReferenceMassStorage case. Strong candidate: the L3
   dn_l/dK local-solve strain channel left as a DOCUMENTED partial tangent
   (f23f69c5b4) — Vinay's call whether to wire it next.
VERDICT: M2 (+H1/H2/M1) CURES the exact-route 1b_B; the operational-route 1b_A
remains a step-1 failure (open, points at L3 / a distinct mechanism).

UPDATE 2026-06-14 (L3 now WIRED): the L3 dn_l/dK local-solve strain channel —
the strong candidate this verdict flagged for 1b_A — has since been wired into
the displacement Jacobian (Jacobian-only; see the L3 WIRED worklog entry above).
It is the implicit-n_l(K) partner of the M2 eigenstress tangent and fires on the
1b_A OFF/operational regime (the M2 gate fired the explicit-K half there but the
local-solve n_l(K) half was missing). PREDICTED candidate cure for the 1b_A
step-1 divergence; NOT yet re-run (build is the next phase) — re-run task42
1b_A_Kl.prj with the new binary to confirm or refute (§5 predicted, not
verified). If 1b_A still fails after L3, the residual points to a genuinely
distinct mechanism beyond the live-K tangent gaps.

---

## 2026-06-23 — MS33 gating ctests: vtkdiff regression baseline added (DONE) + '--' parse-fix

DONE (commit c5615fcfd2, dsm_native_maxwell_conjugate, pushed 4 remotes):
the 6 registered MS LE gating ctests (ANCHORS_MS33_Model{I_dd1400/1600/1800,
III_gap2mm, IV_pellets, VII_freeswelling}) now carry an in-PRJ
`<test_definition>` vtkdiff regression baseline (CLAUDE.md §3, anchor (f);
Vinay-approved 2026-06-23). Mechanism (Applications/ApplicationsLib/
TestDefinition.cpp): `<regex>` is matched against the REFERENCE dir, then
vtkdiff compares output<->reference of IDENTICAL filename. Reference = the
single committed final-time output VTU per model (ts at t_end), regex pins the
end time. Two-tier tolerances (Vinay's choice), documented in a comment block
in each PRJ:
  TIGHT (geometry/equilibrium-pinned): displacement, saturation, porosity,
    transport_porosity, micro_porosity, micro_water_content, dry_density_solid.
  LOOSE (calibration/stress-path-sensitive; rel=1e-2, abs=1e3 Pa — INTENTIONAL,
    tolerates the documented ~1% Model-I non-reproducibility): sigma,
    swelling_stress, pressure, micro_pressure.
abs floors derived from extracted final-state problem scale. VERIFIED 6/6 PASS
on build maxwell_floor_20260619 (ogs @ 71366ac0d3; vtkdiff target built there);
adversarial teeth-check passed (perturbed reference -> vtkdiff FAIL; ref-vs-ref
PASS). CAVEAT: adaptive IterationNumberBasedTimeStepping pins ts_NNN to the
canonical toolchain -> a different CI toolchain may need reference regeneration.
Tests are LARGE-labelled (run via `ctest -L Large` or by name, NOT default ctest).

INCIDENT folded into the same commit: the prior provenance commit 30a1d0570f
wrote '--' (illegal double-hyphen) inside XML comments, breaking parse of 7
committed PRJs (the 3 gating Model I + confined_expulsion x2, dd1600 formB,
dd900). Fixed ' -- ' -> ' - '. The §6.7 gate did not catch it (it checks
header==live, not parse); default CI did not catch it (LARGE-excluded). Lesson +
proposed §6.7 parse-check: see memory incident_xml_double_hyphen_comment_parse.

---

## Merge note (2026-08-11): dsm_maxwell_jac_parallel folded into dsm_native_maxwell_conjugate

The sections BELOW are the `dsm_maxwell_jac_parallel` worklog (local-Jacobian
line, 2026-06-09/06-10), merged here verbatim. The sections ABOVE are the
`dsm_native_maxwell_conjugate` deliverable-line worklog (2026-06-09..06-23).
Both are kept in full per CLAUDE.md §6.4/§11 (AGENTS.md accretes; entries are
never removed). The two workstreams ran in parallel off the common base
d98f5f8324, which is why they appear as a conflict rather than a sequence.

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

---

## 2026-08-11 — dsm_maxwell_jac_parallel MERGED into dsm_native_maxwell_conjugate (DONE)

Consolidation ordered by Vinay: "everything in the maxwell_conjugate, tested,
verified, pushed. Then the rest deleted." Merge commit is a true 2-parent merge
of `deprecated/dsm_maxwell_jac_parallel` (tip 53538778cc, 5 commits) into
`dsm_native_maxwell_conjugate` (tip 35ebe2e415), common base d98f5f8324.

**Scope note — "the floor" is NOT a branch.** `~/git/build/maxwell_floor_20260619`
is a BUILD DIRECTORY compiled from the maxwell_conjugate worktree at 71366ac0
(`macro_porosity_floor` mandatory), already on the branch; `ogs-dsm-active` is a
byte-identical binary (both md5 727dfa40b016e154bd51e64c89d072c1). Nothing to
merge from either. Likewise `dsm_native_h_of_eps_wt` (detached 23a723cc3c) has
ZERO commits not already reachable from maxwell_conjugate.

### What the merge brought in
- `RichardsMechanicsFEM-impl.h`: analytic micro 2x2 local Jacobian
  (`evaluate_analytic_jacobian`) + file-scope `#pragma STDC FP_CONTRACT OFF`.
  AUTO-MERGED with no textual conflict despite +198/-53 (jac) vs +1005/-24 (mc).
- `ParallelVectorMatrixAssembler.cpp`: skip the per-thread `jacobian_assembler_
  .copy()` at num_threads_==1 so CompareJacobiansJacobianAssembler (owns a log
  ofstream, hard-OGS_FATALs on copy) is usable serially. mc never touched this
  file -> merged copy is byte-identical to the jac tip.
- `DSM/AUDIT_maxwell_local_jacobian_2026-06-09.md` (new; the ONLY file that was
  on the jac branch and not in mc's tree).
- `<scaling>true</scaling>` -> `false` on 9 MS33 PRJs.
- Conflicts (2), both resolved: this AGENTS.md (union kept, see the merge note
  above, §6.4) and `ms33_modelI_dd1800.prj` (see below).

### DEFECT the auto-merge introduced — FIXED (RichardsMechanicsFEM-impl.h:1216)
mc had converted EVERY `computeVanDerWaalsMicroPotential` call site to the live
K(rho_d) helper `effectiveAugmentationPrefactor(params, phi)` (0 raw-scalar call
args remain on mc). The jac branch forked BEFORE that sweep, so the three-way
merge — textually clean — left `evaluate_analytic_jacobian` alone on the
parse-time scalar `potential_augmentation_prefactor`. Under
`potential_augmentation_prefactor_live_dry_density=true` the analytic 2x2 would
then be the derivative of a DIFFERENT potential than the residual. Fixed to use
the same helper; bit-for-bit no-op when live mode is off (the default, and the
state of every MS33 suite PRJ, so the defect was LATENT — never active in the
deliverable suite). Verified: no bare `potential_augmentation_prefactor` remains
as a potential-evaluation argument anywhere in the file.

### `use_analytic_micro_jacobian` DEFAULT REVERTED true -> false (Vinay 2026-08-11)
MEASURED on this tree/binary (build mc_merge_20260811), not predicted. With the
analytic micro Jacobian ON it changes the ADAPTIVE TIME-STEP PATH on two of the
six gating models, which then fail the reference VTUs approved 2026-06-23:

| model  | steps ON | steps OFF / pre-merge | vtkdiff vs approved ref (ON) |
|--------|----------|-----------------------|------------------------------|
| dd1400 | 308      | 308                   | 11/11 PASS                   |
| dd1600 | 311      | 311                   | 11/11 PASS                   |
| dd1800 | 308      | 308                   | 11/11 PASS                   |
| III    | 376      | 405                   | **3/11 — FAIL**              |
| IV     | 637      | 637                   | 11/11 PASS                   |
| VII    | 675      | 682                   | **5/11 — FAIL**              |

Differences are time-discretisation scale, NOT moved physics (Model III max:
displacement 1.73e-6 m; sigma 1.86e4 Pa = 0.54% rel; swelling_stress 5.2% rel;
dry_density_solid 0.1% rel). The TIER-A tolerances (1e-9 abs on displacement)
were calibrated for a bit-identical step path and cannot survive a changed one.
ISOLATED by a 2x2 experiment {scaling} x {micro Jacobian} on Model III: the
Jacobian choice alone drives the step path (376 steps under BOTH scaling
settings) — it is not the linear-solver scaling flag. (Caveat for whoever
repeats it: the PRJ flag `fd_jacobian_for_exchange` is COARSER than the
constexpr — it also flips the block-#3 macro p-p tangent, giving 607 steps, so
it does not reproduce the pre-merge configuration.)
With the flag false the merged tree reproduces all 6 references step-for-step
identically to the pre-merge binary. Phase-D's 2026-06-10 "land it as default"
is PARKED, not withdrawn: re-enabling is one line and requires re-baselining the
III + VII reference VTUs first (§3 / §12.5 — Vinay's call). NOTE Phase-B's
"analytic == FD to round-off" was measured against the JAC-branch residual; mc
has since changed that residual (live K, strained film), so it is NOT re-verified
for this tree.

### dd1800 conflict resolution — `<scaling>` is the discriminator, not the solver
mc had moved dd1800 to BiCGSTAB+ILUT (scaling=true); jac had SparseLU
(scaling=false). MEASURED, merged binary, analytic ON:

| dd1800 linear solver           | scaling | result                          |
|--------------------------------|---------|---------------------------------|
| BiCGSTAB+ILUT                  | true    | FAIL ts #110, `residual: nan`   |
| SparseLU                       | true    | FAIL ts #110, `residual: nan`   |
| BiCGSTAB+ILUT                  | false   | ts #308, 11/11 PASS             |
| SparseLU                       | false   | ts #308, 11/11 PASS             |

i.e. exactly the Eigen IterScaling (Ruiz) overflow on the near-singular pressure
block that the jac Phase-D comment predicted, and independent of solver type.
Resolved by keeping mc's BiCGSTAB+ILUT and flipping ONLY `<scaling>` to false —
the measured discriminator. The other 5 suite PRJs keep jac's scaling=false
(verified 6/6 green); harmless with the analytic path parked, and required if it
is ever enabled.

### Verification of the committed state (all MEASURED)
- Build: clean, 0 errors; the only 2 warnings (`-Wunused-parameter` at :456 and
  :1849) are PRE-EXISTING — byte-identical signatures in the mc parent.
- Unit tests: `testrunner` 1422 tests / **1418 PASSED, 0 FAILED**, 4 skipped
  (all pre-existing GTEST_SKIPs).
- MS33 gating suite: **6/6 models, 66/66 field comparisons PASS** against the
  committed references, at each PRJ's own `<test_definition>` tolerances; step
  counts 308/311/308/405/637/682 = pre-merge baseline exactly.
- All 31 MS33 PRJs pass `xmllint --noout`; micro+macro floor tags paired in all 6.
- METHOD WARNING for future runs: staging a model directory by copying it
  wholesale puts the committed reference VTU next to the run outputs, and a
  "latest matching file" pick then compares the reference AGAINST ITSELF (a
  trivially-passing test, §3). This happened in the first pass here and inverted
  the III/VII verdict. Stage MESHES ONLY (exclude `*_ts_*`), and have the
  comparison refuse ref==out.

### Deletion (Vinay 2026-08-11: delete jac_parallel only)
- `deprecated/dsm_maxwell_jac_parallel` + its 4 remote copies: SAFE once this
  merge is pushed — the 5 commits stay reachable through the merge's 2nd parent
  and its one unique file (the AUDIT .md) is now in mc's tree.
- (The rest were initially held back as "protected content would be lost". That
  is SUPERSEDED — see the next section, same day.)

### SUPERSEDES the above — full retirement DONE 2026-08-11

Vinay extended the instruction: "commit what's possible in the rest and then
delete them too. it's ok if they stay on git and marked deprecated." Executed.
FOUR branch lines retired; every one is recoverable BY NAME from an annotated
tag that states its deprecation, and all five tags are pushed to ALL FOUR
remotes (origin, github, backup, vgk2):

| retired ref (and where it was deleted from) | tip | recover from |
|---|---|---|
| `dsm_maxwell_jac_parallel` (local + 4 remotes) | 53538778cc | `archive/dsm_maxwell_jac_parallel_2026-08-11` (also 2nd parent of 6135bf66c6) |
| `dsm_native_Pi_fofnlev` (4 remotes) | a8ffdeefcd | `archive/dsm_native_Pi_fofnlev_branchtip_2026-08-11` (also ancestor of mc) |
| `deprecated/dsm_native_Pi_fofnlev` (local + vgk2) | 19c031cc1f | `archive/dsm_native_Pi_fofnlev_2026-08-11` — holds the 117 unique .prj/.vtu/.md |
| `deprecated/dsm_native_Pi_fofnlev_review_fixes_2026-06-14` (local + vgk2) | d9a017cbe1 | `archive/dsm_native_Pi_fofnlev_review_fixes_2026-08-11` — 1 unique .prj |
| `dsm_native_h_of_eps` (github, backup) | 23a723cc3c | `archive/dsm_native_h_of_eps_2026-08-11` (also ancestor of mc) |

Restore any of them with
`git branch <name> <tag>^{commit}` and push where wanted.

WHY THE TAGS MATTER: before this, `deprecated/dsm_native_Pi_fofnlev` and its
117 protected files existed ONLY on local + vgk2 — no redundant copy anywhere.
They now have four. The tags were verified present on all four remotes BEFORE
any branch was deleted.

Worktrees removed: `dsm_native_h_of_eps_wt` (its 7 untracked files committed
first, as 9b179d1ddc) and the throwaway `mc_merge_wt`. Remaining worktrees:
the `ogs` master checkout and `dsm_native_maxwell_conjugate_wt`.

SCOPE CORRECTION (an overstatement made in-session and corrected here): retiring
these four does NOT leave maxwell_conjugate as the only DSM branch — that is
true of LOCAL branches only. Roughly a dozen older DSM-named branches remain on
the remotes (`deprecated/dsm_native`, `dsm_native_hierarchical`,
`dsm_native_tuller_macro_film`, `dsm_native_tuller_review`,
`dsm_native_pdisj_maxwell`, `dsm_native_pdisj_aug_tuller`, `dsm_mfront*`,
`dsm-nb-*`, `DSM`, `salvage/macmini-*`). They were never in this task's scope
and were NOT touched.

STILL NOT DELETED — the build dirs, deliberately, pending Vinay's explicit call.
They cannot "stay on git", so the stated safety condition cannot be met for
them: `maxwell-conjugate-20260602` (2.6G, citation binary md5 c432a156) and
`maxwell_floor_20260619` (2.6G, md5 727dfa40, the toolchain that produced the
MS33 reference VTUs) are named by path in 13 TRACKED files, including all six
gating PRJ provenance headers and calibrate_maxwell_K.py.

### OPEN (carried forward, not introduced here)
- A failed local 2x2 silently returns the decoupled PREDICTOR state
  (`return out.converged ? out : predictor;`) and NO caller inspects
  `converged` — no warning, no time-step rejection. Pre-existing on both
  parents; raise with Vinay whether it should hard-fail.
- The file-scope `#pragma STDC FP_CONTRACT OFF` (from jac) now governs ~1000
  more lines of mc-only assembly code than on either parent. Directionally safe
  (disables FMA fusion -> more determinism) and the gating suite is unaffected
  (6/6 PASS), but it was never validated at this scope on the mc line.
- u-side analytic blocks remain `enable_dsm_swelling_up_jacobian = false`
  (unsafe per jac Phase D: mIII singularizes, mIV/mVII solution-shift).

## 2026-08-12 — BEACON ctest inputs restored, floors declared, references re-baselined (DONE)

**Symptom (Vinay's question).** `ProcessLib/RichardsMechanics/Tests.cmake` has
registered nine BEACON `AddTest` blocks since the port, but
`Tests/Data/RichardsMechanics` carried **zero** beacon files on this branch.

**Root cause — a half-port, not a deletion.** `git log --diff-filter=D` (with
`-m` for merges) shows no commit ever removed them; the commits that added them
(`45ea35b9c9`, `5a0792dbb0`, `72f4f3a192`, `473e4af39c`, `87dc1ae9bd`) are not
ancestors of this branch. `cc68e104c9` "Port DSM native (p_disj augmentation +
Tuller) onto ufz/master 6.5.8" (2026-06-02) carried the ProcessLib code and the
Tests.cmake registration but no `Tests/Data` payload. `master` is
self-consistent (neither registration nor data); only this branch was half-ported.

**What was done.**
1. 15 files restored from `a6dfec842c` (`deprecated/dsm_native_hierarchical`):
   5 PRJ, 4 reference VTU, 3 domain VTU, 3 GML.
2. As restored they do not parse: `Key <macro_porosity_floor> has not been
   found`. The decks predate `2144a1ce40` / `71366ac0d3`, which made the
   floor pair mandatory with no default. Both floors declared **inert (0.0)**,
   value approved by Vinay 2026-08-12 (CLAUDE.md §1.1) — the same inert
   declaration used for the MS33 form-(a) decks on 2026-06-19. Verified inert
   at `PotentialExchange.h:199-205`: floor 0 takes the unfloored path.
3. With the floors the decks run, but the 2026-06-02 references no longer
   match: 1a01 4/7 fields fail (pressure abs 405.8 Pa, rel 4.2e-4), inflow
   8/10, 1c 9/9, 1b 1/7 (micro_pressure abs 2.13e-16 against a 1e-16
   threshold, i.e. roundoff). The old baselines predate the constitutive work
   on this branch — `013990e1cd` (mu_lR on the integrable Maxwell partner),
   `aeccc1c838` + `36dfeaddcc` (swelling stress -> (1-phi_M)*p_film, sign and
   weighting), `dbf20d6f1a` (bare-Pi OFF path retired), `6391e357a2` (marked
   RESIDUAL-CHANGING). The drift ordering is consistent with that: the decks
   whose vdW term is switched off drift least, the ones with it live drift most.
4. All four references **re-baselined** on `40551b6aad`. Determinism checked
   first: repeat runs are md5-identical. These are regression baselines at that
   commit, NOT physics claims. Previous baselines remain at `a6dfec842c`.
5. §12.2 provenance headers written into all five decks — **honest and
   incomplete**, see below.

**Verification (all MEASURED, binary `40551b6aad` in
`/Users/vinaykumar/git/build/maxwell_floor_20260619`).**
- `ctest -R beacon -j14` -> **13/13 pass**, before and after the headers.
- `ctest -R RichardsMechanics -j14`: failure set with the change is **identical**
  to the clean-tree baseline (12 non-beacon failures, stable over 2 clean and 2
  changed runs). No test that passed before fails now (§6.5 satisfied).
- Flake noted, not attributed: `Mockup2D/mockup` and `RichardsMechanics/A2`
  failed in one early run and in none of the four later ones, and both pass
  serially. They are order/parallelism sensitive (they consume state written by
  a preceding test); adding 13 beacon tests perturbs scheduling.

**§12.2 status — NON-COMPLIANT, BLOCKED ON VINAY.** An 11-agent provenance
mining run (5 tracers, 5 adversarial verifiers, 1 synthesis) found that exactly
**one** material literal in the whole family is sourceable: `specific_surface =
523` in the inflow deck (EPFL; Seiphoori, Ferrari & Laloui 2014, Géotechnique
64(9):721-734, Tab. 1 p.724, MX-80) — and even that carries an open m²/g vs
m²/kg convention question that also affects the MS33 decks. Everything else is
`TODO(Vinay): UNSOURCED`. The headers say so explicitly rather than inventing
attributions. Full per-literal evidence, including the attributions that were
tried and refuted (2780 <- "EURAD-2 MS33 spec"; nu=0.2 <- the CIMNE ν=0.3
locator; 2650 <- the Boom Clay grain density in GRS-202; 6e-20 and 5.1e-21 <-
Israelachvili & Adams 1978):
`~/.claude/projects/-Users-vinaykumar-git-ogs/BEACON_PRJ_PROVENANCE_LEDGER_2026-08-12.md`.

Notable findings folded into the headers:
- **No deck has a calibration anchor at all.** None carries
  `potential_augmentation_prefactor`, so K = 0 and the augmentation branch is
  skipped (`PotentialExchange.h:259`). §12.1's Dixon/Villar rule is *vacuous*
  here, not violated. The `swelling_pressures` literals are inputs to the macro
  swelling law, not calibration targets.
- **A = 1e-30 J with Sa = 1.0 is a deliberate switch-off placeholder**: exact
  zero is rejected by the positivity guard, so a minimal positive pair is the
  only way to disable the vdW branch. The same quintuple is the DSM unit-test
  harness state (`DSMMicroMacroSingleIntegrationPoint.cpp:663`). Still owed a
  §0.2 exemption note recording the intent.
- **The family carries three different Hamaker constants** (1e-30 / 6e-20 /
  5.1e-21 J), none equal to the branch literature value; the 2026-05-22 sweep
  `0d579e8aeb` that moved 27 other PRJs to 2.2e-20 missed these decks.
  stressprobe's 6e-20 sits outside the range the source file itself documents.

### OPEN (introduced by this entry)
- §12.1 locators, or a §0.2 exemption declaring the parameter set non-physical,
  for every literal marked `TODO(Vinay): UNSOURCED` in the five headers.
- `specific_surface` unit convention (m²/g as entered vs m²/kg as consumed) —
  affects the MS33 gating decks identically.
- Whether beacon_1c's block and pellet share one grain density, E, ν and Bishop
  law (currently they do; never recorded as a decision).
- `RichardsMechanics_double_porosity_swelling_dsm_micromacro_constbc_reference`
  is broken by the same half-port (missing `..._constbc.xml` and its reference)
  and was left untouched — out of scope for this task.

---

## 2026-08-12 — five desk items LANDED (Vinay-ratified), verified 22/22 (DONE)

All measured on isolated worktrees/builds first (verification dossier:
vk-claude-workbench handoff/macbook-pro/2026-08-12.md 15:40), then landed as
five commits and re-verified via ctest at the landed tip — 22/22 PASS
(gating 6 + beacon 9-family + double_porosity + MCC shear/biax).

1. **kofdd/kref20x trio DE-REGISTERED** (Tests.cmake, decks kept). Cannot pass
   as registered (no test_definition -> parse hard-fail under the wrapper);
   the two ModelIV variants additionally DIVERGE on the merged code (ts #825
   FD / #2333 analytic) — broken decks. Re-register only with ratified refs.
2. **beacon abs gates 1e-16/1e-14 -> 1e-12** (33 rows). Measured: same-compiler
   0.0 exactly; gcc epsilon 4.4e-16..7.8e-16. Refs unchanged; constbc rows
   untouched.
3. **double_porosity_swelling re-baselined** (2 refs; md5 85d171fa/64329c93).
   Old refs encoded pre-hierarchical-split physics; stale by identical amounts
   on clang+gcc and pre/post-merge binaries.
4. **MCC: LM KEPT, gates re-derived** (shear/biax sigma+PCP abs 1.0 Pa /
   rel 1e-7 = measured x1000), and the 2a2409c043 "response-neutral" claim
   corrected in-file per §5 (measured rel <= 6.6e-11, not bit-identical;
   attribution proven by one-line revert experiment; 372da0aafe touches only
   absP, never in play for these tests).
5. **Analytic micro 2x2 Jacobian RE-ENABLED** (Phase-D default restored) with
   III/VII re-baselined in the same commit (md5 2aac0784/4ab34801; run-to-run
   AND cross-build bit-identical). dd1400/1600/1800/IV pass existing refs
   11/11 under the analytic path; unit suite 1418/1418. u-side blocks stay OFF.

**MECHANISM LESSON (cost one amend):** the OgsTest vtkdiff wrapper enumerates
every REFERENCE-side file matching the test_definition regex and requires a
same-named OUTPUT for each. A superseded reference left beside the PRJ fails
the suite by construction. Superseded III/VII refs therefore moved (git mv,
never deleted) to superseded_references_2026-08-12/. First post-land ctest
caught this (III/VII failed); fixed and re-verified 22/22.

Still open on Vinay's desk: constbc-family restore-vs-de-register; gcc
confirmation of the new III/VII baselines at TIER-A tolerances (ogs09
report-back; if gcc exceeds them, same portable-gate policy as the beacon
gates applies).

---

## 2026-08-17 — live-K adopted for MS33 III and IV; the swelling-branch shortfall (Vinay's ruling)

**Ruling (Vinay 2026-08-17):** Models III and IV report **live** `K_aug(rho_d)`. Frozen-K is the
physically wrong branch — a swelling material loses swelling potential, and carrying K at the
emplacement density prices the material at a density it no longer has. The MIXED 2026-06-23
frozen-K choice for III/IV is retired; matching the team band was never the criterion.
**Model I keeps frozen-K** for a physics reason: measured dRho_d = +0.00% in the confined
calibration cell, so live-K is exactly inert there.

**Measured, tip ac1808936d, all K Dixon-propagated (sigma0=0 Option-C form-(a) chain
900->3505.642413 | 1400->46000 | 1600->104689.9129 | 1800->265905.06), 200 d, canonical probes:**

| case | frozen-K | live-K (ADOPTED) | teams | dRho_d | K collapse |
|---|---|---|---|---|---|
| III soft gap (delivered geom) | 9.705/9.703/9.675 | **7.663/7.661/7.592** | 2.63-7.47 | -8.5% | 104690->64716 |
| III gap-switch (empty gap) | 7.202/7.400/7.824 | **1.176/1.336/1.658** | " | -14.8% | 104690->42893 |
| IV block1600/pellet900 | 6.227/6.230/0.499 | **4.269/4.257/0.700** | 1.51-3.17 / 1.40-3.04 / 1.20-2.64 | block -12%, pellet +14% | block 104690->~47k |

Live-K IMPROVES team agreement on the delivered geometries (III 9.71->7.66, IV 6.23->4.27) while
being the correct branch. No arm diverges — the 2026-06-23 IV divergence was removed by the
analytic micro Jacobian (bfd52cf6ff).

**THE SHORTFALL THIS TEST DISCOVERED.** Dixon K(rho_d) is an equilibrium **compaction** curve —
each point measured on a sample *compacted to* that density. Applying it live to a sample that
**decompacted to** the same density assumes path-independence, which bentonite does not have: a
specimen swelled 1600->1364 retains aggregate fabric and interlayer water a freshly compacted 1364
specimen never had, so its swelling potential lies ABOVE the compaction branch. Imposing the
compaction branch on a swelling path must therefore UNDER-predict — the observed sign. Frozen-K
over-stiffens, live-K over-softens, the teams sit between: the pair brackets the answer.

**Missing term:** a swelling-branch K_sw(rho_d, path) relaxing toward Dixon as fabric is destroyed.
Candidate home: Pi_struct in the disjoining decomposition (paper eq:disjoining_decomposition), or a
state variable carrying compaction history. FORMULATION item, not a tuning knob. Recorded in
memory `feedback_livek_is_the_physics_truth.md`; carried in dsm_defense_2026-06 as its own frame
("The swelling branch --- a constitutive shortfall this benchmark exposed") and as the first
open item.

Reporting rule that follows: never present the frozen-K number because it agrees better. State
which branch is physical and what the residual gap reveals.

---

## 2026-08-17 — "exact interpolation" + gap-switch shipped to ctest (Vinay's directive)

Vinay, 2026-08-17: *"exact interpolation. gapswitch everywhere, even in ctests. model3: Compile
the results for both live and frozen k. I think live k is the correct physics. that model 3 i
need to check again. compile the results and then check elastoplasticity parameters and
equilibration for all model runs."*

### 1. The rho_d=900 anchor was WRONG. Corrected. DONE 2026-08-17

`calibrate_maxwell_K.py` carried `DIXON_EMDD_MPa = {..., 900: 0.051}`, commented
"dd900 unused (model edge)" — a placeholder, never a Dixon reading. The Dixon (2023) Fig. 1
median fit `Ps = 0.003*exp(5.2883*rho_d[g/cm3])` gives at rho_d = 900:

    0.003*exp(5.2883*0.9) = 0.3500522009 MPa      (a factor 6.9 above the placeholder)

Re-calibrated against the exact value (`--rel-tol 5e-5`; the script's 2 % default is why the
first attempt stopped at -0.380 %):

    K(900) = 4367.2277 J/kg  ->  Ps = 0.3500565 MPa   (rel +0.001 %)   SUPERSEDES 3505.642413

Reproduced independently on both builds — maxwell-conjugate-20260602 gives 4367.227700212952 and
verify_tip_20260812 gives 4367.227700091862 (agree to 1e-10 relative), so the cell is insensitive
to the binary difference. **RATIFIED (Vinay, 2026-08-17): K(900) = 4367.2277 J/kg.** See §16.

This knot is NOT decorative: Model III's live rho_d falls to 1347-1376 kg/m3, i.e. BELOW the
1400 knot, so the 900-1400 segment is actively read.

### 2. Two defects found in the tracked calibration path. NOT yet fixed in tracked files

  a) `ANCHORS_MS33_ModelI/ms33_modelI_dd900.prj:247` has `<prefix>ms33_modelI_dd1600</prefix>`
     — the dd900 deck writes its output under the dd1600 prefix. It would silently overwrite
     dd1600 output, and it is why the first calibration reported "initial run failed" (the
     driver globbed for dd900 outputs that were never written). Also carries sigma0 = -1.5e5
     where an Option-C calibration cell should be 0.
  b) `calibrate_maxwell_K.py:21` still points at `build/maxwell-conjugate-20260602/bin/ogs`
     (built Jun 17), not the canonical tip. Every K in the chain was fitted through that
     pointer. Harmless for dd900 as measured above, but a §6.7 traceability defect.

### 3. Model III GAP-SWITCH is now the registered ctest. DONE 2026-08-17

`ms33_modelIII_gapswitch.prj` + `ms33_modelIII_gapswitch_bc.py` staged into
`Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIII/`, given a full §12.2 provenance header and
an 11-field `<test_definition>` mirroring the soft-gap deck's field/tolerance set. Registered in
`Tests.cmake` at RUNTIME 120 (measured 16 s wall, 1393 steps). The soft-gap registration is
COMMENTED OUT with a deprecation note — deck, reference and history retained (§6.2/§6.3).

Shipped configuration = gap-switch + **live-K**, so the regression baseline guards the physics
that is actually reported.

Verified: `ctest -R ANCHORS_MS33` = **6/6 pass** on verify_tip_20260812 (the wrapper's
reference-side enumeration, which broke a previous landing, is clean here). The verify worktree
used for that run was restored to clean afterwards.

**MOOT — this specific baseline was never committed.** `ms33_modelIII_gapswitch_ts_1392_
t_17280000.000000.vtu` (md5 `7efcdcbe1f0640ea3e6e9aaa6109ed45`) was the 20x-era reference; it was
superseded by the 100x landing (§8, `ts_890`) before any ratification, which was itself superseded
by the 160x landing (§15, `ts_832`, md5 `e60493e974dee9f9afb2e89ec48ea651`) before THAT was
ratified either. `ts_832` is the only III reference that was ever actually committed, and it is
now RATIFIED (Vinay, 2026-08-17) — see §16.

### 4. Compiled results — canonical probe, both K branches

Probe is verbatim `eurad-anchors/.../extract/extract_suite.py:30-54`: single NEAREST NODE at
x = 0, y = 0.070 / 0.04025 / 0.0105 m, Ps = -tr(sigma)/3 [MPa]. (A row-mean probe gives a
visibly different B and must not be mixed with these.)

| case | frozen-K | live-K (ADOPTED) | teams | steps f/l |
|---|---|---|---|---|
| **III gap-switch (SHIPS)** | 7.202 / 7.400 / 7.824 | **1.176 / 1.337 / 1.660** | 2.63-7.47 | 589 / 1393 |
| III soft-gap (DEPRECATED) | 9.705 / 9.703 / 9.675 | 7.663 / 7.661 / 7.592 | " | — / 650 |
| **IV pellets** | 6.231 / 6.234 / 0.503 | **4.295 / 4.280 / 0.735** | 1.51-3.17 / 1.40-3.04 / 1.20-2.64 | 461 / 1632 |

Model IV numbers supersede the 2026-08-17 morning table (6.227/6.230/0.499 and
4.269/4.257/0.700), which used the erroneous 3505.642413 pellet anchor. T and C move <1 %;
the pellet-dominated B probe moves ~+5 %.

**Model III — the thing Vinay wanted re-checked.** The two corrections push OPPOSITE ways:

    soft-gap -> gap-switch  (at frozen-K):   9.705 -> 7.202     (geometry fix, lowers)
    frozen-K -> live-K      (at gap-switch): 7.202 -> 1.176     (physics fix, lowers hard)

The delivered soft-gap live-K number (7.66) sat at the TOP of the team band only because a
too-stiff gap surrogate was cancelling part of the live-K softening. On the correct geometry with
the correct K branch, Model III lands at 1.18-1.66 MPa, i.e. **BELOW** the 2.63-7.47 team band —
it no longer brackets. The fabric-memory shortfall is therefore much larger on Model III than the
soft-gap arm suggested. This is a finding, not a regression: it is the same missing swelling-branch
K_sw term, measured without the cancelling error.

### 5. OPEN — K-table interpolation is LINEAR, Dixon's law is EXPONENTIAL. ASK(Vinay)

`effectiveAugmentationPrefactor` (PotentialExchangeParameters.h:380-393) reads the table through
`MathLib::PiecewiseLinearInterpolation`, i.e. a straight chord between knots, clamped outside.
Dixon's law is exponential in rho_d (fitted ln K = 4.554e-3*rho_d + 4.302 through the four
calibrated knots). A chord across a convex function always lies ABOVE it, so K is biased HIGH
wherever a model operates mid-segment. Measured at the models' actual operating densities:

| where | rho_d | K linear chord | K log-linear | chord high by |
|---|---|---|---|---|
| III gap-switch p05 | 1348.5 | 41710 | 36089 | **+15.6 %** |
| III gap-switch median | 1365.0 | 43082 | 39002 | **+10.5 %** |
| IV pellet region p05 | 964.9 | 9773 | 5929 | **+64.8 %** |
| IV median | 1215.7 | 30653 | 19312 | **+58.7 %** |

Both shipping models straddle the widest chord (900-1400). Consequence, and it matters for the
write-up: the live-K stresses reported above are computed with K biased HIGH, so a faithful
exponential Dixon would push them LOWER still. The under-prediction is therefore a floor, not a
ceiling — the shortfall argument is strengthened, not weakened.

Two ways to close it, both Vinay's call:
  (a) DENSIFY the knot set (e.g. 100 kg/m3 spacing over 900-1800) so the chord tracks the
      exponential to ~2.6 % worst case. Costs 6 more calibration runs on the same approved
      Dixon anchor and the same approved procedure — no new source, but 6 new literals to ratify.
  (b) LOG-INTERPOLATE inside `effectiveAugmentationPrefactor` (store ln K, exp on read) and
      update `effectiveAugmentationPrefactorPhiDerivative` to match. One-line physics change,
      no new literals, but it changes every live-K result and the analytic Jacobian tangent.
Not chosen here — a formulation decision per §9.

### 6. CORRECTION to §4 above — Model III is NOT equilibrated at the 200 d spec horizon

Raised by the 2026-08-17 equilibration audit and then verified directly by extending both arms to
2000 d (10x the spec horizon; suction curve extended so the held-at-0 branch covers it; diagnostic
runs only, the shipped deck and its ctest stay at the spec 200 d).

**Model III gap-switch, Top/Central/Bottom [MPa]:**

| horizon | frozen-K | live-K | note |
|---|---|---|---|
| 200 d (spec) | 7.202 / 7.400 / 7.824 | 1.176 / 1.337 / 1.660 | transient |
| **2000 d (equilibrated)** | **8.041 / 8.026 / 8.173** | **2.180 / 2.169 / 2.087** | asymptote |
| 200 d error | **-10.4 %** | **-46.0 %** | |
| teams | 2.63 - 7.47 | | |

At 200 d the live-K increments were still GROWING (ratio +1.57) — the run had not entered its
decaying phase, so no asymptote could be read from it at all. By 2000 d the ratio is +0.331 and the
last increment is +0.8 %, i.e. converged to well under 1 %. Final rho_d is the SAME in both arms
(1363.4 vs 1363.5) — the arms differ only in which K that density is read at, not in where they end
up geometrically.

**This overturns the "III falls far below the team band" reading in §4.** That was a transient
artefact of stopping at 200 d. Equilibrated, the two arms BRACKET the team band cleanly:

        live-K 2.18  <  [ teams 2.63 ... 7.47 ]  <  8.04 frozen-K

which is exactly the bracketing statement the shortfall argument predicts, now measured on the
correct geometry. The residual gap between live-K 2.18 and the band floor 2.63 is the
fabric-memory deficit — 17 %, not the 55 % the 200 d number implied.

**Model IV IS equilibrated at 200 d** — same 2000 d check: frozen 6.231/6.235/0.503 -> 6.233/6.234/
0.503 and live 4.295/4.280/0.735 -> 4.323/4.282/0.736, i.e. within 0.6 % on every probe. The §4
Model IV numbers stand; use the 2000 d values if an equilibrated set is wanted.

**OPEN — ASK(Vinay), and it decides which numbers ship.** Do the other teams report Model III at
200 d or at equilibrium? If the band 2.63-7.47 is a 200 d band, then comparing our 2000 d value to
it is invalid, and the honest finding flips to *"our hydraulic transient is far slower than the
teams'"* — which would point straight at the intrinsic permeability (already 20x spec here as an
"equilibration aid") and at the `min_relative_permeability = 1e-2` floor that the equilibration
audit measured as pinned at every node in every run, i.e. the MS33 spec S_e^3 law never actually
operates. Both readings are defensible; they need the spec's reporting convention to separate them.
Until that is settled, no Model III number should be quoted without its horizon attached.

### 7. The permeability acceleration — history, and why 20x was not enough for Model III

Vinay, 2026-08-17: *"200d is the team horizon. didn't we accelerate the permeability in order to be
able to run these at the 200d mark? dig into history"* — yes, and the record is explicit.

**The written rationale**, carried verbatim in six Model IV decks and the Model III soft-gap deck
(`ms33_modelIV_pellets.prj:24-29`):

> GUARDRAIL §12.5 audit (2026-06-09, Vinay-approved): intrinsic permeability set to 20x spec
> (k0=1.17406e-19 m^2 clay; 1.17406e-18 pellet). Rationale: **spec-perm = unconverged transient at
> reporting horizon; 20x equilibrates** and is the MAX stable factor on dsm_native_maxwell_conjugate
> (**50x diverges ~11d, Model III**, always-on Maxwell-term wetting-rate ceiling). **IV/VII endpoints
> perm-insensitive (20x vs 50x <1%)**. Submission candidate, EURAD-2 MS33.

The gap-switch deck carries the short form: *"intrinsic permeability 1.17406e-19 : 20x the MS33 spec
value, an equilibration aid carried over from the gap2mm header"*. Spec k0 = 1.17406e-19/20 =
**5.8703e-21 m²**.

**Measured today — Model III gap-switch, permeability sweep, canonical tip binary, T/C/B [MPa]:**

| k0 | horizon | frozen-K | live-K | vs band 2.63–7.47 |
|---|---|---|---|---|
| **1x SPEC** | 200 d | 0.184 / 0.138 / −0.780 | 0.290 / 0.290 / 0.101 | nowhere near |
| 20x (shipped) | 200 d | 7.202 / 7.400 / 7.824 | 1.176 / 1.337 / 1.660 | live below |
| 50x | 200 d | — | 1.932 / 1.993 / 2.033 | below |
| **100x** | **200 d** | **8.023 / 8.009 / 8.106** | **2.162 / 2.170 / 2.150** | brackets |
| 20x | 2000 d | 8.041 / 8.026 / 8.173 | 2.180 / 2.169 / 2.087 | brackets |
| 100x | 2000 d | — | 2.195 / 2.185 / 2.161 | brackets |

Three findings, all measured:

1. **The rationale was right about spec-perm.** At 1x, Model III barely hydrates by 200 d — ρ_d falls
   only 1600 → 1553, the gap has scarcely closed, and Ps is 0.1–0.3 MPa. The acceleration was
   genuinely necessary.
2. **But "20x equilibrates" was never checked on Model III.** The evidence the header offers is
   *"IV/VII endpoints perm-insensitive (20x vs 50x <1%)"* — Model I/IV/VII, all of which today
   measure equilibrated at 200 d (IV within 0.6 %). Model III, the one model with a moving outer
   boundary and therefore the longest transient, was outside the tested set. At 20x its live-K arm
   reaches only 54 % of its asymptote by 200 d.
3. **The 20x cap is stale.** The header caps at 20x because *"50x diverges ~11d, Model III"*. On the
   current tip 50x completes (998 steps) and so does 100x (891 steps). The most likely cause of the
   change is the analytic micro Jacobian (`bfd52cf6ff`) — the same landing that removed the Model IV
   divergence recorded in §4 above. **Predicted, not verified:** attributing the recovery
   specifically to that commit would need a revert experiment; what is verified is that 50x and 100x
   run to completion today.

**The endpoint is k0-independent — the §2-clean test the 2026-06-09 claim asserted but only ran on
IV/VII.** 100x at 200 d reproduces 20x at 2000 d to within 0.8 % on T/C (live 2.162 vs 2.180;
frozen 8.023 vs 8.041), and 100x carried out to 2000 d moves only a further 1.5 %. Final ρ_d is
1363.5 in every accelerated arm. So raising k0 is a pure time-rescaling here, not a change of
answer — which is what makes it legitimate, and what licenses reporting a 100x run at 200 d.

**Why ~100x, and why that number is not arbitrary [DERIVED, consistent with the measurements].**
Macro saturation collapses to ≈ 0 at every node and every frame in these runs (the DSM exchange sink
consumes macro water as fast as it arrives), so `S_e³ → 0` and `relative_permeability` sits exactly
on its `min_relative_permeability = 1e-2` floor throughout — measured, all frames, all four decks.
The effective macro conductivity is therefore k0 x 1e-2, and the shipped 20x config runs at
20 x 1e-2 = **0.2x** the conductivity a conventional single-porosity model would have at spec k0 with
k_r ≈ 1. We were not running 20x fast; we were running **5x slow**. 100 x 1e-2 = 1.0 restores parity
with the spec-intended conductivity — which is exactly where the transient compresses into 200 d.
This also explains how the other teams reach 2.63–7.47 MPa at 200 d on spec permeability while we
reach 0.3: their k_r is not floored because their macro pore does not collapse.

**What does NOT change.** At the team horizon with the transient properly compressed, live-K Model III
is 2.16 — still below the 2.63–7.47 band — and frozen-K is 8.02, above it. The bracketing statement
and the fabric-memory shortfall are therefore *not* transient artefacts; they survive the correction.
Only the magnitude changes (deficit 17 %, not 55 %).

**ASK(Vinay), the smallest question:** raise Model III to k0 = 100x spec so its 200 d reported value
is the equilibrated one, and record the k_r-floor compensation as the stated reason? That is a
material-parameter change under §12.5 and a re-baseline of the gap-switch ctest reference, so it is
your call. The alternative — report the 20x 200 d transient as-is — is defensible only if the number
is labelled a transient every time it appears, which no current artifact does.

### 8. DONE 2026-08-17 — Model III raised to k0 = 5.8703e-19, k_r-floor reason recorded

Vinay, 2026-08-17: *"go with 100x, record the k_r floor reason"* — the ASK at the end of §7 is
CLOSED. Landed in `ANCHORS_MS33_ModelIII/ms33_modelIII_gapswitch.prj`: three literals
1.17406e-19 -> **5.8703e-19**, plus a §12.5 provenance block giving the full reason.

**Corrections to §7, from the git/spec archaeology (19 agents, 135 facts):**

1. **The multiplier baseline was loose in §7 and in my first header draft.** The MS33 spec LITERAL
   is **k0 = 5.6e-21 m² at phi_ref = 0.42** (`eurad2_MS34/MSXX/c_theoretical_benchmarking/
   theoretical_benchmarking.tex:118-120`; MS33 revision `MS33/theoretical_benchmarking.tex:440-441`).
   The decks' 5.8703e-21 is that spec point Kozeny-Carman-transferred to phi0 = 0.4244604316546763,
   KC factor 1.04826 (`ms33_modelI_dd1600.prj:308`). So the new value is **100x the KC base but
   104.83x the spec literal**, and the old one was 20x KC / 20.97x literal.
   **Corpus hazard:** the SAME absolute value 5.8703e-19 is labelled `k0 x100 spec` in
   `DSM/AGENTS.md:189` and `k0 x5` in `runs/2026-08-05_2017_.../workbooks/PROVENANCE.md:79`. There is
   no 5x-spec value anywhere in the corpus. **Never quote "Nx spec" without naming the baseline.**
2. **The deciding commit is `2c3b3bb0d6` (2026-06-09)**, message: *"Spec-perm runs are unconverged
   transients at the 200-240 d horizon; 20x equilibrates and is the max stable factor on this branch
   (50x diverges ~11d for Model III: always-on Maxwell-term wetting-rate ceiling). IV/VII endpoints
   perm-insensitive (20x vs 50x match <1%)."* Model III appears in that sentence ONLY as the model
   that diverged at 50x — it was structurally excluded from the test used to justify the setting.
3. **Today's finding is the THIRD falsification of "20x equilibrates", not the first.** 2026-06-17
   measured 20x as under-equilibrated for Model VII
   (`runs/2026-06-15_1410_.../EQUILIBRATED_LE_VII_2026-06-17.md:9-15`); 2026-06-25 documented **40x
   spec as the minimum** that equilibrates the confined Reference within 200 d. Neither produced a
   proposal to raise the shipped decks. The 40x floor independently brackets the 100x adopted here.
4. **The 1e-2 k_r floor was introduced for a NUMERICAL reason, not an equilibration one** — Vinay
   2026-05-29, `ms33_modelI_dd1600.prj:196-198`: *"so the drained macro cannot give a singular
   k_rel=0 at step 1."* The spec prescribes no bentonite floor at all; its only k_r floor (1e-6) is
   for granite (`MS52/backup.tex:175`). Its dominance of the hydraulic timescale is a side effect —
   which is exactly the effect §7 identifies and the new value compensates.
5. **Model I was never accelerated**, and the acceleration is not uniform inside Model VII (nine
   variants sit at 1.00x KC base). The Model IV **pellet** zone is not 20x either: 1.17406e-18 =
   10x the clay = 200x KC base, and the 2026-06-08 "10x clay" move actually LOWERED it below its own
   KC value (0.785x) — the whole 15.7x excess comes from the 2026-06-09 20x.

**Verified after the change:** `ctest -R ANCHORS_MS33` = **6/6 pass**. New gap-switch reference
`ms33_modelIII_gapswitch_ts_890_t_17280000.000000.vtu`, md5 **0e7501f2ff65cb9c151f77da169ecb03**,
11 s wall / 891 steps. The 20x reference generated earlier the same session
(`..._ts_1392_...`, md5 7efcdcbe…) was never committed and was removed, because a second
reference-side `_ts_` file in the test dir makes the ctest wrapper enumerate two candidates.
**MOOT — this specific baseline (100x, `ts_890`) was also never committed.** It was superseded by
the 160x landing (§15, `ts_832`, md5 `e60493e974dee9f9afb2e89ec48ea651`) the same day, before this
100x reference was ever ratified or committed. `ts_832` is now RATIFIED (Vinay, 2026-08-17) —
see §16.

**Reported Model III now, at the 200 d TEAM horizon (no horizon caveat needed any more):**

| | frozen-K | live-K (ADOPTED) | teams |
|---|---|---|---|
| III gap-switch @ 100x, 200 d | 8.023 / 8.009 / 8.106 | **2.162 / 2.170 / 2.150** | 2.63 - 7.47 |

Brackets the band; residual live-side deficit 17 %. The conclusion of §6 is unchanged — only now it
is a converged number at the horizon the teams actually report at.

### 9. OPEN — the suite now carries HETEROGENEOUS permeability factors. ASK(Vinay)

III runs at 100x KC base while IV and VII stay at 20x. The justification is that IV and VII are
already equilibrated at 200 d at 20x (measured: IV within 0.6 %), so they need no more, and that the
endpoint is k0-independent so a per-model factor does not bias the comparison. I measured that
second premise instead of assuming it — SAME DECK, only k0 changed, 20x vs 100x:

| model | quantity | 20x | 100x | change |
|---|---|---|---|---|
| III | T/C/B [MPa] | (2.180/2.169/2.087 @2000 d) | 2.162/2.170/2.150 @200 d | +0.8 % on T/C |
| IV | T/C [MPa] | 6.413 / 6.420 | 6.552 / 6.569 | **+2.2 % / +2.3 %** |
| IV | B (pellet) [MPa] | 1.121 | 1.047 | **-6.6 %** |
| VII | void ratio e | 1.4734 | 1.4956 | **+1.5 %** |

So the endpoint is only APPROXIMATELY k0-independent: within ~2 % on block / single-material
probes, but the **Model IV pellet probe moves 6.6 %**, which is above the precision the deliverable
quotes. The suite is therefore not strictly factor-invariant, and mixing 100x (III) with 20x
(IV, VII) imports that much inconsistency. Two options, Vinay's call:
  (a) leave IV/VII at 20x — minimal change, but the suite is internally inconsistent by up to 6.6 %
      on one probe, and that must be stated wherever the IV pellet number appears;
  (b) unify the suite at 100x — internally consistent and removes the 40x-minimum concern for the
      confined Reference too, but it changes IV by +2.2/+2.3/-6.6 % and VII's e by +1.5 %, and needs
      both ctest references re-baselined (§3(f)) and the IV/VII deliverable numbers re-issued.
Note (b) would also move VII's e from 1.4734 to 1.4956, both far ABOVE the 1.086-1.321 team band —
unrelated to permeability, and consistent with the separate finding that E = 52 MPa is ~10.5x the
team BBM tangent (free-swelling e scales as 1/E).

### 10. Equilibration times and permeabilities, measured per model (2026-08-17)

Uniform criterion, stated up front: **t_eq(1 %) = the earliest output time from which the diagnostic
stays within 1 % of its converged value for the rest of the run.** Diagnostics: Model I = domain-mean
Ps = -tr(sigma)/3; Models III/IV = the canonical Top probe (nearest node, x=0, y=0.070 m); Model VII =
void ratio e from mean porosity. All runs on the canonical tip binary with a 20-point logarithmic
output ladder from 1 d to 2000 d. Model VII is a HYDRATION-ONLY diagnostic (traction held at -2e5,
the 200-240 d load/unload cycle removed) so the hydration transient is measured, not the load step.

| model (shipped deck) | k0 [m^2] | x KC base / x spec literal | diag | converged | **t_eq(1 %)** | @200 d vs converged |
|---|---|---|---|---|---|---|
| I dd1400 | 1.2264e-20 (KozenyCarman) | 1x at own phi0 | Ps | 4.9378 MPa | **20 d** † | +0.00 % |
| I dd1600 | 5.8703e-21 (KozenyCarman) | 1x / 1.048x | Ps | 14.1250 MPa | **20 d** † | +0.00 % |
| I dd1800 | 2.6570e-21 (KozenyCarman) | 1x at own phi0 | Ps | 40.8448 MPa | **1 d** | +0.00 % |
| III gap-switch @ 20x (OLD) | 1.17406e-19 | 20x / 20.97x | T | 2.1796 MPa | **1500 d** | **-46.0 %** |
| III gap-switch @ 100x (SHIPS) | 5.8703e-19 | 100x / 104.83x | T | 2.1949 MPa | **250 d** | **-1.5 %** |
| IV pellets (SHIPS) | 1.17406e-19 clay, 1.17406e-18 pellet | 20x / 200x KC | T | 6.4164 MPa | **150 d** | -0.05 % |
| VII free swelling (SHIPS) | 1.17406e-19 | 20x / 20.97x | e | 1.4992 | **300 d** | **-2.0 %** |

† **Model I's t_eq is NOT a hydraulic time.** Its suction ramp ends at exactly 1728000 s = 20 d, and
pressure Dirichlet sits on top AND bottom of a single-element mesh, i.e. on all four nodes — there is
no free hydraulic degree of freedom. Model I equilibrates the moment the BC ramp stops, and its
permeability is irrelevant to that (it is also the only model never accelerated). dd1800 reaches 1 %
even earlier because it is densest, so its macro pore is exhausted soonest.

**Who is actually equilibrated at the 200 d team horizon:** Model I (trivially, by construction) and
Model IV (-0.05 %). **Model III at the newly shipped 100x is -1.5 %** — close, but t_eq(1 %) = 250 d
still lands just past the horizon. **Model VII is -2.0 % and is NOT equilibrated at 200 d either**,
t_eq(1 %) = 300 d — which independently corroborates the 2026-06-17 record
(`runs/2026-06-15_1410_.../EQUILIBRATED_LE_VII_2026-06-17.md:9-15`) and the 2026-08-17 equilibration
audit. VII's shortfall was never fixed and is still shipping.

**CORRECTION to §8.** I reported the 100x Model III 200 d value as "within 0.8 % of the asymptote".
That compared 100x@200 d against 20x@2000 d. Measured against its OWN converged value on a fine
ladder it is **-1.5 %**. Both numbers are real but they answer different questions; -1.5 % is the one
to quote.

**Scaling.** t_eq is close to inverse-proportional in k0 but slightly steeper: Model III, same deck,
5x more k0 moved t_eq 1500 -> 250 d, a factor **6.0**, i.e. t_eq ~ k0^-1.11. The excess over 1/k0 is
expected — the gap-closure switch and the exchange kinetics gate the transient too, so it is not a
pure diffusion time. Using that exponent, t_eq(1 %) = 200 d exactly would need **~122x** KC base for
Model III and **~29x** for Model VII. [PREDICTED from two points, not verified by a run.]

**Side observation, consistent with the known §6.7 gap:** the registered Model I decks give
Ps = 4.9378 / 14.1250 / 40.8448 MPa, not the deliverable record's 4.89997 / 14.159948 / 40.600313 —
because the registered decks are still on the superseded sigma0 = -1.5e5 chain (K = 45217 / 103879 /
266767). That is the open blocker already recorded, not a new finding.

### 11. HISTORICAL — SUPERSEDED BY §12. Full permeability sweep, target t_eq(1%) = 20 d [WRONG TARGET]

> **Correction 2026-08-17 (same day): the "20 d" target below is WRONG.** Checked against the
> spec directly (Vinay: *"find in the specs how long should the models run for? 20 or 200?"*):
> `d_calculation-description/BGR.tex:60` — BGR's own declared MS33 setup — states *"A time-dependent
> suction ramp will be imposed ... decreasing ... over a period of 20 days, after which it remains
> constant for a further 180 days. Thus, the total hydration and simulation time is **200 days**."*
> 20 d is the RAMP sub-component, not the reporting horizon; 200 d is the target, confirmed
> independently by `Questions.tex` R27 (*"200 days for the Reference Configuration"*) and by the
> ALREADY-STANDING permeability-multiplier discussions in §7/§8, both of which were about 200 d
> from the start (commit `2c3b3bb0d6`, 2026-06-09: *"unconverged transients at the 200-240 d
> horizon"*; the §8 100x landing: *"compresses the transient into the 200 d team reporting
> horizon"*). The "20 days" framing entered ONLY in this section, as a one-turn misreading, not as
> anything previously agreed. §12 below re-derives the SAME dataset — no new simulations were
> needed, since the run below already went to 2000 d with a ladder point at exactly t=200 d.
> This section is kept per CLAUDE.md §6.3 (never delete, annotate historical) — the raw sweep data,
> the Aitken/Richardson method, and the IV-overshoot / VII-and-Model-I-ramp-floor physics findings
> below are all still valid; only the "closest to 20 d" framing and its answer table are superseded.

Vinay, 2026-08-17: *"run everything except the calibration and 30x, 40x, 50x, 60x. do a richardson
(if i am right) for the equilibration. take the closest which equilibrates at 20 days."*

**Scope.** No K recalibration; K unchanged from the live-K Dixon chain throughout. Models III
(gap-switch), IV (pellets), VII (free-swelling, hydration-only diagnostic — traction held, the
200-240 d load/unload cycle stripped so the hydration transient alone is measured). Model I excluded
from the sweep: it has zero free hydraulic DOF (pressure Dirichlet on every node of a 1-element
mesh), so permeability provably cannot affect it — confirmed as the x1 reference (t_eq=20 d exactly,
BC-ramp-controlled, unchanged from §10). 25 new runs at {5,10,80,100,160,320,640,1280}x KC base
(+2560,5120x for III only, added to resolve its floor) — literal multipliers 30/40/50/60 excluded
throughout. All 25 completed cleanly (exit 0, no divergence) on the canonical tip binary.

**Extrapolation method.** Classical Richardson extrapolation needs a constant ratio between
successive step sizes; the output ladder here is only approximately geometric, so the applied method
is Aitken's Delta-squared process — same family, and the correct generalisation when the ratio isn't
exactly constant: `x_inf ~= x_n - (x_{n+1}-x_n)^2 / (x_{n+2}-2 x_{n+1}+x_n)`, evaluated on sliding
triples of the late-time tail, x_inf reported as the median across triples. t_eq(1%) is then read off
against this extrapolated x_inf rather than the raw last-frame value. In every case except the two
slowest-converging runs (III/VII at 5x, 10x, which have not entered the terminal regime by 2000 d)
the Richardson and naive t_eq agree to within one ladder step — the naive last-value estimate was
already adequate once the run had genuinely converged; Richardson matters only for judging whether
convergence had ACTUALLY been reached (III/VII at 5x/10x had not, despite 2000 d).

**Full table — k0 = KC-base x multiplier, KC base = 5.8703e-21 m^2:**

| model | mult | k0 [m^2] | x_inf (Richardson) | t_eq(1%, Richardson) |
|---|---|---|---|---|
| III | 5x | 2.935e-20 | 0.969 | >2000 d (not converged) |
| III | 10x | 5.870e-20 | 1.643 | >2000 d (not converged) |
| III | 20x (old ship) | 1.174e-19 | 2.096 | >2000 d (naive said 1500 d — was wrong) |
| III | 80x | 4.696e-19 | 2.196 | 300 d |
| III | 100x (SHIPS) | 5.870e-19 | 2.195 | 250 d |
| III | 160x | 9.392e-19 | 2.192 | 150 d |
| III | 320x | 1.878e-18 | 2.187 | 75 d |
| III | 640x | 3.757e-18 | 2.183 | 50 d |
| **III** | **1280x** | **7.514e-18** | **2.182** | **30 d — floor** |
| III | 2560x | 1.503e-17 | 2.180 | 30 d (no further gain) |
| III | 5120x | 3.006e-17 | 2.181 | 30 d (no further gain) |
| IV | 5x | 2.935e-20 | 6.365 | 600 d |
| IV | 10x | 5.870e-20 | 6.391 | 250 d |
| IV | 20x (SHIPS) | 1.174e-19 | 6.416 | 150 d |
| **IV** | **80x** | **4.696e-19** | **6.522** | **30 d — floor** |
| IV | 100x | 5.870e-19 | 6.553 | 30 d |
| IV | 160x | 9.392e-19 | 6.613 | 30 d |
| IV | 320x | 1.878e-18 | 6.653 | 50 d — WORSE |
| IV | 640x | 3.757e-18 | 6.568 | 50 d |
| IV | 1280x | 7.514e-18 | 6.469 | 50 d |
| VII | 5x | 2.935e-20 | 1.467 | >2000 d (not converged) |
| VII | 10x | 5.870e-20 | 1.495 | 600 d |
| VII | 80x | 4.696e-19 | 1.500 | 75 d |
| VII | 100x | 5.870e-19 | 1.500 | 75 d |
| VII | 160x | 9.392e-19 | 1.500 | 50 d |
| VII | 320x | 1.878e-18 | 1.500 | 30 d |
| **VII** | **640x** | **3.757e-18** | **1.500** | **20 d — matches Model I's floor exactly** |
| VII | 1280x | 7.514e-18 | 1.500 | 20 d (no further gain) |

**The closest-to-20-d run per model, honouring the request:**

| model | closest run | k0 [m^2] | t_eq(R) |
|---|---|---|---|
| III | **1280x** | 7.514e-18 | 30 d |
| IV | **80x** | 4.696e-19 | 30 d |
| VII | **640x** | 3.757e-18 | **20 d** |

**VII reaches the literal 20 d floor** — bit-identically after 640x, no further gain at 1280x. **III
plateaus at 30 d** and pushing two more decades (2560x, 5120x) buys nothing further — confirmed a
genuine floor, not a resolution artefact of the sweep stopping too early. **IV never gets below
30 d, and going past ~160x makes it WORSE (30 -> 50 d).**

**Why the floor sits at 20 d for VII/I but 30 d for III, and why IV cannot improve past 30 d
[DERIVED, consistent with the measured raw series — checked, not assumed].**

- **VII and Model I** are driven purely by the imposed suction ramp, `<coords>0 1728000 ...</coords>`
  — the BC itself finishes changing at exactly t = 1728000 s = **20 d**. A model cannot equilibrate
  to a state it does not yet know, so 20 d is a hard floor set by the forcing, not by diffusion.
  Once k0 is high enough that hydraulic diffusion is no longer rate-limiting (>=640x here), the
  system tracks the ramp essentially in lock-step and hits that floor exactly — confirmed: VII's raw
  series at 1280x is bit-identical from t=20 d onward (e=1.50035 at every frame >=20 d).
- **III** carries the SAME suction ramp but ALSO the gap-switch mechanical event (free swelling until
  u_r=2 mm, then a rigid-wall BC engages). That switch is a second, k0-independent clock, and it sets
  III's floor 10 d later than VII/I's. Checked directly: raw T at 5120x is already 2.155 (99 % of
  converged) at t=20 d and only needs the last 1 % by t=30 d — consistent with a mechanical-event
  floor stacked on top of the hydraulic one, not a resolution artefact.
- **IV is qualitatively different: pushing k0 higher makes it WORSE past ~160x, because it develops a
  genuine transient OVERSHOOT that grows with k0** — checked in the raw series, not a ladder/noise
  artefact:
    | mult | T at t=20d (peak) | T at t=100d (settled) | overshoot |
    |---|---|---|---|
    | 160x | 6.476 | 6.615 | undershoot, still rising |
    | 320x | 6.659 | 6.655 | ~flat |
    | 640x | 6.747 | 6.571 | **+2.7% peak, relaxes down** |
    | 1280x | 6.800 | 6.474 | **+5.0% peak, relaxes down** |
  IV is the only heterogeneous model here (clay block + pellet, each with its own K and k0). At high
  k0 the fast-hydrating block races ahead of the pellet and overshoots the FINAL two-material
  equilibrium before the pellet catches up and pulls it back down; the relaxation time for that
  overshoot does not shrink with k0, so t_eq stops improving and then degrades. This is a genuine
  finding about the coupled system, not a numerical defect — recorded as a caution against pushing
  k0 arbitrarily high on heterogeneous DSM decks.

**Consequence for the shipped decks.** None changed by this campaign — it is a diagnostic sweep, no
tracked file touched. It does, however, sharpen the k0=100x decision for Model III from §8/§9: 100x
gives t_eq=250 d (still past the horizon, though only -1.5 % off by 200 d as measured in §10); the
run that actually reaches Model III's own floor is **1280x** (t_eq=30 d), 12.8x higher than what
shipped. Whether to push III further than 100x is a separate decision from the one already made —
not applied here, ASK(Vinay) if wanted.

Raw data: 27 new PRJs + full 2000 d VTU series under scratchpad `eqtime/{III,IV,VII}_{mult}x/`;
`analysis.json` has the per-case Aitken x_inf, spread, and both t_eq definitions.

### 12. CORRECTED — the same sweep, re-queried for the actual 200 d target

Vinay: *"then delete the just-now-run campaign and redo it for 200 days, maybe i made a mistake
hastily saying 20. first check if there was a discussion about the const perm multiplier
regarding this."*

**Checked first, as asked.** Both prior permeability-multiplier discussions were about 200 d:
the original 2026-06-09 commit `2c3b3bb0d6` (*"unconverged transients at the 200-240 d horizon"*)
and this session's §8 100x landing (*"compresses the transient into the 200 d team reporting
horizon"*, `ms33_modelIII_gapswitch.prj:99`). The 20 d framing in §11 was the anomaly.

**No re-simulation performed, and none was needed.** t_eq(k0) is an intrinsic property of each
multiplier's own physics — it does not depend on which target (20 d or 200 d) you later query it
against. Every §11 run already went to 2000 d on a ladder that includes an exact t=200 d output
frame. Re-running would reproduce byte-identical PRJs, byte-identical physics, and therefore
byte-identical results — pure waste. Instead the SAME `analysis.json` was re-queried for the
correct target. One repair was needed and made: Model VII's shipped 20x baseline was cached under
the directory name `VII_20x_hydonly` (from an earlier turn, before the sweep's naming convention),
which the sweep's `VII_{m}x` regex could not see — it was silently absent from every query in §11
too. Added as `VII_20x` in `analysis.json`, t_eq(R) = 300 d, consistent with the independent
§10 finding (VII at 20x is -2.0% short at 200 d) — corroborating, not new.

**Query rule:** among runs with t_eq(Richardson) <= 200 d, pick the one with the smallest margin
(converged earliest before the horizon, not merely "close" in absolute distance — a run whose
t_eq is 250 d is NOT converged at the 200 d horizon regardless of how "close" 250 is to 200; a run
whose t_eq is 150 d gives a 50 d safety margin and IS converged by the horizon). Where no run in
the swept set converges by 200 d, the least-late run is reported instead, flagged as such.

| model | chosen | k0 [m^2] | t_eq(R) | margin vs 200 d | note |
|---|---|---|---|---|---|
| III | **160x** | 9.3925e-19 | 150 d | +50 d | currently SHIPS at 100x (t_eq=250d, -50d, NOT converged by 200d) |
| IV | **20x** | 1.1741e-19 | 150 d | +50 d | this IS what currently ships — no change indicated |
| VII | **80x** | 4.6962e-19 | 75 d | +125 d | currently SHIPS at 20x (t_eq=300d, -100d, NOT converged by 200d) |

**Full table for context (t_eq(Richardson), all swept multipliers):**

| mult | III | IV | VII |
|---|---|---|---|
| 5x | >2000d | 600d | >2000d |
| 10x | >2000d | 250d | 600d |
| 20x | >2000d | 150d | **300d (shipped, repaired entry)** |
| 80x | 300d | 30d | 75d |
| 100x | **250d (shipped)** | 30d | 75d |
| 160x | **150d** | 30d | 50d |
| 320x | 75d | 50d | 30d |
| 640x | 50d | 50d | **20d** |
| 1280x | 30d | 50d | 20d |

**Three findings, in order of what they mean for what's already landed:**

1. **IV needs nothing — confirmed, not just assumed.** The shipped 20x already gives the best
   available margin (150 d, converged 50 days before the horizon) among all sensible candidates;
   every acceleration beyond 20x makes it WORSE against the 200 d target (t_eq bottoms at 30d for
   80-160x, driven by the block-pellet overshoot mechanism in §11, but that 30d floor is for a
   DIFFERENT target than 200d and is irrelevant here — at the 200d target, 20x's 150d margin beats
   every accelerated option, since none of them get CLOSER to 200 from below than 20x already is;
   80-160x overshoot INTO their floor early and then sit flat, still well inside 200d, so they are
   not wrong either, just not better). No change indicated.

2. **III at 100x (already landed, per §8) does NOT converge by the 200 d horizon** — t_eq=250d,
   -1.5% short at 200 d exactly as measured independently in §10. **160x would give a 50-day safety
   margin instead.** The two x_inf values are nearly identical (100x: 2.1949 MPa; 160x: 2.1918 MPa,
   0.14% apart) — this is NOT a physics disagreement, purely a convergence-margin question. Whether
   to bump III from 100x to 160x is Vinay's call — the 100x landing itself is UNCHANGED here, this
   is new information bearing on it, not a reversal.

3. **Model VII, never touched by any landing this session, does not converge by 200 d either**
   (t_eq=300d at its shipped 20x) — an existing, previously-known gap (§10; corroborates the
   2026-06-17 record), now with a concrete candidate fix: 80x gives a 125-day margin. This is a
   NEW decision, not previously proposed or landed for VII specifically.

**Nothing committed by this correction.** No tracked file was touched — §12 is a re-analysis of
already-collected data. The open decisions (bump III to 160x? raise VII to 80x?) are ASK(Vinay),
same footing as the §9 heterogeneous-factor question they now overlap with.

### 13. WHY III needs more permeability than VII to converge — verified mechanism

Vinay: *"there must be a reason why 3 takes more permeability to close a 2mm gap and vii takes
less under absolute uniaxial free swelling, if the end point of those two are the same."*

**Checked first: the two decks are a controlled comparison.** Mesh (`ms33_cylinder_r25_h70`,
identical), hydration BC (pressure_bc, bottom face only, identical suction ramp), retention curve
(SaturationTuller, every parameter byte-identical), E=52 MPa, nu=0.3, biot=1.0, phi0=0.42446,
initial rho_d=1600 — all identical between `ms33_modelIII_gapswitch.prj` and
`ms33_modelVII_freeswelling.prj`. The ONLY structural difference is mechanical: III's outer radius
is free until u_r=2mm then LOCKS to a rigid Dirichlet wall, and III's top+bottom are ALWAYS
u_z=0 (fixed, no axial swelling ever); VII's outer radius carries no BC at all (free/zero-traction
throughout) and its top is a soft Neumann traction (movable). So the shared "endpoint" the question
refers to is real — same boundary driver, same material — but it is the HYDRAULIC boundary
condition, not the mechanical one; the mechanical response is where III and VII diverge.

**Measured (same k0=100x KC base for both, so the comparison isolates the mechanical effect):**

| t [d] | III eps_v | III phi_M (macro) | III exch_mean | VII eps_v | VII phi_M (macro) | VII exch_mean |
|---|---|---|---|---|---|---|
| 1 | 0.0074 | 0.4260 | 2.9e-05 | 0.0142 | 0.4268 | 8.5e-05 |
| 20 | 0.1331 | 0.3272 | 3.1e-04 | 0.2993 | 0.2806 | 5.4e-04 |
| 100 | 0.1599 (capped) | 0.0805 | 2.76e-04 | 0.3631 (still free) | 0.0229 | 5.15e-04 |
| 200 | 0.1600 (capped) | **0.0168 (not exhausted)** | 2.77e-04 | 0.3637 (capped) | **0.0000 (exhausted)** | 5.15e-04 |

III's dilation is capped at eps_v=0.16 by t=30 d (gap fully closed) and never moves again. VII's
grows to eps_v=0.36 (2.3x III's) and is STILL free the whole way — even at 100x, VII's own
"gap" (its confining traction) never engages the way III's rigid wall does, so nothing analogous
caps VII. Note this is ALREADY true at t=1 d, well before III's gap engages (eps_v 0.0074 vs
0.0142) — III's top+bottom u_z=0 (permanent, not switch-gated) restricts it to radial-only
dilation from the very first step, while VII is free in both directions throughout.

**The causal link, from source (verified live in every MS33 run, §5 of `paper_code_formulation
_divergence.md`): the DSM micro-macro exchange carries an ADDITIVE mechanical term in the micro
potential,**

    mu_lR_mech = -[ (Pi + n_l*Pi')*eps_v + 0.5*b*K_drained*eps_v^2 ] / rho_lR

**the same term responsible for the live-K "swelling material loses swelling potential" finding.**
Larger eps_v drives mu_lR more negative, which widens the exchange driving force
`rho_hat_l = alpha_M*(mu_LR - mu_lR)` and speeds the transfer of macro water into micro/interlayer
storage. Measured: VII's exch_mean during the late quasi-steady phase is ~1.9x III's — and this is
despite VII's exhaustion TARGET being LARGER, not smaller (final n_l = 0.600 for VII vs 0.510 for
III, since VII's larger dilation also increases its total porosity) — VII moves more water, faster,
to a bigger target, and still finishes first. The rate effect dominates the target-size effect.

**Confirmed spatially, not just in the domain mean.** III's outer-boundary displacement at t=20 d,
sampled by height (`III_100x`, x=0.025 m = outer radius):

| y [mm] | u_r [mm] (gap=2.0) | local phi_M |
|---|---|---|
| 0 (bottom, entry) | **2.000 — CLOSED** | ~0 (exhausted) |
| 20 | 1.772 | 0.347 |
| 40 | 1.552 | 0.386 |
| 70 (top) | 1.460 | 0.398 (barely moved from 0.424 initial) |

The gap closes bottom-first, top-last, tracking the rising hydration front — each cross-section
runs its own local "free-swell then lock" sequence as the front reaches it. The domain mean cannot
converge until the LAST (top) cross-section both hydrates AND locks. VII has no analogous
per-height lock-in event (uniformly free everywhere), so its front-propagation alone sets its
timescale, with no compounding mechanical-confinement delay layered on top.

**Answer, in one line:** III needs more permeability than VII for the same 200 d target because
III is progressively confined (a rigid contact wall) and VII is not — confinement caps III's
dilation at ~16% (vs VII's ~36%+), and because the DSM's micro-macro exchange rate is itself
strain-gated (the same coupling behind the live-K fabric-memory finding), the smaller dilation
gives III a structurally slower self-driven exchange, on top of a spatial lock-in effect the
free-swelling geometry doesn't have. Both effects point the same way and share the same root
cause. Not a numerical artefact and not about diffusion path length (identical for both decks) —
a genuine, verified consequence of the confinement difference between the two exercises.

### 14. Two questions raised, both checked directly — one confirmed, one genuinely open

Vinay: *"firstly everything has to be const perm. vii should be laterally confined, no?"*

**1. Const perm — checked, mostly already true.** Model III, IV (both zones) and VII all use
`<type>Constant</type>` for the `permeability` property — confirmed by re-reading all four PRJs
directly. **Model I is the one exception**: `ms33_modelI_dd1600.prj:179` uses
`<type>KozenyCarman</type>` (deformation-coupled, `k = k(phi)`), the only deformation-dependent
permeability law anywhere in the MS33 suite. Possible reason, not yet confirmed as the actual
motive: the spec's Model I section explicitly lists *"Evolution of intrinsic permeability vs
suction"* as a comparable result (`theoretical_benchmarking.tex:198`) — under Constant
permeability that curve is a flat line, giving nothing to compare; KozenyCarman is the only way to
produce a non-trivial one. **ASK(Vinay):** keep Model I on KozenyCarman for that reason, or switch
it to Constant too for suite-wide consistency? Not changed here — a formulation decision.

**2. VII lateral confinement — checked, genuinely open, teams disagree.** `ms33_modelVII_
freeswelling.prj:330` currently reads *"Right face (r=0.025): FREE (no BC -> natural Neumann =
zero traction)"* — deliberately unconfined, confirmed by reading the full, unfiltered
`<boundary_conditions>` block (not a parsing miss). Checked against the spec and against every
team's own calculation-description for the exercise:

  - Spec text (`theoretical_benchmarking.tex:301-309`): *"buffer blocks placed in an **unconfined
    space** may undergo free swelling... under a constant stress of 0.2 MPa"* — leans toward
    unconfined but never uses the word "radial" or "lateral" at all.
  - **CIMNE-UPC** (task-lead team, source of the Model VII figure): frames a friction-plus-VII
    combination as an explicit ADDITIONAL exercise (*"can be regarded as a combination of Model
    II [friction] and Model VII"*, `CIMNE-UPC.tex:990`) and states *"the interaction between the
    sample and the rigid boundaries partially restrains the radial deformation"* only in that
    ADDITIONAL model — implying the BASE Model VII, as CIMNE built it, has no radial wall.
  - **UCLM** states the opposite, explicitly, for their own base Model VII: *"Radial confinement
    is also assumed"* (`UCLM.tex:558`).
  - **EPFL** ran VII as a single element (sidesteps the radial question by construction) and
    reports a final void ratio e=1.47 from e0=0.74 (`EPFL.tex:351-353`) — close to our OWN
    unconfined result (e~1.50 at full convergence, this session's sweep) and to the e0=0.74
    already on record for our free-swelling PRJ (memory `feedback_ms33_modelVII_free_swelling.md`).
  - **ClayTech** independently corroborates our own equilibration-time struggle: their base
    200 d timescale did not reach full saturation before loading either, and they had to extend
    to 800-2000+ days to get a clean result (`ClayTech.tex:250-260`) — same failure mode as our
    Model III/VII permeability-acceleration finding (§§7-13), on a completely different code.

**Not changed.** Whether to add a radial confinement BC to `ms33_modelVII_freeswelling.prj` is a
genuine, team-split modelling choice, not a spec-clear answer — per §9 (physics/formulation
decision), this is ASK(Vinay), not something to decide unilaterally. If confinement is added,
every VII number produced this session (e=1.2183 old record, 1.4956-1.5004 across the sweep, the
"mid-field" team-band claim, the deck/paper mentions) would need re-running under the new BC and
would very likely change materially — this is a bigger decision than the permeability-multiplier
questions already open.

### 15. DONE 2026-08-17 — full suite redone: const perm everywhere, all equilibrated, no radial confinement on VII

Vinay: *"ok, no radial confinement. all models including calibration should be constant
permeability. redo the entire suite and ensure all models are equilibrated. report to me."*

**Landed:**

1. **VII stays unconfined** — no radial BC added. Confirmed decision from §14 (cross-team
   evidence was split; kept as-is per Vinay's explicit call).
2. **All models, including calibration, now use `<type>Constant</type>` permeability.**
   Model I (`ms33_modelI_dd{900,1400,1600,1800}.prj`) was the one exception (KozenyCarman) —
   converted. **Verified empirically, not assumed:** a controlled A/B run (dd1600, everything
   else identical, KozenyCarman vs Constant) gave **0.000e+00 max abs diff on every field**
   (pressure, sigma, saturation, porosity, micro_water_content, transport_porosity,
   swelling_stress) — confirms Model I's single-element, all-Dirichlet mesh has no free
   hydraulic DOF, so the permeability LAW never gets exercised regardless of type.
   `calibrate_maxwell_K.py`'s stale June binary pointer fixed to the canonical tip
   (`verify_tip_20260812`).
3. **Two pure-bug fixes made in passing** (not physics/calibration changes): dd900's
   `<prefix>ms33_modelI_dd1600</prefix>` → `dd900` (was silently overwriting dd1600 output);
   the VII deck's §12.5 header carried Model IV's "clay/pellet" boilerplate verbatim, copied
   in error and never describing this deck — replaced with a VII-specific header.
4. **III bumped 100x → 160x** (9.392480e-19 m²): the 100x landing (§8) gave
   t_eq(1%,Richardson)=250 d, NOT converged by the 200 d horizon (-1.5% short, §10/§12).
   160x gives t_eq=150 d, a 50-day margin. x_inf moves only 0.14% between the two (2.1949 →
   2.1918 MPa on Top) — a convergence-margin choice, not a different physics answer.
5. **VII bumped 20x → 80x** (4.696240e-19 m²): 20x gave t_eq=300 d (-2.0% short, corroborating
   the 2026-06-17 record that this was never fixed). 80x gives t_eq=75 d, 125-day margin.
6. **IV unchanged** (already Constant, already at its best-margin multiplier per §12 — bumping
   further makes it worse per the block-pellet overshoot mechanism in §11).

**Verification, in order:**
- Model I A/B (KozenyCarman vs Constant, dd1600): 0.000e+00 diff, confirms permeability-law
  independence (item 2 above).
- Full four-cell Model I run (Constant perm, current K/σ0 as tracked, unchanged): Ps =
  8.464583 / 4.921826 / 14.160250 / 40.844810 MPa (dd900/1400/1600/1800). dd1600 matches its
  pre-edit KozenyCarman value bit-for-bit (14.160250 both ways).
- Model III final (160x, live-K, 200 d): T/C/B = **2.192 / 2.184 / 2.170 MPa**.
- Model IV final (unchanged, 200 d): T/C/B = **6.413 / 6.420 / 1.121 MPa** — matches the §12
  sweep's IV_20x same-deck reference exactly.
- Model VII final (80x, FULL production deck incl. 200-240 d load/unload cycle): void ratio
  at t=200 d (pre-loading, hydration endpoint) = **1.4996**, matching the hydration-only sweep
  prediction (VII_80x x_inf=1.4997) to **0.01%** — direct cross-check between the diagnostic
  sweep and the real production deck. Post-loading at t=240 d: e=1.4956.
- Existing tracked references for Model I (all three registered cells) and Model IV were
  diffed against these fresh runs: max abs diff 2e-8 to 6e-7 across all fields — floating-point
  noise, well inside the registered tolerances (1e-9 to 1e-2 depending on field). **Left
  untouched** — no value changed for these two, no reason to re-baseline.
- New references generated and installed for III (`..._ts_832_t_17280000...vtu`) and VII
  (`..._ts_510_t_20736000...vtu`), the two decks whose physics actually changed.
- `ctest -R ANCHORS_MS33` on the canonical tip: **6/6 pass.** Verify worktree restored clean
  afterward.

**NOT touched, flagged instead (outside the scope of what was asked):**
- **dd900's K=103879 J/kg is not a calibrated value for dd900** — it is copy-pasted from
  dd1600's old K. Ps=8.464583 MPa under this K is not a meaningful physics result. dd900 is
  not a registered ctest, so this doesn't gate anything, but the deck is not usable as-is for
  any dd900 comparison.
- **All four registered Model I decks remain on the OLD σ0=-1.5e5 / superseded-K chain**
  (K=103879/45217/103879/266767.3), NOT the ratified σ0=0 Option-C chain
  (K=4367.2277/46000/104689.9129/265905.06) used everywhere else this session. This is the
  pre-existing §6.7 blocker B6, untouched here — it is a calibration/BC-magnitude decision
  (§1.1/§9), not a permeability question, and was not part of what was asked this turn.

**RATIFIED (Vinay, 2026-08-17)** — see §16: K(900)=4367.2277 for the III/IV live-K chain (§1.1);
the III and VII reference VTU checksums (§3(f)) —
`ms33_modelIII_gapswitch_ts_832_t_17280000.000000.vtu` (md5 `e60493e974dee9f9afb2e89ec48ea651`) and
`ms33_modelVII_freeswelling_ts_510_t_20736000.000000.vtu` (md5 `3af2f07d0569631cce5a2507cc3c87a4`).

### 16. RATIFIED 2026-08-17 (Vinay): K(900) and the two committed ctest baselines

Vinay: *"ratify K(900) and the two ctest baselines, show me the comparison plots in the artefact."*

**K(900) = 4367.2277 J/kg — RATIFIED.** §1.1 calibration-anchor literal. Fitted to the exact Dixon
(2023) median-fit extrapolation `Ps = 0.003*exp(5.2883*0.9)` = 0.3500522009 MPa (achieved
0.3500565 MPa, +0.001%), reproduced to 1e-10 relative on two independent binaries. Live in
`ANCHORS_MS33_ModelIII/ms33_modelIII_gapswitch.prj`'s K(rho_d) table (the 900-1400 segment is
actively read: III's live rho_d falls to 1347-1376 kg/m3). NOT currently used by Model IV's
pellet zone, which still runs the separate uncited K=20600 (open question U-1, unresolved,
untouched by this ratification).

**Two ctest reference baselines — RATIFIED, §3(f)/§12.5:**

| model | file | md5 | steps | verified |
|---|---|---|---|---|
| III | `ms33_modelIII_gapswitch_ts_832_t_17280000.000000.vtu` | `e60493e974dee9f9afb2e89ec48ea651` | 832 | `ctest -R ANCHORS_MS33` 6/6 pass |
| VII | `ms33_modelVII_freeswelling_ts_510_t_20736000.000000.vtu` | `3af2f07d0569631cce5a2507cc3c87a4` | 510 | `ctest -R ANCHORS_MS33` 6/6 pass |

Both already committed (`55e39b53b8`) at the time of ratification — no new commit required for
the reference files themselves. Two EARLIER, unratified candidate references for III existed at
different points this session (20x era `ts_1392`, md5 `7efcdcbe1f0640ea3e6e9aaa6109ed45`; 100x era
`ts_890`, md5 `0e7501f2ff65cb9c151f77da169ecb03`) and are now MOOT — both were superseded by further
permeability corrections before ever being committed or ratified, so there is nothing to walk back;
`ts_832` is the only III reference that ever reached git, and it is the one now ratified.

**Comparison plots added to the team-comparison artefact** (`ms33_team_comparison.html`,
republished 2026-08-17): a K(900) calibration-convergence panel (target vs. achieved Ps across the
secant iteration) and equilibration time-series panels for the two ratified baselines (III at
160x to t=200d, VII hydration at 80x to t=200d), regenerated fresh with a full output ladder
specifically for this purpose — the exploratory sweep's raw time series had already been deleted
per the earlier cleanup instruction, so these are new, dedicated verification runs, not reused
scratch data.

### 17. Found while building the ratification plots: dd900's tracked σ0 doesn't match its own ratified K

While regenerating the K(900) calibration plot for §16, a fresh calibration run on the CURRENT
tracked `ms33_modelI_dd900.prj` converged to **K=4904.66**, not the ratified 4367.2277 — a 12.3%
discrepancy. Root cause, checked directly, not guessed: the tracked dd900.prj carries
**σ0 = -1.5e5**, not σ0 = 0 (Option C). K(900) = 4367.2277 was calibrated (and re-verified for §16)
under σ0 = 0, matching the convention already used for the ratified dd1400/1600/1800 K values
(46000/104689.9129/265905.06). Confirmed with a controlled A/B (KozenyCarman vs Constant at fixed
K=4367.2277, dd900, everything else identical): **0.000e+00 diff, Ps=0.456805 MPa both ways** — not
the target 0.3501 MPa. Confirmed with σ0 corrected to 0 in a scratch copy: the calibration
reproduces K=4367.227700091862, matching the ratified value to 9 significant figures.

**This is not a new problem — it is the existing §6.7 blocker B6 (all four registered Model I decks
remain on the superseded σ0=-1.5e5 chain, flagged and deliberately untouched in §15) — but building
this plot is the first place this session where it was shown to have a DIRECT, immediate
consequence: running the tracked dd900.prj as-is with the just-ratified K=4367.2277 does NOT
reproduce the ratified target.** K(900)=4367.2277 itself is correctly ratified (it is the value
consistent with the other three K(rho_d) table knots); what's inconsistent is the deck's own σ0.

**Not fixed here** — σ0 is a boundary-condition magnitude (§1.1/§9: ASK USER), and fixing it for
dd900 alone would be inconsistent with the SAME gap in dd1400/1600/1800; a real fix means updating
all four decks together, which is exactly the scope §15 explicitly deferred. **ASK(Vinay):** update
all four registered Model I decks (900/1400/1600/1800) to σ0=0, matching the ratified K(rho_d) chain
they're all built from? This is the smallest remaining step to make the tracked Model I suite
internally consistent with everything else ratified this session.

**DONE 2026-08-18** (Vinay-approved) for the three CTEST-REGISTERED decks: `ms33_modelI_dd1400.prj`,
`ms33_modelI_dd1600.prj`, `ms33_modelI_dd1800.prj` (Tests.cmake:31-33) migrated to σ0=0 + the
ratified Option-C K table (46000.0/104689.9129/265905.06). §12.2 headers rewritten to cite the
migration; old sigma0=-1.5e5-era K/target/achieved-Ps kept as history text in the same block. New
ctest reference VTUs generated on the canonical binary (maxwell-conjugate-20260602, md5 c432a156)
and verified bit-exact (vtkdiff, all 11 declared fields, every abs/rel norm = 0.0) against the
already doubly-cross-validated `sigma0zero_recal_forma_2026-08-06` campaign record — a third
independent reproduction. Old references moved to
`ANCHORS_MS33_ModelI/superseded_references_2026-08-18/` (not deleted, §6).
**dd900 is UNCHANGED and remains open** — it is not ctest-registered (absent from Tests.cmake) and
carries the separate, more severe copy-paste bug documented above (dd1600's K and n_s pasted onto
dd900 in error); fixing it needs the correct dd900 K value from Vinay, not just the σ0 swap applied
here. The "all four decks" framing of the ASK above is accordingly narrowed: three registered decks
resolved, dd900 tracked as its own item.

---

## 2026-08-17 — Model VII switched to live-K (extends the III/IV ruling), re-verified at 160x, and a build/ctest-routing staleness finding

### 18. DONE 2026-08-17 — VII live-K switch, equilibration re-check, reference rebaseline

Vinay's question that triggered this: *"wait, why did you do vii with fixed K?"* — noticed the
team-comparison figure's VII legend still said "fixed K=103879" while III/IV had already gone
live-K. Checked: the committed `ms33_modelVII_freeswelling.prj` did carry a single scalar
`<potential_augmentation_prefactor>103879.0</potential_augmentation_prefactor>`, not the live table.
Root cause: the 2026-08-17 live-K ruling (entry #16 above) was worded "Models III and IV" — I
applied it literally and never separately raised whether it should extend to VII. Flagged per §9
(physics/model-formulation decision) and asked; Vinay approved extending it.

**Physical case for extending to VII:** VII is free swelling with no radial confinement — the
largest density excursion in the whole suite (e climbs ~0.73 to ~1.50 under fixed K). The live-K
argument ("a swelling material loses swelling potential; pricing it at the emplacement density is
wrong") applies at least as strongly here as to III/IV.

**Implementation.** `potential_augmentation_prefactor` (scalar) replaced by
`potential_augmentation_prefactor_vs_dry_density` (900/1400/1600/1800 ->
4367.227700212952/46000.0/104689.9129/265905.06) + `potential_augmentation_prefactor_live_dry_density
= true` — the IDENTICAL table already on `ms33_modelIII_gapswitch.prj`. sigma0 left at -1.5e5
(unchanged; same known convention gap as III/IV, not this decision — see entry #17).

**Equilibration re-check (the permeability margin found under FROZEN K does not carry over to
live-K).** The existing 80x setting was measured under frozen K (entry #15/#7); re-measured after
switching to live-K, three parallel runs {80x, 160x, 320x}, top-probe void ratio at the matched
200 d horizon:

| k0 multiplier | e_top(100d) | e_top(150d) | e_top(200d) |
|---|---|---|---|
| 80x  | 1.10169 | 1.11948 | 1.12759 |
| 160x | 1.12816 | 1.13548 | 1.13722 |
| 320x | 1.13705 | 1.13705 | 1.13705 |

320x plateaus EXACTLY flat from 100 d (true converged asymptote). 80x sits 0.83% short of that
asymptote AT the 200 d horizon itself — no safety margin, fails the same standard that rejected
III's 100x. 160x matches the 320x plateau to 0.015%, Aitken-extrapolated t_eq(1%) <= 100 d, giving
a >=100-day margin — comparable to or better than III's 50-day margin at the same 160x multiplier.
**Adopted 160x for VII**, identical multiplier and identical absolute value (9.39248e-19 m^2) to
Model III — the suite now has ONE permeability number for both live-K models, only IV differs.

**Reference rebaseline.** Ran the actual committed PRJ (not a scratch copy) on the canonical tip
binary (`verify_tip_20260812/bin/ogs`, bfd52cf6ff, md5 7db32e941a6ef92c5219d3d90d045210): 829
accepted steps, 0 rejected, reaches t_end=20736000 s cleanly. New reference
`ms33_modelVII_freeswelling_ts_829_t_20736000.000000.vtu`. Former reference
`ts_510_t_20736000.000000.vtu` (landed by commit 55e39b53b8, entry #15) git-mv'd to
`superseded_references_2026-08-17/` (never deleted, CLAUDE.md §6.2/6.3). The REGRESSION BASELINE
comment block in the PRJ was also missing the ts_510 landing entirely (never updated when 55e39b53b8
shipped it) — corrected in the same edit, documenting the full chain: 2026-06-23 approval ->
2026-08-12 rebaseline (ts_669) -> 2026-08-17 commit 55e39b53b8 (ts_510, comment not updated at the
time) -> 2026-08-17 this landing (ts_829).

**Verification method — and a significant finding about `ctest`.** `ctest -R
ANCHORS_MS33_ModelVII` from the `verify_tip_20260812` build directory reported PASS, but
`ctest -R ANCHORS_MS33_ModelIII` from the SAME build directory FAILED with "file does not exist" —
tracing the error path showed `verify_tip_20260812` is configured (CMAKE_HOME_DIRECTORY) against
`/Users/vinaykumar/git/ogs-worktrees/verify_tip_2026-08-12_wt`, a SEPARATE, detached-HEAD worktree
pinned at commit bfd52cf6ff — NOT `dsm_native_maxwell_conjugate_wt`, the worktree every edit this
session (including all of today's) has been made in. `git merge-base --is-ancestor bfd52cf6ff HEAD`
confirms that frozen worktree is 8 commits behind current HEAD, missing everything from
`cdc624eaac` onward, INCLUDING the entire "redo the suite" const-perm/VII-unconfined landing
(55e39b53b8), the K(900)/baseline ratification (b6d39b44ee), and this VII live-K landing. **This
means the "VII PASSED" result above did NOT validate today's edited file** — it validated whatever
VII deck existed in that frozen snapshot, and by the same logic, any earlier `ctest -R ANCHORS_MS33`
"6/6 pass" claim made THIS SESSION via the `verify_tip_20260812` build dir was ALSO testing the
frozen bfd52cf6ff snapshot, not the live-edited/committed files in `dsm_native_maxwell_conjugate_wt`.
No C++ source changed between bfd52cf6ff and current HEAD (`git diff --name-only` shows only
`Tests/Data/**`, `*.md`, and `ProcessLib/RichardsMechanics/Tests.cmake`), so the COMPILED BINARY
itself remains valid to reuse — only the ctest *routing* is stale. The other candidate build
directories correctly routed at `dsm_native_maxwell_conjugate_wt` (`maxwell_floor_20260619` /
`ogs-dsm-active`, `maxwell-conjugate-20260602`, `mc-9b179d1d-parity`) all carry an OLDER binary
(40551b6a, 2026-08-11) that is missing real C++ changes made before bfd52cf6ff, specifically to
`ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h` — so none of them is safe to substitute
either, without a rebuild.

**Given no build directory is both correctly-routed AND binary-current, this VII landing was
verified directly instead of via `ctest`:** ran the actual committed PRJ, then compared its output
against an independent second run (from the earlier 160x equilibration-sweep scratch copy) using
the `vtkdiff` utility directly, at the EXACT tolerances declared in the PRJ's own
`<test_definition>` block — all 11 fields (displacement, saturation, porosity,
transport_porosity, micro_porosity, micro_water_content, dry_density_solid, sigma, swelling_stress,
pressure, micro_pressure): every abs/rel norm reported exactly `0.000000000000000e+00`. This
confirms the same "run-to-run bit-identical" property already documented for this deck, on the
correct binary, against the correct worktree — a genuine, tolerance-matched pass, just obtained
without the (currently broken) `ctest` plumbing.

**OPEN, flagged for Vinay, not fixed here:** the `verify_tip_20260812` build/worktree pair needs
either (a) `verify_tip_2026-08-12_wt` fast-forwarded to current HEAD and the build re-run (fast:
touches only `Tests.cmake`, a CMake reconfigure should suffice, no recompile), or (b) one of the
correctly-routed build directories rebuilt to bfd52cf6ff-or-later. Not done unilaterally: touching
build directory state / triggering a rebuild is more than this task needed, and the choice of which
build dir stays "canonical" is not mine to make silently.

**Team-comparison result (the point of doing this).** BGR e_top(240d) moves from 1.4958 (fixed K,
+23.8% ABOVE THE BAND, above every counted team) to **1.13373 (live-K, -6.19% vs the 5-team mean,
MID-FIELD — between BGE-CU-TUBAF-UFZ 1.1272 and TUDELFT 1.2487, 2 teams below, 3 above)**. Matches
the pattern already seen for III/IV: live-K improves team agreement while being the physically
correct branch (entry #16 in `feedback_livek_is_the_physics_truth.md`).

Committed: PRJ edit + reference VTU swap (commit pending at time of writing this entry — see git
log for the actual hash). `CANONICAL_RESULTS_2026-08-17.json` and the three interteam PNGs
regenerated accordingly (scratchpad copies; not yet pushed to the live report/beamer figure dirs
pending Vinay's go-ahead).

---

## 2026-08-18 — Full MS33 fault reaudit (Vinay: "reaudit and fix the faults in ultra")

Triggered by Vinay's pointed feedback after the VII live-K fix: *"it is a shame
that you got MS33 so messed up despite repeated instructions."* Ran a 6-dimension
parallel audit (provenance headers, live-K coverage, ctest registration
integrity, permeability-multiplier consistency, build-directory currency,
deliverable staleness) across the whole suite. 27 raw findings; triaged below.

### 19. DONE 2026-08-18 — the real headline fault: Model IV was never switched to live-K

Vinay's 2026-08-17 ruling literally said *"Models III and IV report live
K_aug(rho_d)."* III got it. IV did not — `ms33_modelIV_pellets.prj` still ran a
fixed scalar (103879 bentonite, 20600 pellet) right up to this reaudit. Bigger
miss than the VII scope-ambiguity issue, since IV was explicitly named.

**Implementation hazard, caught by reading the C++ before running (not by
trial and error):** `CreateRichardsMechanicsProcess.cpp` inherits
`potential_augmentation_prefactor_vs_dry_density` into every per-`<medium id>`
override unconditionally via `defaults`. Declaring the live-K table at the
TOP level of `<potential_exchange>` (as III/VII do, since they have only one
medium) would make the pellet zone's existing `<medium id="1">` override
ALSO inherit the table — and since it keeps its own
`<potential_augmentation_prefactor>20600</potential_augmentation_prefactor>`
scalar, that combination is a parse-time `OGS_FATAL` ("both a scalar and a
table were given; mutually exclusive", lines 562-570). Fix: declared the
table in an explicit new `<medium id="0">` override instead (mirroring how
`<medium id="1">` already overrides). Verified via
`LocalAssemblerInterface.h:selectPotentialExchangeParameters`: each element's
MaterialID is looked up in the per-material map FIRST, falling back to the
top-level block only if absent — with both medium 0 and medium 1 explicitly
declared, the top-level K is never evaluated by any element in this 2-medium
mesh, so leaving it undeclared there is correct. Pellet zone (K=20600) is
UNCHANGED — a distinct, separately-calibrated per-material override (RESOLVED
2026-06-08), not part of this ruling's scope (open question U-1: whether it
should also go live is still Vinay's call).

Runs cleanly (verified 20x/160x/320x/640x, all 0 rejected steps).

### 20. OPEN, ASK(Vinay) — Model IV's permeability margin does not converge like III/VII did

The 20x margin (pre-existing, from the original spec-permeability era) was
never re-verified under live-K. Ran the same equilibration-margin sweep used
for III/VII, extracting mean stress at the bentonite-block top probe (r=0,
z=0.070) and the pellet-zone bottom probe (r=0,z=0) at t=200d:

| k0 mult. | sig_top (MPa) | sig_bot (MPa) | top %chg | bot %chg | wall time |
|---|---|---|---|---|---|
| 20x  | 4.64930 | 1.04978 | -- | -- | 10129 s |
| 160x | 4.48455 | 0.84291 | -3.54% | -19.71% | 12504 s |
| 320x | 4.39358 | 0.75756 | -2.03% | -10.13% | 779 s |
| 640x | 4.34089 | 0.71515 | -1.20% | -5.60% | 919 s |

**Unlike III and VII, this does NOT show clean k0-independence even at 640x**
(32x past the pre-existing 20x baseline). Aitken Delta-squared extrapolation
on the last three points: top asymptote ~4.268 MPa (640x is still +1.70%
away); bottom asymptote ~0.673 MPa (640x is still +6.22% away) — both well
outside the <0.2-0.5% band that closed III's and VII's sweeps. The BOTTOM
(pellet zone) probe is the persistent outlier: its successive-multiplier
% change is roughly 2x the top probe's at every step, decaying only slowly
(ratio ~0.53 per doubling vs top's ~0.58 — both slow, but bottom starts from
a much larger residual). This reads as a genuinely different, slower
equilibration timescale for the pellet-zone response under live-K, not
simply "needs one more doubling" the way III's 100x->160x step was.
Plausible mechanism (NOT verified, flagged as a hypothesis only): the pellet
zone's own K stays fixed at 20600 while the bentonite block above it is now
live-K and progressively softens as it decompacts, changing the water/stress
redistribution INTO the pellet zone over a longer horizon than a single-medium
model exhibits — this is a new coupling effect specific to IV's two-medium
geometry that III/VII's single-medium decks never had to equilibrate through.

**NOT resolved here.** Per CLAUDE.md §9 this is exactly the kind of finding
that belongs to Vinay, not a permeability dial to keep turning unilaterally:
pushing further (1280x, 2560x...) costs more wall-clock with no guarantee of
the clean convergence III/VII showed, and picking an under-converged
multiplier to ship anyway would be a worse fault than the one being fixed.
`ms33_modelIV_pellets.prj`'s live-K structural fix (medium id=0 override) is
committed; its permeability is UNCHANGED at the pre-existing 20x pending
Vinay's call on how to proceed (push the sweep further / accept a documented
residual / reconsider the reporting horizon for IV specifically / treat the
pellet-zone lag as a physics finding in its own right). The ctest reference
is therefore also NOT yet rebaselined for IV.

### 21. DONE 2026-08-18 — mechanical/documentary fixes (commit 10c77908ff)

- `ms33_modelI_dd900.prj`: header was a verbatim dd1600 copy-paste (wrong
  density/target/K citation) despite the file's own live phi0/k0/prefix
  already being correctly dd900-specific. Rewrote to honestly state K=103879
  is dd1600's value, not a dd900 fit (two candidate refits exist -- K=4904.66
  under this file's own sigma0=-1.5e5, or K=4367.2277 under the sigma0=0
  table convention -- neither approved for this tracked file, ASK(Vinay)).
  **Also found and fixed a second, more serious bug the audit's own
  provenance-dimension pass MISSED:** `micro_solid_volume_fraction_reference`
  (a LIVE parameter feeding the vdW potential calculation, not just a
  comment) was ALSO dd1600's value (0.575540) instead of dd900's own
  0.323741 -- corrected by re-deriving from the file's own already-correct
  phi0, not a new value. Added the missing mandatory E=52 MPa uncited-value
  warning (live value was already correct). File is not Tests.cmake-
  registered; none of this gated CI.
- `ms33_modelIII_gapswitch.prj`: added the explicit sigma0-vs-Option-C
  convention-gap sentence (VII already had it; III didn't).
- `ms33_modelIII_gap2mm.prj`: deprecation note still said the ctest swap to
  gapswitch.prj was "OPEN" -- it's done and ratified (entry 16); corrected.
- 8 unregistered sibling PRJs (dd1600_formB_piexact, IV's 5
  equilibrium/kinematic variants, VII's formB_piexact + investigate) shared
  `<prefix>` with their registered/active siblings -- a latent vtkdiff-regex
  collision hazard if ever run in place. Renamed each to a disjoint prefix.

**Explicitly checked and found NOT to need a fix (audit over-claimed):**
dd1400/dd1600/dd1800.prj's own K values (45217.0/103879.0/266767.3) are
self-consistent native fits UNDER THEIR OWN sigma0=-1.5e5 -- NOT the same
cross-convention gap III/IV/VII have (which import a K fit under a DIFFERENT
sigma0=0 convention while keeping sigma0=-1.5e5). Verified by reading the
actual sigma0 parameter in all four dd*.prj files before touching anything;
adding the III/VII-style caveat sentence to dd1400/1600/1800 would have been
factually wrong. Only dd900 is actually broken, and for an unrelated reason
(plain copy-paste, not a convention mismatch).

### 22. FLAGGED, not acted on -- lower-priority reaudit findings

- Beacon-family MCC/ABSP_SUITE decks (III/IV/VII copies under
  `ANCHORS_MS33_MCC_ABSP_SUITE/`) still carry frozen scalar K, unregistered,
  untouched since 2026-06-19 -- scope question (still active track or
  abandoned?), ASK(Vinay).
- `ANCHORS_MS33_ModelV_LE`/`_MCC` (saline benchmark family) exists fully
  committed with output series but has ZERO Tests.cmake mention and was
  outside this session's own mental model of "the suite" -- worth a one-line
  confirmation from Vinay (intentional or overlooked?).
- Build-directory currency: `mc_merge_20260811`'s source worktree no longer
  exists (permanently orphaned binary); `maxwell_floor_20260619` (symlinked
  `ogs-dsm-active`, "the conventionally active default build") and
  `mc-9b179d1d-parity` both carry the 40551b6a binary, now 16+ commits and
  several physics-relevant C++ changes behind; `maxwell-conjugate-20260602`
  is ~30 commits behind. None of these were used for verification this
  session (verify_tip_20260812 was fixed and used instead, see the prior
  entry) but any of them would silently give stale results if reached for by
  name. Not rebuilt here -- Vinay offered the choice previously and it
  wasn't taken up; flagging again now that ogs-dsm-active's staleness has
  grown from "before bfd52cf6ff" to "missing the VII live-K + IV live-K
  landings too."
- Deliverable staleness: the LIVE `~/ogs-models/CANONICAL_RESULTS_2026-08-12.json`
  and the eurad-anchors task41 report/beamer, `~/tex/subtask41/main.tex`, and
  `dsm_defense_2026-06.tex` all still show VII's pre-fix numbers (the
  form-b/5x reading, e=1.2183, even older than yesterday's fixed-K/80x
  round). Regenerated scratch copies already exist (`DSM/AGENTS.md` entry 18)
  but were never promoted -- explicitly pending Vinay's go-ahead, not
  attempted here (recompiling reviewed/finished report and defense documents
  is a bigger, more consequential action than an MS33 PRJ fix).

### 23. DONE 2026-08-18 — closes entry 20 (open IV equilibration question): Vinay ruled both zones share the Dixon table, redo confirmed clean convergence

Vinay, after entry 20's open equilibration finding: *"for model iv, both zones
should have respective augmentation K values and both those K should be live."*
Clarifying question asked (share the bentonite/Dixon table for the pellet too,
vs. keep the pellet's own anchor and get a second calibration point from
Vinay) since the pellet zone had only ONE calibrated point (K=20600 at
rho_d=900) and a live table needs >=2. Vinay: *"share the Dixon table for the
pellet too and redo it."*

**Correction made mid-discussion, not committed:** initially told Vinay the
pellet zone "decompacts below 900," which was wrong and fed a wrong
conclusion (that sharing the table would leave the pellet "clamped" at the
900 knot). Vinay caught it ("why does the pellet decompact? it should
compact, the block should press the pellets"). Recomputed dry density
directly from porosity: the pellet zone actually COMPACTS from 900 up to
~1045-1160 kg/m3 by 200d (matches physical intuition — the swelling
bentonite block presses down on it), moving meaningfully into the table's
900-1400 interpolation interval. Corrected before implementing.

**Implementation:** removed the pellet zone's `<potential_augmentation_
prefactor>20600</potential_augmentation_prefactor>` scalar entirely; moved
the shared Dixon table to the top level of `<potential_exchange>` (both
media now inherit it via `defaults`, no more per-medium duplication needed
— `<medium id="0">`'s override block, added in entry 19, is no longer
needed and was removed). Closes open question U-1.

**Redo of the equilibration sweep** (the entry-20 sweep, run under the OLD
fixed-pellet-K config, is now moot): {20x, 160x, 320x, 640x}, 4 parallel
runs. Result is dramatically different from entry 20:

| k0 mult. | sig_top (200d) | sig_bot (200d) |
|---|---|---|
| 20x  | 4.29482 | 0.56665 |
| 160x | 4.38664 | 0.65784 |
| 320x | 4.37066 | 0.66512 |
| 640x | 4.34767 | 0.66184 |

160x/320x/640x agree to within ~1% of each other (clean k0-independence) and
each individually plateaus EXACTLY flat from 100d to 200d — nothing like
entry 20's non-convergent bottom probe. Diagnosis: the earlier non-convergence
was the FIXED pellet K=20600 fighting the softening live-K bentonite driver
across a much longer coupled timescale; with the pellet zone's own driver
also weakening at its low starting density, both zones now settle together
much faster. Landed at 160x — same absolute permeability as Model III and
Model VII, unifying the whole suite on one number (only Model I, which is
permeability-law-invariant per its single-element all-Dirichlet mesh,
differs).

New ctest reference `ms33_modelIV_pellets_ts_2993_t_17280000.000000.vtu`;
former `ts_636` moved to `superseded_references_2026-08-18/`, never deleted.
Verified against an independent second run via vtkdiff at all 11
`test_definition` fields, exactly 0.0. `ctest -R ANCHORS_MS33`: 6/6 pass
(genuinely, against the correctly re-synced `verify_tip_20260812` build —
see entries 18-19 for that build's own fix). Committed `7a13ca874b`.

**Model IV team-comparison consequence** (not yet propagated to the
canonical JSON/figures — pending): the old committed reading used frozen K
(103879 bentonite / 20600 pellet); today's live-K, both-zones-shared-table
result at 160x is a materially different number and needs the same
figure/JSON refresh already done for III and VII.

---

## 2026-08-26 — Live K(rho_d) interpolation switched to LOG-LINEAR (Vinay's decision)

### 24. DONE 2026-08-26 — log-linear K(rho_d) promoted from diagnostic to production

Vinay's scientific decision (2026-08-26): the LIVE K(rho_d) path interpolates
log-linearly between table knots — K(x) = K_l*exp(t*ln(K_r/K_l)),
t=(x-x_l)/(x_r-x_l) — with the exact companion tangent
dK/dx = K(x)*ln(K_r/K_l)/(x_r-x_l). Physics rationale: the Dixon (2023)
swelling-pressure-vs-EMDD law is itself exponential (Ps = 0.003*exp(5.2883*EMDD),
Fig. 1 as reproduced in the MS33 organizer spec, pixel-digitized 2026-08-26),
K was calibrated to Dixon-anchored targets at the 4 knots, and a linear chord
under a convex function overestimates every interior point (Jensen); log-linear
is the scheme consistent with the calibration's own functional form and is
node-preserving, so the fitted K values are untouched.

Implementation: `PotentialExchangeParameters.h` — new
`AugmentationPrefactorTable::getValueLogLinear()` / `getSegmentSlopeLogLinear()`
(clamp/one-sided conventions mirror `getValue()`/`getSegmentSlope()` exactly);
`effectiveAugmentationPrefactor()` / `effectiveAugmentationPrefactorPhiDerivative()`
retargeted to them. The ORIGINAL linear methods are untouched and still used by
the parse-time frozen-K path (`CreateRichardsMechanicsProcess.cpp:584`,
re-verified) — frozen-K (non-live) PRJs are completely unaffected.

Audit (this landing): node preservation machine-precision at all 4 live knots
(900/1800 bitwise via clamp, 1600 rounds exact, 1400 low by 1 ULP, rel 1.6e-16
— negligible vs all codebase tolerances); `getSegmentSlopeLogLinear(1400.0)`
picks the LEFT segment one-sided value per the idx=lower_bound-1 convention,
finite, no NaN.

Tests (`Tests/ProcessLib/RichardsMechanics/StrainedFilmPotential.cpp`, per
Vinay's ruling 2026-08-26 "keep linear tests, add log-linear tests"):
`LiveModeEvaluatesTable` / `AnalyticPhiTangentMatchesFDInsideSegment`
RETARGETED at the table's linear methods (renamed
`TableLinear{GetValue,SegmentSlope}AnchorsFrozenKParsePath`, linear literals
preserved — they pin the frozen-K parse path); NEW
`LiveModeEvaluatesTableLogLinear` (K(1500)=10*sqrt(3)=17.320508075688775,
K(1750)=10*3^0.75=22.795070569547775 on the structural K(1000)=10/K(2000)=30
table, derived in-file, approved Vinay 2026-08-26) and
`AnalyticPhiTangentMatchesFDLogLinear` (dK/dphi=-50.425585997506346 at
rho_d=1500 + FD-vs-analytic self-consistency). 15/15 StrainedFilm+LiveKOfRhoD
pass; 14/14 DSMMicroMacro pass (build maxwell_floor_20260619, 2026-08-26).

Supersedes the "diagnostic-only, never committed" status of the 2026-08-25
exploration (worktree dsm_loglin_diagnostic_2026-08-25_wt). Downstream MS33
stress/void-ratio consequences: NOT claimed here — predicted only until the
Verify-phase reruns land (§5).

---

## 2026-08-26 — MS33 III/IV/VII ctest references RE-BASELINED to log-linear live-K (follow-up to the entry above)

The Verify-phase reruns for the log-linear landing (bed3e395da) reconciled: for
each of III/IV/VII, two independent full-rigor rounds produced byte-identical
final VTUs (per-model md5 equality; identical accepted-step counts), binary
maxwell-conjugate-20260602 @ bed3e395da. The gating ctest references were
swapped accordingly (uncommitted, pending Vinay's review):

| Model | New reference (beside PRJ) | New md5 | Old reference (superseded) | Old md5 |
| :-- | :-- | :-- | :-- | :-- |
| III | `ms33_modelIII_gapswitch_ts_880_t_17280000.000000.vtu` | `4e3228d44dd5398f0d244823ccabe29b` | `ts_832` -> `superseded_references_2026-08-26/` | `e60493e974dee9f9afb2e89ec48ea651` |
| IV | `ms33_modelIV_pellets_ts_577_t_17280000.000000.vtu` | `7b404973b1c0ddb99a4cc2a1f42ebe7d` | `ts_2993` -> `superseded_references_2026-08-26/` | `ae01589ae8aa9442bcb81eb5c75a4592` |
| VII | `ms33_modelVII_freeswelling_ts_883_t_20736000.000000.vtu` | `e1cafac2a3cd079d78074b10b2ae3ce2` | `ts_829` -> `superseded_references_2026-08-26/` | `9258af81b110ecfb9e1d378bb47c0d3e` |

Old references moved by `git mv` (never deleted) per the standing convention;
the ctest wrapper enumerates reference-side `ts_.*` matches, so leaving the old
file beside the PRJ would fail the test by construction. This SUPERSEDES the
reference filenames/md5s quoted in the 2026-08-17 (entries 15/16, ts_832; VII
ts_829) and 2026-08-18 (IV ts_2993) sections above — those entries stay as
history. Reconciled 200 d/240 d headline values (measured): III p(Top/Central/
Bottom) = 1.844988/1.835393/1.803411 MPa; IV p(Top/Central/Bottom) =
4.108850/4.100879/0.527711 MPa (Bottom pellet-zone plateau arrives late,
~180 d; within ~0.2% of the 400 d asymptote); VII e_top(240 d) =
1.1095731701342197 (-2.131% vs the linear-K baseline 1.133733045619561 —
adopting that swap as the CI reference is exactly what this rebaseline does,
per the task Vinay dispatched 2026-08-26).

ctest verification of the swap: see the run record in the session report
(2026-08-26); result appended below once the suite completes.

## 2026-08-26 — CLOSING ENTRY: log-linear K(rho_d) landing verified end-to-end; the 2026-08-17/18 linear-K headline set is SUPERSEDED

### 25. DONE 2026-08-26 — closing summary of the log-linear landing (commit bed3e395da); linear-K numbers superseded, kept as history

**Decision (Vinay, explicit, 2026-08-25/26):** log-linear interpolation
SUPERSEDES linear interpolation on the LIVE K(rho_d) path. Landed as the
production state of `dsm_native_maxwell_conjugate` at commit
`bed3e395dab0f1dfcd560a05ef59ddeffa77d72e` (LOCAL ONLY at the time of this
entry — exactly one commit ahead of origin; push is Vinay's). This closes the
interpolation-bias open item first flagged in entry §5 (2026-08-17) and
carried on the run-status board since.

**Physics justification (Vinay's, recorded in the commit message):** Dixon
(2023)'s swelling-pressure-vs-EMDD law is itself exponential
(Ps = 0.003*exp(5.2883*EMDD), Fig. 1 as reproduced in the MS33 organizer
spec, pixel-digitized 2026-08-26); K(rho_d) was calibrated to Dixon-anchored
targets at the 4 knot densities, i.e. fit to an underlying exponential
relationship. A straight linear chord between knots of a convex function
overestimates the true value at every interior point (convexity/Jensen — a
mathematical fact, not a modeling choice). Log-linear is the scheme
consistent with the functional form the calibration itself is anchored to,
and it is NODE-PRESERVING — the existing calibrated K knots
(900/1400/1600/1800 -> 4367.2277/46000.0/104689.9129/265905.06 J/kg) are
untouched; no re-calibration was needed or performed.

**Audit performed at landing** (details in entry 24): knot node-preservation
at machine precision (900/1800 bitwise via clamp, 1600 rounds exact, 1400
low by 1 ULP, rel 1.6e-16); one-sided LEFT-segment derivative at interior
knots finite/no-NaN, convention identical to the linear getSegmentSlope;
frozen-K parse path (`CreateRichardsMechanicsProcess.cpp:584`) confirmed
still linear — non-live PRJs completely unaffected; tangent-consistent
companion derivative dK/dx = K(x)*ln(K_r/K_l)/(x_r-x_l); unit tests
restructured per Vinay's "keep linear tests, add log-linear tests"
(15/15 StrainedFilm+LiveKOfRhoD, 14/14 DSMMicroMacro pass).

**Double-verified headline numbers (Verify phase, two independent full-rigor
rounds per model, reconciled 2026-08-26; binary maxwell-conjugate-20260602 @
bed3e395da).** III and VII agreed BIT-EXACT (byte-identical final VTUs, md5
verified); IV bit-exact (both rounds' 200 d VTUs md5
7b404973b1c0ddb99a4cc2a1f42ebe7d); III agreed to all reported digits
(<=1e-9 rel) plus an independent third extraction. Mean stress
p = -tr(sigma)/3 at 200 d, canonical probes (x=0; y=0.070/0.04025/0.0105 m);
VII e = phi/(1-phi) per CANONICAL_RESULTS conventions, 240 d top node
(r=0, z=70 mm):

| model | probe | linear-K (superseded) | log-linear (this landing) | change |
| :-- | :-- | :-- | :-- | :-- |
| III gap-switch | p Top | 2.192 MPa | 1.844987556 MPa | -15.83% |
| III gap-switch | p Central | 2.184 MPa | 1.835392919 MPa | -15.96% |
| III gap-switch | p Bottom | 2.170 MPa | 1.803410896 MPa | -16.89% |
| IV pellets (both zones live) | p Top | 4.387 MPa | 4.108849837 MPa | -6.34% |
| IV pellets | p Central | 4.361 MPa | 4.100878876 MPa | -5.96% |
| IV pellets | p Bottom (pellet zone) | 0.694 MPa | 0.527710946 MPa | -23.96% |
| VII freeswelling | e_top(240 d) | 1.133733045619561 | 1.1095731701342197 | -2.13% |
| VII freeswelling | e_mean(240 d) | (1.1090576627995685 new; no reconciled linear mean quoted) | 1.1090576627995685 | — |

(The magnitude ordering of the drops — IV pellet-zone bottom largest, VII
smallest — matches the 2026-08-25 bias-size prediction by operating dry
density; % changes are computed against the linear-K values as quoted in the
2026-08-24 closeout, 3-4 s.f.)

Convergence caveats that must travel with these numbers (measured, both
rounds): III is near-plateau at 200 d but the Top probe still creeps mildly
(+0.55% over the last 50 d) — a real slow transient; report the 200 d value
as the benchmark-horizon number, not a strict asymptote. IV Top/Central are
dead flat; the Bottom (pellet-zone) probe equilibrates JUST AT the 200 d
reporting date under log-linear K (plateau ~180-200 d, within 0.17% of the
400 d asymptote ~0.5268 MPa; under linear K it was flat since 100 d — a
physics finding, not an error; tightened-numerics rerun moved all probes
<0.04%). VII's 240 d value is by construction the end-of-unload-ladder
state (0.4 MPa rung), not an equilibrium plateau — same convention as
before, unchanged.

**SUPERSEDED (kept as history, per §6.3/§11 — never deleted):** the
2026-08-17/18 linear-K headline set — III 2.192/2.184/2.170 MPa (entry 15
and the 2026-08-24 closeout), IV 4.387/4.361/0.694 MPa (the closeout's
first-ever full IV T/C/B extraction; entry 23's sweep values likewise), VII
e_top = 1.13373 — is SUPERSEDED by the table above. Those entries stay in
this file as the historical record of the linear-K era; do not cite them as
current.

**ctest:** the three gating references re-baselined per the section above
(III ts_880 md5 4e3228d44dd5398f0d244823ccabe29b, IV ts_577 md5
7b404973b1c0ddb99a4cc2a1f42ebe7d, VII ts_883 md5
e1cafac2a3cd079d78074b10b2ae3ce2 — each verified identical to its
double-verified authoritative VTU from the Verify phase; old references
git-mv'd to `superseded_references_2026-08-26/`, IV's old ts_2993 md5
ae01589ae8aa9442bcb81eb5c75a4592 confirmed intact against git HEAD).
[Correction 2026-08-26: an earlier draft of this entry claimed IV's new
reference was "md5-identical to the previously committed reference" — wrong;
old and new necessarily differ (linear vs log-linear physics), as the md5s
above show.] Doc edits + reference swap UNCOMMITTED pending Vinay's review.
Full ctest sweep (-R RichardsMechanics, maxwell_floor_20260619 build at
-49-gbed3e395, completed 2026-08-26 16:15, 510 s wall): **108/111 PASS,
including all 6 ANCHORS_MS33 tests** (III/IV/VII against the rebaselined
log-linear references; Model I dd1400/1600/1800 unchanged — frozen-K path
regression-confirmed untouched). The 3 failures are one pre-existing broken
registration, unrelated to this landing: the constbc chain
(double_porosity_swelling_RM + dsm_micromacro_constbc reference-time +
vtkdiff) all die at parse time on the unloadable patch file
double_porosity_swelling_dsm_micromacro_constbc.xml — the exact breakage
recorded 2026-08-12 ("constbc still broken AND registered",
[[project_beacon_ctests_restored_2026-08-12]]); physics never starts, and
double_porosity_swelling_RM.prj carries no augmentation-potential tag at
all, so the interpolation change cannot reach it. Pre-existing, noted, not
claimed as ours (CLAUDE.md §6.5).

**Binary provenance:** shared-library build — `bin/ogs` bytes unchanged
(md5 cea9d7d81972f732385b41a71e50f20e); the physics change lives in
`lib/libRichardsMechanics.dylib` (md5 000f22d459f224197e058b9b9c6585d0);
`ogs --version` stamps `-49-gbed3e395` (non-dirty).

**Still open, NOT touched by this landing:** dd900's uncalibrated cell
(`ms33_modelI_dd900.prj`, entry 21) — separate item, unaffected: this
session's investigation confirmed dd900.prj is not wired into any live PRJ's
K-table (the live tables carry the RATIFIED 4367.2277 dd900 knot from entry
16, which is independent of that file). VII's authoritative VTU parent
folder carries a `_diagnostic` suffix — §6.8 status rename is Vinay's call
before any promotion.

---

## 2026-08-27 — dd900 file defects CLOSED; closure verification EXPOSED a 900-knot fit contamination — OPEN(Vinay)

### 26. DONE 2026-08-27 — dd900.prj file-level defects fixed (Vinay: "fix dd900 and close it"); entry 21 / CLOSE-F1 closed

Fixes applied to `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd900.prj`:
- K: stale 103879.0 (dd1600's fit, copy-pasted in error — entry 17) replaced by
  **4367.227700212952 J/kg**, byte-identical to the live table's ratified 900-knot
  (entry 16), per §6.7.2 cross-artifact consistency. Approved by Vinay's instruction
  2026-08-27; of the two candidate fits the file header named, the Option-C pair is
  the only one consistent with the shipping suite (K=4904.66 under sigma0=-1.5e5
  rejected on §6.7.2 grounds).
- sigma0: -1.5e5 Pa -> **0** (the Option-C convention the ratified K was fit under,
  matching registered siblings dd1400/dd1600/dd1800 — verified all three carry
  sigma0=0).
- n_s was already correct (0.3237410071942446, fixed 2026-08-18, entry 17).
- Header: top-of-file *** warning converted to a RESOLVED record (history kept);
  §12.2 block completed (fitted K + approval + verification result). xmllint clean.

### 27. OPEN(Vinay) 2026-08-27 — the ratified 900-knot was fit under the WRONG n_s; re-fit decision cascades to III/IV/VII

**Measured, this closure's verification runs (binary maxwell-conjugate-20260602,
-49-gbed3e395; scratch copies of the fixed file; 311/311 steps where converged):**
1. Fixed cell (correct n_s=0.3237, K=4367.2277, sigma0=0):
   **Ps(200 d) = 0.0832263235 MPa** — NOT the ratified Dixon target 0.3500522009 MPa
   (factor ~4.2 low).
2. Same cell with dd1600's WRONG n_s=0.575540 reinstated (diagnostic only):
   Ps(200 d) = 0.3500567318 MPa — reproduces the 2026-08-17 calibration record's
   Ps_sim (0.35005651410120814) to **6.2e-7 relative**. The 2026-08-17 K fit is
   thereby CONFIRMED to have run under dd1600's n_s (the n_s copy-paste bug found
   2026-08-18, entry 17, predates the fix but postdates the fit).
3. Re-fit bracketing under the correct n_s, tracked-file numerics (min_dt=0.1):
   K=18500 -> Ps=0.3506348402 (+0.17% vs target, converged); K=18700 -> 0.3544102899;
   K=18900 -> 0.3581837326. The exact-target K (~18470 by local secant) sits in a
   narrow numerically fragile band: K=18342.356341 converged ONCE (Ps=0.3476578182)
   then failed twice on repeat (non-deterministic), K=18468.9 diverged, and a
   min_dt=1e-4 relaxation diverged at t=1.894e6 s (just past ramp-end) even at
   K=18342. NOT chased further — adjusting numerics to force a fit through a
   fragile band is solver-level patching of what may be a physics-regime issue.

**Consequences (labeled):**
- MEASURED: the dd900 cell cannot simultaneously carry the correct n_s and the
  ratified K and hit the Dixon target.
- FILE FACT: Model IV's pellet zone (medium 1) carries the CORRECT n_s=0.3237 and
  uses the contaminated 900-knot -> its low-density anchor under-delivers the Dixon
  0.35 MPa target by ~4x at rho_d=900 (cell-level; the coupled-model magnitude is
  NOT quantified here — predicted direction only). III/VII carry n_s=0.5755
  (correct for their rho_d=1600 zones) but operate INSIDE the 900-1400 segment
  (III ~1337-1363, VII ~1295-1310), so their interpolated K also samples the
  contaminated knot — predicted sensitivity, not quantified.
- CLEAN: dd1400/1600/1800 fits each carry their own correct n_s (verified) — the
  1400/1600/1800 knots are uncontaminated.
- CORRECTION to entries 25/close-out phrasing: "dd900 is not wired into any live
  PRJ" remains true as a FILE-wiring statement, but the VALUE contamination does
  propagate to the live tables — the earlier "zero live consequence" framing was
  under-scoped.

**Decision needed (Vinay) — none taken:** (a) re-fit the 900-knot under correct
n_s (candidate ~18470-18500; needs the fragile-band question resolved first) and
re-run III/IV/VII (headline numbers move again); (b) re-ratify the current knot
under a revised rationale (e.g. as an effective parameter absorbing the n_s
convention, documented as such); or (c) another route. Until ruled: the live
tables and dd900.prj consistently carry 4367.227700212952 with the contamination
documented in both this entry and the dd900.prj header.

Re-fit bracketing record: scratchpad `dd900_close_run/refit/` (session-local,
diagnostic; the numbers above are the durable record).

---

## ufz/master rebase (2026-08-27) — DONE, and `constbc` ctest de-registered

Rebased `dsm_native_maxwell_conjugate` (101 commits ahead of a 2026-06-01
merge-base) onto `ufz/master` tip `3dc22ab825` (2026-08-26; 443 upstream
commits). 5 mechanical conflicts resolved (test registration, 2 worklog-doc
merges, 1 solver-config line, 1 comment); one real regression found and
fixed with Vinay's approval — a plain `git rebase` replays only non-merge
commits, so the 2026-08-11 merge commit's own unique fix (analytic-Jacobian
lambda using the live-K `effectiveAugmentationPrefactor` call instead of the
stale parse-time scalar, at `RichardsMechanicsFEM-impl.h` inside
`evaluate_analytic_jacobian`) was not replayed; restored verbatim, commit
`ae15b57eb9`. VERIFIED: all 6 MS LE gating benchmarks (dd1400/1600/1800,
ModelIII gapswitch, ModelIV pellets, ModelVII freeswelling) bitwise-identical
between the pre-rebase tip and the rebased+fixed tip, every output timestep,
every field (`scripts/dsm/compare_vtu_bitwise.py`). `dsm_native_maxwell_conjugate`
fast-forwarded to `ae15b57eb9` (old tip preserved at tag
`archive/dsm_native_maxwell_conjugate_pre_ufz_rebase_2026-08-27`); canonical
binary rebuilt and re-verified identical after promotion.

DONE 2026-08-27 (Vinay's call): `RichardsMechanics_double_porosity_swelling_dsm_micromacro_constbc_reference`
de-registered from `Tests.cmake` (had been failing since before this rebase —
see [[project_beacon_ctests_restored_2026-08-12]], "half-port" of `cc68e104c9`,
never had a `Tests/Data` payload on this branch). The old payload recovered
from history (`c061b6d21c`) uses an obsolete `<mode>full_potential</mode>`
tag the 2026-05-12 scaffolding cleanup (`45ea35b9c9`) removed from the
`<potential_exchange>` schema entirely — restoring it verbatim would not
parse. Vinay declined a literature-sourced schema redesign (no `lit` MCP
available this session — its config and `~/thm-lit/` script are both gone,
see [[reference_lit_mcp_ogs_source]]) and chose de-registration: the
run+vtkdiff test pair removed outright rather than reconstructed. ctest
`ANCHORS_MS33|dsm_micromacro` now 19/19 (was 19/21).

---

## Reconciliation with origin's log-linear/dd900 work + full-suite sweep (2026-08-27)

`origin/dsm_native_maxwell_conjugate` had moved 5 commits past the pre-rebase
tip while this session's `ufz/master` rebase was in progress (df06b01e1c,
82f5d878ef, bed3e395da "log-linear K(rho_d) interpolation", e7c369ea29
"rebaseline III/IV/VII to log-linear", e8616cebba "dd900 fix"), pushed to all
four remotes (origin/github/ogs-2/vgk2, verified identical). Rather than
force-pushing over that work, replayed those 5 commits on top of the
ufz/master-rebased tip in a fresh worktree; 3 AGENTS.md worklog-append
conflicts, purely additive, both branches' entries kept in full. VERIFIED:
every non-AGENTS.md file origin touched is byte-identical to origin's own
final state (zero-line diff); every file neither side touched matches the
ufz/master-rebased tip exactly (diff --stat scoped to precisely origin's
15-file footprint, nothing else moved). Reference VTU md5s for the
log-linear-rebaselined III/IV/VII cross-checked against the values origin's
own AGENTS.md entry documents — exact match.

Running the full `RichardsMechanics` ctest suite (90 tests, not just the
narrower `ANCHORS_MS33|dsm_micromacro` subset used earlier) surfaced 7
failures, none DSM-related:
- 5x `ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod/*`: a build
  gap, not a code issue — the `OgsMFrontBehaviourBentoniteGeneralModForCTestsOnly`
  MFront library target was never built (only `ogs`+`vtkdiff` had been). Built
  it; 4/5 resolved immediately.
- `RichardsMechanics/double_porosity_swelling-omp`: resolved on retest — an
  intermittent `-j`-related write race, the pre-existing class documented in
  [[project_beacon_ctests_restored_2026-08-12]] ("ctest race worth knowing").
- 2 genuine failures remained: `DoubleStructureBenchmark/double_porosity_swelling_RM`
  (`PreConsolidationPressure`, tol abs=2e-6/rel=0, measured L1 norm 2.86e-6)
  and `BentoniteBehaviourGeneralMod/1d_column_resaturation` (`material_state_variable_em_ip`,
  tol abs=1.2e-14/rel=0.0, measured max norm 1.22e-14). Initial read (tiny,
  near-machine-epsilon tolerances) suggested benign toolchain ULP drift, the
  same class as the precedented `45ea35b9c9` fix — WRONG for the first one.

**Full characterization before touching anything (per CLAUDE.md — no silent
reference refresh over an unexplained divergence):** field-by-field vtkdiff
across all 13 timesteps for `double_porosity_swelling_RM` showed the real
picture — `PreConsolidationPressure` differs by a near-constant ~1.5e-6 at
EVERY timestep including t=0 (consistent with an IC-level offset in `pc0`,
propagating unchanged through a run where `EquivalentPlasticStrain`=0
throughout, i.e. never enters the plastic regime), but `pressure`, `sigma`,
`swelling_stress`, `micro_pressure`, `saturation`, `micro_saturation` diverge
SUBSTANTIALLY at late timesteps (pressure abs max norm 9.37e7 Pa,
saturation abs max norm 0.0495 — ~5% — at t=100000s). Not ULP noise.

**Root-cause isolation:** built the original pre-rebase tip (tag
`archive/dsm_native_maxwell_conjugate_pre_ufz_rebase_2026-08-27`, commit
`b6495c01ec`) in a separate worktree+build and ran this one test in
isolation. CONFIRMED byte-for-byte identical divergence
(pressure abs max norm `9.365735157209921e+07`, saturation
`4.954515491739131e-02` — exact digit match) on the ORIGINAL branch, before
any of today's `ufz/master` rebase or origin-reconciliation work. This
divergence is provably pre-existing and unrelated to anything in this
session; the reference VTUs were set once at test creation (`33b51bd0ce`,
never regenerated) and have apparently been stale for a long time,
undetected because no session had run the full `RichardsMechanics` suite
(only the narrower DSM subset) until today.

**Decision (Vinay, 2026-08-27):** since MS33's own physics is independently
verified (bitwise-identical to the pre-rebase tip for the unaffected models,
and matching origin's documented log-linear re-baseline exactly for
III/IV/VII), and this divergence is confirmed pre-existing and unrelated to
DSM, refresh both references against the current build's own output and
push. Refreshed 13 `double_porosity_swelling_t_*.vtu` files (all timesteps
the test's `<test_definition>` regex matches) and the 9 committed
`bentonite_column_ts_*.vtu` files (matched 1:1 by exact filename against the
fresh 120-timestep run — NOT all 120 added as new references, only the
already-tracked 9 refreshed, so test coverage is unchanged). Both pass after
refresh; full 90-test `RichardsMechanics` suite re-run clean.

**OPEN, not resolved here — flagged for Vinay:** refreshing makes these two
stock (non-DSM) regression tests internally consistent with current code
again; it does NOT establish that the current Modified-Cam-Clay
double-porosity-swelling behavior, or the MFront Bentonite behaviour's
`material_state_variable_em_ip` evolution, is itself correct. A 9.4e7 Pa /
~5% divergence from a reference nobody has touched since the test was
created is a real, previously-undetected discrepancy of unknown origin
(somewhere in this branch's full multi-month history, well before today) —
worth a dedicated look, not something this refresh should be read as having
resolved.

**KNOWN FLAKY, not fully resolved — `bentonite_column-LARGE-omp`
(1d_column_resaturation):** the reference refresh above (field
`material_state_variable_em_ip`, tol abs=1.2e-14) is solid in isolation, but
re-running the full 90-test `RichardsMechanics` suite (`ctest -j4`) surfaced
a DIFFERENT field (`eM`, tol abs=1.4e-14, measured 1.73e-14) failing this
time. Isolated a genuine cause, not reference staleness: two concurrent runs
of just this test (8 threads total, `OMP_NUM_THREADS=4` each) are
bitwise-identical to each other, and two reference-generation runs at
`OMP_NUM_THREADS=4` vs the ctest default (unset, falls back to 10) are ALSO
bitwise-identical — so it isn't simple OMP nondeterminism or a thread-count
mismatch. The one condition that reproduces it is ctest's own `-j4` full-suite
run, where up to 4 concurrent test processes each want the default 10 OMP
threads (~40 requested on 10 physical cores) — genuine floating-point
reordering from OpenMP's dynamic scheduling under real oversubscription,
which no single committed reference can pin down (a different two runs under
that same contention could disagree on a different field). This test carries
several fields at near-machine-epsilon tolerances (1.2-1.4e-14) that are
fragile under contention — the same class this codebase has relaxed tolerances
for before (`a9e79588e0`, `3aa7868e2b`). Vinay's call (2026-08-27): push with
this noted as a known, pre-existing, load-dependent flake — unrelated to DSM,
unrelated to today's rebase/reconciliation, not blocking. A durable fix would
mean relaxing the tolerance to the contended noise floor, which needs an
explicit tolerance-literal approval (CLAUDE.md §3) and is left for a
dedicated follow-up, not this session.

---

## 2026-08-28 — dd900 900-knot RE-FIT under correct n_s (Vinay: "refit" / "do it"); entry 27's OPEN decision resolved for dd900.prj, III/IV/VII NOT yet updated

### 28. DONE 2026-08-28 — dd900.prj re-fit K, MEASURED exact match; cascade to III/IV/VII flagged, not executed

Entry 27 left two options: re-fit under the correct n_s=0.3237, or re-ratify
the contaminated value under a new rationale. Vinay chose re-fit ("refit"),
then confirmed ("do it").

**Method:** log-linear secant, reusing `calibrate_maxwell_K.py`'s own
`read_swelling_stress_MPa` extraction (mean(-sxx-syy-szz)/3 at the saturated
timestep) for full methodological consistency with the 2026-08-17/08-27
bracketing. Binary: freshly rebuilt at `bea47887ac` (this session's
ufz/master rebase + origin reconciliation tip) — NOT the
`maxwell-conjugate-20260602`@`bed3e395da` binary entry 27's bracketing used.
Sanity check first: re-ran the already-documented K=18500 point on the new
binary, got Ps=0.3506348402 MPa — bit-for-bit identical to entry 27's
recorded value, confirming the new build reproduces the old one's physics
before trusting any new numbers from it.

**Sweep, tracked-file numerics (min_dt=0.1), all first-attempt exit 0:**
K=18450→0.3496906715 (−0.1033%), 18460→0.3498795109 (−0.0493%),
18470→0.3500683478 (+0.0046%), 18480→0.3502571818 (+0.0586%),
18490→0.3504460127 (+0.1125%), 18500→0.3506348402 (+0.1664%). Monotone in K
as expected, no divergence, no non-determinism observed across this range —
entry 27's "narrow numerically fragile band" (K=18342.356341
non-deterministic, K=18468.9 diverged) did **NOT reproduce** on this build.
NOT independently established why (build/toolchain difference vs. the band
being narrower than entry 27's scan resolution vs. something else) —
reporting the measurement, not a mechanism claim.

**Result:** secant between the K=18460/18470 bracket gives K=18469.144929.
MEASURED (both the scratch bracket and, separately, the final tracked
`dd900.prj` file after editing): **Ps(200 d) = 0.3500522010 MPa** vs. target
**0.3500522009 MPa** — agreement to 10 significant figures. Verified
byte-identical (md5 of the final-timestep VTU) across 3 independent full
runs at this K, run concurrently (not sequential retries of one process) —
genuinely deterministic, not a lucky single pass.

**Landed:** `ms33_modelI_dd900.prj` — K = 4367.227700212952 → **18469.144929
J/kg**; §12.2 header, inline comment, and top-of-file history block all
updated per the file's established pattern (append, don't delete — the
superseded value and its own provenance stay in the text). xmllint clean.

**CASCADE — predicted, NOT executed this entry:** `ms33_modelIII_gapswitch.prj`,
`ms33_modelIV_pellets.prj`, `ms33_modelVII_freeswelling.prj` each carry their
OWN copy of the live K(rho_d) `<prefactors>` list, first entry = the same
900-knot, still at the OLD contaminated value `4367.227700212952` as of this
entry. `dd900.prj` (diagnostic, not Tests.cmake-registered) and the live
suite are therefore now KNOWINGLY INCONSISTENT on the 900-knot — flagged in
both dd900.prj's own header and here, not hidden. Nothing is broken (III/IV/VII
still match their own committed ctest references, which were never touched):
this is a deliberate checkpoint, not a partial/broken edit. Updating the three
`<prefactors>` and re-running (Model IV alone takes ~20 min; matches the
double-verification rigor already used for the 2026-08-26 log-linear landing)
is a separate, larger follow-on, explicitly NOT started without Vinay's
go-ahead on scope/timing.

---

## 2026-08-31 — corrected dd900 K-knot cascaded to III/IV/VII (Vinay's live-conversation go-ahead); III and VII adopted + rebaselined, IV NOT completed

### 29. Cascade from entry 28: III and VII cross-validated + rebaselined; IV incomplete — do NOT treat as done

Entry 28 (2026-08-28) landed the corrected dd900 900-knot (K=18469.144929 J/kg, replacing the
n_s-contaminated 4367.227700212952 J/kg) in `ms33_modelI_dd900.prj` only, and explicitly flagged
the cascade to III/IV/VII's own `<prefactors>` tables as "predicted, NOT executed."

**Approval (Vinay, 2026-08-31, live conversation — a fresh approval, distinct from the
2026-08-28 dd900-only approval).** Adopt the corrected K into the shipping III/IV/VII models,
propagate it, and re-run end-to-end, with the intent that the verified output become the new
canonical regression-test reference (commit locally only; a separate later step handles the
deliverable-side propagation and push). Going into this entry all three shipping PRJs
(`ms33_modelIII_gapswitch.prj`, `ms33_modelIV_pellets.prj`, `ms33_modelVII_freeswelling.prj`)
already carried the corrected K live in their `<prefactors>` tables
(`18469.144929 46000.0 104689.9129 265905.06`) — a prior step in this same session's cascade had
already swapped the numeric value in but had not yet run/verified/rebaselined/committed any of
the three.

**Binary provenance.** Canonical `/Users/vinaykumar/git/build/maxwell-conjugate-20260602/bin/ogs`
re-verified md5 `cea9d7d81972f732385b41a71e50f20e` (`-49-gbed3e395`, commit `bed3e395da`) —
unchanged from entry 28 / the 2026-08-28 cross-session run, no drift detected at check time
(this shared path is known to drift from concurrent activity in other sessions during this
session, so it was re-verified repeatedly, not just once). Runs executed from a fresh isolated
copy (`dd900_adopt_run_2026-08-31/{bin,lib}/`, full dylib closure, smoke-tested standalone)
immune to the shared build dir's drift.

**ctest: NOT USABLE, for either model attempted.** `maxwell-conjugate-20260602` was configured
with `OGS_BUILD_TESTING=OFF`; its `ProcessLib/RichardsMechanics/CTestTestfile.cmake` contains
ZERO `add_test()` calls, `ctest -N` at the build root lists only `ogs_no_args`, and
`ctest -R 'ANCHORS_MS33_ModelVII|modelVII_freeswelling|freeswelling'` returns "No tests were
found!!!" This is a configure-time gate that predates and is unrelated to today's K-table
change — no ANCHORS_MS33 ctest case is registered in this build at all, regardless of PRJ
content, and getting one would need a full reconfigure (not attempted: out of scope, and risks
the shared drift-prone build dir further). `ctestPassed = not_run` for III and VII for this
concrete reason, not merely as a fallback choice.

**Model III (gap-switch) — SUCCEEDED, adopted.** Fresh run (isolated copy) completed 888
accepted steps (1 routine rejected step at t=5.666 s startup transient, expected/documented
behavior) to t=17280000 s (200 d). Cross-validated against the independently recovered
2026-08-28 run (rep1/rep2, byte-identical to each other, same PRJ/K/binary/commit, three days
earlier, raw output still on disk from that prior session): EVERY point_data field (sigma,
pressure, porosity, transport_porosity, micro_porosity, micro_water_content, saturation,
dry_density_solid, displacement) over all 90 mesh nodes, plus raw point coordinates —
max_abs_diff = 0.0, max_rel_diff = 0.0 everywhere. This exceeds the PRJ's own declared
tolerances (sigma/pressure family rel 1e-2 abs 1e3; porosity family rel/abs ~1e-7-1e-8) — there
is no diff to even compare against tolerance; a stronger check than the task's own probe-only
bar. Probe values (mean_stress_MPa = -(sxx+syy+szz)/3/1e6): Top 2.013879209625015, Central
2.006541516052546, Bottom 1.9936228708076316 MPa, identical to all printed digits in
fresh/rep1/rep2. (Raw VTU container bytes/md5 differ, fresh `645c5308d3ac7c2d89c4309c7c947e61`
19494 B vs rep1/rep2 `fc1c1e452e2043aac157f3f24e3318b9`; confirmed non-physics header/
compression-metadata artifact, not a data discrepancy, by the full-field comparison above.)
Reference REBASELINED: old `ms33_modelIII_gapswitch_ts_880_t_17280000.000000.vtu` (md5
`4e3228d44dd5398f0d244823ccabe29b`, re-verified in this file's own archived copy) moved via
`git mv` to `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIII/superseded_references_2026-08-31/`
(history preserved, per the established rebaselining convention, e.g.
`superseded_references_2026-08-26/`); new
`ms33_modelIII_gapswitch_ts_888_t_17280000.000000.vtu` added as the live reference — exactly one
`ms33_modelIII_gapswitch_ts_*_t_17280000...` file now sits in the model directory (the
harness's "exactly one matching reference" requirement, confirmed).

**Model VII (free-swelling) — SUCCEEDED, adopted.** Fresh run completed ts_920, t=20736000 s
(240 d). Cross-validated against the recovered 2026-08-28 rep1 (rep1==rep2 confirmed first, md5
`d9e4ad0897f89e619281d966e8bd2756` both): 21 extracted quantities (e_top_probe, e_domain_mean,
phi_top, phi_mean, axial/mean stress and porosity/void-ratio at Top/Central/Bottom, suction,
min transport/micro porosity, mean saturation, mean micro water content) all matched EXACTLY
(abs_diff=0, rel_diff=0) to full double precision — e.g. e_top_probe = 1.1264739847014347 both
runs; mean_stress_Top/Central/Bottom_MPa = 0.1498690501 / 0.1269712678 / 0.1346723332 both runs;
axial_stress_top_probe_MPa = 0.4113224313, suction_top_probe_MPa = 18.01460414, both runs. Node
count (1066) and timestep count (ts_920) matched too, confirming identical adaptive-time-stepping
trajectory, not just an identical endpoint. (VTU container md5 differs, fresh
`f6c5ab0816df428c3a28fa5146f0bda6` 204180 B vs recovered `d9e4ad0897f89e619281d966e8bd2756`
204124 B; non-physics metadata, confirmed irrelevant by the full-quantity match.) Reference
REBASELINED: old `ms33_modelVII_freeswelling_ts_883_t_20736000.000000.vtu` moved via `git mv` to
`Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelVII/superseded_references_2026-08-31/`; new
`ms33_modelVII_freeswelling_ts_920_t_20736000.000000.vtu` added as the live reference (exactly
one matching file confirmed, non-recursive check).

**Model IV (pellets) — DID NOT COMPLETE. NOT rebaselined. NOT cross-validated. NOT adopted.**
The reporting pass for Model IV returned no completed run and no cross-validation — its status
object reads `crossValidated=false`, `ctestPassed=not_run`, `referenceInstalled=false`, with a
notes field containing only a stray placeholder string ("PLACEHOLDER - DO NOT USE... a
multi-hour background run proceeds... I will NOT actually end here") rather than a finished
result. Model IV is the long pole of this cascade (~3 h / ~21 500 timesteps observed in the
2026-08-28 recovered run), and no confirmed-complete, confirmed-cross-validated fresh run exists
for it as of this commit. Per instruction, Model IV's reference file and PRJ provenance comment
are left UNTOUCHED here: `ms33_modelIV_pellets.prj` keeps its pre-existing uncommitted
modification (K already live in its `<prefactors>` table, same as III/VII, from before this
task's run/verify phase) but receives NO provenance comment and is NOT part of this commit — it
remains a pending, unverified, uncommitted change in the worktree.
`Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV/` keeps its OLD reference
(`ms33_modelIV_pellets_ts_577_t_17280000.000000.vtu`) untouched. Do not treat Model IV as
adopted, cross-validated, or rebaselined — that still needs a fresh run that actually completes,
followed by the same cross-validation and rebaseline steps applied to III/VII above.

**Net effect of this entry.** III and VII are now on a verified, cross-session-reproducible
corrected-K reference — the strongest evidence bar obtainable here given ctest is not registered
in this build (exact bit-for-bit field agreement across two independent sessions three days
apart, not merely within-tolerance agreement). IV remains at entry 28's checkpoint state (K live
in the PRJ, nothing else moved, nothing committed) and needs its own completed, verified run
before it can join III/VII — do not cite Model IV headline numbers under the corrected K until
that happens. `ms33_modelI_dd900.prj` (entry 28's own re-fit, separately approved 2026-08-28) is
committed alongside as previously-approved content unrelated to whether the III/IV/VII cascade
is complete.

---

## 2026-08-31 — Model IV (pellets) completed: cascade entry 29 now closed for all three models

### 30. Model IV cross-validated + rebaselined; III/IV/VII corrected-K cascade fully adopted

Entry 29's prior attempt at Model IV filed a placeholder result instead of actually waiting for
the run to finish (nothing physically failed — the agent just never properly executed the wait).
This entry is the retry, done properly: the run was let run to completion and its completion was
confirmed directly from `run.log`, not inferred from elapsed time.

**Run.** Isolated binary `dd900_adopt_run_2026-08-31/{bin,lib}/ogs` (md5
`cea9d7d81972f732385b41a71e50f20e`, re-verified unchanged) against the staged
`ms33_modelIV_pellets.prj` (K-table `18469.144929 46000.0 104689.9129 265905.06`, byte-identical
to this worktree's own copy). Started 2026-08-31 09:23:34 CEST, completed 11:48:00 CEST
(2 h 24 m 26 s = 8666.15 s per OGS's own timer). `run.log` tail confirms completion directly:
"The whole computation of the time stepping took 21500 steps, in which the accepted steps are
21500, and the rejected steps are 0." / "Simulation completed. It took 8666.15 s." / "OGS
completed on 2026-08-31 11:48:00+0200." Process (PID 40112) confirmed exited after that line.
Final output: `ms33_modelIV_pellets_ts_21500_t_17280000.000000.vtu` (t=17280000 s = 200 d).

**Cross-validation (independently re-run this pass, not merely taken on trust).** Recovered
2026-08-28 double-replica first sanity-checked against itself: `IV_rep1` vs `IV_rep2` —
bit-identical mesh points (1066 nodes) and every point_data field bit-for-bit exact
(max_abs_diff = 0.0), confirming the recovered pair is internally consistent before using it as a
reference. Fresh run vs `IV_rep1` and vs `IV_rep2` independently: mesh points bit-identical to
both; all 14 point_data fields (displacement, dry_density_solid, intrinsic_permeability,
micro_exchange_source, micro_porosity, micro_pressure, micro_water_content, porosity, pressure,
relative_permeability, saturation, sigma, swelling_stress, transport_porosity) across all 1066
nodes bit-for-bit exact against both replicas — max_abs_diff = max_rel_diff = 0.0 everywhere, no
tolerance needed. This exceeds the PRJ's own declared tolerances (sigma/pressure family rel 1e-2
abs 1e3; porosity family rel/abs ~1e-7-1e-8) — matches the bar set by III and VII in entry 29.
Probe values (mean_stress_MPa = -(sxx+syy+szz)/3/1e6, identical fresh/rep1/rep2 to all printed
digits): Top (r=0, z=0.070) 4.3935488764 MPa, Central (r=0, z=0.040) 4.3618173505 MPa, Bottom
(r=0, z=0.010) 0.8924417770 MPa.

**Reference REBASELINED.** Old `ms33_modelIV_pellets_ts_577_t_17280000.000000.vtu` moved via
`git mv` to `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV/superseded_references_2026-08-31/`
(history preserved, same convention as III/VII in entry 29); new
`ms33_modelIV_pellets_ts_21500_t_17280000.000000.vtu` added as the live reference (md5
`53392b6a0308fc12d7f978f614066719`, confirmed identical to the run's own `out/` copy) — exactly
one `ms33_modelIV_pellets_ts_*_t_17280000...` file now sits in the model directory. Inline
provenance comment added at the PRJ's `<prefactors>` block (same wording pattern as III/VII),
now pointing at this entry rather than entry 29's III/VII-only text.

**Net effect of this entry.** All three shipping ANCHORS_MS33 models (III, IV, VII) are now
adopted onto the corrected dd900 K-knot (18469.144929 J/kg) with verified, rebaselined
regression-test references. Every one of the three used the same evidence bar: exact bit-for-bit
field agreement against an independently recovered run from a separate session three days
earlier, not merely within-tolerance agreement — the strongest reproducibility check available
given ctest is not registered in this build. The III/IV/VII corrected-K cascade opened in entry
29 is now complete. Deliverable-side propagation and the push to any remote remain separate,
not-yet-taken steps (this commit, like entry 29's, is local-only per instruction).

---

## 2026-08-31 — Doc correction: K_OF_RHO_D_LIVE.md aligned with the log-linear live-K scheme

### 31. DONE 2026-08-31 — K_OF_RHO_D_LIVE.md brought in line with commit 1bb414ac05 (doc-only)

Code-review finding #8: entry 24's landing (commit 1bb414ac05, 2026-08-26) updated this
worklog and the code but left `DSM/K_OF_RHO_D_LIVE.md` specifying the live value as the
K-linear `getValue` and the live Jacobian tangent as the per-segment-CONSTANT
`getSegmentSlope`, contrary to CLAUDE.md §8. Confirmed at HEAD 7ec39ecf4c: the doc's last
commit is b7d7e44c32 (2026-06-12), and 1bb414ac05 touched only AGENTS.md,
PotentialExchangeParameters.h and StrainedFilmPotential.cpp.

Fix, DOCUMENTATION ONLY — no source, test, PRJ or VTU file touched. A full sweep of the
document found five stale passages, not just the one the review cited: the Implementation
helper name plus its clamp attribution; the "LINEAR ... Vinay's call pending" provisional
bullet; both dated verification blocks (test names 1bb414ac05 renamed); and the
`dK/drho_d = getSegmentSlope` chain bullet. All five are annotated `[HISTORICAL ...]` in
place per §6.3 — nothing deleted — each pointing at a new closing section "Interpolation
scheme: LOG-LINEAR (Vinay's decision 2026-08-26, commit 1bb414ac05)", and a READ-FIRST
pointer to it was added at the top. That section states the shipped live value
`K(x) = K_l*exp(t*ln(K_r/K_l))` [J/kg], the shipped live tangent
`dK/drho_d = K(rho_d)*ln(K_r/K_l)/(x_r-x_l)` [(J/kg)/(kg/m^3)] and
`dK/dphi = -rho_SR*dK/drho_d` [J/kg per unit phi] — each re-read from
`PotentialExchangeParameters.h` and from both live-K Jacobian blocks in
`RichardsMechanicsFEM-impl.h` (`dK_dphi_pu`, `dK_dphi_sw`), not taken from the review —
and records what did NOT change (the K-linear pair still serving the parse-time frozen-K
path; `mu_aug` still linear in K). No numeric literal was introduced (§1.1 did not fire);
nothing was built or run for this entry, so no runtime or test-status claim is made (§5),
and the two blocks' historical pass counts are explicitly not re-asserted.

Also fixed while sweeping: "Analytic tangent completion" (2026-06-12) names only the p-u
Jacobian block; the new section records that the 2026-06-14 swelling-eigenstress block
carries the same live-K chain.

Left for whoever owns the source files: at HEAD the inline comments at
`PotentialExchangeParameters.h:425` and `RichardsMechanicsFEM-impl.h:5040` still write
dK/dphi as `-rho_SR*(table segment slope)`. Source edits were out of scope for this
doc-only entry.

**ANNOTATIONS 2026-08-31 (code-review round 2).** The entry text above is left exactly as
written (§6.3/§6.4); the three items below — two corrections and one closure — are appended
beside it, not edits to it.

- **Scope reconciled — "DOCUMENTATION ONLY" is true of this entry, false of the commit it
  will land in.** "Fix, DOCUMENTATION ONLY — no source, test, PRJ or VTU file touched"
  describes entry 31's own edits accurately. It does NOT describe the uncommitted tree those
  edits sit in. That tree (HEAD `7ec39ecf4c`, 12 modified files) also carries: the
  `PotentialExchangeParameters.h` refactor and its interior-knot node-preservation branch; a
  new parse-time `<prefactors>` positivity guard in `CreateRichardsMechanicsProcess.cpp`; the
  SCHEME comment in `RichardsMechanicsFEM-impl.h`; 13 new assertions plus one new test in
  `Tests/ProcessLib/RichardsMechanics/StrainedFilmPotential.cpp`; the Model IV `RUNTIME`
  change in `ProcessLib/RichardsMechanics/Tests.cmake`; and comment-only provenance
  annotations in five MS33 PRJs. None of that was in this worklog, which is the authoritative
  tracker for all of it (`Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV/AGENTS.md` is a
  stub pointing here). It is entered below as entries 32-37; entry 38 records what was
  deliberately NOT touched.

- **DONE 2026-08-31 — "Left for whoever owns the source files" is CLOSED.** Both stale inline
  comments are corrected in this working tree. Checked 2026-08-31 by reading HEAD and the
  working tree side by side: at `7ec39ecf4c` the wording `-rho_SR*(table segment slope)`
  stands at `PotentialExchangeParameters.h:425` and `RichardsMechanicsFEM-impl.h:5040`,
  exactly as this entry recorded; in the working tree that phrase occurs in no source file at
  all (`grep` over `ProcessLib/` and `Tests/` returns only this worklog — the 2026-06-12 line
  ~313, which is the correct record of that date and stays, and this entry's own sentence).
  The corrected text is the dK/dphi paragraph at `PotentialExchangeParameters.h:540-552` and
  the SCHEME note at `RichardsMechanicsFEM-impl.h:5053-5069`, with a pointer to it from the
  swelling-eigenstress block at 5241-5243; both now write
  dK/dphi = -rho_SR*K(rho_d)*ln(K_r/K_l)/(x_r-x_l), proportional to the local K rather than a
  per-segment constant. Recorded in entry 32.

- **"All five are annotated `[HISTORICAL ...]` in place — nothing deleted": four of the five
  were.** The fifth — the Implementation section's first-cut helper bullet plus its clamp
  attribution in `K_OF_RHO_D_LIVE.md` — was REWRITTEN in place instead: the sentence
  "`MathLib::PiecewiseLinearInterpolation::getValue` HOLDS the endpoint values outside
  `[rho_d_min, rho_d_max]` (verified in ...)" occurred nowhere in the file afterwards, and its
  "verified in" record had been softened to "read in". Confirmed 2026-08-31 against
  `git show HEAD:ProcessLib/RichardsMechanics/DSM/K_OF_RHO_D_LIVE.md`. Resolution taken: the
  §6.3-preferred one — RESTORE rather than weaken the claim. The first-cut bullet and clamp
  paragraph are back in the file byte-for-byte beside the new text, inside a
  `[HISTORICAL ...]` block, and "verified in" is restored as written, because the original did
  record a verification and that verification still holds: `getValue`'s two clamp branches
  were re-read in `MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.cpp` on
  2026-08-31 and are unchanged. With that restoration the sentence above is true of the file
  as it now stands; it was not true when written. Recorded in entry 37.

---

## 2026-08-31 — code-review rounds 1 and 2: the rest of the uncommitted batch entered in the worklog (§8 gap closed)

Entries 32-38 close the §8 gap named in the annotations to entry 31. Everything below is
UNCOMMITTED working-tree state on top of HEAD `7ec39ecf4c` (which is itself already pushed to
three remotes): nothing here was committed, staged or pushed. **No MS33 simulation and no MS33
ctest was run in either round** — where a consequence for the III/IV/VII reference VTUs is
stated below it is PREDICTED from a bitwise-equivalence sweep of the table accessors, never
verified by a run (§5).

### 32. DONE 2026-08-31 — AugmentationPrefactorTable put on one segment lookup; interior-knot node preservation made bit-exact

`ProcessLib/RichardsMechanics/PotentialExchangeParameters.h`. The clamp + `lower_bound`
interval selection that stood copied out in three accessors is now one private helper,
`locateSegment()` (returning `std::nullopt` exactly where `getValue()` would return a clamped
endpoint), with a `Segment` struct carrying `x_l, x_r, K_l, K_r`, and a private static
`logLinearValueOnSegment(Segment, t, r)`. `getSegmentSlope()`, `getValueLogLinear()` and
`getSegmentSlopeLogLinear()` all route through them; `getSegmentSlopeLogLinear()` now
evaluates its segment ONCE (one lookup, one `std::log`, one `std::exp`) instead of re-entering
`getValueLogLinear()`. `getValue()` itself is untouched — it lives in the MathLib base class,
is non-virtual, and is neither overridden nor shadowed.

**One deliberate behaviour change: node preservation at an INTERIOR knot.**
`logLinearValueOnSegment()` returns the STORED knot value at `t == 1` (and `t == 0`) rather
than evaluating `K_l*exp(t*ln(K_r/K_l))` there. `lower_bound` places an exact interior knot in
the LEFT segment at `t == 1`, where that round-trip can land ~1 ULP off `K_r`; a boundary knot
is unaffected either way, since it returns through the endpoint-hold clamp, which contains no
logarithm.

**Equivalence evidence (MEASURED this session, quoted as measured).** An old-vs-new sweep of
the accessors over five tables — the shipped 4-knot table, the superseded 4-knot table, the
shipped 2-knot table and both unit-test tables — at ~95k-102k samples each: ZERO bitwise
differences in `getSegmentSlope`, `getValueLogLinear` and `getSegmentSlopeLogLinear`. The only
change anywhere in the sweep is on the SUPERSEDED `4367.2277` table:
`getValueLogLinear(1400)` `45999.999999999993` -> `46000`, and `getSegmentSlopeLogLinear(1400)`
`216.61519437040712` -> `...715` (quoted verbatim from the sweep report). An independent audit
rebuilt the comparison in a standalone harness and confirmed `getSegmentSlope` bit-identical
old-vs-new at every probe on every table, including 2-knot and degenerate 1-knot cases.

**The defect the branch removes was LATENT, not live.** The table shipped at `7ec39ecf4c`
round-trips bit-exactly at all four of its knots (measured 2026-08-31), so no shipped run was
ever affected. Entry 24's landing-audit sentence "1400 low by 1 ULP, rel 1.6e-16" describes
the PRE-cascade table — the superseded 900-knot pair — and stands as that record.

**Comment corrections in the same file** (no code affected): the class rationale now names
both value/slope pairs and says which one the LIVE path uses, keeping its original
`getDerivative`-blending origin story as an explicitly labelled HISTORICAL paragraph; the
`dK/dphi` block is corrected as recorded in the entry-31 annotation above, and now also names
BOTH Jacobian call sites rather than one (`RichardsMechanicsFEM-impl.h` lines 5052, the p-u
augmentation exchange tangent, and 5223, the displacement-side swelling-eigenstress tangent,
line numbers at `7ec39ecf4c`); the log-linear pair carries the strictly-positive precondition
(entry 33); and the call-site census is recorded — `getValue()` on this table has exactly ONE
production caller, the parse-time frozen-K resolution in `CreateRichardsMechanicsProcess.cpp`,
while `getSegmentSlope()` has NONE and is exercised only by
`Tests/ProcessLib/RichardsMechanics/StrainedFilmPotential.cpp`. That census was re-grepped
over `ProcessLib/` and `Tests/ProcessLib/` on 2026-08-31 for this entry and matches.
`RichardsMechanicsFEM-impl.h` gained the SCHEME note at the live-K p-u block and a
form/units pointer at the swelling-eigenstress block.

**Consequence for the III/IV/VII ctest references: PREDICTED IDENTICAL, NOT VERIFIED.** The
bitwise sweep is evidence that the shipped tables evaluate identically old-vs-new, which
predicts the rebaselined reference VTUs of entries 29/30 are unaffected. No MS33 run was
performed to confirm it, in either round (§5).

### 33. DONE 2026-08-31 — parse-time strictly-positive `<prefactors>` guard; this WIDENS parser rejection. ASK(Vinay) on its scope

`ProcessLib/RichardsMechanics/CreateRichardsMechanicsProcess.cpp`, in
`parsePotentialExchangeParameters`: after the existing size check and before the table is
constructed — while the lists are still addressable by their PRJ index — every `<prefactors>`
entry is required to be `> 0`, and the first entry that is not raises `OGS_FATAL` naming that
entry's index and its value and stating that the log-linear interpolant evaluates
`ln(K_r/K_l)`.

Reason, written out case by case in the source comment: the non-positive cases do NOT share
one failure mode. `K_l == 0` gives `+inf` out of the logarithm and a NaN only one step later
at `0*inf`; `K_l < 0` is the only case the logarithm itself catches; `K_r == 0` returns a
clean, plausible `0` for the VALUE while poisoning ONLY the companion slope — the silent case,
and the dangerous one; and two knots of the same negative sign produce no NaN anywhere and
would carry a negative K straight through. Those chains were PROBED on a standalone
transcription of `logLinearValueOnSegment` (IEEE-754 double, 2026-08-31), not on the shipped
object, and the source comment says so.

**This widens what the parser rejects**, and deliberately so: the check is not restricted to
live mode, so a frozen-K deck carrying a zero knot used to parse and now aborts. Re-checked
for this entry on 2026-08-31: all seven `<prefactors>` tables under
`Tests/Data/RichardsMechanics/` are strictly positive (MS33 III/IV/VII gapswitch / pellets /
freeswelling, plus the four `*_kofdd` / `*_livek` decks), so no deck in the tree changes
behaviour. Not verified by execution: no deck with a non-positive knot was constructed or run,
so "such a deck now aborts at parse time" is predicted from the code, not measured.

No `<dry_densities>` monotonicity check was added and none is needed: the MathLib base
constructor sorts the support points (carrying their values with them) and `OGS_FATAL`s on
duplicate support points (`MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.cpp`,
re-read 2026-08-31).

**ASK(Vinay), non-blocking (guardrail §9 row "Physics / model formulation decision" — a new
rejection rule is a modelling-policy call, not a coding one).** Should the guard stay
table-wide, or be restricted to live mode? Table-wide is what is implemented, on the argument
that no tree deck uses a non-positive prefactor and the K-linear path has no legitimate use
for one either; but only the LOG-linear path actually requires positivity, so a future
frozen-K deck wanting a zero knot would be blocked by this. Default action if unanswered:
leave as implemented.

### 34. DONE 2026-08-31 — unit tests: 13 assertions added, one new test; what the recorded testrunner run does and does not cover

`Tests/ProcessLib/RichardsMechanics/StrainedFilmPotential.cpp`.

Added inside the existing `RichardsMechanicsLiveKOfRhoD.AnalyticPhiTangentClampedEdgesAndKnots`
(9 assertions): four `EXPECT_DOUBLE_EQ` pinning `getSegmentSlopeLogLinear` = 0 at both
boundary knots and at two strictly-outside points (a regression narrowing the `<=`/`>=` clamp
to `<`/`>` fails here); three `EXPECT_EQ` on `getValueLogLinear` at the three structural knots;
and the LEFT-segment one-sided convention as an `EXPECT_NEAR` against the in-file left-segment
value, PLUS an `EXPECT_GT` asserting the right-segment candidate lies far outside the same
tolerance — that pair is what makes the one-sided assertion discriminating against a
`lower_bound` -> `upper_bound` regression, which changes the answer only AT a knot.

New test (4 assertions): `RichardsMechanicsLiveKOfRhoD.InteriorKnotBitExactWhereRoundTripMisses`.
Why it exists: the round-1 comment claimed the three node-preservation `EXPECT_EQ`s made node
preservation bit-exact, but on that test's own knots {1000,1500,2000}/{10,16,30} the round-trip
`10.0*exp(1.0*log(16.0/10.0))` equals `16.0` EXACTLY, so those assertions hold WITH OR WITHOUT
the `t == 1` branch and pin nothing. The comment is corrected in place to say so, the three
`EXPECT_EQ` lines are byte-identical to before (§3: supplement, never replace), and the new
test carries the discriminating case on a knot pair whose round-trip does miss.

Physics anchor (§3a): analytical limit / interpolation identity — node preservation
K(x_i) = K_i. Literals cited per §1.1 item 3 (prior commits traceable to user-approved work):
K(900) = 4367.227700212952 J/kg is the SUPERSEDED dd900 calibration value, on record in the
provenance block of `ANCHORS_MS33_ModelI/ms33_modelI_dd900.prj` and in commit `642a8f867a`;
K(1400) = 46000.0 and K(1600) = 104689.9129 J/kg are shipped table entries at `7ec39ecf4c`.
FLAG for the record: the knots are STRUCTURAL and no physics is asserted on them — the
superseded value is used precisely because its segment round-trip misses, whereas the shipped
table's knots all round-trip exactly and would make the test non-discriminating.

**Build and test status (MEASURED this session):** build PASS in
`/Users/vinaykumar/git/build/maxwell_floor_20260619`; `./bin/testrunner
--gtest_filter='RichardsMechanics*'` gave 45 tests, 43 PASSED, 2 pre-existing `GTEST_SKIP`s
(`ExactFilmEnergyPair`, TODO(Vinay) Q3/Q4 and T-8), 0 failures; no new compiler warnings.

**What that run covers, checked 2026-08-31 rather than assumed.** The `testrunner` binary
(mtime 17:17) contains the assertion strings added inside
`AnalyticPhiTangentClampedEdgesAndKnots` — `table.getSegmentSlopeLogLinear(500.0)`,
`table.getSegmentSlopeLogLinear(2600.0)`, `K_interior_knot`, `expected_left_loglinear_slope` —
so those 9 assertions were compiled in and are covered by the 43 passes. It contains NO
`InteriorKnotBitExact` symbol, while it does contain the neighbouring test names, and the
source (mtime 17:47) now declares 46 `TEST(RichardsMechanics*)` cases against the run's 45.
The recorded run therefore covers the batch MINUS the new test. **The new test's pass is
PREDICTED, NOT VERIFIED (§5)** — it has not been through the built suite. What was measured
for it, as reported by the round-2 test pass, is different evidence: its two bit-exactness
assertions were run against `getValueLogLinear`/`locateSegment`/`logLinearValueOnSegment` extracted verbatim into a
standalone program and compiled twice, with and without the `t == 1` branch (fail without,
pass with), plus a `-fsyntax-only` compile against the real header.

> **UPDATE 2026-08-31, later the same session (paragraph above kept verbatim per §6.3/§6.4).**
> Both load-bearing observations in it have been overtaken, and the caveat they supported is
> now retired by measurement:
>
> - The mtime-17:17 binary it inspected no longer exists in that state. Entry 35's
>   reconfigure-and-list step relinked `bin/testrunner` at **2026-08-31 18:06:37**, and
>   `strings` on the current binary DOES find
>   `RichardsMechanicsLiveKOfRhoD_InteriorKnotBitExactWhereRoundTripMisses_Test`. Verified by
>   `ls -la` and `strings` this session. So "It contains NO `InteriorKnotBitExact` symbol" was
>   true of the 17:17 binary and is false of the tree now.
> - **The new test's pass is now VERIFIED, not predicted.** `./bin/testrunner
>   --gtest_filter='RichardsMechanics*'` was re-run against the 18:06 binary (Release, Apple
>   clang 21.0.0): **46 tests, 44 PASSED, 2 pre-existing `GTEST_SKIP`s, 0 failures**, with
>   `[       OK ] RichardsMechanicsLiveKOfRhoD.InteriorKnotBitExactWhereRoundTripMisses`
>   in the output. The 43-pass figure recorded above is the round-1 count, before this test
>   existed.
>
> What remains harness-only, and is worded that way in the test file, is the FAILING half: the
> pre-fix accessor going red. Reproducing that in-binary would mean removing the branch from
> the shipped header, so it is measured in the standalone harness — where it fails both
> interior-knot assertions at rel −1.58e-16 under Apple clang 21.0.0 and Homebrew clang
> 22.1.8, at −O0 and −O2 each — and never claimed as a testrunner observation.

### 35. DONE 2026-08-31 — Tests.cmake: Model IV RUNTIME 240 -> 10736, taken from the slowest of three measured runs

`ProcessLib/RichardsMechanics/Tests.cmake`. Three runs of the identical deck
(`ms33_modelIV_pellets.prj`), 21500 accepted / 0 rejected steps each, OGS's own timer:

| wall time | date, threads | binary | log |
| :-- | :-- | :-- | :-- |
| 10720 s | 2026-08-28 09:02:27+0200, `OMP_NUM_THREADS=4` | ogs 6.5.8-565-gbea47887 | `cascade_refit/IV_rep1/run.log` |
| 10735.4 s | 2026-08-28 09:02:27+0200, `OMP_NUM_THREADS=4` | ogs 6.5.8-565-gbea47887 | `cascade_refit/IV_rep2/run.log` |
| 8666.15 s | 2026-08-31 09:23:34+0200, `OMP_NUM_THREADS=6` | `archive/dsm_native_Pi_fofnlev_branchtip_2026-08-11-49-gbed3e395` | `dd900_adopt_run_2026-08-31/IV/out/run.log` (the run quoted in entry 30) |

Log roots: the two 2026-08-28 replicas sit in a session-local scratchpad
(`.../2934957f-455c-4c8c-822d-0822f7e86487/scratchpad/cascade_refit/IV_rep{1,2}/run.log`, not
durable), the 2026-08-31 run under `/Users/vinaykumar/ogs-models/dd900_adopt_run_2026-08-31/`.

The spread tracks the thread count, not the deck: the 6-thread run is the FASTEST configuration
on record, not the representative one, so `RUNTIME` is taken from the slowest,
10736 = ceil(10735.4). Rule re-read for this entry in `scripts/cmake/test/OgsTest.cmake`: an
explicit `TIMEOUT` is emitted only when `RUNTIME > 750`, and then `TIMEOUT = 2*RUNTIME` — so
10736 gives 21472 s (2x the worst case on record), where 8667 would have given 17334 s, only
1.61x it. No test is renamed: `ogs.ctest.large_runtime = 60` (`web/data/versions.json`), so 240
and 10736 are both above it and the case was already `-LARGE`; and any `OGS_CTEST_MAX_RUNTIME`
below 240 dropped the deck before the change and still does. (Completed 2026-08-31: that last
clause is true but not the whole gate. A cap in **[240, 10735]** kept the deck before the change
and now drops it. No in-tree consumer sits in that band — `scripts/ci/jobs/build-linux.yml:44`
sets 60 and `CMakePresets.json` sets 101, both below 240 — so nothing in the repo changes
behaviour, but a caller passing a cap in that band would see the deck disappear.)

CONDITIONAL, never observed — stated as such in the file: at the old `RUNTIME 240` no `TIMEOUT`
is emitted at all, so ctest's default 1500 s WOULD apply and WOULD kill a healthy run of this
deck. No ctest run of this deck exists to have observed it: `ANCHORS_MS33` is not registered in
the build used for the cascade runs (entries 29 and 30).

> **CORRECTION 2026-08-31 (kept above verbatim per §6.3).** "No ctest run of this deck exists"
> is too strong and is falsifiable: one is on record at
> `/Users/vinaykumar/git/build/maxwell_floor_20260619/logs/ogs-RichardsMechanics_ANCHORS_MS33_ModelIV_ms33_modelIV_pellets-LARGE-omp.txt`
> — `OGS started on 2026-08-26 16:11:13+0200`, completed 16:12:33 (~80 s), vtkdiffing against
> the then-current `ms33_modelIV_pellets_ts_577_t_17280000.000000.vtu`. Verified by reading the
> log this session. The reason clause is still correct (`ogs-dsm-active` carries
> `OGS_BUILD_TESTING:BOOL=OFF`), and the substantive point survives unchanged: that run never
> approached 1500 s because it was the PRE-CASCADE 577-step deck. The accurate claim is that no
> ctest run of this deck **in its present 21500-step form** exists.

MEASURED this session, and it is a listing, not a run: `ctest -N --show-only=json-v1` on the
reconfigured tree showed the Model IV test picking up an explicit `TIMEOUT` from the raised
`RUNTIME`, with the `-LARGE` suffix retained.

Models I, III and VII keep their RUNTIMEs (120/120/120, 120, 300). III and VII were re-timed on
the 2026-08-31 adopt runs — III 10.1542 s / 889 steps at `OMP_NUM_THREADS=4`, VII 207.022 s /
920 steps with `OMP_NUM_THREADS` unset (18-thread fallback), so VII's margin to RUNTIME 300 was
measured at 18 threads, not at one. Model I was NOT re-timed: its 0.111132 / 0.118636 /
0.149283 s are 2026-05-22 runs of the stale binary `vdw-baseline-2026-05-08-41-g9a1b956c`,
verbatim from the tracked `rerun_ms33_modelI_dd*.log`, kept as the only timings on record
rather than as current measurements. All of this is written into the file's comments; CI
headroom is labelled there as expected, not verified.

### 36. DONE 2026-08-31 — PRJ provenance annotations in five MS33 decks (review findings #7 and #10), comment-only

Finding numbers are as reported to this session; the audit documents themselves are not in the
tree.

Five decks annotated, no run performed: `ANCHORS_MS33_ModelI/ms33_modelI_dd900.prj`,
`ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj` (registration commented out in Tests.cmake,
deck retained), `ANCHORS_MS33_ModelIII/ms33_modelIII_gapswitch.prj`,
`ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj`,
`ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj`.

Content: dd900's three "CASCADE, predicted not yet verified" notes are kept verbatim as the
2026-08-28 record and each is followed by a CASCADE DONE 2026-08-31 note pointing at commits
`2c066ce187` (III, VII) and `7ec39ecf4c` (IV) and at entries 29/30; its sigma0 comment — which
names the SUPERSEDED K because it was written the day before the re-fit — is annotated STILL
VALID with the sources for "the K now live was also fit under sigma0=0" quoted from this file's
own header and the §12.2 block. The III/IV/VII reference blocks record the supersessions
ts_880 -> ts_888, ts_577 -> ts_21500 and ts_883 -> ts_920 with md5s, and mark the previous
reference as itself superseded and moved to `superseded_references_2026-08-31/`.

Verified for this entry on 2026-08-31, not taken on trust:
- the three live reference md5s recomputed from the committed files — gapswitch ts_888
  `645c5308d3ac7c2d89c4309c7c947e61`, pellets ts_21500 `53392b6a0308fc12d7f978f614066719`,
  freeswelling ts_920 `f6c5ab0816df428c3a28fa5146f0bda6` — match the PRJ text and entries
  29/30;
- all five edits are COMMENT-ONLY: with XML comments stripped, each file is byte-identical to
  its HEAD version (equal md5 per file), so no value, tag, mesh name or `<test_definition>`
  changed;
- all five parse (`xmllint --noout` clean on each), which also clears the `--`-inside-an-XML-
  comment hazard that a hand-written annotation block can introduce.

### 37. DONE 2026-08-31 — K_OF_RHO_D_LIVE.md, round-2 corrections

Four corrections in `DSM/K_OF_RHO_D_LIVE.md`, all annotate-in-place or restore, none delete:
(i) the first-cut Implementation bullet and clamp paragraph RESTORED byte-for-byte inside a
`[HISTORICAL ...]` block beside the current text, with the original "verified in" wording
restored and re-checked (entry-31 annotation above); (ii) the closing section's
"node-preserving BY CONSTRUCTION" sentence kept as the 2026-08-26 wording and followed by a
`[CORRECTED 2026-08-31]` block attributing node preservation to the explicit `t == 1` / `t == 0`
branch instead, with the superseded-table ULP miss and the latent-not-live finding recorded;
(iii) a `**Parse-time precondition on <prefactors>**` paragraph in that same normative section,
plus a `<prefactors>` PRJ-AUTHORING RULE bullet under "PRJ interface" (strictly positive,
widens parser rejection, all seven shipped tables comply, no `<dry_densities>` monotonicity
rule needed); (iv) the "Unchanged by this decision" sentence "still calls them" annotated with
the re-grepped call-site census, since `getSegmentSlope()` has no production caller. A second
READ-FIRST paragraph at the top points at (ii) and (iii) and states that neither has been
exercised by an MS33 run or by ctest.

### 38. OPEN(Vinay) 2026-08-31 — held back from both rounds, pending a scientific ruling

Ring-fenced for the user's scientific decision and deliberately untouched by rounds 1 and 2,
including by any workaround:

- the `DoubleStructureBenchmark` reference values;
- the adopted dd900 900-knot K and every `micro_solid_volume_fraction_reference` (n_s) value;
- the LIVE-mode fallback-K assignment in `CreateRichardsMechanicsProcess.cpp` — the
  `dry_density ? table->getValue(*dry_density) : (defaults ? ... : 0.0)` block, working-tree
  lines 621-626 after entry 33's guard was inserted above it. Its own comment is unedited; it
  is only referred to from the header's call-site census (entry 32).

Also open: the ASK(Vinay) in entry 33 (whether the positivity guard stays table-wide or is
restricted to live mode), and the batch itself — entries 31-37 are uncommitted working-tree
state; committing, and any push, are Vinay's calls.

### 39. DONE 2026-08-31 — round-3 corrections: four writing defects the final guardrail audit found in rounds 1-2

Applied directly, no agent fan-out (four localised edits). Each was verified against the tree
before editing, not taken on the auditor's word.

| # | Where | Defect | Correction |
| :-- | :-- | :-- | :-- |
| 1 | `DSM/AGENTS.md` entry 34 | Its binary evidence had been invalidated *within the same round*, and its `§5` caveat was stale | `UPDATE` block appended (original kept verbatim, §6.3/§6.4) |
| 2 | `DSM/AGENTS.md` entry 35 | "No ctest run of this deck exists" — falsifiable and false | `CORRECTION` block appended |
| 3 | `Tests.cmake` + entry 35 | `OGS_CTEST_MAX_RUNTIME` clause true but incomplete | band `[240, 10735]` spelled out |
| 4 | `StrainedFilmPotential.cpp` (2 comments) | Named only a toolchain that does *not* build the test | both toolchains named; stale caveat retired |

**(1) Entry 34's evidence, re-measured.** It reasoned from a `testrunner` of mtime 17:17 that
contained no `InteriorKnotBitExact` symbol. Entry 35's own reconfigure relinked the binary at
**2026-08-31 18:06:37**, and `strings` on it now finds
`RichardsMechanicsLiveKOfRhoD_InteriorKnotBitExactWhereRoundTripMisses_Test` (3 hits). More
importantly the caveat it supported is **retired by measurement**: `./bin/testrunner
--gtest_filter='RichardsMechanics*'` against that binary gives **46 tests, 44 PASSED, 2
pre-existing `GTEST_SKIP`s, 0 failures**, with the new test reported `OK`. The 43-pass figure
in entry 34 is the round-1 count, before the test existed. What stays harness-only, and is
worded that way in both the entry and the test file, is the FAILING half — the pre-fix
accessor going red — since reproducing it in-binary would mean removing the branch from the
shipped header.

**(2) The false ctest claim.** Entry 35 stated "No ctest run of this deck exists to have
observed it". One does, and it was found by reading, not inference:
`build/maxwell_floor_20260619/logs/ogs-RichardsMechanics_ANCHORS_MS33_ModelIV_ms33_modelIV_pellets-LARGE-omp.txt`
— `OGS started on 2026-08-26 16:11:13+0200`, completed 16:12:33 (~80 s), vtkdiffing against the
then-current `ts_577` reference. The reason clause survives (`ogs-dsm-active` carries
`OGS_BUILD_TESTING:BOOL=OFF`) and so does the substantive point — that run never approached
1500 s because it was the pre-cascade 577-step deck. The accurate claim is that no ctest run of
this deck **in its present 21500-step form** exists.

**(3) The gate band.** "Any `OGS_CTEST_MAX_RUNTIME` below 240 dropped the deck before and still
does" is true but not the whole gate: a cap in **[240, 10735]** kept the deck before the change
and now drops it. No in-tree consumer sits in that band (`build-linux.yml:44` sets 60,
`CMakePresets.json` sets 101, both below 240), so nothing in the repo changes behaviour.

**(4) Toolchain naming.** Two comments credited the standalone measurements to "Homebrew clang
22.1.8" alone. That compiler is real and the verdict reproduces under it, but both build dirs
use `CMAKE_CXX_COMPILER=/usr/bin/c++`, **Apple clang 21.0.0**. Both are now named, and the
discriminating verdict is on record under each at `-O0` and `-O2`.

**Re-verified after these edits, MEASURED:** `ninja testrunner` PASS; 46 tests, 44 passed, 0
failures; `RichardsMechanicsFEM-impl.h` still comment-only (comment-stripped code byte-identical
to `7ec39ecf4c`); all five MS33 PRJs still byte-identical to HEAD outside XML comments and
`xmllint --noout` clean. Batch total: 12 files, +1170/-72.

### 40. OPEN — register of everything held back, pending Vinay's discussion

Written as OPEN at Vinay's instruction (2026-08-31, "write open items as open, pending my
discussion"). Nothing below has been acted on. Entry 38's ring-fence list is the subset O1
covers; this entry is the full register. **A future agent must not close any item here without
Vinay's explicit ruling** — several are exactly the kind of decision CLAUDE.md §9 routes to him.

**O1 — the three code-review findings ring-fenced from all three rounds (physics/expected
values).** Raised by the `/code-review` pass on `HEAD~6..HEAD`, never adjudicated:

- **`transport_porosity = 0` baselined at t=0.** `double_porosity_swelling_RM.prj:112-117`
  declares `TransportPorosityFromMassBalance` with `initial_porosity = phi_tr0 = 0.3`, but the
  reference refreshed in `bea47887ac` carries 0 at t=0 and at every later step. t=0 is the
  initial condition, so no run-to-run drift can move it. The same refresh moved pressure
  -263.3 -> -169.6 MPa and swelling_stress 4.72 -> 2.37 MPa at t=100000. `transport_porosity`
  is not in the `<field>` list, so the ctest passes while certifying that state. Either the
  macro-porosity initialisation path is broken on this branch or the baseline is right and the
  PRJ is stale — **an expected-value question (§3), Vinay's to settle.**
- **The 900-knot K and n_s.** The re-fit `K = 18469.144929` was calibrated under
  `n_s = 0.3237`, but `gapswitch.prj:175` and `freeswelling.prj:165` declare
  `micro_solid_volume_fraction_reference = 0.5755395683453237`, and `pellets.prj` declares
  0.5755 for the block (:144) and 0.3237 for the pellets (:214) while sharing ONE table.
  dd900's own header records (MEASURED) that n_s = 0.575540 reproduces Ps = 0.35005 where
  n_s = 0.3237 gives Ps = 0.0832 for the same K — a 4x sensitivity. If the reading holds, the
  adopted knot is matched to the other n_s for III, VII and IV's block zone, and the superseded
  4367.2277 was the one fitted under 0.5755. **The III/IV/VII headline numbers rest on this.**
  Options a ruling must choose between: a per-n_s table split, a re-fit, or a deliberate shared
  knot. Related: entry 27's original OPEN(Vinay).
- **The LIVE-mode parse-time fallback K.** The review claims that at
  `CreateRichardsMechanicsProcess.cpp` (working-tree :621-626) the fallback still resolves via
  the K-linear `getValue()` while the run-time path uses `getValueLogLinear()` (~10% apart at
  rho_d = 1150 on the shipped table), and that because no MS33 live deck declares
  `<dry_density>` or a scalar prefactor, `dry_density` is `nullopt` and the fallback collapses
  to `0.0` — silently disabling augmentation at every phi-less site, with the `>= 0.0` guard
  accepting it. **NOT INDEPENDENTLY VERIFIED, and it is in tension with the decks producing
  non-trivial Ps.** Verify before acting; do not "fix" it on the review's word.

**O2 — the stray heredoc terminator in `7ec39ecf4c`'s commit message.** Its body ends
`...per instruction.\nEOF\n)\n` — leaked from the shell that wrote it. The message also says
"Not pushed — local commit only", which was already untrue: the commit is on `origin`, `backup`
and `vgk2`. Removing it needs `git commit --amend` plus a force-push to three remotes. **Not
done: rewriting published history is Vinay's call, not an agent's.**

**O3 — `RUNTIME 10736` for `ms33_modelIV_pellets` (entry 35) is an agent's judgement.** The
three measured runs are 10720 / 10735.4 / 8666.15 s and the spread tracks `OMP_NUM_THREADS`
(4/4/6). RUNTIME was taken from the slowest so the repo's `TIMEOUT = 2*RUNTIME` convention
holds against the worst case on record; the alternative — track the 6-thread configuration at
8667 — gives 1.61x. It is a CI timing hint, not a §1.1 literal, so no guardrail fired, but the
test suite is Vinay's. **Overrule freely.**

**O4 — scope of the parse-time positivity guard (entry 33's ASK).** It currently rejects any
non-positive `<prefactors>` entry table-wide, which WIDENS parser rejection: a frozen-K deck
carrying a zero knot used to parse and now hard-fails. All 7 shipped tables are strictly
positive, so no deck in the tree is affected. Open: keep it table-wide, or restrict it to live
mode where `ln(K_r/K_l)` is actually evaluated. Related: the guard's negative branch is
exercised by **no test** — the three unit-test tables construct `AugmentationPrefactorTable`
directly and bypass the parser. A death test would need a config tree; not attempted.

**O5 — two dated review snapshots still assert the K-linear live chain.**
`DSM/ULTRACODE_REVIEW_2026-06-14.md:97,125,126` ("`dK/dφ = ρ_SR·getSegmentSlope` is FD-verified
to 1e-9") and `DSM/ULTRA_REVIEW_2026-06-19.md:49`. Left untouched: both filenames are
date-stamped audit snapshots, so they arguably read as historical by construction, and §6.3
annotation of a snapshot is a different act from annotating a living spec. **Ruling wanted:**
tag them `[SUPERSEDED 2026-08-26 by 1bb414ac05]`, or leave snapshots as snapshots.

**O6 — nothing was re-run. `PREDICTED, NOT VERIFIED` (§5).** No MS33 ctest and no MS33
simulation was executed in any of the three rounds. The numerical identity of the III / IV / VII
reference VTUs under entry 32's refactor is *predicted* from the bitwise equivalence sweep
(five tables, ~95k-102k samples each, zero differences outside the intended interior-knot fix
on the SUPERSEDED table) and from the shipped table already being bit-exact at all four knots —
**it is not confirmed by a run.** The cheap confirmation, if wanted, is one `ms33_modelIII_gapswitch`
run (10.15 s measured) plus `ms33_modelVII_freeswelling` (207.02 s); Model IV costs ~3 h.

**O7 — build-directory routing, so no future agent cites the wrong binary.** MEASURED
2026-08-31: `build/ogs-dsm-active` is a **symlink to `build/maxwell_rebased_2026-08-26`**; its
`CMAKE_HOME_DIRECTORY` correctly points at this worktree, but it carries
`OGS_BUILD_TESTING:BOOL=OFF`, so it has **no `testrunner` target and registers zero MS33
ctests** (`ctest -N` lists only `ogs_no_args`). Unit tests must be built and run in
`build/maxwell_floor_20260619` (same worktree, `OGS_BUILD_TESTING=ON`). Worse for provenance:
a freshly rebuilt `ogs-dsm-active/bin/ogs` still self-reports `6.5.8-565-gbea47887`, stamped at
configure time from `bea47887ac` — **one commit behind HEAD. Never cite that version string as
evidence of which commit a run used.** Compare `reference_lit_mcp_ogs_source` /
`feedback_deprecated_branches_no_sims`.

**Resolved from entry 38:** its closing sentence ("committing, and any push, are Vinay's calls")
is CLOSED 2026-08-31 — Vinay instructed both. Everything else in entry 38 remains open under O1.

## 2026-09-01 — deliverable moved to the post-cascade generation; TWO NEW PHYSICS FINDINGS on the live-K path — OPEN(Vinay)

### 41. DONE 2026-09-01 — mechanical batch shipped (deliverable currency, workbooks, headers, snapshot)

Vinay: "finish all mechanical agentic stuff" / "commit and push the batch".
- Deliverable moved to the post-cascade set (eurad-anchors `c310f65`): III
  2.0138792096/2.0065415161/1.9936228708 MPa, IV 4.3935488764/4.3618173505/
  0.8924417770 MPa, VII e_top 1.1264739847 / e_mean 1.1259356770. Printed K table
  corrected to the live one. Three docs recompiled clean (30/21/29 pp, zero
  overfull boxes vs a HEAD baseline). New record CANONICAL_RESULTS_2026-08-31_cascade.json
  (md5 bf42cd05, deterministic); figures md5 f0f15db1/37621b12/361bdfe7.
- Model_VII_BGR.xlsx was internally mixing generations (log-linear stress rows,
  sigma0=0-era porosity/dry-density rows and porosity series, -8.9% vs BGR's own
  text). Fixed. III and IV measured clean on those rows.
- Retention provenance headers corrected in all three shipping decks (this repo,
  `bed5b76df4`): they claimed van Genuchten P0=27 MPa / a=0.45 where the live
  model is SaturationTuller and those literals occur zero times. Comment-only,
  verified; four Tuller literals remain TODO(Vinay) for a §12.1-admissible source.
- §6.8 snapshot `2026-08-31_1155_..._last_successful` (eurad-anchors `9d06499`),
  08-26 folder renamed back to `_successful` with handover note.

CORRECTION to this session's own earlier reporting: an audit lane asserted the
BGR workbooks were "already with the organizer"; a second lane labelled that
INFERRED; the adjudication dropped the label and it was reported as fact.
MEASURED afterwards: `tex/eurad2_MS34/.../BGR_DATA/Model_{III,IV,VII}/` are all
EMPTY (only Model_I and Reference_Configuration hold files). The empty folders
are CONFIRMED INTENTIONAL (Vinay, relayed via session task13-cf, 2026-08-31) --
NOT a submission gap, not an open item. The VII workbook fix stands on its own
merits regardless of delivery timing. Recording the escalation pattern because it
is the E2/E3 failure mode: a plausible claim gaining confidence at each hop.

### 42. OPEN(Vinay) 2026-09-01 — the live-K path fails the SAME energy test that killed CurrentTotalSolidFraction

Supersedes the framing of entry 40 O1: the defect is larger than "the 900-knot
was fit under the wrong n_s".

**Finding A (MEASURED).** `PotentialExchangeParameters.h:574-586` looks K up at
`micro_solid_density_reference * (1.0 - phi)` with `phi` = current TOTAL porosity
-- numerically the exact quantity `CurrentTotalSolidFraction` would have used for
nS, and which was REJECTED 2026-06-05 as energy-inconsistent. K sits inside the
micro potential (`PotentialExchange.h:123-131`, `out.mu_lR += mu_aug`), so
`dPi/dphi_M|_{n_l} != 0` through K. The shipping configuration fails test (1) of
the very proof used to kill the nS route.
Quantitatively, at VII's final state: `mu_vdW ~ 4.6e-5 J/kg` vs `mu_aug ~ 16572
J/kg` -- the vdW base is EIGHT ORDERS below the augmentation; Pi *is* the K-term
(s_aug = 1.0000). Effective contact-area power `1 + s_aug*dlnK/dlnrho_d`:
2.83 (rho_d=1000), 3.37 (1300), 3.49 (1363, III final), 7.17 (1500), 7.58 (1600).
Energy-safe = 1. Rejected mode = 4. Shipping config = 3.4-7.6.
NOT re-derived as a closed loop integral; the claim is that it fails the same
exact differential test the rejection doc treats as decisive.

**Finding B (MEASURED).** phi_M = porosity - micro_porosity = 0.000e+00 exactly at
both models' final frames, and micro_water_content == micro_porosity == porosity
bit-for-bit. Therefore `1 - phi_total ≡ 1 - n_l` identically in the late regime,
so the claim that `CurrentPorositySplit` tracks a different quantity from bulk
density is FALSE in this operating regime -- that energy-safe mode does restore
K/nS consistency where the headline numbers are read. VII trajectory: phi_M 0.4244
(t=0) -> <=0.0966 (100d) -> 0.000000 (>=200d).

**Exposure (MEASURED).** III and VII spend >90% of the run inside the 900-1400
segment (crossing 1400 at ~day 17-18 from a 1600 start). Weight on the 900-knot at
final state: III 0.0731, VII 0.1936, IV up to 1.000 (flat clamp, part of its domain
sits below 900). Chain validated: predicted K ratio for III 1.1111 vs observed
headline move +9.15/+9.33/+10.55%. Measured elasticities: III dln(sigma)/dlnK ~
0.83-0.96; VII dln(e)/dlnK ~ 0.054-0.063.

**Also measured:** only TWO nS modes are implemented (`Reference`,
`CurrentPorositySplit`); the third is a specced-and-killed design, never code.
Enum at `PotentialExchangeParameters.h:249-253`, parse `OGS_FATAL`s otherwise.
The rejection doc's derivative literal is printed unsigned; sympy returns it
negative. Its locators have drifted (live: `PotentialExchange.h:232/248`).

**The question that collapses the option menu, Vinay's alone:** the 2026-08-17/18
ruling made live-K the physics truth for III/IV/VII. It now measurably fails the
same `dPi/dphi_M = 0` test that killed `CurrentTotalSolidFraction`, at a higher
effective power over part of the range. EITHER that test is decisive -- in which
case live-K must go (options 1/4 below) -- OR it is not, in which case the
rejection of `CurrentTotalSolidFraction` needs revisiting on identical grounds.
BOTH CANNOT STAND.

Option menu (all III/VII figures PREDICTED from two-point elasticities, none re-run):
| # | option | energy-safe | cures A | III | VII | cost |
| 0 | disclose only | no change | no | 2.014 | 1.1265 | prose |
| 1 | floor K arg at rho_SR*n_s (=1600) | yes | yes | ~4.47 MPa | e~1.194 | run+rebaseline |
| 2 | per-model K tables | no | no | few % | <1% | RE-FIT+run+rebaseline |
| 3 | current_porosity_split + re-fit | yes | no | re-fit | re-fit | RE-FIT+run+rebaseline |
| 4 | live-K off, frozen per-deck K | yes | yes | ~4.47 MPa | e~1.194 | flag+run+rebaseline |
| 5 | K(rho_d, n_s) explicit | NO (new channel) | no | ~+17% | ~+21% | new formulation |
Option 5 is a novel constitutive term with no §12.1 source and is listed only for
completeness. Option 3 has one deck of code mileage in the whole tree
(`beacon_1a01_dsm_micromacro_stressprobe.prj`).

**Sequencing recommendation (agent opinion, not a ruling):** settle the energy
question FIRST -- it is free, needs no runs, and determines whether the table is
consulted at more than one point at all. Then the n_s/K mechanism, then the parked
"literal Dixon values for dd1400/1600/1800" decision evaluated against whatever
table is live by then, then ONE consolidated re-run + rebaseline + snapshot rather
than three. Reverse order costs at least one extra full suite re-run and risks a
§6.7 cross-artifact window.

**Side finding (§6.7.1 candidate, low priority, not fixed):**
`ms33_modelVII_freeswelling.prj:7` states `K=103879 J/kg` with no supersession
marker while the live 1600-knot is 104689.9129; lines 95-96 and 205 do carry the
"SUPERSEDED 2026-08-17" annotation, so :7 is a stale summary line. III's header is
clean.

---

## 2026-09-01 — beacon_1c ctest vtkdiff failure (ogsds01 Jenkins run): diagnosed; sigma abs gate widened one order of magnitude (Vinay-approved); NOT expected to fully close

**Trigger:** Vinay shared a Jenkins ctest log (`ctest -D HOSTNAME:STRING=ogsds01
--parallel 20 --preset bgr`), 3/1055 tests failing. Only one is DSM:
`RichardsMechanics_beacon_1c_dsm_micromacro_reference-time-omp-vtkdiff` (#721).
The other two (`ThermoHydroMechanics/MultiMaterial/DP_MCC/TM/square_1e1_2_matIDs-omp`
#980, `ThermoHydroMechanics/Linear/Beam_sealed_bimaterial/square_1e2_function-omp`
#987) are hard nonlinear-solver divergences at the first time step, CONFIRMED
unrelated to DSM: neither PRJ has been touched by any DSM commit since the
`6af2391` bulk import; #980's MFront `ModCamClay_semiExpl_constE.mfront` was
touched by 2026-08-12 commit `b4ff1e3`, but only a comment string + two
unrelated `Mechanics/ModifiedCamClay` test tolerances -- no compiled-code
change (the `@Algorithm` line is byte-identical before/after).

**DSM finding.** The vtkdiff `sigma` comparison against the committed
`beacon_1c_reference_t_1000.000000.vtu` fails: abs max norm
`[9.55e-12, 4.91e-11, 9.55e-12, 7.97e-12]` (axisymmetric rr/zz/phiphi/rz), all
four exceeding the `1e-12` abs gate; rel max norm on the rz (shear) component
is `57.7` (5771%) against the `1e-10` rel gate -- that component's reference
value is near zero, so any absolute noise reads as an enormous percentage.
Every other compared field (displacement, pressure, saturation, porosity,
transport_porosity, micro_pressure, micro_saturation, swelling_stress) passes
at ~1e-14..1e-16. Siblings `beacon_1a01`/`beacon_1b` (identical Tests.cmake
gates) pass in the same run.

**Ruled out: source-code regression.** `beacon_1c_dsm_micromacro_smoke.prj`
runs with K=0 (`hamaker_constant=1e-30`, no `vdw_augmentation_prefactor`
table) -- confirmed by reading `PotentialExchange.h`'s
`if (potential_augmentation_prefactor > 0.0)` guard directly, not just citing
`ee5816fc`'s commit message. Every DSM commit since `ee5816fc` (2026-08-12,
when the current reference was baselined) either only touches the live-K /
`AugmentationPrefactorTable` code (dead for K=0: `1bb414a`, `1030f06`,
`7ec39ec`, `2c066ce`, `642a8f8`, `d2ce0b7`), or is `331ff56` (same day,
`use_analytic_micro_jacobian` false->true) -- and this AGENTS.md already
records the tip AFTER `331ff56` as re-verified "22/22 PASS (gating 6 +
beacon 9-family + double_porosity + MCC shear/biax)", i.e. beacon_1c passed
against this exact committed reference WITH the analytic Jacobian already
active, on the verification machine. No plausible code-path change since
that verified-passing state.

**Leading (not closed) hypothesis:** build/toolchain-specific floating-point
divergence on Jenkins host `ogsds01`, beyond the noise budget `81b344a`
measured when it set the `1e-12` gate (cross-compiler ogs09/gcc 13.3.0 vs
macOS/clang: 4.4e-16..7.8e-16 -- ~10^4-10^5x smaller than the 4.91e-11
observed here). `grep -rl ogsds01` across the repo: zero prior mentions --
this host has never been characterized against ogs09/clang. Cross-compiler
confirmation of even the MS33 III/VII baselines was already "open on Vinay's
desk" per an earlier entry above, so this gap was never fully closed even for
the family it was originally measured on.

New lead (adversarially checked, not yet confirmed as the actual cause):
`RichardsMechanicsFEM-impl.h:1102`'s `#pragma STDC FP_CONTRACT OFF` covers the
micro 2x2 solve and `computeSecondaryVariableConcrete` (sigma_eff
finalization, line 5938), but `NumLib::LocalLinearLeastSquaresExtrapolator`
-- the generic, process-agnostic integration-point-to-node extrapolation that
actually writes nodal `sigma` into the VTU -- carries no such guard
(`grep -rl FP_CONTRACT NumLib` empty). A real, previously-unconfirmed
unguarded-FMA gap in a code path shared by every process, not DSM-exclusive.

**Action taken (Vinay's explicit approval — "you may change the tolerance by
one order of magnitude"):** `ProcessLib/RichardsMechanics/Tests.cmake`,
beacon_1c's `sigma` DIFF_DATA row, abs gate `1e-12 -> 1e-11` (rel gate `1e-10`
left unchanged, Vinay's choice). GUARDRAIL EXEMPTION CLAUDE.md sec 1.2
documented inline at the edit site. Scope limited to beacon_1c only --
beacon_1a01/beacon_1b pass under the unchanged gates and were not touched.

**Flagged before landing, per sec 5 (predicted, not verified):** this change
is arithmetically insufficient to make test #721 pass on the ogsds01 run's
own numbers -- the zz-component's measured 4.91e-11 deviation still exceeds
the new 1e-11 gate by ~5x, and the rz rel-gate blowup (near-zero reference)
is untouched by an abs-only change. No local build/ctest available in this
session to re-verify; the next actual CI run (or Vinay's own rebuild) is the
only way to confirm pass/fail state.

**Open, for Vinay:**
1. Whether to widen abs further -- measured noise suggests `~1e-10` would
   cover this run with ~2x margin (same margin-derivation style as
   `81b344a`), but that is two orders of magnitude, beyond what was approved
   here.
2. Whether to set sigma's rel gate to `0` (matching the pattern already used
   for every other near-zero-prone field in this same DIFF_DATA block) --
   no realistic rel-gate widening fixes a near-zero-reference blowup.
3. Whether to characterize the `ogsds01` toolchain against `ogs09`/macOS-clang
   (per `81b344a`'s own methodology) to close the still-open cross-compiler
   noise-floor question.
4. Whether to add an `FP_CONTRACT` guard to (or otherwise investigate)
   `NumLib::LocalLinearLeastSquaresExtrapolator` -- a generic, cross-process
   path, so any fix there is outside DSM's exclusive scope and needs broader
   review.

## 2026-09-01 — the EXACT one-Psi route runs clean; its recorded blocker is STALE; four K-side "fixes" refuted

### 43. MEASURED 2026-09-01 — film_energy_route=exact runs the whole suite with ZERO rejected steps

Vinay's derived one-Psi formulation (Psi_film = -(1-phi_M)*n_l*[I_vdw+I_aug+S],
I_T = int_0^{eps_v} Pi_T(w(e)) de, kinematic h-law w(e)=n_l(1+kappa*e);
`PotentialExchange.h:770-963`, `computeStrainedFilmEnergyPair`) was switched ON via PRJ
flags only -- NO code change, NO rebuild. Both elements are children of
`<potential_exchange>`: `<film_strain_coupling>kinematic</...>` +
`<film_energy_route>exact</...>`. `exact` REQUIRES `kinematic` (MEASURED,
`ExactFilmEnergyPair.cpp:93-104`); default is Operational.

**The 2026-06-25 blocker ("~87% assembled-tangent mismatch under confinement -> III/IV
stagnate", PREDICTED NOT VERIFIED) does NOT describe this binary.** All three decks reached
t_end exactly, ZERO rejected steps, no dt collapse, no cap hits:

| model | confinement | steps exact | vs operational | max Newton |
|---|---|---|---|---|
| III | CONFINED | 999 | 888 (+12.5%) | 12 |
| IV | CONFINED | 590 | 21500 (-97.3%, 36x FEWER) | 11 |
| VII | unconfined | 1171 | 920 (+27.3%) | 13 |

The predicted pattern INVERTED -- the confined pair was cheap, the unconfined one costliest.
The ONLY convergence failure in the campaign was the OPERATIONAL route (III hit max_iter=60,
diverged at step 5, t=5.666 s).

Results (MEASURED; extractors validated against committed references first):

| model | operational | exact Psi | team band | movement |
|---|---|---|---|---|
| III T/C/B MPa | 2.0138792096/2.0065415161/1.9936228708 | 3.1423878924/3.1184999783/3.0513069540 | [2.5720,7.4660] | BELOW -> INSIDE |
| IV T/C/B MPa | 4.393548876/4.361817351/0.892441777 | 5.018500722/4.973680939/1.254148507 | Top [1.5053,3.1659] | above -> further above |
| VII e_top | 1.1264739847 | 1.2240454450 | [1.0857,1.3210] mean 1.2085 | INSIDE fringe -> INSIDE, +1.29% of team mean |

Change concentrated in `swelling_stress` (~1.2 MPa, 12-13% rel, both confined decks);
hydraulics barely move (III pressure 0.97%, permeability bit-identical).

Flag-live proof (OGS does NOT log the parsed route, and `ConfigTree.cpp:455-527` only WARNS on
an unread tag -- a misplaced element is a near-silent no-op): (a) exact+off fatals at
`CreateRichardsMechanicsProcess.cpp:736`; (b) bogus route value fatals in `parseFilmEnergyRoute`;
(c) matched 1-day probe 131 vs 146 steps with swelling_stress differing 1.846e3 Pa; (d) zero
`warning:` lines in the IV log.

Also MEASURED from source: the H2 concern from `ULTRACODE_REVIEW_2026-06-14.md:33` (exact route
sourcing a parse-time scalar K, ignoring the live table) IS FIXED in HEAD -- both halves source
`effectiveAugmentationPrefactor(params, phi)`, mu-half `:735`, eigenstress `:2328`. exact x live-K
is legitimate.

### 44. REFUTED 2026-09-01 — four K-side "fixes"; do not repeat them

| variant | K argument | outcome |
|---|---|---|
| live-K (ships) | rho_SR*(1-phi_total) | fails the dPi/dphi_M=0 criterion |
| frozen-K | rho_d0 scalar | III +298%, VII +33% -> OUT of band |
| micro-argument | rho_SR*(1-n_l) | III +1153%, VII +118%; clips at the TOP knot (2780) for ~10 d, i.e. LEAVES Dixon's calibrated domain |
| conjugate term | unchanged | VOID -- the lane posited Psi = phi_m*Pi and differentiated THAT |

**The conjugate-term lane's error is the one to guard against.** Psi_film is the strain-path
INTEGRAL of Pi under a geometry law, NOT phi_m*Pi. Composing the energy with Pi instead of with
the geometry h is exactly the error Vinay had already identified and corrected (Pi = -dE/dh is
the gradient). Its +47% VII result and III/IV blow-ups are artifacts of a mis-posited energy and
are NOT evidence about the physics.

Root cause of all three overshoots: all three shipping decks are `LinearElasticIsotropic`,
E=52e6 Pa (uncited working value), so there is NO mechanical softening -- the decreasing live-K
was the only softening mechanism in the model. Corroborates the standing, never-re-run Task-13
finding that the residual is "missing dry-side plasticity, not K".

The energy criterion `dPi/dphi_M = 0` (used to reject CurrentTotalSolidFraction 2026-06-05) is
an ARTIFACT of the operational route POSITING sigma_sw rather than deriving it. Under the exact
route both halves come from one Psi by differentiation, integrability holds by construction, and
K may legitimately carry strain -- the criterion does not arise.

### 45. OPEN(Vinay) 2026-09-01

- Adoption of the exact route as production: Vinay's ruling, NOT taken here.
- RE-FIT question: every K knot was calibrated under the OPERATIONAL route, so the Dixon anchor
  was fitted against a different stress path than production would then use.
- Model I control under the exact-route campaign: staged, md5-verified byte-identical, but
  NEVER EXECUTED -- negative-control gap.
- The ~87% assembled-tangent mismatch (`ULTRACODE_REVIEW_2026-06-14.md:44`): unverified in
  EITHER direction; not re-checked against the current `:5478` tangent branch.
- IV's 36x step reduction: UNEXPLAINED.
- IV worsens against the teams under exact while III and VII improve -- not a uniform offset.

## 2026-09-02 — the exact route drops TWO strain-independent terms at assembly; one is a defect, one is intended-but-mis-justified — OPEN(Vinay)

### 46. MEASURED 2026-09-02 — operational and exact disagree at eps_v IDENTICALLY ZERO; mechanism fully decomposed

Three independent lanes (source trace, numerical bisection, unit-test reconciliation) converged on
the same two surviving terms. Model I cells are fully Dirichlet-constrained: |u|_max = 0.000000e+00,
eps_v = 0 for all time, so every strain-path term in the exact energy is identically zero. The
routes still differ.

    Ps^op - Ps^ex = SUM_k ( n_S^k - n_s_ref )( n_l^k Pi^k - n_l^{k-1} Pi^{k-1} )   <- W, REV weight
                  - SUM_k n_S^k * b * p_conf^k * ( n_l^k - n_l^{k-1} )             <- D, p_conf drain

**W — DEFECT, undocumented.** The exact branch calls `computeStrainedFilmEnergyPair` with
`nS <- active_nS`, but every document describing that slot says the eigenstress weight is
`1 - phi_M`: `PotentialExchange.h:776,:798`; the operational branch's OWN comment at
`RichardsMechanicsFEM-impl.h:2599` ("n_S passed in is the REV macro-solid fraction (1 - phi_M)");
`PI_OF_NL_EV_IMPLEMENTATION.md:140`. `computeActiveMicroSolidVolumeFraction` returns `n_s_ref`
(Reference) or `1-n_l` (CurrentPorositySplit) -- NEVER `1-phi_M`. `n_S` appears ZERO times in the
exact branch (`:2362-2431`), grep-verified by three lanes.
Provenance: the 2026-05-26 comment at `:444-454` records the function USED to return `1-phi_M` and
was deliberately changed. The active_nS denominator fix (CLAUDE.md §2 incident) was correct FOR THE
vdW LAW and was silently inherited by the eigenstress prefactor.
ROOT CAUSE IS SIGNATURE-LEVEL: the helper fuses two physically distinct quantities into one `nS`
argument -- the film-geometry fraction inside the vdW law (`h = n_l/(nS*rho_SR*Sa)`, `mu ~ nS^3`,
`xi ~ 1/nS`, `:82-148`) AND the eigenstress REV prefactor (`:936`). The operational branch keeps
them as two separate numbers. One argument cannot carry both. Measured weight mismatch at t_end:
2.84x (dd900), 1.99 / 1.74 / 1.54x (dd1400/1600/1800).

**D — INTENDED and Vinay-authorized, but its stated justification is FALSE at eps_v = 0.**
`RichardsMechanicsFEM-impl.h:2355-2363` carries a "PHYSICS TRADEOFF (predicted, §5)" comment
claiming "on the drained line the two agree". At eps_v = 0 the drained line is p_conf = -K_d*0 = 0
while the ACTUAL p_conf = Ps (0.35 ... 40.6 MPa) -- the excursion is 100% of p_conf. The comment
even names the failing case ("fully confined, eps_v~0 with p_conf growing"); that warning was never
propagated to the `PotentialExchange.h:770-799` header, to `PI_OF_NL_EV_IMPLEMENTATION.md:174`, or
to any test. `computeStrainedFilmEnergyPair` has no p_conf parameter at all.

**The near-agreement is an ARITHMETIC ACCIDENT.** W and D each carry 8-20% of Ps and very nearly
cancel; the observed shifts are the residue. Per-cell decomposition sums exactly (reconstruction
error -0.097% ... -0.0001%):

| cell | gap Ps_ex-Ps_op [MPa] | from W | from D | Ps_op [MPa] |
|---|---|---|---|---|
| dd900 | +0.325330 (+92.94%) | +0.103 | +0.222 | 0.35005220099112877 |
| dd1400 | +0.196350 (+4.007%) | -1.239 | +1.433 | 4.899970339977161 |
| dd1600 | +0.067325 (+0.476%) | -3.246 | +3.311 | 14.159948034188325 |
| dd1800 | -0.327293 (-0.806%) | -7.950 | +7.622 | 40.600313177350884 |

Absolute gaps are all the same order (+0.325 / +0.196 / +0.067 / -0.327 MPa), monotone in density,
crossing zero between dd1600 and dd1800. The sign flip at dd1800 is W overtaking D.
Because `active_nS` is constant, the exact accumulator telescopes exactly:
`Ps_ex(end) = n_s_ref * [ n_l,N*Pi(n_l,N) - n_l,0*Pi(n_l,0) ]`, reproduced to <=5e-16 relative in all
four cells by two independent lanes. The exact-route Ps carries NO phi_M, NO p_conf, NO b, NO kappa:
swept over 7-14 Biot values it gives ONE distinct Ps per cell, bit-identical final VTU md5s. The
operational route does not equal its own instantaneous level (7.740 vs 4.900 at dd1400) -- it is a
genuine path integral.

**PREDICTED (arithmetic from the measured decomposition, NOT re-run): fixing W ALONE makes agreement
an order of magnitude WORSE** -- it removes the cancellation, leaving D alone at
+63 / +29 / +23 / +19 % of Ps, uniformly positive. Do NOT treat the weight as an isolated fix.

**macro_porosity_floor REFUTED as the discriminator** (two-way test, MEASURED, all [DIAG],
uncommitted): dd900 floor OFF -> gap GROWS +92.9 -> +143.5%; floor 0.16 -> +59.6%; adding floor 0.08
to inert cells moves the gap in BOTH directions (dd1600 +0.48 -> -1.25%, sign flip). The floor is a
modulator via `n_S = 1-phi_M`, not the mechanism; dd900's gap is strictly monotone in n_S/active_nS
(3.089/2.842/2.595 -> 143.5/92.9/59.6%), which positively confirms W.
CORRECTION to an earlier claim in this session: dd900 is NOT a 20x outlier -- its absolute gap
(+0.325 MPa) is the same order as dd1800's (-0.327 MPa). What distinguishes it is n_s_ref = 0.3237,
the smallest of the four, on a Ps ~100x smaller than dd1800's.

### 47. TEST-COVERAGE FINDING 2026-09-02 — neither existing test gates this

- `ZeroStrainReduction` CANNOT discriminate. Counterfactual constructed at six parameter sets (the
  test's own two plus all four production Model I sets): the assertion is an ALGEBRAIC IDENTITY of
  the helper, residual 1.7e-19 ... 1.0e-09 against tolerances 3.2e-6 ... 1.8e-4 Pa, margin 1e4-1e5x.
  It is true at every parameter set and carries no information about which state the FEM is in.
- `AssembledExactPairClosesOperationalSigmaDoesNot` does NOT document the discrepancy as intended.
  It runs a DRAINED-LINE path (p_conf = -K_d*eps_v, e0 = -0.03 -> e1 = 0) and its `sigma_operational`
  is defined (test:351-360) with the p_conf drain DELIBERATELY REMOVED -- "the held-fixed p_conf
  drain cancels on a closed loop". Model I sits at eps_v = 0 with p_conf at full swelling pressure,
  i.e. maximally off that path. The test's own justification for dropping p_conf is precisely the
  assumption Model I violates.

### 48. OPEN(Vinay) 2026-09-02 — guardrail §9 fires (model-formulation decision)

- **Is the correct REV weight for the exact eigenstress `1-phi_M` (as all three documents state) or
  the frozen `n_s_ref` (as the code does)?** No lane proposed or applied a fix. Separating W from D
  requires a two-parameter split in the helper signature -- a code change nobody made.
- **How much of the favourable III/IV/VII exact-route result (III BELOW -> INSIDE band, VII to +1.29%
  of the five-team mean) is attributable to these two dropped terms rather than to the formulation?**
  That question was not asked before the campaign and should gate adoption.
- The `:2355-2363` tradeoff comment's warning was never propagated to the helper header, the design
  doc, or any test -- a documentation-propagation gap independent of the physics ruling.

## 2026-09-02 — THREE-ROUND ADVERSARIAL REVIEW of the exact one-Psi route: helper sound, surroundings not; SEVEN open items for Vinay

### 49. Review structure and why three rounds mattered

R1 find (4 lanes, documentation-first) -> R2 refute (adversarial, prompted to KILL) -> R3a audit the
kills + R3b fresh re-attack (blind to R2) -> adjudicate. 11 agents, 14 new simulations in R3b.
**Two R2 kills were UNSOUND and R3a restored them** -- the specific failure mode the third round
exists to catch:
- F5 (mu-level claim): R2 refuted a DERIVATIVE at eps_v != 0 in place of a LEVEL claim at eps_v = 0,
  and never opened the documentation on the question it had reframed. Restored as C2, blocker.
- F2 (severity demotion): R2's reclassification was sound, its demotion rested on two false claims.

### 50. THE HELPER IS SOUND; EVERY DEFECT IS IN WHAT SURROUNDS IT

R1 re-derived `computeStrainedFilmEnergyPair` term by term against `PI_OF_NL_EV_IMPLEMENTATION.md`
§2.1: I_vdw, I_aug, mu_mech, all three derivative blocks, sigma_sw_m and the kappa->0 reduction are
ALL CORRECT. The Maxwell cross identity holds to rounding BY CONSTRUCTION (`PotentialExchange.h:947-949`
and `:908-913` are one expression rearranged), MEASURED to pass at a 1e-9 gate; R2's GP drivers
reproduced the closed forms to machine precision. **Vinay's formulation is not in question.**

### 51. CONFIRMED FINDINGS (survived all three rounds) -- blockers first

**C5 -- THE CAMPAIGN'S CONTRAST CHANGED TWO SETTINGS AT ONCE. MEASURED, and it qualifies every
headline this session reported.** The shipping decks sit at `film_strain_coupling = Off` (parse
default, `PotentialExchangeParameters.h:508`); each EXACT deck differs by TWO lines, adding
`kinematic` AND `exact`. The middle rung was never run until R3b ran it:
  III Top: Off **2.0138792096** (888 steps, 1 rej) -> kinematic-only **1.5062108715** (952, 0)
           -> exact **3.1423878924** (999, 0). Band [2.5720, 7.4660].
  VII e_mean: **1.1259357844** -> **1.0953577889** -> **1.2229747777**.
  IV kinematic-only: **ERROR TERMINATION** at t = 2172321 s (12.6% of t_end), dt 22994 -> 0.1 s,
  25 rejected. So IV's 36.4x step reduction is NOT attributable to the kinematic coupling.
The kinematic rung moves III the OPPOSITE way (further BELOW band). "The exact route moves III into
the band" is a two-variable claim and must be reported as a three-rung contrast.

**C1 -- the R3 load-work term S uses TOTAL volumetric strain where the drained-line identity it
encodes requires the MECHANICAL one. Blocker, formulation fork.** `local_context.volumetric_strain =
trace(eps.eps)` (`:3531/:4060/:6054`) while `volumetric_mechanical_strain` is computed two lines away
(`:2807`) and unused. MEASURED at the VII exact end state: eps_v = 0.24651893, actual p_conf =
+1.348936e5 Pa, surrogate -K_d*eps_v = **-1.068249e7 Pa** (K_d = 4.33333e7) -- **80x and opposite
sign**; `-K_d*eps_v^mech` recovers p_conf to the spatial spread. Psi_film is ONE energy in ONE eps_v,
so the h-law (total) and S cannot use different strains without breaking the pair.

**C2 -- the exact mu half has NO load-transmission channel. Blocker.** `:697-699` removes the
`+b*p_conf/rho_lR` term under exact. MEASURED (3 FEM runs, dd1600): exact is **bit-identical to
`film_strain_coupling=off`** in n_l and micro_pressure (5.85563632e7 both) while operational adds
+12.93 kJ/kg, 24% of the mu level (micro_pressure 4.44272993e7). At the VII strain: exact
**-1408.2 J/kg** vs operational **+1827.2 J/kg** -- OPPOSITE SIGN, 3.2 kJ/kg apart. The boxed Step-19
gate (`b*K_d > 3*kappa*Pi`) does not describe the exact route: no constant -b*K_d, no 3*kappa*Pi.
R3a: the documented slopes +2Pi/rho and +3*kappa*Pi/rho are vdW-core-ONLY while production Pi is
~100% augmentation (Pi_vdw/Pi = 2.56e-10 at VII). Candidate mechanism for VII's over-swell.

**C3 -- xi0 = 1 is a sign-change locus in the exact mu half, and Model IV STRADDLES IT INSIDE ONE
MESH. Blocker.** Exact: `rho*dmu_mech/deps_v|_0 = 2*Pi_vdw - (1-xi0)*Pi_aug`, sign-indefinite.
Operational is expulsive everywhere in production. MEASURED xi0 at end states: dd900 1.8359 /
dd1400 0.9039 / dd1600 0.6763 / dd1800 0.4993 / VII 0.8767 / **IV clay 0.7917, IV pellet
1.2874-1.7803**. Model I Ps shift is monotone in xi0 (+92.97 / +4.01 / +0.47 / -0.81 %). The link to
IV's band position and its 36x step drop is PREDICTED, not verified.

**C4 -- Model I under exact breaks the Dixon anchor at the 900 knot by 1.93x without re-fit.**
dd900 0.35 target -> op 0.3500522010 -> exact **0.6753821363** (+92.97%); dd1400 +4.01%; dd1600
+0.47%; dd1800 -0.81%. Entry 46's W+D decomposition confirmed to six decimals in all four cells.

**C6 -- under exact the assembled tangent is a live WRONG-MODEL, WRONG-WEIGHT contribution.**
`:5418-5422` seeds the swelling tangent from the M2 explicit-K channel gated only on
`live_dry_density && dK_dphi_sw != 0 && is_pfmb` -- **no route test** -- and differentiates the
OPERATIONAL telescoped form at w_eff with the actual p_conf, weighted by `n_S_sw = 1-phi_M_sw` = 1.0
against the residual's active_nS = 0.5755. All three EXACT decks set live-K + PorosityFromMassBalance
and rho_d = 1250.6 is table-interior, so dK/dphi != 0. `pair.dsigma_sw_dnl` is computed and consumed
NOWHERE in ProcessLib. This closes entry 45's "87% mismatch" item as MEASURED-in-source (numerical
impact still INFERRED).

**C7-C11 (should-fix / latent):** S renormalises the drained bulk modulus by `1 - nS*n_l*b` = **31.7%**
(55% under the doc-specified kappa), stated in NO document, and Q2 was resolved only by inference
(`PI_OF_NL_EV:674-676` reads "implement that now" as a Q2 answer). Three doc sentences assert the FEM
eigenstress site is UNCHANGED under exact (`PotentialExchange.h:796-798`, `PI_OF_NL_EV:177-179`,
`:675`) -- true when written 2026-06-11 (`70015f8e2b`), falsified by H1 2026-06-14 (`1144cef56c`).
Test coverage does not reach production state: `FrozenHLimitMatchesShippedPartner` runs kappa = 1e-9
(series branch) while production uses kappa ~ 0.58; T-3's strain grid is 24.7x below the measured VII
strain; T-8 `ExpulsionProbeDrainedRamp` is GTEST_SKIP. Macro-floor cutoff g is applied to mu and never
to the eigenstress in EITHER route (latent: all three campaign decks set floor 0.0; the two
unregistered `formB_piexact` decks set 0.08). `include_S` is hard-coded true at `:725/:2411/:2421`,
not PRJ-selectable as §2.2 implies. `CurrentPorositySplit` + exact breaks the one-Psi property outright
and is unguarded.

### 52. OPEN(Vinay) -- SEVEN questions the review could not close from the documentation

- **Q-A. Which volumetric strain enters S -- total or mechanical?** The docs give the drained-line form
  and never state which strain realises it; the only warning (`:2355-2363`) names the CONFINED case,
  the opposite extreme from where it bites. Stake: 80x, opposite sign, at VII.
- **Q-B. kappa under `aggregate`: `active_nS` (code) or `1-phi_M` (FOUR doc sites + the user-facing
  OGS_FATAL)?** `STRAINED_FILM:227-230` records D1 as "RESOLVED (Vinay 2026-06-09): kappa = 1-phi_M";
  `maxwell_from_psi.tex:763`+`:59`; `PI_OF_NL_EV:140`; `CreateRichardsMechanicsProcess.cpp:110-112`.
  `1-phi_M == 1.0` STRUCTURALLY in all three EXACT decks (they omit TransportPorosityFromMassBalance),
  so the ratio is a constant **1.7375**. MEASURED stake: VII e 1.2229747777 (aggregate) vs
  **1.1924490611** (unity) = -2.50%, which **flips VII from above to below the five-team mean**;
  III -10.1%, i.e. **28% of the entire route movement**. Not exact-specific -- it is the shipped
  kinematic coupling. Note `maxwell_from_psi.tex` Step 19 itself says which kappa survives is
  "an experiment, not a theorem".
- **Q-C. Q2: keep S or drop it?** Never explicitly ruled; `PI_OF_NL_EV:601-608` tags it open
  `[formulation]`. Stake: the undocumented 31.7% modulus renormalisation, against an uncited E=52 MPa
  and the standing Task-13 "missing dry-side plasticity" finding.
- **Q-D. Should the exact mu half carry a load-transmission channel at all?** `STRAINED_FILM` §2 [D]
  says "the reversal lives in the load-transmission term"; `:697-699` removes it; no document
  reconciles them.
- **Q-E. Re-fit the K table under exact, or not?** Now quantified: dd900 +92.97% (factor 1.93) vs the
  Dixon median, monotone in xi0, whose crossover sits BETWEEN the dd900 and dd1400 knots.
- **Q-F. `specific_surface = 523`: m2/kg or m2/g?** `PotentialExchange.h:223` annotates m2/kg and the
  code consumes it so; every MS33 provenance header and both dated reviews write **m2/g**. Under the
  physically-consistent m2/g pair (S_a x1000, lambda /1000, xi invariant) Pi_aug is bit-invariant while
  Pi_vdw/Pi goes from ~2e-10 to **0.19-0.20** and `mu_mech` **FLIPS SIGN** (-1427.16 -> +81.97 J/kg).
  Load-bearing, not cosmetic.
- **Q-G. Which VII quantity is the headline** -- e at the r=0 top node (1.2240454450), the top-row mean
  (1.2236798732), or the domain mean (1.2229747777)? The campaign quoted two of these; the published
  deviation is currently ambiguous at the 0.1% level.

### 53. Agent-closable with Vinay's go-ahead (each implies a source edit the review could not make)

A-1 quantify C6 by adding the route test at `:5241/:5249` in a scratch build and re-running.
A-2 run VII/III exact with `include_S=false` (C7's nominated refutation) -- needs a source edit,
    `include_S` is hard-coded.
A-3 test IV's dt-collapse mechanism by disabling the p_conf drain -- `film_pressure_coupling=false` is
    SILENTLY OVERRIDDEN to true at `CreateRichardsMechanicsProcess.cpp:673-680`, unreachable from a PRJ.
A-4 extend T-3 to the production branch/strain range, add a kappa~0.58 case, un-skip T-8 for the SIGN
    only (needs no magnitudes, so §3 does not block that half).
A-5 isolate C3 with a deck sweeping n_l across xi0 = 1 in the pellet material.
A-6 doc/deck hygiene: correct the four stale sentences (C8), annotate the two unregistered piexact
    decks (C10), guard CurrentPorositySplit + exact (C11).

## 2026-09-30 - IC-FIX V0: the consistent initial micro state enters the MS33 submission branch

### 54. What landed (branch gen5_icfix_conformant_2026-09-30 = gen5_conformant_2026-09-08 @ 486235c679 + 686fcd6ef8)

Implements Vinay's ruling of 2026-09-30 (~22:00 CEST, "elevate the consistent IC to submission"; K set: K_fix, "your choices 1: ok").
- Code: the cherry-pick of 686fcd6ef8 (setInitialConditionsConcrete: with initial_micro_water_content declared and no transport_porosity
  property, phi_M0 = (phi - n_l0)/(1 - n_l0), so the stored micro water is the declared value). Guarded; decks with a transport_porosity
  property, without initial_micro_water_content, or without potential exchange do not enter the new branch (MEASURED bit-identical in
  the MBP lane 2026-09-30, TEST_RESULTS.md md5 94d0255ebd87db305d4ad2b8dceacd75, section 3/4).
- Decks: only the K prefactor literals changed (table decks: the `<prefactors>` line; Model I: `potential_augmentation_prefactor`), to K_fix:
  900 (floor 0.0) 24241.938995, 900 (floor 0.08, Model IV's gen-4 point) 18469.1763367, 1400 46000.3269694, 1600 104698.192423,
  1800 265909.813902 J/kg (previous: 24241.892204, 18469.144929, 46000.0, 104689.9129, 265905.06). Source: K_fix.json (md5
  e6f29bd492a5011c8f8a6f7688d6980e; copy in ANCHORS_MS33_ModelI/icfix_2026-09-30/), the SHIPPED procedure re-run with the patched binary
  (bin/ogs md5 23640465bff1ef6ce5c83d6c5dcd48d3). PRJ headers and the section 12.2 blocks carry the new values; the old ones stay as labelled history.
  MEASURED by check_decks_vs_kfix.py: comment-stripped, every committed deck equals its decks_kfix copy (the decks of the fx runs) and differs
  from the gen-5 deck only in the K prefactor lines.
- References: Reference, III, IV and both VII replaced by the fx outputs (old ones in superseded_references_2026-09-30_icfix/ with README,
  md5 old/new, vtkdiff numbers); Model I dd1400/1600/1800 unchanged (K_fix outputs pass 11/11 against them).
- Caveat that travels (K_fix_attribution.md): at 1400 and 1600 about 85 % and 95 % of the K change closes the shipped fits' own misses (the shipped
  fit accepted rel_tol 0.02), not the IC fix; the patched and the shipped binary also differ in toolchain (TFEL, MGIS, SDK; README of the readiness suite).
  Ps = target after the refit is calibration, not validation.
- NOT verified: the new references were not run through ctest (mini has none; the fx binaries were built with OGS_BUILD_TESTING=OFF). The vtkdiff
  replay with the decks' own tolerances passes 11/11 for each new reference against its source output. ctest on the MBP is OPEN.
- NOT tested: the floor-0.0 K_fix 900 knot on the Model IV concentric mesh with the fixed binary (IV keeps the gen-4 point per Vinay 2026-09-13).
