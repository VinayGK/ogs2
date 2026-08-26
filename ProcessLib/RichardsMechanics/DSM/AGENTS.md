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
