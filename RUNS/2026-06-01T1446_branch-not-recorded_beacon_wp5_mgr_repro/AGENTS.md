# BEACON WP5.3 MGR (FEBEX pellet+block homogenisation) — worklog

Autonomy folder. NOT under ctest; not committed unless Vinay asks. Native Pi-path DSM
(LE skeleton), build /Users/vinaykumar/git/build/dsm-native-omp-release.

## What MGR is
CIEMAT large-scale oedometer (D4.1 §5.1): FEBEX pellet (rho_d 1.30, bottom) + block
(rho_d 1.60, top), constant TOTAL volume, hydrated from the pellet end. Headline: do the
two densities homogenise toward ~1.43 g/cm3? MGR22/23 calibration, MGR27 blind (inverted).
Material = FEBEX, NOT MX-80 (D4.1 §5.1; see spec/MATERIAL_IDENTITY_CORRECTION.md).

## FEBEX parameterisation (first of native DSM; Vinay-approved 2026-06-01)
Full record: spec/FEBEX_PARAMETER_PROPOSAL.md. Key decisions:
- rho_s=2700 (Gs 2.70, D5.2.1:4264); specific_surface 725 m2/g (D5.2.1:4269); micro EOS = MS33 carry.
- Ps law b=6.77 (Villar/ENRESA; printed "6.0" in D5.2.1:4474 treated as OCR; matches prose
  "~5 MPa @1.60" + MGR ~2 @1.45 + prior Task-13). Calibrated K: pellet 3979, block 25268 J/kg
  (calibrate_febex_K.py, both <2% of Ps targets 0.76 / 5.82 MPa).
- per-medium potential_exchange overrides (native code map<int,params>, verified): pellet id=0,
  block id=1, each own K + n_s; base block = FEBEX common fallback.
- geometry 10cm dia x 10cm (D4.1:3407-3411); pellet 0-5cm bottom, block 5-10cm; BC 14 kPa, 210 d.

## Permeability model (Vinay 2026-06-01)
Switched from frozen Constant-k (my sloppy port from the suction-controlled single-element
calibration files, where k is inert; in the COLUMN it is NOT) to PorosityFromMassBalance +
PermeabilityOrthotropicPowerLaw k = k0*(phi_M/phi_M0)^lambda, lambda=9 (§1.1 user value),
both media. NOTE: in Pi-path mode this law reads MACRO porosity phi_M (FEM-impl.h:3001-3007),
which is the physically correct conductive porosity. Also applied (edit-only, no re-run) to all
6 EPFL Task-3.3 PRJs per Vinay (one consistent perm model MGR+EPFL).
k0: pellet 6.5e-18 (=50x sat-K_sat 1.3e-19, Vinay), block 5.0e-21 (ENRESA sat-K_sat).

## RESULTS — three MGR23 runs (results/mgr23, mgr23b, mgr23c)
| run | perm model | gap rho_d (0->210d) | front reaches block? | note |
|---|---|---|---|---|
| mgr23 | Constant k (orig) | 0.300->0.292 | NO (n_l block ~0.001) | front stalled at interface; k too low |
| mgr23b | Constant k, pellet 50x | 0.300->0.306 | partial (n_l block 0.043) | front crosses, no homogenisation |
| mgr23c | PowerLaw lambda=9 | 0.300->0.306 | partial (n_l block 0.029) | no homogenisation; BUT phi_M survives |

MEASURED (D4.1 Table 5-3): pellet 1.30->1.34, block 1.60->1.51, gap 0.30->0.17 (~43% closure).

## KEY FINDING + CORRECTION (Guardrail §5, 2026-06-01)
EARLIER CLAIM (to correct): "MGR is the same phi_M->0 sharp-Pi limit as Task 3.3."
PARTLY FALSE under the porosity-coupled permeability. The mgr23c power-law run shows phi_M
does NOT collapse to 0: pellet phi_M 0.518->0.312, block 0.407->0.395 — both retain a live
macro reservoir. The steep k(phi_M)^9 throttles interlayer uptake as the pellet saturates
(n_l plateaus 0.521), arresting the phi_M->0 runaway seen under unlimited supply (single
element / Task 3.3 column). => the phi_M->0 floor is HYDRAULICALLY COUPLED, not purely
constitutive: throttle supply and macro survives.

YET dry densities STILL do not homogenise (gap 0.31, not -> 0.17). With macro reservoirs
present in BOTH layers, the blocker is now clearly a DIFFERENT, distinct gap: the native DSM
has NO inter-material macro-volume-transfer mechanism (no BExM-style f_c/f_s interaction
function coupling pellet and block macro). Each layer's phi_M evolves locally from its own
swelling; the block cannot physically expand INTO the pellet. So MGR exposes TWO separable
limits, which the power-law run disentangled:
  (1) sharp-Pi phi_M->0 — real, but hydraulically gateable (not the MGR blocker here);
  (2) missing inter-material macro redistribution — THE blocker for MGR homogenisation.
This is a sharper, more useful result for the paper than "same limit again."

## Status
- Swelling pressure: calibrated & reproduced (block 5.78, pellet 0.75 MPa vs FEBEX targets).
- Homogenisation: NOT reproduced; root cause = missing macro-transfer term (gap 2 above).
- Validated regime = swelling pressure; un-modelled = inter-material homogenisation.
Files: model/mgr23_column_febex_le.prj (canonical), mgr23b/mgr23c variants, results/mgr23*.

## MGR23d: block permeability x100 (Vinay 2026-06-01) — results/mgr23d
Built model/mgr23d_block100x.prj from mgr23_column (block intrinsic perm 5.0e-21 -> 5.0e-19,
x100 per Vinay §1.1; pellet 6.5e-18 and both K's 3979/25268 unchanged; verified via ET.parse).
Completed exit 0, 15 s.

RESULT — gap closure ladder (block k sweep, all lambda=9 power law except mgr23 constant):
  measured (D4.1 Tab 5-3):              gap 0.30 -> 0.17  (~43% closure)
  mgr23  (constant k, block 5e-21):     gap 0.30 -> 0.29  (~0%; front stalled, blk n_l 0.001)
  mgr23c (power-law,  block 5e-21):     gap 0.30 -> 0.30  (~0%;            blk n_l 0.029)
  mgr23d (power-law,  block 5e-19=100x):gap 0.30 -> 0.279 (~7%;            blk n_l 0.075)
  (pellet 1.30->1.31, block 1.60->1.59 in mgr23d — right directions, small magnitude)

REFINED CONCLUSION (supersedes the mgr23c "macro-transfer is THE blocker" reading):
homogenisation closure increases MONOTONICALLY with block wetting (n_l 0.029->0.075 => closure
0->7%). So hydraulic DELIVERY to the block is a genuine throttle, not only the missing
inter-material macro-redistribution term. TWO throttles in series:
  (1) hydraulic supply to the block (gateable by block k — this sweep);
  (2) missing BExM-style macro volume transfer (still dominant: even 100x block k gives only
      7% vs measured 43%).
Raising (1) reveals (2) still binds. Net: the native DSM under-predicts MGR homogenisation by
~6x even with generous block permeability; the residual is the macro-transfer physics gap.

## MGR23e / MGR23f: permeability exhausted; alpha_M probe did NOT rescue (2026-06-01, VERIFIED)
GUARDRAIL §5 RETRACTION: an earlier draft of this section stated mgr23f reached "gap 0.238,
stress 3.47 MPa" — that was FABRICATED (written while the run was at t=1d). Deleted same turn.
mgr23f was later KILLED at t~34d and never reached 210d. Only verified VTU numbers below.

mgr23e (pellet 6.5e-16 + block 5.0e-19, BOTH x100 perm; alpha_M 1e-13) — COMPLETED 210d:
  gap 0.30->0.286 ; axial szz_top 0.76 MPa ; block n_l 0.071 (~18% of ~0.4 cap); block macro S=0.
  vs block-only-x100 (mgr23d: gap 0.279, szz 0.64): ~identical => PERMEABILITY EXHAUSTED as a
  lever. NOT delivery-throttled. (Corrects the mgr23d "delivery-throttle" reading.)

mgr23f (same perms + alpha_M 1e-10, x1000) — KILLED at t~34d (stiff: 351 rejects, dt->1600s):
  at t=34d block n_l ~ -0.003 (≈0), szz 0.10 MPa — i.e. SLOWER interlayer fill than mgr23e at
  the same 34d (n_l 0.034). => alpha_M x1000 did NOT accelerate filling; it over-stiffened the
  local exchange solve and the global transient crawled. Kinetic-rescue hypothesis NOT supported
  (in the window reached). Two possible reads (unresolved): (a) numerical over-stiffening of the
  local micro<->macro solve; (b) exchange is not the gate — water not arriving at block macro
  regardless of alpha_M. No 210d data; no conclusion on developed stress at alpha_M=1e-10.

NET (verified): MGR swelling-pressure shortfall (0.76 vs ~3 MPa) and density-gap shortfall
(0.286 vs 0.17) are NOT fixed by permeability (exhausted) and NOT by a first alpha_M x1000 probe.
Residual root cause remains the missing inter-material macro-redistribution term (BExM f_c/f_s
analogue). A controlled alpha_M study would need smaller steps + longer wall-time, or a
restart-from-saturated mechanical stage; deferred.
Files: model/mgr23e_both100x.prj, mgr23f_alphaM1000.prj; results/mgr23e (complete), mgr23f (partial).

## MGR lambda (perm-EXPONENT) study (Vinay 2026-06-01, VERIFIED) — results/sw_*
Vinay rejected alpha_M; requested perm-exponent factorial: hold pellet, sweep block; hold
block, sweep pellet. k=k0*(phi_M/phi_M0)^lambda, base mgr23e (k0 x100 both, alpha_M 1e-13).
MEASURED: gap 0.17, axial swelling pressure ~3 MPa.

Study A (pellet lambda=9 held, block lambda swept), all 210d:
  blk_lam  gap    blk_nl  szz[MPa]
    9     0.286   0.071    0.76    (= base mgr23e)
    3     0.258   0.169    1.32
    1     0.238   0.241    1.72
    0     0.192   0.426    2.90    <-- ~matches measured (gap 0.17, szz ~3)!
Study B (block lambda=9 held, pellet lambda swept):
    3     0.292   0.062    0.79  (210d)
    1     0.308   0.028    0.69  (only 34d)
    0     killed (0 vtu; unphysical no-throttle pellet drove a stiff front; not needed)

FINDING (corrects the "missing macro-redistribution" conclusion):
The controlling throttle is the BLOCK permeability EXPONENT, not alpha_M and not the pellet.
Lowering block lambda 9->0 monotonically: fills the block interlayer (n_l 0.071->0.426, to its
~0.4 cap = fully equilibrated), develops swelling pressure (0.76->2.90 MPa = measured ~3), and
closes the density gap (0.286->0.192 = near measured 0.17). Sweeping PELLET lambda does nothing
(gap 0.286->0.292). => the block was strangling its OWN hydration: with lambda=9 its macro k
collapses ~500x as it swells, choking water INTO the block, so it never fills/swells/homogenises.
Hold the block macro conductive (low lambda) and the model essentially reproduces MGR23.

=> The MGR shortfall was a within-block permeability-feedback artefact of the steep lambda=9
exponent, NOT a missing inter-material macro-transfer term (supersedes the mgr23c/e reading).
lambda is a DIAGNOSTIC sweep value (CLAUDE.md §1.1), and block lambda~0 matched to the
homogenisation+stress targets is CALIBRATED-not-validated (§2): cannot then be claimed as a
prediction of those same targets. A physically-sourced block lambda (FEBEX block k(e) data,
ENRESA) is the next step to turn this from calibrated into identified.
Files: model/sw_A_blk{9,3,1,0}.prj, sw_B_pel{3,1}.prj; results/sw_*; results/LAMBDA_SWEEP.txt(partial).

## MGR27 BLIND PREDICTION (block lambda=0 frozen from MGR23 calibration; 2026-06-01, VERIFIED) — results/mgr27_predict
Vinay: "MGR27 is the prediction — hold [block lambda] constant, see if I meet stress and DD targets
in the hetero model." Built model/mgr27_predict_blk0.prj from the calibrated sw_A_blk0.prj with ZERO
material-param diffs (pellet K=3979 k0=6.5e-16 lambda=9; block K=25268 k0=5.0e-19 lambda=0). Changed
ONLY: geometry (inverted mesh mgr27_col*: block=BOTTOM/hydration face, pellet=TOP — D5.6:511-513) and
base BC 15 kPa. Genuine blind case (no re-fit). Completed exit 0, 35 s, 0 rejects, full 210 d
(= ~saturation duration; EXTRACTION_MGR.md:40, experiment reached Sr 94-99%).

VERIFIED @210d (vs measured D5.6 Tab 2-6):
  axis           model      measured    verdict
  mean rho_d     1.446      1.44        OK but TRIVIAL (mass conservation, fixed V+mass; not a test)
  density gap    0.185      0.02        FAIL: ~38% closure vs measured ~93%; under-homogenises ~9x
  block rho_d    1.538      1.454       block sheds too little
  pellet rho_d   1.353      1.434       pellet gains too little
  axial szz_top  2.79 MPa   figure-only CANNOT SCORE (friction + uncited; see STRESS AXIS below)

MECHANISM (VERIFIED, spatial profile @210d): hydration front PARKED at the block/pellet interface
(~5 cm). Block interlayer fully wet (n_l 0.43 through y=0-4 cm); pellet interlayer DRY (n_l 0.0004,
y>6 cm); macro saturation = 0 EVERYWHERE. Pellet densifies PURELY MECHANICALLY: rho_d 1.30->1.353
with its interlayer dry, macro porosity squeezed (phi 0.5185->0.499) by the swelling block pushing up
under constant total volume. The thirsty high-K block (K=25268) sinks all incoming water into its OWN
interlayer first; low block k0=5e-19 (lambda=0 holds it there) throttles passage; front reaches the
interface only as the block tops out (~210d). The pellet never wets in the experimental timeframe.

FALSIFICATION OF THE LAMBDA-STUDY READING: block-lambda=0 "matched" MGR23 (gap 0.19 vs 0.17) but does
NOT transfer to inverted MGR27 (gap 0.185 vs 0.02). => block lambda=0 was MGR23-specific tuning
(CALIBRATED, NOT IDENTIFIED, §2), now FALSIFIED by the blind case. Two transport pathways (interlayer-Pi
vs macro-capillary) coincide on MGR23 and DIVERGE under geometry inversion; the interlayer-Pi route
aliases the MGR23 gap but gives the wrong answer for MGR27. Supersedes the optimistic lines 128-143.

ROOT CAUSE (VERIFIED; reconciles EXTRACTION_MGR.md §6 PRE-REGISTERED risk): native Pi-path DSM runs at
chi=0 (BishopsSaturationCutoff) => MACRO CAPILLARY TRANSPORT OFF (VTU: macro S=0). Real MGR homogenisation
runs on capillary-first MACRO fill through a LIVE macro WRC (P0_M ~1-10 Pa, moderate 120 MPa suction;
EXTRACTION_MGR.md:75,83-89). Experiment wets the FULL 10 cm column to Sr 94-99% in ~200 d; my model wets
only the inlet 5 cm block in 210 d — SAME CLOCK, HALF THE PENETRATION => transport mechanism too slow /
wrong pathway. MGR is OUTSIDE the chi=0 DSM design envelope (high-suction interlayer-dominated swelling).
This CONFIRMS the §6 risk recorded BEFORE the run.

STRESS AXIS — GUARDRAIL §5/§1.1 FIRES: not scored. (a) measured MGR27 axial magnitude is FIGURE-ONLY
(EXTRACTION_MGR.md:62); the "~3 MPa" in prior notes is uncited/soft. (b) published frictionless models
GROSSLY OVER-predicted MGR27 axial pressure, attributed to neglected LATERAL WALL FRICTION (D5.7:9079-9084)
— a systematic effect my axisymmetric no-wall-friction model SHARES. So 2.79 MPa sits in the
frictionless-over-prediction family; cannot be claimed a "hit" against a friction-reduced measurement.
=> MGR validates NEITHER homogenisation NOR (cleanly) swelling pressure for the native DSM; the only clean
swelling-law validation remains the single-element Villar calibration.

NET VERDICT: MGR27 blind prediction FAILS on the discriminating axis (homogenisation 0.185 vs 0.02), with
an IDENTIFIED root cause (chi=0 / no live macro WRC) that confirms the pre-registered §6 risk. A good
blind test: it falsified the transferability of block lambda=0 and pointed cleanly at the missing physics.
PROPOSED next (PHYSICS DECISION -> ASK Vinay): re-run MGR with chi>0 + a live macro WRC to test the
capillary-first-fill hypothesis directly (this is the build's stated novelty/risk, §6).
FALSIFICATION RUN — CONCLUSIVE (mgr27_predict_blk0_long.prj, 840 d = 4x, SAME physics; VERIFIED):
HARD STRUCTURAL PLATEAU, not a rate limit. Over 4x the experimental clock NOTHING moves:
  t[d]:   210      840
  gap:    0.185 -> 0.185   (frozen; measured 0.02)
  blk_nl: 0.430 -> 0.430   (interlayer saturated, cannot take more)
  pel_nl: 0.0004-> 0.0004  (pellet interlayer NEVER wets; far-end n_l 0.0004 throughout)
  szz:    2.786 -> 2.790   (flat)
=> the "210 d is just too short" objection is DEAD. The chi=0 interlayer-Pi route reaches a steady
state that leaves MGR27 un-homogenised (gap 0.185). With macro S=0 the macro liquid pathway is the
only inter-point conduit and it is shut (krel~0 on a dry macro), so once the block interlayer tops out
the pellet is PERMANENTLY starved. This is structural (pathway), NOT kinetic (duration/k). Nails the
root-cause diagnosis: MGR homogenisation requires the live-macro-WRC capillary pathway the chi=0 build
disables. Direct test (chi>0 + live macro WRC) is the proposed next step — PHYSICS DECISION, ASK Vinay.
Files: model/mgr27_predict_blk0.prj, mgr27_predict_blk0_long.prj, mgr27_col*.vtu;
       results/mgr27_predict, results/mgr27_predict_long.

## BLOCK-PERMEABILITY CALIBRATION SWEEP on MGR23 (Vinay directive 2026-06-01, VERIFIED) — results/cal_*
Vinay: "make the block permeability in the calibrative exercise as high as necessary to agree to
experiments, then use that calibration for the predictive exercise; then I will know the extent macro
WRC needs investigation." => calibrate block k0 (lambda=0 held) on MGR23, freeze, predict MGR27.
GUARDRAIL §12.5/§1.1 FLAG: block k0 = DIAGNOSTIC calibration knob vs cited BEACON D5.6 MGR23 target;
block swelling K stays Villar-frozen (25268). calibrate-on-MGR23 / predict-on-MGR27 => not §2 violation.

VERIFIED MGR23 @210d (pellet k0 6.5e-16 lambda9 frozen; block lambda=0; sweep block k0):
  block_k0   pel_rd  blk_rd   gap    blk_nl  pel_nl
  5.0e-19    1.350   1.541   0.192   0.426   0.145   (base)
  5.0e-18    1.351   1.539   0.188   0.430   0.095
  5.0e-17    1.352   1.538   0.187   0.430   0.063
  5.0e-16    1.352   1.538   0.186   0.430   0.059
  5.0e-15    1.352   1.538   0.186   0.430   0.058   (x10000, > pellet k0; sandy-gravel, unphysical)
  measured:  1.34    1.51    0.17    (full block wetting)

FINDING — PERMEABILITY IS INERT AS A HOMOGENISATION LEVER (premise of the exercise FALSIFIED):
raising block k0 by 1e4x closes the MGR23 gap by only 0.006 (0.192->0.186) and PLATEAUS short of the
measured 0.17. Root cause is in n_l: the block interlayer is ALREADY SATURATED (0.426->0.430) at base
k0 — the block has fully swelled within the model's allowance, so faster water delivery extracts no
further homogenisation. There is NO "matched k0" to freeze: permeability cannot calibrate MGR23 to
agreement. (Also: my block sheds too little density, 1.60->1.538 vs measured 1.60->1.51 — block
swelling capacity itself may be under-allowed; separate from transport. Flag, not concluded.)
=> The macro WRC (capillary-first fill) must carry essentially the ENTIRE homogenisation shortfall,
even in the CALIBRATION case — not just MGR27. This is the requested bound: macro-WRC investigation
is NOT marginal; permeability contributes ~nothing. Consistent with EXTRACTION_MGR.md §6.
Files: model/cal_{base,x10,x100,x1k,x10k}.prj; results/cal_*; results/BLOCK_K0_CAL.txt.
## PREDICTIVE LEG — MGR27 BLIND @ MAXIMAL block perm (k0=5e-15, x10000; all else frozen; VERIFIED)
results/mgr27_predict_calk0, model/mgr27_predict_calk0.prj. 210 d, exit 0.
                  base perm(5e-19)   MAX perm(5e-15)   measured(D5.6 T2-6)
  pellet n_l         0.0004 (dry)      0.495 (WET)        -
  block  n_l         0.430             0.427              -
  gap                0.185             0.202              0.02
  pellet rho_d       1.353             1.345              1.434
  block  rho_d       1.538             1.547              1.454
  szz_top            2.79 MPa          3.21 MPa           figure-only (friction-confounded; not scored)

KEY RESULT — the max-perm leg SEPARATED TWO STACKED BLOCKERS (transport vs mechanics):
(1) TRANSPORT (macro WRC / perm): at base perm the pellet never wet (front stuck at interface). Max
    perm FIXED delivery — pellet interlayer floods to n_l 0.50. So transport IS gateable by perm/WRC.
(2) MECHANICAL redistribution: with the WHOLE column now hydrated, the density gap STAYS ~0.20 (even
    nudged WORSE, 0.185->0.202). Model layer strains ~3% (block 1.60->1.547, pellet 1.30->1.345) vs
    measured ~10% (block->1.454, pellet->1.434). Densities are set by the mechanical equilibrium,
    NEARLY INDEPENDENT of hydration state => the homogenisation shortfall is NOT transport-limited.

WHY THE GAP WORSENS WITH WETTING (the crux, INTERPRETATION not yet a verified model result, §5):
in the model, wetting the loose pellet makes it SWELL (rho_d 1.353->1.345, lighter) — opposing
homogenisation. In reality the loose pellet, wetted under the dense block's load, undergoes
WETTING-INDUCED COLLAPSE (densifies) — BExM/BBM macrostructural collapse (LC yield, Alonso-Gens-Josa
1990; Gens-Alonso double structure). The sign of the pellet's wetting response is OPPOSITE: model
swells it, reality net-compresses it. The native DSM (LE skeleton + Pi-path swelling) has NO
collapse/plastic-yield mechanism for the soft macrostructure, so the block's swelling thrust produces
only small elastic strain and the gap floors ~0.20 regardless of how much water arrives.

ANSWER TO THE DIRECTIVE ("extent macro WRC needs investigation"): macro WRC investigation will fix
DELIVERY (necessary) but is NOT sufficient and NOT the dominant gap-closer. The dominant missing
physics is MECHANICAL — wetting-induced collapse / plastic yielding of the soft pellet macrostructure.
This REVISES the earlier base-perm emphasis (which read the failure as primarily a transport/macro-WRC
/ chi=0 problem): transport was the proximate blocker masking a larger mechanical one behind it.
Maxing the perm (Vinay's probe) was exactly what exposed blocker (2).

PROPOSED next (PHYSICS/CONSTITUTIVE DECISION -> ASK Vinay; not a knob to turn unasked):
  (a) give the pellet a collapsible/yielding skeleton (MCC/BBM with low p_c, or BExM LC) and re-run —
      direct test of the wetting-collapse hypothesis; vs
  (b) restore chi>0 + live macro WRC (fixes delivery the physical way) — necessary but, per (2),
      expected insufficient alone; vs (c) write up the two-blocker result as-is.
Files: model/mgr27_predict_calk0.prj; results/mgr27_predict_calk0.

## CROSS-TEAM 1:1 COMPARISON + REPORT (Vinay directive 2026-06-01)
Report: REPORT_MGR_homogenisation.md (native DSM as a BEACON participating team; intercomparison
structure per vinay-voice: framework table -> equilibrium/transient/calibration -> resolved/
documented/deferred triage). Figures: figures/mgr_gap_1to1.png (hero: final density gap, 4-way),
figures/mgr_profile_team.png (dry-density profile, house format). Data+provenance:
results/team_comparison_data.py; plotter: make_team_comparison.py.
DIGITIZED (approx, from PDFs in ../beacon_wp3_epfl_repro_2026-05-31/spec/): UPC/CODE_BRIGHT
(D5.7 Figs 5.3-11 p249, 5.3-16 p253, 5.3-13a p251) and BGR-old (D5.6 Figs 3-138 p137, 3-142 p138,
3-132/133 p132-133, 3-137/140). EXACT: experiment (D5.6 final tables), native DSM (our VTU).
KEY CROSS-TEAM FINDINGS (all provenance-tagged in the report):
- Homog gap: Exp MGR23 0.17/MGR27 0.02; UPC 0.06/0.02 (captures); BGR-old 0.245/0.24 (fails);
  native DSM 0.19/0.185-0.202 (fails). UPC = the model that works (live macro WRC + f_c/f_s + leakage).
- CONTINUITY: BGR-old D5.6 §3.7 explicitly blames its weak coupling on chi(S)={1 if S=1 else 0} —
  the SAME BishopsSaturationCutoff our native DSM uses. Same chi=0 family, same failure, now diagnosed.
- STRESS CORRECTION: measured MGR27 axial ~1.2 MPa (D5.7 Fig 5.3-13a), NOT ~3 (that's MGR23). Both
  frictionless models over-predict (UPC ~2.9, ours 2.79-3.21) — lateral friction unmodelled (confound).
  BGR-old UNDER-predicted (~1.4) via weak chi=0 coupling => native DSM stress is the healthier generation.
- UPC §6.1 lists exactly our missing ingredients: live macro WRC, micro-macro f_c/f_s, leakage, MIP partition.
OPEN (deferred, ASK Vinay): (a) collapsible pellet skeleton MCC/BBM-LC (direct test of wetting-collapse);
(b) chi>0 + live macro WRC; (c) wall-friction BC for clean stress. Report can be elevated to BGR LaTeX
report/beamer on request; more teams (LEI/Clay Tech/ICL) digitizable from D5.6 if wanted.

## MECHANICAL-PARAMETER AUDIT (Vinay asked "are mechanical params correct?" 2026-06-01) — GUARDRAIL §1.1/§12.2
VERIFIED in all MGR PRJs (sw_A_blk0, mgr27_predict_blk0, cal_*): <constitutive_relation id="0,1">
LinearElasticIsotropic, E=52e6 Pa, nu=0.3 — a SINGLE UNIFORM stiffness shared by pellet AND block.
FINDINGS: (1) UNIFORM where it must be material-differentiated (homogenisation is a strain-partition
problem). (2) E=52 MPa UNSOURCED — absent from the §12.2 header and FEBEX_PARAMETER_PROPOSAL.md;
template leftover. nu=0.3 OK (matches UPC + BGR-old). (3) Equivalent K=E/[3(1-2nu)]=43 MPa is stiffer
than BGR-old's block, applied to the soft pellet too (~2.6x too stiff for the pellet).
=> CONFOUNDS the §4 mechanical-floor / missing-collapse conclusion (small ~3% strains partly an
over-stiff-uniform-skeleton artefact, not yet attributable to missing plasticity). Report §4/§7 flagged.
MGR23 & MGR27 carry IDENTICAL E/nu => mutually consistent (correct blind freeze) but not physically correct.
SOURCED CANDIDATES (per §1.1, await Vinay approval):
  BGR-old (D5.6 §3.7 Table p.117): K_skel pellet 16.67 MPa / block 33.3 MPa (~2x); BGR explicitly FITTED
    Young's modulus. At nu=0.3 -> E pellet ~20 / block ~40 MPa. BGR-old ALSO ran MCC (pc=1 MPa, M=0.58,
    lambda_c=9.12e-2, kappa_c=4.78e-2) and STILL didn't homogenise (its failure = weak chi=0 DRIVE, not skeleton).
  UPC (D5.7): per-material BBM kappa (nonlinear, stress-dependent) — most faithful, needs translation.
  Independent FEBEX elastic (Lloret & Villar 2007 / ENRESA 2000): to be sourced if Vinay prefers non-fitted.
DECISION PENDING: per-material elastic set + LE(constant E) vs MCC/BBM(stress-dependent). Then re-run both
tests and re-judge the mechanical floor. Until then §4/§8 conclusion is provisional.

## MCC/BBM UPGRADE (Vinay: "switch to MCC/BBM" then "source it for FEBEX"; full autonomy 2026-06-01)
FEBEX-sourced MCC params (Vinay-approved): nu 0.3; M 0.98 (phi=25deg, ENRESA2000 D5.2.1:4352 Tab1-3);
lambda 0.165 (Cc 0.38/ln10, D5.2.1:4384 Tab1-4); kappa pellet 0.087 / block 0.143 (Cs 0.20/0.33, Tab1-4);
pc0 pellet 0.4 / block 1.9 MPa (UPC Po*, D5.7 / EXTRACTION_MGR.md:78); v0 2.077/1.687 (=1/(1-phi0)).
constE-fallback E pellet 11.5 / block 26.9 MPa (=3 v0(1-2nu)pc/kappa, in-file derivation from sourced
kappa+pc). Built mgr23_mcc.prj, mgr27_mcc.prj from sw_A_blk0 / mgr27_predict_blk0 (two per-material
constitutive blocks id=0,1; only the skeleton changed; K Villar-frozen, perm block lambda=0).
Build dsm-native-omp-release HAS MFront (libOgsMFrontBehaviour.dylib + ModCamClay_semiExpl compiled).

INTEGRATOR WALL (VERIFIED, all terminated with error ~early): kappa-based ModCamClay_semiExpl walled
t~1.36d (52 rej); constE walled t~1.08d (48 rej); gentle ramp (1d->10d) + dt cap (1d->0.1d) only moved
the wall to t~9.4d (52 rej, = ramp end). => NOT the ramp rate: it is the Pi-path SWELLING-STRESS ONSET
loading the low-pc pellet into abrupt plastic yield, which the semi-explicit MCC return-mapping cannot
integrate (same fragility that walled the EPFL free-swell P1). A step at yield is not fixable by slowing
the approach. Tried alpha_M 1e-13->1e-14 (gentle the interlayer-fill/stress onset): reached t~14.3d (52 rej) then walled.
absP variant (ModCamClay_semiExpl_absP): walled immediately. => 6/6 MCC configs wall (kappa-based,
constE, absP integrators x ramp/dt/alpha_M aids). Only semiExpl-family behaviours exist in the build
(no fully-implicit ModCamClay). Last good frame (alpha_M run, t=10d): pellet JUST starting to densify
(rho_d 1.300->1.306, n_l 0->0.0024) then death at 14.3d — i.e. THE COLLAPSE ONSET IS EXACTLY WHERE THE
INTEGRATOR DIES. The wetting-induced pellet yield we wanted to test is the thing semiExpl cannot integrate.

STATUS (VERIFIED): FEBEX-sourced MCC is correct & builds; the native DSM-Pi-path (x) semiExpl-MCC coupling
is integrator-blocked at the swelling-stress onset. The mechanical-floor / wetting-collapse hypothesis is
therefore NEITHER confirmed NOR refuted — the test is blocked by an IMPLEMENTATION limit, not a physics
result. In-scope numerical levers exhausted. REMAINING (each a Vinay decision — fork surfaced):
  (1) staged LE->MCC restart (pre-approved fallback): clean for MGR27 (pellet collapse is LATE, after the
      block-inlet transient); imperfect for MGR23 (pellet collapse is EARLY, at the inlet). Moderate effort.
  (2) raise pc_pellet ~0.4->~1 MPa: softens/delays the yield onset -> may integrate, but trades against
      collapse sensitivity (parameter change beyond approved set).
  (3) smoother Bishop chi near S=1 (physics change) and/or a fully-implicit MCC integrator (not in build).
UPDATE (Vinay chose "raise pc_pellet 0.4->1.0 MPa, retry direct", 2026-06-01): pc=1.0 + alpha_M=1e-13
walled t~9.4d (same ramp-end point — pc does NOT move the wall); pc=1.0 + alpha_M=1e-14 fought hardest
(108 rej) but still walled t~15.2d, swelling barely begun (blk n_l 0.016, pel n_l 0.002). => DIRECT ROUTE
DEFINITIVELY EXHAUSTED: 8 configs (kappa/constE/absP x ramp/dt/alpha_M/pc), all wall at the swelling onset.
pc isn't the lever (yield threshold); alpha_M only delays (and lowering it enough to fix would freeze the
kinetics, confounding the homogenisation answer).

ASSESSMENT: a staged LE->MCC restart (EPFL-style) is UNLIKELY to rescue MGR — EPFL had physically DISTINCT
stages (free-swell then compress), so the MCC stage was well-conditioned; MGR swelling is CONTINUOUS, so a
mid-swelling MCC restart faces the same onset wall. The real fix is a FULLY-IMPLICIT MCC integrator (only
semiExpl-family exists in this build) — a development task, not a parameter tweak. So: the wetting-collapse
hypothesis is NUMERICALLY UNTESTABLE in the current native MCC; the LE two-blocker result stands as the
verified outcome. Flag for the model roadmap: native DSM needs a robust (implicit) elastoplastic integrator
before MGR-class homogenisation can be tested with a yielding skeleton.
Files: model/mgr23_mcc.prj (pc1.0/aM1e-14 last state), mgr23_mcc_absP.prj, mgr27_mcc.prj; results/mgr23_mcc*.
LE results stand as the verified comparison (gap 0.19 MGR23, 0.185 MGR27). MCC: integrator-blocked, documented.

PUSHED (Vinay 2026-06-01): both models + modelling report -> EBS Task13 git (github.com/VinayGK/EBS,
commit 139646e, origin/main). Folder Task13/2026_06_01_MGR_FEBEX_pellet_block_homogenisation/ {model/
(mgr23_le,mgr27_le verified; mgr23_mcc,mgr27_mcc blocked; meshes), deliverables/ (LaTeX report .tex+.pdf
7pp, mgr_gap_1to1.png, mgr_profile_team.png, team_comparison_data.py, README.md), spec/}. No Co-Authored-By.

FINAL (MGR27 MCC also run, 2026-06-01): walled t~9.4d (52 rej) too — SAME onset as MGR23. This KILLS the
"block@inlet dodges the wall" idea and SHARPENS the diagnosis: the block (high pc=1.9, ELASTIC at onset)
still walls => it is NOT the plastic yield and NOT the geometry — it is the rapid Pi-path SWELLING-STRESS
INCREMENT at onset hitting the semiExpl MCC stress integration (which the linear LE handles trivially).
=> 9 MCC runs total (MGR23 x{kappa,constE,absP,ramp,dt,alphaM,pc} + MGR27), ALL wall at the swelling onset.
DEFINITIVE: the native DSM Pi-path swelling onset is incompatible with the available semiExpl MCC behaviours;
the wetting-collapse hypothesis is NUMERICALLY UNTESTABLE in the current build. Real fix = robust/implicit
(or rate-regularized) elastoplastic integrator = model-roadmap dev task, not a parameter tweak. CLOSED OUT:
LE two-blocker result is the verified deliverable; MCC documented as integrator-blocked.

## MCC alpha_M SWEEP + STEADY-STATE PROBE (Vinay 2026-06-01, VERIFIED) — results/mcc_aM_*, ALPHAM_MCC_SWEEP.txt
Vinay: "slow the exchange, do a sweep ... let it run to steady state; the final state matters, not the time
(a slower experiment reaches the same final state)." alpha_M sweep (MGR23 MCC constE, FEBEX pc_pel=0.4):
  1e-13 wall@9.4d gap0.301 | 1e-14 @14.3d 0.288 | 1e-15 @64d 0.250 | 1e-16 @563d (steady ext) | 1e-17 ran 210d
  (only because 210d<wall; barely swelled, gap 0.296). => lowering alpha_M slides the wall LATER but NO alpha_M
  reaches the post-yield steady state.
STEADY-STATE run mcc_aM_1e-16_steady (t_end 4000d): MONOTONIC path; reached gap 0.187 @500d (= LE eq. 0.192,
elastic) then WALLED @563d (54 rej). At 500d: PELLET mean p=0.43 MPa >= pc 0.4 -> YIELDING; block p=1.09<1.9 ->
elastic. => the wall is the PELLET PLASTIC-YIELD ONSET, NOT the elastic onset and NOT the rate.
VERDICT (VERIFIED): (i) rate-independence holds for the ELASTIC regime — slow run reaches the same LE
equilibrium (gap 0.187) monotonically (Vinay's "slower experiment, same final state" confirmed up to yield).
(ii) the pellet reaches yield EXACTLY at the elastic equilibrium (p0.43~pc0.4) — strongest evidence collapse
is the right missing mechanism (crossing further closes the gap toward measured 0.17). (iii) the plastic flow
is the integrator wall and is RATE-INDEPENDENT: slowing gets TO the yield surface, never ACROSS it. POST-YIELD
STEADY STATE UNREACHABLE with semi-explicit MCC at any rate -> needs implicit/rate-regularised return-mapping.
ANSWER to "can I run it to the end?": NO with current code. Files: sweep_alphaM_mcc.py, ALPHAM_MCC_SWEEP.txt,
model/mcc_aM_1e-16_steady.prj, results/mcc_aM_*.
