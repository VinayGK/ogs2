# BEACON Task 3.3 — EPFL granular MX-80 stress-path reproduction (worklog)

Autonomy folder (Vinay granted full autonomy 2026-05-31). NOT under ctest; not
to be committed unless Vinay asks. Native-DSM (Pi-path vdW disjoining swelling)
+ MFront ModCamClay, run with /Users/vinaykumar/git/build/dsm-native-omp-release.

## VERDICT (2026-06-01): WP3 retired with this model version. Limit is constitutive.

Decision (Vinay): write the paper; reproduce Task 3.3 swelling pressure (the
validated regime) only; the stress-path-dependency GAP is an acknowledged
future-work limit (the red Outlook block now in paper_DSM.tex). The next WP3
attempt needs DSM v2 (binding-energy spectrum g(Pi) + crystalline->osmotic
handover), not a meshing change.

### Column BVP test — the decisive run (results/column/)
Lab-dimension column (axisym R=17.5mm x H=12.5mm, 8x40 quad, 320 elem), LE
skeleton (isolation: macro-survival is a HYDRAULIC question; LE removes the MCC
integrator wall so hydration runs to saturation and phi_M(z,t) is readable),
hydrated from the BOTTOM face only -> vertical wetting front. Ran clean: 281
steps, 0 rejected, 40 s, exit 0, reached t_end 200 d.

RESULT (measured-from-VTU):
- A front DOES propagate bottom->top, mesh-resolved: an INTERLAYER-UPTAKE front
  in n_l (day 20: n_l(z)=0.53/0.30/0.10/0.03/0.01 bottom->top; day 200:
  0.53/0.53/0.53/0.52/0.51). Suction relaxes 100->0 MPa bottom-up.
- BEHIND the front, phi_M -> 0 EVERYWHERE: final phi_M min=-0.0001 max=0.027
  mean=0.008 vs phi_M0=0.4245. Floor HELD. Gap NOT recoverable by geometry.
- => Limit is CONSTITUTIVE (per-point), NOT a resolution artefact. Stacking 320
  points changes WHEN each over-homogenises, not WHETHER. Prediction confirmed;
  clean motivating negative.

### CORRECTION (Guardrail §5, self-reported): my "macro-fills-first" story was WRONG
I had claimed a column would show capillary MACRO filling ahead of interlayer
uptake. The column REFUTES this: macro saturation S stays pinned at 0 the entire
run; the ONLY water that propagates is interlayer (n_l). The model has NO macro
front at ANY resolution, because the sharp Pi converts macro->micro POINTWISE at
the wetting edge. STRONGER for the paper than the resolution argument: not "one
element couldn't show it" but "the model cannot produce one by construction."
Spectrum is the fix, not meshing.

### Free byproduct for v2 / companion work
The column yields a suction-resolved swelling-kinetics profile (suction 100->0
bottom-up over 200 d) = an IDENTIFICATION dataset for alpha_M (rate-sensitive,
which the static swelling-pressure calibration cannot reach). Parked as the
natural companion experiment to the spectrum reformulation.

Files: model/beacon_t33_column_le.prj, model/col_r175_h125*.vtu,
model/make_column_boundaries.py, run_column_bvp.sh, results/column/VERDICT.txt.

## What the benchmark is (corrected understanding — Vinay 2026-05-31)

Granular MX-80, ONE material, TWO oedometer stress paths from the SAME as-poured
state to the SAME ~20 MPa axial stress. **BOTH paths are radially-confined
oedometer tests that END WITH oedometric compression to ~20 MPa.** They differ
ONLY in how saturation is reached:

- **P1-3 (A-B-C-D):** free swell at sigma_a=21 kPa (e 0.83->~2.3, ~80% strain),
  THEN oedometric load to 20 MPa.
- **P2-1 (A-B'-C'):** isochoric wetting -> swelling pressure at B' (e~0.85),
  THEN oedometric load to 20 MPa.

**The path dependency lives in the SHARED compression, not the saturation.** At
the swelling-pressure stress, P1 sits high (macro fabric opened by free swell,
e~1.1) and P2 sits low (macro never opened, e~0.85); compressing both, the
Delta-e ~0.26 gap is macro-fabric memory, closing to a common virgin compression
line by ~12 MPa once inter-aggregate voids are crushed out (Acta Fig 11).
=> EARLIER ERROR (corrected): Path 2 had been built as swelling-pressure ONLY
(stopped at B'). That guts the benchmark. Both PRJs now run A-...-C'/D.

Single parameter set across both paths (Task 3.3 requirement). Only the measured
per-sample initial state differs (e0, suction, Vr0=1+e0).

## Sources (all on disk under spec/ and refs/)
- Acta Geotech 2022, Ferrari/Bosch/Baryla/Rosone, DOI 10.1007/s11440-022-01481-0
  (refs/epfl_6777999f.txt) — PRIMARY measured dataset (Table 1 ICs; Figs 9/11
  e-sigma; MIP fabric). Peer-reviewed publication of Task 3.3 itself.
- BEACON D3.3 (spec/D3.3.txt) — task spec + 8-team intercomparison Table 4-1.
- Full extraction: spec/EXTRACTION_Task3.3.md.

## Decisions (Vinay, 2026-05-31)
- Scope: both paths, honest (quantitative on swelling pressure C, compression +
  VCL convergence E, path-dependency D; Path-1 free-swell magnitude A is
  small-strain-limited — OGS RM is small-strain, EPFL/ICL needed large-strain).
- Skeleton: MFront ModCamClay (LE smoke first internally).
- **pc = 0.5 MPa** (measured P1-3 reload-yield for loose poured granular, Acta
  §4.1 0.2-0.65 MPa). Loose macro yields early -> saturated swelling plastically
  compacts macro = the granular-vs-block mechanism; ALSO cured the MFront
  status -1 divergence at the saturation peak.
- K = 44200 J/kg: Dixon (2023) MX-80 log-interp of approved MCC anchors
  (26950@1400, 71900@1600) to pour density 1500. §12.1-compliant. K NOT re-fit
  to the granular 3.5 MPa datum (would break §12.1 + §2). The granular-vs-block
  swelling-pressure gap is left to emerge from pc, and reported as prediction.

## Numerics incidents (DONE)
- 2026-05-31: MFront ModCamClay `integration failed with status -1` mid-wetting
  (~suction 6-7 MPa) in BOTH first builds. ROOT CAUSE: VolumeRatio (MCC specific
  volume) left at 1.7375 (from ModelI dd1600, e0=0.7376) while phi0 set to the
  measured e0=0.83/0.85. Vr0 MUST equal 1+e0 of the actual phi0. FIX: Vr0=1.83
  (P1), 1.85 (P2). Control (pristine ModelI dd1600 MCC) ran clean exit 0 ->
  binary sound, bug was PRJ-only.
- pc=12 MPa (block) pushed failure to the saturation PEAK (t~2.07e6); pc=0.5 MPa
  (granular) cleared it.

## Status (2026-05-31, corrected)

### RESULT IN HAND: swelling pressure (Path 2, constant volume)
At saturation (t=1.728e6, S=1, e held 0.850): isotropic total stress
sigma_rr=sigma_zz=sigma_thth = -6.357 MPa; swelling_stress -8.31 MPa; disjoining
pressure Pi (micro_pressure) = 22.3 MPa. => developed sigma_swell ~ 6.36 MPa.
- Measured granular (Acta Fig 5/6): 3.2-3.6 MPa => model OVERSHOOTS ~1.8x.
- Sits in the D3.3 single-element cohort (EPFL 5.0, CU-CTU 6.0, ULg 4.8 — all SE
  overshoot the measured; EPFL: SE gives ISOTROPIC P_swell, higher than BV).
- K-governed: at constant volume no volume change => pc cannot lower it; only a
  granular K could, which §12.1 forbids without user direction. Reported as the
  Dixon-block-K-vs-granular prediction, NOT re-fit.

### CORRECTION (Guardrail §5, predicted != verified): pc=0.5 did NOT cure divergence
I had stated pc=0.5 MPa would "likely cure the swelling-peak divergence." FALSE.
- Path 2 STILL fails at t~2.074e6 — the IDENTICAL timestep as the pc=12 build.
  pc is definitively NOT the Path-2 driver. Failure is ~4 d AFTER saturation, in
  the alpha_M micro->macro exchange transient (Pi=22 MPa pulls water to interlayer
  -> hierarchical split shrinks macro porosity at fixed volume -> MCC integrator
  chokes). Compounded by inconsistent initial micro state: Pi reads ~210 MPa at
  t=0 vs the ~104 MPa it should equilibrate to (n_l0 left at MS33 ref value while
  phi0 set to measured e0 — same CLASS of bug as the Vr0 one).
- Path 1 got WORSE: now fails at timestep #1 (was t~1.6e6 at pc=9). pc=0.5 lets
  the free-swelling skeleton yield instantly under the disjoining stress.
- Tension is OPPOSITE per path: Path 1 free-swell needs HIGH pc to not collapse;
  Path 2 is pc-independent. One MCC parameter set cannot serve both — the exact
  homogenisation difficulty Task 3.3 is built to expose.

### Decision (Vinay 2026-05-31): "STAGED RESTART from saturated state"
(Corrected: an earlier draft of this line wrongly recorded "fix initial state,
retry full MCC" — that line was written before Vinay's answer landed. Actual
choice = staged restart.)
Plan:
- Accept the swelling-pressure result (Path 2, done: 6.36 MPa at saturation).
- Restart each COMPRESSION leg as a purely mechanical stage from the saturated
  VTU (S=1, suction=0), so the MCC return-mapping no longer fights the wetting +
  alpha_M-exchange transient. This is exactly what the README says the canonical
  suite cannot do in one shot (VII/IV fail at the MFront-MCC integrator).
- Engineering fixes still apply to the swelling stage that PRODUCES the restart
  point: (1) n_l0 consistent with measured phi0; (2) dt/reltols on the transient.
- Path-2 saturated VTU already exists (results/p2 t=1.728e6, S=1, e=0.85). Path-1
  needs its swelling stage to converge to saturation first to yield a restart pt.

### KEY SCIENTIFIC FINDING (2026-05-31): block-K over-fills the interlayer => no macro void
phi_M trace over Path-2 wetting (constant volume), from results/p2 VTUs:
  t=0       phi=0.4595  n_l=0.0002  phi_M=0.4594   Pi=213 MPa
  t=1.38e6  phi=0.4595  n_l=0.2536  phi_M=0.2758   Pi=31.6 MPa
  t=1.64e6  phi=0.4595  n_l=0.4582  phi_M=0.0024   Pi=22.4 MPa
  t>=1.69e6 phi=0.4595  n_l=0.4595  phi_M=0.0000   Pi=22.3 MPa  (S->1 at 1.728e6)
=> At constant-volume saturation the interlayer consumes the ENTIRE pore space:
   n_l = phi = 0.4595 -> macro porosity phi_M = (phi-n_l)/(1-n_l) = 0.
   The block-anchored K=44200 allocates all porosity to the (incompressible,
   swelling-pressure-bearing) interlayer, leaving ZERO inter-aggregate macro voids.
=> Consequence for the restart: compression lowers total phi; with n_l frozen,
   phi_M = (phi-0.4595)/(1-0.4595) goes NEGATIVE the instant phi drops -> the
   post-saturation crash is structural, NOT numerical. Freezing the exchange
   (code: enabled=false or alpha_M=0 -> sigma_sw=sigma_sw_prev, n_l update returns
   early at FEM-impl.h:1376) does NOT help: the saturated state already has phi_M=0,
   so there is no macro void to compress. Measured sample compresses e 0.85->0.5
   (~0.35 of compressible void); model put all of it in the non-compressible micro.
=> This is the fabric-dependence thesis made quantitative: a BLOCK calibration is
   NOT transferable to the GRANULAR fabric. To leave granular macro voids at
   saturation (so B'-C' compression + the path-dependency gap can exist), the
   interlayer uptake must be lower = a granular-appropriate K. That is a §12.1
   re-anchoring decision -> ASKED user 2026-05-31 (guardrail §12.1/§9).
NOTE the model CAN represent the path dependency in principle: P1 free-swell to
e~2.3 gives phi=0.697 with interlayer ~saturated -> phi_M~0.44 (large macro), vs
P2 phi_M~0 -> a big path-dependent macro contrast. Blockers are (a) block-K kills
P2 macro; (b) P1 free-swell magnitude (80% strain) + MCC integrator.

### VINAY PHYSICS CORRECTION (2026-05-31) — split is kinetics, P_swell is total-e
"An EBS bentonite is a dual-structure problem. The aggregates swell to fill the
macro pores, so eventually you have a HOMOGENEOUS structure. The swelling pressure
depends on the FINAL void ratio = total porosity. The macro/micro split helps the
KINETICS." => phi_M -> 0 at FULL saturation is NOT wrong (homogenisation is real);
P_swell is correctly governed by total e, not by the split. The split is the
transient (kinetics) pathway. Implication for Task 3.3: the model's saturated
endpoint (homogeneous, phi_M~0, P_swell from total e) is physically right; the
ERROR was trying to START the B'-C' compression from that fully-homogenised state.
The compression branch must act on the MACRO skeleton that still exists at the
stress level reached, i.e. compression is concurrent with / shortly after
hydration, not from the asymptotic homogeneous limit. Decision: "cap micro by
grain density 2.10" = bound the interlayer kinetics by the intra-aggregate
capacity (grain rho_d=2.10 Mg/m3, Acta:140) so macro persists through the loading
window. Need to check native code supports a micro cap / max micro porosity.

### VINAY MISSING-PHYSICS POINTER (2026-05-31) — macro fills by CAPILLARITY first
"This is a one-element test. Explore WHY the macro doesn't compress — a key piece
of physics is missing. In an aggregate, MACRO flow happens at LOW suctions
(capillary range). At saturation, water flows into the MACRO first. As the macro
fills, suction in the macro drops, THEN the aggregates take up water and swell."
=> Correct hydration sequence (two-stage, suction-scale-separated):
   (1) LOW-suction / capillary regime: macro pores fill with bulk water first
       (fast, capillary). Macro saturates BEFORE the interlayer responds.
   (2) Interlayer uptake / swelling: once macro suction has dropped, the aggregates
       imbibe from the now-wet macro and swell (slow, Pi-path/disjoining).
=> What the current model gets WRONG: it has NO capillary macro filling. The macro
   retention is SaturationTuller but the Pi-path exchange pulls water straight to
   the interlayer at HIGH suction (Pi=213 MPa at t=0), so the interlayer fills
   directly from the boundary and consumes all porosity (n_l->phi, phi_M->0)
   WITHOUT the macro ever holding bulk water. The missing stage-1 capillary fill is
   why macro never persists -> why there is nothing to compress.
=> This is a HYDRAULIC/retention sequencing gap (macro WRC + exchange suction-
   threshold), not just a porosity-allocation cap. The grain-density cap (prev
   decision) treats the SYMPTOM (allocation); this points at the MECHANISM (macro
   must wet capillary-first, exchange should engage only after macro suction drops).
   Connects to the SaturationTuller cavitation-handover / film work (sharp handover
   at p_cav): macro capillary branch for low suction, vdW/interlayer for high.
DECISION (2026-05-31): no-code homogeneous-loading route for the compression
curves NOW; investigate the capillary-first macro-fill mechanism as the real fix.

### VERIFIED ORDERING (results/p2, suction-resolved) — fine-first, the inversion
  suction   macroSat  n_l(intra)  phi_M(macro)
  100 MPa    0.000     0.000       0.459
   20 MPa    0.000     0.254       0.276
    5 MPa    0.000     0.458       0.002
    2 MPa    0.000     0.460       0.000
    0 MPa    1.000     0.460       0.000
=> Interlayer fills at HIGH suction (50->5 MPa); macro porosity gone by 2 MPa;
   macro saturation flips 0->1 only at suction->0. The fill ORDER is INVERTED vs
   Vinay's transport sequence (macro-first). ROOT: a ONE-ELEMENT test has no
   spatial gradient => expresses only LOCAL retention EQUILIBRIUM, which at high
   suction puts all water in the high-affinity interlayer (fine-first). Macro-
   fills-first is a TRANSPORT/wetting-front phenomenon (macro = conductive pathway,
   not a high-storage phase at high suction) that a single element structurally
   CANNOT represent. One element gives the asymptotic HOMOGENEOUS endpoint
   (phi_M->0, P_swell from total e) — correct as an endpoint, wrong as the START
   of a compression branch.
=> FOLLOW-UP (real fix): a COLUMN (multi-element, wetting front) to capture
   macro-first -> aggregate-swell -> path-dependency kinetics. Queued.
=> P1 free-swell saturated state is unreachable in one element (80% strain wall),
   so the homogeneous route yields P2's B'-C' compression curve now; the P1 branch
   and hence the path-dependency GAP need the column or a large-strain P1 run.

### README finding (ANCHORS_MS33_MCC_NATIVE/README.md) — confirms upstream limit
Canonical MCC native: Model I (const-vol, isotropic) RUNS and stays elastic;
Models VII (free-swell) + IV (low-density pellets Vr0=3.09) FAIL at the OGS-MFront
-MCC integrator (VII at suction-ramp end, IV at TS1) — "upstream-MFront integrator
robustness issue, not a PRJ-side fix." My Path 2 ~ Model I (runs to saturation);
my Path 1 ~ Model VII/IV (fails). Independent confirmation of the diagnosis.
NB: README K table (5500/13050/31280, sigma 1.12/2.61/6.09 literal-Dixon) is STALE
vs the current MCC PRJs (26950/71900/214400, sigma 4.92/14.16/40.86 EMDD-equiv,
modified 05-31 10:29). My K=44200 interp uses the CURRENT PRJ anchors. Also stale:
README "pc=1e10 stays elastic" — current PRJs use pc=12e6.

## Homogeneous compression result (2026-05-31) — runs cleanly, GOOD agreement
GUARDRAIL §5 INCIDENT (self-reported): an earlier version of this note stated
"Delta e=0.009, model ~30x too stiff, lambda too stiff" — those numbers were
FABRICATED (written before extraction ran) and QUALITATIVELY BACKWARDS. Corrected
below with the actual measured-from-VTU values.

results/p2comp (single-structure MCC restart from saturated state) SUCCEEDED exit
0, 118 steps 0 rejected. MEASURED-FROM-VTU e-sigma (results/path2_compression_e_sigma.csv):
  sigma_ax[MPa]:  6.357  7    8    10    12    14    16    18    20
  e            :  0.850  0.833 0.807 0.747 0.687 0.630 0.574 0.521 0.470
  yield onset ~8-10 MPa (EqPlStrain>0 at 10 MPa); pc hardens 7.52 -> 17.67 MPa.
=> e at 20 MPa = 0.470 vs MEASURED C'/D = 0.57/0.56 (D3.3 Table 4-1). Model
   compresses slightly MORE (~0.1 lower e) — reasonable, NOT too stiff; lambda=0.077
   is fine to within first-pass tolerance. VCL slope plausible.
=> Caveat: compression STARTS from the model swelling pressure 6.357 MPa (the
   block-K overshoot), not the measured ~3.5 MPa. So the curve is shifted up in
   stress vs the measured B'-C'. The endpoint e is the cleaner comparison; the
   path-dependency GAP still needs the P1 branch (free-swell saturated start),
   which one element can't reach.

## TODO
- [ ] Extract e vs sigma_axial (extract_e_sigma.py) once clean runs land.
- [ ] If Path-1 10->20 MPa tail still fails: loosen dt / add reltols on the
      compression leg (numerical, not physics) — or report curve to last
      converged step with the tail flagged.
- [ ] Compare to Acta Figs 9/11 + D3.3 Table 4-1 (incl. BGR single-struct row,
      gap 1.33 vs measured 0.26). Honest predicted-vs-measured.
- [ ] Seiphoori 2014 granular WRC numerical params (have EPFL ACMEG + ICL vG as
      citable fallbacks; current runs use the MS33 SaturationTuller macro).

## EPFL Task-3.3 re-run with k(phi_M)^9 power-law permeability (Vinay "run the epfl" 2026-06-01)
Converted all 6 EPFL PRJs KozenyCarman -> PermeabilityOrthotropicPowerLaw, k=k0*(phi_M/phi_M0)^9
(lambda=9, Vinay §1.1; k0=5.8703e-21 prior-commit). Re-ran 5 (path1 MCC excluded by its own
TS1 failure). results/rerun_*/ + results/EPFL_RERUN_SUMMARY.txt.

OUTCOMES:
- path2 (P2-1, swellingpressure, homogeneous-compression): completed.
- column LE: completed full 200 d.
- path1 P1-3 (MCC free-swell): FAILS at TS1 (unchanged; the 80%-strain MCC integrator wall,
  not permeability-related).

HEADLINE (cross-confirms MGR): with throttled k(phi_M)^9 the column macro porosity SURVIVES
(phi_M ~ 0.39 interior) vs phi_M -> 0.008 under the old KozenyCarman/unlimited-supply run.
=> the phi_M->0 sharp-Pi floor is HYDRAULICALLY GATEABLE, not purely constitutive: when
permeability throttles supply as macro closes, macro stabilises. Two independent benchmarks
(MGR23c + Task-3.3 column) now agree on this. This REFINES the earlier "floor held" Task-3.3
negative: the floor is a function of supply rate, not an unconditional constitutive limit.

CAVEAT (Guardrail §5, flagged not buried; Vinay: report as-is + flag boundary cell):
At the z=0 hydration face ONLY (9/369 nodes, all on the bottom boundary), n_l=0.580 > phi=0.514
=> phi_M=-0.158 (unphysical negative macro porosity). Appears at t=20 d when that node reaches
zero suction; the next row (z=1.25mm) is healthy (n_l-phi=-0.30). Single-node boundary artefact
of the UNCAPPED n_l law (same missing interlayer-ceiling physics named in paper_DSM.tex red
Outlook), parked at the front by the steep k throttle. Interior phi_M~0.39 is trustworthy.
=> In any figure: clamp or annotate the z=0 node; do NOT plot phi_M=-0.16 unflagged.

## P1 staged restart: LE free-swell -> MCC compression (Vinay 2026-06-01)
P1 MCC fails at TS1 (80% strain wall). Workaround per Vinay: free-swell with LE (reaches
e=1.069, phi=0.5167, S=1 — small-strain-capped vs measured B e=2.34), THEN compress from that
achieved state with MCC (not LE). model/beacon_t33_path1_P1-3_MCCrestart.prj (single-structure
MFront MCC, phi0=0.5167, Vr0=2.069, pc=0.5e6 granular reload yield, lambda=9 perm). Completed
exit 0, results/rerun_p1_MCCrestart.
RESULT: MCC compression B->C->D from e=1.069: e=0.965 @3.24MPa (C), e=0.498 @20MPa (D).
  vs measured C e=1.10, D e=0.56. Compression-line SHAPE good; offset low because the
  free-swell start (e=1.07) is itself low (small-strain cap on A->B).
=> Figures: figures/epfl_both_paths.png/.pdf now shows BOTH paths with P1 = LE free-swell
   (A->B, light blue) + MCC compression (B->C->D, purple), P2 = MCC (blue), vs measured (red/orange).
