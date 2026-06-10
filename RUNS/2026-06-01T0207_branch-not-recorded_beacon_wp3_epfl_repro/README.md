<!-- RUN BUNDLE (CLAUDE.md §13) -->
**Bundle:** 2026-06-01T0207_branch-not-recorded_beacon_wp3_epfl_repro
**Run finished:** 2026-06-01 02:07 (+0200) — last `OGS completed on 2026-06-01 02:07:28+0200` in the run logs
**Branch / commit:** branch not recorded. Commit: the run logs report version string `vdw-baseline-2026-05-08-74-g62bc63ba` (i.e. `62bc63ba`, 74 commits past tag `vdw-baseline-2026-05-08`); no branch name appears anywhere in the folder's own records.
**Binary:** `/Users/vinaykumar/git/build/dsm-native-omp-release` (per the folder's AGENTS.md)
**Full raw outputs (untracked, local):** `/Users/vinaykumar/git/ogs/beacon_wp3_epfl_repro_2026-05-31/`
**Report / beamer:** TODO (no formal report/beamer in the folder). Working analysis = `AGENTS.md` (copied here). The campaign feeds the red Outlook block of `paper_DSM.tex` (per AGENTS.md VERDICT) — pointer only, not in this bundle.

---

# BEACON Task 3.3 — EPFL granular MX-80 stress-path reproduction (2026-05-31 → 06-01)

Snapshot bundle consolidated 2026-06-10 from the campaign folder above. All
content below is distilled from the folder's own `AGENTS.md` and
`results/EPFL_RERUN_SUMMARY.txt`; no new numbers were computed.

## What was run, and why

Native Pi-path DSM (vdW disjoining swelling) + MFront ModCamClay on the BEACON
Task 3.3 benchmark: granular MX-80, one material, two oedometer stress paths
(P1-3 free-swell-then-load, P2-1 isochoric-wetting-then-load) to the same
~20 MPa axial stress. Single parameter set across paths; K = 44200 J/kg by
Dixon (2023) log-interpolation of approved MCC anchors to pour density 1500
(NOT re-fit to the granular datum — §12.1/§2). Also a lab-dimension column BVP
(axisym 8x40, hydrated bottom-only) as the decisive macro-front test, and a
final re-run of all PRJs with k(phi_M)^9 power-law permeability (lambda=9,
Vinay §1.1; k0 = 5.8703e-21 prior-commit).

## Headline results (quoted from AGENTS.md / EPFL_RERUN_SUMMARY.txt)

- **Swelling pressure (Path 2, constant volume):** developed sigma_swell
  ~6.36 MPa at saturation vs measured granular 3.2-3.6 MPa — model overshoots
  ~1.8x, sitting in the D3.3 single-element cohort (EPFL 5.0, CU-CTU 6.0,
  ULg 4.8). K-governed; reported as the Dixon-block-K-vs-granular prediction.
- **Key scientific finding:** block-anchored K over-fills the interlayer at
  constant-volume saturation (n_l -> phi = 0.4595, phi_M -> 0) — no macro void
  left to compress; a block calibration is NOT transferable to granular fabric.
- **Vinay physics correction:** the split is kinetics; P_swell is total-e.
  phi_M -> 0 at full saturation is the correct homogeneous endpoint; the error
  was starting the B'-C' compression from it.
- **Homogeneous compression (single-structure MCC restart from saturated
  state):** ran clean (118 steps, 0 rejected); e at 20 MPa = 0.470 vs measured
  C'/D = 0.57/0.56 — VCL slope plausible, lambda=0.077 fine to first-pass
  tolerance. Curve shifted up in stress (starts at the 6.357 MPa overshoot).
- **P1 staged restart (LE free-swell -> MCC compression):** e = 0.965 @3.24 MPa
  (C), 0.498 @20 MPa (D) vs measured 1.10 / 0.56; compression-line shape good,
  offset low because the LE free-swell start (e=1.069) is small-strain-capped
  vs measured B e=2.34.
- **k(phi_M)^9 re-run headline (cross-confirms MGR):** with throttled
  permeability the column macro porosity SURVIVES (phi_M ~0.39 interior) vs
  phi_M -> 0.008 under KozenyCarman — the phi_M->0 sharp-Pi floor is
  hydraulically gateable, not purely constitutive.
- **VERDICT (Vinay, 2026-06-01):** WP3 retired with this model version; limit
  is constitutive. Next attempt needs DSM v2 (binding-energy spectrum g(Pi) +
  crystalline->osmotic handover), not a meshing change.

## Caveats (as flagged in the folder's own records)

- path1 P1-3 full MCC fails at TS1 (80%-strain MFront-MCC integrator wall),
  unchanged by permeability; independently confirmed by ANCHORS_MS33_MCC_NATIVE.
- Boundary artefact in the k^9 column: z=0 hydration-face nodes show
  n_l > phi => phi_M = -0.158 (uncapped n_l law); interior trustworthy —
  clamp/annotate that node in any figure.
- Two Guardrail §5 self-reports are recorded in AGENTS.md (a fabricated
  earlier Delta-e claim, corrected from VTU; and the wrong "macro-fills-first"
  story, refuted by the column run).

## Bundle contents / omissions

Copied: AGENTS.md, run scripts, extract/plot scripts, `model/` (PRJs + meshes
as run), `figures/`, results CSVs + summary/verdict TXTs, run logs (logs >1 MB
tailed to `*.log.tail`, 200 lines), `final_state/` = largest-t VTU per model
variant per run dir (15 files).
Omitted (size policy, in the original folder only): full VTU time series,
`spec/` BEACON deliverable PDFs/txt (~105 MB; spec extraction .md kept),
`refs/` Acta Geotech 2022 PDF+txt (~1.6 MB), `data/` (empty).
