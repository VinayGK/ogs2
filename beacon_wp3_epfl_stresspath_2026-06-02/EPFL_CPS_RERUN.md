# EPFL Task-3.3 re-run under current_porosity_split (2026-06-05)

Re-ran the stress paths with the disjoining switched to the energy-consistent
micro-only mode `micro_solid_volume_fraction_mode = current_porosity_split`
(nS = 1 - n_l), and K re-derived from the recalibrated MX-80 anchors.

## Setup
- Binary: `build/pdisj_maxwell_revref_20260605/bin/ogs` (branch `dsm_cps_mx80`).
- Mode added to the DSM prjs: `path2_P2-1_dsm_mcc`, `path1_P1-3_LE`
  (`path1_P1-3_MCCrestart` is single-structure MCC, no potential_exchange → unchanged).
- **K = 56183 J/kg @ rho_d=1500** — log-linear interpolation of the recalibrated
  current_porosity_split MX-80 Dixon anchors (36796 @1400, 85776 @1600); was 44200
  (interp of the older 26950/71900). DERIVED, not measured (block locus, not granular),
  same guardrail as the original.
- `initial_dt` 10 -> 0.5 (the split mode is stiff at the S_L->1 saturation transition;
  the dnS/dn_l chain is not yet in the local-solve Jacobian — dt is the interim).

## Result (all three completed to t_end)
e vs sigma_v vs Ferrari et al. (2022) — `figures/plot_epfl_cps.py`,
`figures/epfl_cps_both_paths.png`:
- P2 (const-volume): model endpoint e=0.498 @16.5 MPa vs measured C' 0.57 — close.
- P1 free-swell peak: e=1.02 vs measured B 2.34 — small-strain LE limit (was ~1.5
  under reference mode; the switch lowered it).
- P1 compression (MCC, B->D): e=0.308 @20 MPa vs measured D 0.56 — over-compresses.

Verification, not calibration: K was fit to the Dixon swelling-pressure anchors
(ANCHORS MS33 Model I); these EPFL e-sigma_v loci are the independent check.
The full single-shot `path1_P1-3_dsm_mcc` still diverges at the wetting front
(gate-independent, hydraulic-side limit) — hence the LE-swell + MCC-restart split.
