<!-- RUN BUNDLE (CLAUDE.md §13) -->
**Bundle:** 2026-06-02T2330_branch-not-recorded_beacon_wp3_epfl_stresspath
**Run finished:** 2026-06-02 23:30 (+0200) — last `OGS completed on 2026-06-02 23:30:15+0200` in the run logs
**Branch / commit:** not recorded. The run logs report version string `maxwell-conjugate-pre` (no hash); the control runs report `vdw-baseline-2026-05-08-74-g62bc63ba.dirty`. No branch name or commit hash for the main binary appears in the folder's own records. (The successor folder calls this source "pre-`d542d6c08c`", but that is not recorded here.)
**Binary:** `/Users/vinaykumar/git/build/maxwell-conjugate-20260602/bin/ogs`; control = `/Users/vinaykumar/git/build/dsm_native_hier_wt-release/bin/ogs` (per FINDINGS.md). OMP_NUM_THREADS=4.
**Full raw outputs (untracked, local):** `/Users/vinaykumar/git/ogs/beacon_wp3_epfl_stresspath_2026-06-02/`
**Report / beamer:** TODO (no formal report/beamer). Findings note = `FINDINGS.md` (copied here).

---

# EPFL Task-3.3 stress-path suite on the Maxwell-conjugate binary (2026-06-02)

Snapshot bundle consolidated 2026-06-10 from the campaign folder above. All
content below is distilled from the folder's own `FINDINGS.md` and
`results/EPFL_MAXWELL_SUMMARY.txt`; no new numbers were computed.

## What was run, and why

The EPFL Task-3.3 stress-path PRJ suite (ported unchanged from
`beacon_wp3_epfl_repro_2026-05-31/model`, only the retired
`<accumulate_swelling_contributions>` tag stripped; K = 44200 J/kg Dixon 2023
carried over) on the **Maxwell-conjugate binary** (full-Pi closure + the
Maxwell-conjugate micro-potential term, residual-only; analytic Jacobian block
not yet wired). Purpose: see whether the coded gate p' >= Pi fires on a real
wetting+loading stress path.

## Headline results (quoted from FINDINGS.md)

- Run outcomes 4/5 complete: column_le (200 d), path2_P2-1, path2_P2-1
  swellingpressure (240 d each), path2 compression homogeneous (40 d) all
  completed; **path1_P1-3 diverged** at 918 685 s (~10.6 d), step 62.
- **The coded gate as written is inert on the real EPFL path:** on path2_P2-1,
  max p' = 10.38 MPa vs Pi = 22-28 MPa — the gate `p' >= Pi` **never fires**;
  the REV-consistent gate `p' >= phi_m*Pi` (9.2-10.2 MPa) **fires near peak
  load**. Third independent confirmation that the coded gate mixes scales;
  strong empirical support for gate-scale decision #4, option 1. (Observed,
  not yet a re-run with the rescaled gate.)
- **path1 divergence attributed, NOT the Maxwell term:** the control re-run on
  `dsm_native_hier_wt-release` (same full-Pi closure, no Maxwell term)
  diverges at the bit-identical time 918684.609 s — intrinsic to the full-Pi
  closure + MFront-MCC on the P1-3 path; consistent with the term being inert
  there (Pi ~ 203 MPa at the dry initial state).

## Caveats

- No parameter changed vs the 05-31 baseline; this is a binary/term probe,
  not a recalibration.
- Gate firing is observed via post-processing of phi_m*Pi, not via a rescaled
  gate in code — that re-run is the successor folder
  (`beacon_wp3_epfl_stresspath_termON_phimPi`).

## Bundle contents / omissions

Copied: FINDINGS.md, run_epfl_maxwell.sh, `model/` (PRJs + meshes as run),
`figures/` (plot + script), results summary TXT, run logs (>1 MB tailed to
`*.log.tail`, 200 lines), `final_state/` = largest-t VTU per variant (8 files).
Omitted: full VTU time series (in the original folder only).
