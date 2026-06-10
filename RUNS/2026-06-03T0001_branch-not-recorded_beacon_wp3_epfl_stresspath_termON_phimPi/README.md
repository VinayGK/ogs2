<!-- RUN BUNDLE (CLAUDE.md §13) -->
**Bundle:** 2026-06-03T0001_branch-not-recorded_beacon_wp3_epfl_stresspath_termON_phimPi
**Run finished:** 2026-06-03 00:01 (+0200) — last `OGS completed on 2026-06-03 00:01:48+0200` in the run logs
**Branch / commit:** branch not recorded; commit `d542d6c08c` (recorded in the folder's COMPARISON_FINDINGS.md for the rescaled-gate binary). Run logs report version string `maxwell-conjugate-pre`.
**Binary:** `/Users/vinaykumar/git/build/maxwell-conjugate-20260602/bin/ogs`, rebuilt 2026-06-02 23:56 with the gate rescale (per COMPARISON_FINDINGS.md). OMP_NUM_THREADS=4.
**Full raw outputs (untracked, local):** `/Users/vinaykumar/git/ogs/beacon_wp3_epfl_stresspath_termON_phimPi/`
**Report / beamer:** TODO (no formal report/beamer). Findings note = `COMPARISON_FINDINGS.md` (copied here).

---

# EPFL Task-3.3 — term-ON (p' >= phi_m*Pi) vs term-OFF baseline (2026-06-02/03)

Snapshot bundle consolidated 2026-06-10 from the campaign folder above. All
content below is distilled from the folder's own `COMPARISON_FINDINGS.md`; no
new numbers were computed.

## What was run, and why

Decision #4 = option 1: the Maxwell-conjugate gate threshold rescaled from the
intrinsic micro Pi to the REV-consistent partial stress phi_m*Pi = n_S*n_l*Pi
(commit `d542d6c08c`; trigger-only, S1 keeps the full Pi). Same ported EPFL
PRJs as the term-OFF baseline folder (`beacon_wp3_epfl_stresspath_2026-06-02`),
re-run on the rebuilt (rescaled-gate) binary to see whether the term now acts.

## Headline results (quoted from COMPARISON_FINDINGS.md)

- **term-ON ≡ term-OFF, bit-for-bit:** max |Delta| over ALL fields, ALL steps
  = 0.000e+00 for both `path2_P2-1_dsm_mcc` and
  `path2_P2-1_swellingpressure_dsm_mcc`. The e-sigma figure is identical to
  the baseline (`epfl_both_paths_termON_IDENTICAL.png`).
- `path1_P1-3_dsm_mcc` still diverges at the identical time (918 685 s) — the
  gate-independent wetting-front closure limit.
- **Why the rescaled term still does not bite:** p' hugs the phi_m*Pi
  operating point along the whole path (`gate_trajectory_path2.png`); the REV
  gate is open only where the term cannot act (start: eps_v ~ 0; end: only the
  last 1-2 steps, and the explicit gate reads the previous converged stress so
  it lags below threshold). Empirically confirms the note's §6 prediction:
  phi_m*Pi ~ |sigma_sw| ~ the swelling pressure, so p' equilibrates AT the
  gate rather than crossing it. The rescale (option 1) is correct and
  necessary, but the EPFL paths do not exercise the term.

## What it would take to demonstrate the term firing (recorded, not yet done)

1. A path loaded well past the swelling pressure (phi_m*Pi ~ 8-10 MPa here)
   with live strain — BC magnitude is Vinay's call.
2. And/or fix the explicit-lag gate (evaluate on the current iterate).
3. And/or fix path1_P1-3 convergence (the wetting-front limit).

## Bundle contents / omissions

Copied: COMPARISON_FINDINGS.md, `model/` (PRJs + meshes as run), `figures/`
(plots + scripts incl. gate_trajectory_path2), run logs (>1 MB tailed to
`*.log.tail`, 200 lines), `final_state/` = largest-t VTU per variant (7 files).
Omitted: full VTU time series (in the original folder only). The term-OFF
control lives in the sibling bundle
`RUNS/2026-06-02T2330_branch-not-recorded_beacon_wp3_epfl_stresspath/`.
