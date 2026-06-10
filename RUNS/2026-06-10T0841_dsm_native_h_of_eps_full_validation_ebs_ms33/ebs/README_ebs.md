# EBS half of the validation bundle — re-run 2026-06-10

Binary: `/Users/vinaykumar/git/build/h_of_eps_20260609/bin/ogs`
(`--version` reports `archive/dsm_consolidate_on_film_2026-06-08-4-gd98f5f83.dirty`;
branch per task spec `dsm_native_h_of_eps` @ 7ff8861847, `film_strain_coupling`
defaults OFF = maxwell_conjugate baseline). `OMP_NUM_THREADS=3` per case, all four
cases launched concurrently (08:33:10–08:33:57 CEST → **total wall time ≈ 47 s**;
the long pole was mgr27 at 41 s — the other cases blocked/failed fast).

Source cases: `/Users/vinaykumar/git/ogs/validation_2026-06-09/{ebs_task13,epfl_t33,mgr23,mgr27}`.
Full output series left in place under `out_<case>/run_<prj>/` in this directory
(uncommitted). Case runner scripts as executed: `_run_<case>.sh` (only binary path,
output dir, and OMP adapted; PRJs untouched).

## Per-case summary

| Case | What it is (case's own description) | PRJ(s) run | Steps | Wall | Status |
|---|---|---|---|---|---|
| ebs_task13 | "EBS Task 13 (Task-13-Description-Stages-1-2-3; Villar et al., CIEMAT report; digitized expt_stg1_05.csv / expt_stg1_6.csv)" — Stage-1 GAP tests at 0.5 and 6 MPa suction, axisymmetric t-13 mesh, MCC (process-common.xml: `ModCamClay_semiExpl_constE`), python top-gap BC | `stg1/12a_t-13_MCC_0.5MPa.prj`, `stg1/12b_t-13_MCC_6MPa.prj` | 0 (parse) | ~1 s | **BLOCKED** (parse) |
| epfl_t33 | "EPFL/Beacon Task-3.3 stress-path-dependency suite" (run_epfl_t33.sh header) — granular MX-80, paths P1/P2 in e–σv space | the 5 PRJs its own `run_epfl_t33.sh` runs | 4×0 (parse); 118 (homog) | ~1 s | **4/5 BLOCKED**, 1/5 ran |
| mgr23 | "BEACON WP5.3 MGR23 — FEBEX pellet+block homogenisation column (native Pi-path DSM, LE) … does the dense block swell into the soft pellet so the system homogenises toward ~1.43 g/cm3 (D4.1 Table 5-3)?" (PRJ header) | `mgr23_column_calibrated.prj`, `_robust.prj`, `_bicg.prj` (all three present in the case's own 2026-06-09 logs) | fails at t≈3.47e6 s of 1.8144e7 s (≈40 d of 210 d) | 9+9+10 s | **FAILED** (all 3 variants, dt-floor) |
| mgr27 | BEACON WP5.3 MGR27 — FEBEX block(bottom, hydrated)/pellet(top) inverted column, LE, calibrated-K ("calk0"); note: the PRJ header comment is an inherited MGR23 text, the mesh/BC/output (`mgr27_predict_calk0`) are MGR27 per `team_comparison_data.py` ("MGR27: block 0-5 cm (bottom, hydrated), pellets 5-10 cm (top) [inverted]") | `mgr27_le_calk0.prj` | 355 accepted / 6 rejected, t_end=1.8144e7 s reached | 41 s | **RAN — completed** |

### BLOCKED detail (do-not-modify-physics rule applied)

All 6 blocked PRJs (both ebs_task13 PRJs; epfl_t33 column_le, path1_P1-3_dsm_mcc,
path2_P2-1_dsm_mcc, path2_P2-1_swellingpressure_dsm_mcc) fail at parse with the
identical error:

```
error: ConfigTree: ... at path <processes/process/potential_exchange>:
Key <vdw_augmentation_decay_length> has been read 1 time(s) less than it was
present in the configuration tree.
critical: .../BaseLib/ConfigTree.cpp:293 assertNoSwallowedErrors()
error: OGS terminated with error on 2026-06-10 ...
```

The h_of_eps worktree source contains no `vdw_augmentation_decay_length` reader
(grep over `ProcessLib/RichardsMechanics/` is empty) — the binary lacks the
feature these PRJs require. The same PRJs failed the same way on the
mc_20260608 binary in the case's own 2026-06-09 records (`EPFL_T33_SUMMARY.txt`:
4×"terminated with error"; `ebs_task13/stg1/out_test05/run.log`: same swallowed
key). So this is a pre-existing block, not a regression introduced today.

### FAILED detail — mgr23

All three numerical variants (calibrated / robust / bicg) abort with
"The new step size … is the same as that of the previous rejected time step"
at t≈3,470,857 s ≈ 40.2 d (4 output VTUs each: 0, 1, 10, 34 d). This mirrors
the case's own 2026-06-09 logs (`mgr23_run.log`, `mgr23_robust_run.log`,
`mgr23_bicg_run.log` — same failure on the mc binary). **Consistency check:**
the reduced metrics from this run's partial VTUs reproduce the case's own
`mgr23_reduced_partial.json` to ≥4 decimals at every common time (see table
below) — the new binary is behaviourally identical on this case up to the
failure point.

## Comparison figures (regenerated 2026-06-10)

Generator: `figures/make_ebs_comparison_figs.py` — adapted ONLY in input paths
from the cases' own scripts/data (`team_comparison_data.py` for MGR experiment/
team values; measured EPFL P2 locus copied unchanged from
`plot_epfl_both_paths_maxwell.py`, Ferrari et al. 2022 Fig 9 / D3.3 Tab 4-1,
digitized). Numbers in `figures/comparison_numbers.json`.

1. `figures/mgr27_profile_team_20260610.png` — MGR27 final dry-density profile
   (this run) vs experiment (D5.6, exact), UPC BExM and BGR-old (digitized).
2. `figures/mgr23_partial_consistency_20260610.png` — MGR23 partial gap and
   axial-stress evolution (this run vs the case's own 2026-06-09 partial json),
   experiment final values for context. NOT a final-state validation — run failed.
3. `figures/epfl_t33_p2_compression_20260610.png` — Path-2 homogeneous MCC
   compression (the only runnable sub-case) vs the measured P2 locus.
   The case's full both-paths figure can NOT be regenerated: 4/5 PRJs blocked.

## Headline numbers (per the cases' own metrics)

### MGR27 (completed)

| quantity | model (this run) | experiment (case's own data file) | verdict |
|---|---|---|---|
| pellet ρ_d final [g/cm³] | 1.342 | 1.434 (D5.6 tab, exact) | model pellet too light by 0.09 |
| block ρ_d final [g/cm³] | 1.550 | 1.454 (exact) | model block too dense by 0.10 |
| gap block−pellet [g/cm³] | **0.208** | **0.02** (exact) | model retains ~10× the measured gap — does NOT homogenise |
| axial p_sw [MPa] | 2.39 (final; plateau from ~10 d) | 1.2 (D5.7 Fig 5.3-13a, digitized ±0.3) | model over-predicts ~2× |

The under-homogenisation is the documented χ=0-family deficit (same verdict as
the case's own 2026-06-01 compilation: prior calk0 record gap 0.202). The VTU
fields of this run are identical (to printed precision, all 6 output times,
σ_zz top/bottom) to the case's own 2026-06-09 result VTUs — binary parity
confirmed on this case.

Note on p_sw: `team_comparison_data.py` records psw=3.21 for "mgr27_predict_calk0";
this run (and the case's own 2026-06-09 VTUs, which it matches exactly) gives a
2.39 MPa final-time top axial stress. The 3.21 value was compiled 2026-06-01 from
the older beacon_wp5 results directory; the discrepancy is between the two
*records*, not between this run and the 2026-06-09 baseline.

### MGR23 (failed at 40 d / 210 d — partial only)

| t [d] | gap this run | gap case record (json) | σ_zz,top this run [MPa] | record |
|---|---|---|---|---|
| 0 | 0.3000 | 0.30000 | 0.0215 | 0.02148 |
| 1 | 0.2983 | 0.29828 | 0.0666 | 0.06657 |
| 10 | 0.2665 | 0.26651 | 1.2698 | 1.26979 |
| 34 | 0.2389 | 0.23888 | 1.7874 | 1.78738 |

Experiment final (exact, D5.6): gap 0.17, p_sw ≈ 3.0 MPa (digitized ±0.3).
No final-state comparison possible from this run.

### EPFL T3.3 (1/5 sub-cases ran)

Path-2 homogeneous MCC compression: model reaches e = 0.4695 at σ_v = 20 MPa vs
measured P2 locus e = 0.57 at 20 MPa (Ferrari 2022 Fig 9, digitized) — the model
ends ~0.10 below the measured void ratio at the maximum stress. Identical output
(9 VTUs, last t = 3.456e6 s) to the case's own 2026-06-09 run of this PRJ.

### ebs_task13

No model output (both PRJs blocked at parse). Measurements
(`expt_stg1_05.csv`, `expt_stg1_6.csv`: t, WC, DD, Sr, gap; 74 data rows each)
copied to `measurements/ebs_task13/` for the bundle; no comparison generated —
nothing to compare against.

## Bundle layout

```
ebs/
  README_ebs.md                this file
  _run_<case>.sh               as-executed case runners
  run_<case>.log               case-level logs (start/end/rc per PRJ)
  out_<case>/run_<prj>/        full outputs + per-PRJ run.log (left in place, uncommitted)
  prj/<case>/                  PRJs as run (ebs_task13 incl. prj-common includes)
  measurements/<case>/         case-own measurement/comparison data files
  figures/                     regenerated figures + generator + comparison_numbers.json
  final_state/                 final (mgr27, epfl) / last-partial (mgr23) VTUs only
```

Size: ~6.5 MB (< 20 MB target). Full series also remain in
`validation_2026-06-09/<case>/results/` from the prior attempt (untouched).
