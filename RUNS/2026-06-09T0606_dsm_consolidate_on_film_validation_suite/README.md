<!-- RUN BUNDLE (CLAUDE.md §13) -->
**Bundle:** 2026-06-09T0606_dsm_consolidate_on_film_validation_suite
**Run finished:** 2026-06-09 06:06 (+0200) — last `OGS started/terminated` in `mgr23/logs/mgr23_bicg_run.log` (06:06:18); last successful completion `mgr27` at 06:05:59.
**Branch / commit:** `dsm_consolidate_on_film` @ `295dc649` — from the version string printed in every run log: `archive/dsm_consolidate_on_film_2026-06-08-2-g295dc649`. (The mgr27 PRJ schema-port comment also names the binary `mc_20260608`.)
**Binary:** `/Users/vinaykumar/git/build/mc_20260608/bin/ogs` — from `epfl_t33/run_epfl_t33.sh` and `EPFL_T33_SUMMARY.txt`; same version string in the mgr/ebs logs.
**Full raw outputs (untracked, local):** `/Users/vinaykumar/git/ogs/validation_2026-06-09/` (full VTU series under `mgr23/results/{.,robust,bicg}/`, `mgr27/results/`, `epfl_t33/results/run_*/`).
**Report / beamer:** not recorded — the campaign folder carries no report/beamer pointer; the only in-folder summaries are `epfl_t33/EPFL_T33_SUMMARY.txt` and `mgr23/mgr23_reduced_partial.json` (both copied here). TODO.

---

# Validation suite on `dsm_consolidate_on_film` (mc_20260608) — four cases, one campaign (2026-06-09)

Snapshot bundle consolidated 2026-06-10 from the historical campaign folder
above. One campaign, four validation cases, each in its own subfolder. All
four ran on the same binary (version string above) within 06:02–06:06 on
2026-06-09. Everything below is distilled from the folder's own records (run
script headers, PRJ provenance headers, logs, in-folder summaries); numbers
quoted, not recomputed. The campaign folder contains no overall report — the
per-case purpose statements below come from the PRJ headers and script
comments themselves.

## Cases and outcomes (from logs / in-folder summaries)

| case | what (per its own records) | outcome |
|---|---|---|
| `ebs_task13` | EBS Task-13 stage-1 oedometer-type case, MCC, two PRJs (12a 0.5 MPa, 12b 6 MPa) with shared `prj-common/` includes + experimental CSVs | FAILED at PRJ parse: `Key <vdw_augmentation_decay_length> has been read 1 time(s) less than it was present` (schema mismatch on this binary). No output VTUs. |
| `epfl_t33` | EPFL/Beacon Task-3.3 stress-path-dependency suite; PRJs reused verbatim from `beacon_wp3_epfl_stresspath_2026-06-02` (maxwell-conjugate port), "NO config adapted"; K=44200 J/kg (Dixon 2023 MX-80 block locus, ρ_d=1500 log-interp; per-PRJ §12 header) — from `run_epfl_t33.sh` | 1/5 completed: `beacon_t33_path2_compression_homogeneous_mcc` (rc=0, Time: 3456000 s, 9 VTUs). The other 4 PRJs terminated at parse with the same `vdw_augmentation_decay_length` schema error (per `EPFL_T33_SUMMARY.txt` + logs). |
| `mgr23` | BEACON WP5.3 MGR23 — FEBEX pellet+block homogenisation column (native Π-path DSM, LE); first FEBEX parameterisation of the native DSM; K calibrated block 25268 / pellet 3979 J/kg vs Villar/ENRESA Ps law (per-PRJ §12.2 provenance header, quoted in full in `mgr23/model/`) | All 3 variants (calibrated / robust / bicg) TERMINATED with nonlinear-solver divergence at the dt floor (0.1 s). Partial series reach ts_122, t = 2 937 600 s (34 d); partial homogenisation history extracted to `mgr23_reduced_partial.json` (last row, t = 34 d: pellet ρ_d 1.30→1.328, block 1.60→1.567, density gap 0.30→0.239 g/cm³, σ_zz top 1.79 MPa). |
| `mgr27` | BEACON WP5.3 MGR27 — inverted block(bottom, hydrated)/pellet(top) column, LE, calibrated-K prediction (`mgr27_le_calk0.prj`; SCHEMA PORT 2026-06-09 comment in the PRJ); cross-team comparison data with per-value EXACT/DIGITIZED provenance in `team_comparison_data.py` | COMPLETED: ts_355, t = 18 144 000 s (210 d), 37.7 s wall. |

## Caveats (from the folder's own records)

- The recurring `vdw_augmentation_decay_length` ConfigTree error (ebs_task13 +
  4/5 epfl_t33 PRJs) is a PRJ-schema mismatch against this binary, not a
  physics failure — the PRJs were ported from older branches; only
  `mgr27_le_calk0.prj` carries an explicit "SCHEMA PORT 2026-06-09" fix note.
- mgr23 PRJ header flags (Vinay): native Tuller macro retention kept for
  tractability instead of the ENRESA FEBEX vG WRC; equilibrium density is
  WRC-insensitive under BishopsSaturationCutoff=1.
- `team_comparison_data.py` marks each cross-team value EXACT vs DIGITIZED
  (±0.02 g/cm³ density, ±0.3 MPa stress); digitized values are
  comparison/plotting only, never parameters.
- No figures exist in the campaign folder (none were generated or none were
  kept there — not recorded).

## Bundle contents / omissions

- `ebs_task13/` — both PRJs, `prj-common/` includes, meshes, experimental
  CSVs, parse-failure log.
- `epfl_t33/` — `run_epfl_t33.sh`, all 5 PRJs + both meshes (column + single
  element) + `make_column_boundaries.py`, `EPFL_T33_SUMMARY.txt`, all 5 run
  logs.
- `mgr23/` — 3 PRJ variants, meshes, 3 run logs, `mgr23_reduced_partial.json`.
- `mgr27/` — PRJ, meshes, run log, `team_comparison_data.py`.
- `final_state/<case>/` — last available VTU per variant (mgr27 ts_355 final;
  mgr23 ts_122 last-before-abort for each of calibrated/robust/bicg; epfl_t33
  ts_118 final of the one completed run). OMITTED: full VTU series and `.pvd`
  indices (local only, path above). ebs_task13 produced no VTUs.
- All logs < 1 MB, copied in full.

Bundle size ≈ 5 MB (< 30 MB policy).
