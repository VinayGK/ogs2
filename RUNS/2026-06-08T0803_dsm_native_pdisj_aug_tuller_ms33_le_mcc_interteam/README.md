<!-- RUN BUNDLE (CLAUDE.md §13) -->
**Bundle:** 2026-06-08T0803_dsm_native_pdisj_aug_tuller_ms33_le_mcc_interteam
**Run finished:** 2026-06-08 08:03 (+0200) — last `OGS completed on` line in the campaign logs (`LE/ModelIV/out_perK/run.log`). Post-processing figures `ms33_3x3_50x*.png` were regenerated later the same day (file mtimes up to 22:15).
**Branch / commit:** `dsm_native_pdisj_aug_tuller` @ `ae96b70086` (worktree `ogs-worktrees/ufz_integration_wt`) — from RUN_REPORT.md header.
**Binary:** `/Users/vinaykumar/git/build/ufz_integration_20260602/bin/ogs` (OGS `6.5.8-15-g9bd1c32b.dirty`, MFront on, serial/OMP) — from RUN_REPORT.md header.
**Full raw outputs (untracked, local):** `/Users/vinaykumar/git/ogs/ms33_pdisj_aug_tuller_2026-06-08/` (full VTU time series under `LE/Model*/…_out/`, `MCC/Model*/…_out/`, and all `calib_pellet/out_K*/` sweeps).
**Report / beamer:** in-folder report = `RUN_REPORT.md` (copied here). Tex report/beamer: not recorded in the campaign folder (MS33 tex artifacts live under `~/tex/cc2024/VK_SB_EURAD_DSM/` per project memory, but this folder itself carries no pointer — TODO).

---

# MS33 LE+MCC on `dsm_native_pdisj_aug_tuller` — run + inter-team comparison (2026-06-08)

Snapshot bundle consolidated 2026-06-10 from the historical campaign folder
above. Everything below is distilled from the folder's own documents
(`RUN_REPORT.md`, `TEAM_DATA_MAP.md`, `figures/PROVENANCE.md`); numbers are
quoted, not recomputed.

## What was run, why

EURAD-2 MS33 theoretical benchmark Models I (dd1400/1600/1800), III (gap 2 mm),
IV (clay+pellets), VII (free swelling) with the native DSM (Π-path
disjoining-pressure swelling + vdW augmentation) + SaturationTuller macro-WRC
(char. pore size 1e-5 m, p_cav = 1.4e8 Pa), BishopsSaturationCutoff(1), in two
constitutive variants: LE and MCC (native). PRJs copied from the branch (LE) /
ported from the `ANCHORS_MS33_MCC_NATIVE` audit tree (MCC); only mesh paths
fixed, no parameter edits. Purpose: convergence scorecard on this branch,
validation against the Dixon (2023, EMDD≡ρ_d) calibration anchor, and
inter-team comparison figures (team data mapping in `TEAM_DATA_MAP.md`).

## Headline results (quoted from RUN_REPORT.md)

- **Convergence: LE 5/6, MCC 3/3** (MCC deliberately run only for cases known
  to converge: I dd1400/1600, III).
- **Model I LE vs Dixon (2023) anchor:** 5.03 / 14.44 / 40.86 MPa at
  ρ_d = 1.40/1.60/1.80 g/cm³ vs 5.0 / 14.16 / 40.0 MPa experimental
  (+0.6% / +2.0% / +2.2%) — "the calibration survives the port to this branch".
- **I dd1600 LE diverged in the saturated hold** (step 104, t≈22 d,
  singular-Jacobian "Eigen linear solver initialization failed") *after* the
  swelling-pressure plateau (14.44 MPa) was reached; plateau value used,
  annotated. Not rescued by tuning (numerical decision left to Vinay). MCC
  converges at dd1600 where LE diverges.
- **Model IV corrected to per-material K** (clay 85312.6 / pellet 13064 J/kg
  via `<medium id="1">` override; pellet K secant-fit to the Dixon-consistent
  0.350 MPa target on a separate single-element case — no fit-and-verify
  crossing per §2). Center mean stress 9.12 → 5.52 MPa; clay ρ_d 1.60→1.51,
  pellet 0.90→0.96 (partial homogenisation). Still above the team band
  (~1–2 MPa) — genuine model spread per the report.
- **Model III:** center p = 8.78 MPa (LE), 6.94 MPa (MCC), inside the team
  band. **Model VII:** void ratio 0.75→1.0, loading-phase axial-stress spike
  ~5 MPa near 210 d, tracks the team cluster.

## Caveats (quoted)

- Model III "gap closure" for BGR is a proxy (|u_r| at r=23 mm), not a
  contact-mechanics aperture; labeled on the figure.
- Team curves with unreadable columns silently absent (see
  `figures/PROVENANCE.md`).
- BGR multi-element center probes largely dry (Sl≈0) at end-time — hydration
  under Tuller is slow.

## Bundle contents / omissions

- `RUN_REPORT.md`, `TEAM_DATA_MAP.md` — campaign report and team-data mapping.
- `scripts/` — all post-processing/calibration scripts
  (`summarize_runs.py`, `fig_modelI.py`, `fig_III_IV_VII.py`, `fig_3x3.py`,
  `fig_3x3_50x.py`, `gen_provenance.py`, `calibrate_pellet_K.py`,
  `run_with_watchdog.sh`).
- `figures/` — all campaign figures + `PROVENANCE.md`. The `ms33_3x3_50x*.png`
  pair postdates `RUN_REPORT.md` (generated evening 2026-06-08 by
  `fig_3x3_50x.py`); the campaign folder carries no prose describing them.
- `reduced/` — center-probe history CSVs (input to all figures).
- `prj/` — PRJs **as run** + input meshes, preserving the `LE/`, `MCC/` layout.
- `calib_pellet/` — pellet-K calibration PRJ template, `run.prj`, meshes,
  fitted `K_pellet.txt`, the accepted-K (1.306e4) run log, and its final-state
  VTU. OMITTED: the 13 other `out_K*` sweep directories (full series + logs of
  the secant iterations) — local only, see raw-outputs path above.
- `final_state/` — final-time-step VTU per variant (largest `_t_` per output
  dir), prefixed with the variant path. OMITTED: full VTU time series (local
  only).
- `logs/` — full `run.log` where < 1 MB, else `tail -200` (`*.tail200`:
  LE III/IV/IV_perK/VII, MCC III).
- OMITTED: `__pycache__/`, `.DS_Store`.

Bundle size ≈ 7 MB (< 30 MB policy).

## Reproduce (quoted from RUN_REPORT.md §5)

```
cd /Users/vinaykumar/git/ogs/ms33_pdisj_aug_tuller_2026-06-08
# runs (grinder-watched): ./run_with_watchdog.sh <prj> <outdir> <max_wall_s> <stall_s> <omp>
HOME=/tmp python3 summarize_runs.py        # VTU -> reduced/*_history.csv (center probe)
HOME=/tmp python3 fig_modelI.py            # Model I figures
HOME=/tmp python3 fig_III_IV_VII.py        # III/IV/VII figures
```
