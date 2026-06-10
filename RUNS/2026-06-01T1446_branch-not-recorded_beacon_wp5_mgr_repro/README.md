<!-- RUN BUNDLE (CLAUDE.md §13) -->
**Bundle:** 2026-06-01T1446_branch-not-recorded_beacon_wp5_mgr_repro
**Run finished:** 2026-06-01 14:46 (+0200) — last sim event in the logs is `OGS terminated with error on 2026-06-01 14:46:37+0200` (the final MCC attempt, `mcc_aM_1e-16_steady`, integrator-blocked); last successful completion `OGS completed on 2026-06-01 10:09:39+0200`.
**Branch / commit:** branch not recorded. Commit: the run logs report version string `vdw-baseline-2026-05-08-74-g62bc63ba` (i.e. `62bc63ba`); no branch name appears in the folder's own records.
**Binary:** `/Users/vinaykumar/git/build/dsm-native-omp-release` (per the folder's AGENTS.md; build includes libOgsMFrontBehaviour + ModCamClay_semiExpl).
**Full raw outputs (untracked, local):** `/Users/vinaykumar/git/ogs/beacon_wp5_mgr_repro_2026-06-01/`
**Report / beamer:** canonical deliverable (per the REPORT banner) = `ogs-models/EBS/Task13/2026_06_01_MGR_FEBEX_pellet_block_homogenisation/deliverables/MGR_FEBEX_homogenisation_{report,beamer}.{tex,pdf}` on `github.com/VinayGK/EBS`, commit `89b8ab5`. Working analysis note = `REPORT_MGR_homogenisation.md` (copied here, marked HISTORICAL/SUPERSEDED by that deliverable).

---

# BEACON WP5 MGR — FEBEX pellet/block homogenisation with the native DSM (2026-06-01)

Snapshot bundle consolidated 2026-06-10 from the campaign folder above. All
content below is distilled from the folder's own `REPORT_MGR_homogenisation.md`
and `AGENTS.md`; no new numbers were computed.

## What was run, and why

Native Pi-path DSM (LE skeleton, chi=0 BishopsSaturationCutoff) on the BEACON
MGR family: half loose FEBEX pellets (rho_d ~1.30) / half compacted block
(rho_d ~1.60), hydrated at constant total volume — does the density step
homogenise? Two legs kept strictly apart (calibrated ≠ identified ≠
validated): **MGR23 calibration** (block k0 swept x10^4 up to 5e-15 m2,
lambda=0) and **MGR27 blind prediction** (geometry inverted, parameters
frozen). Plus alpha_M and lambda sweeps, a FEBEX-K calibration script, and a
FEBEX-sourced MCC skeleton built to test the wetting-collapse attribution.

## Headline results (quoted from the folder's own report)

- **Calibration leg — permeability is not the lever:** block k0 x10^4 closes
  the final gap by only 0.006 (0.192 -> 0.186, plateau), well short of the
  measured 0.17; block interlayer already at ceiling (n_l ~0.43). No matched
  k0 exists to freeze.
- **Blind leg MGR27 — delivery fixed, gap unmoved:** max-k0 wets the far
  pellet fully (n_l 0.0004 -> 0.495), but the gap goes 0.185 -> 0.202 vs
  measured 0.02; axial sigma_zz 2.79-3.21 MPa vs ~1.2 measured
  (friction-affected). **The sign is backwards:** the model's wet pellet
  swells (rho_d 1.353 -> 1.345) where the experiment's collapses.
- **Verdict:** the native DSM sits in the BGR-old chi=0 lineage — same
  homogenisation failure (gap stalls ~0.19-0.20 vs UPC/CODE_BRIGHT which
  captures it), now diagnosed: a transport floor (gateable) hiding a
  mechanical floor (residual after full hydration). Headline: "MGR is outside
  the chi=0 envelope, and the residual is mechanical collapse, not retention."
- On stress the native DSM is the healthier BGR generation (right order via
  the Pi-path; MGR27 over-prediction is the un-modelled wall friction,
  shared with UPC).

## Caveats (as flagged in the folder's own records)

- **Wetting-collapse attribution remains open/untested:** initially
  confounded by an un-cited uniform E = 52 MPa for both materials (§4
  confound flag); the corrective FEBEX-sourced MCC set (ENRESA 2000 Cc/Cs +
  phi; UPC Po*) is built and correct but **integrator-blocked** — 9 MCC
  configurations wall at the swelling-stress onset (semiExpl family only in
  this build). Real fix = robust/implicit MCC integrator (a development task).
- The markdown report here is HISTORICAL/SUPERSEDED — the canonical numbers
  live in the Task-13 deliverable (banner above).

## Bundle contents / omissions

Copied: AGENTS.md, REPORT_MGR_homogenisation.md, sweep/calibration/comparison
scripts, `spec/` extraction + parameter-proposal .md files, `model/` (PRJs +
meshes as run), `figures/`, results summary TXTs, run logs (>1 MB tailed to
`*.log.tail`, 200 lines — incl. the 125 MB `mgr23f_run.log` and 33 MB
`mcc_aM_1e-16_steady_run.log`), `final_state/` = largest-t VTU per variant per
run dir (31 files).
Omitted (size policy, in the original folder only): full VTU time series,
`spec/` deliverable PDFs/txt (D3.3/D4.1/D5.x, ~15 MB), `__pycache__`.
