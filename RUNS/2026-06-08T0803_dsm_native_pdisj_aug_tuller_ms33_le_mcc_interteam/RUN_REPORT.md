# MS33 LE+MCC on `dsm_native_pdisj_aug_tuller` — run + inter-team comparison

**Date:** 2026-06-08 · **Branch:** `dsm_native_pdisj_aug_tuller` (worktree `ogs-worktrees/ufz_integration_wt`, HEAD `ae96b70086`)
**Binary:** `/Users/vinaykumar/git/build/ufz_integration_20260602/bin/ogs` (OGS `6.5.8-15-g9bd1c32b.dirty`, MFront on, serial/OMP)
**Constitutive:** native DSM (Π-path disjoining-pressure swelling + vdW augmentation) + **SaturationTuller** macro-WRC (char. pore size 1e-5 m, cavitation p_cav=1.4e8 Pa), BishopsSaturationCutoff(1).
All PRJs copied from the branch (LE) / ported from `ANCHORS_MS33_MCC_NATIVE` audit tree (MCC); only mesh paths fixed, **no parameter edits**.

## 1. Convergence scorecard

| Model | variant | result | t reached | saturated $P_s$ / note |
|------|---------|--------|-----------|------------------------|
| I dd1400 | LE | ✓ completed | 200 d | 5.03 MPa |
| I dd1600 | LE | ✗ diverged in hold (step 104) | 22 d | **14.44 MPa plateau captured** (reached Sl=1 before crash) |
| I dd1800 | LE | ✓ completed | 200 d | 40.86 MPa |
| III gap2mm | LE | ✓ completed (slow crawl ~8 min) | 200 d | center p=8.78 MPa |
| IV pellets | LE | ✓ completed (per-material K) | 200 d | center p=5.52 MPa (was 9.12 at single K) |
| VII freeswelling | LE | ✓ completed | 240 d | center p=0.15 MPa (seating) |
| I dd1400 | MCC | ✓ completed | 200 d | 4.92 MPa (elastic) |
| I dd1600 | MCC | ✓ completed | 200 d | 13.24 MPa (plastic — yields) |
| III gap2mm | MCC | ✓ completed | 200 d | center p=6.94 MPa |

**LE 5/6 · MCC 3/3.** Per the agreed scope, MCC was run only for the cases known to converge (I dd1400/1600, III); MCC dd1800/IV/VII were not run (documented tension-apex non-convergence).

### Guardrail / numerical notes (flagged, not silently fixed)
- **dd1600 LE divergence** is a genuine numerical failure on this branch (singular-Jacobian "Eigen linear solver initialization failed", dt collapsed to floor) at t≈22 d **in the saturated hold** — *after* the swelling-pressure plateau was reached. dd1400/dd1800 LE complete; the older `vdw-baseline-2026-05-08` binary completed all three. I did **not** tune the time-stepper/tolerances to rescue it (numerical decision = yours). The plateau value (14.44 MPa) is physically valid and used for the swelling-pressure figure, annotated as plateau-captured.
- **MCC converges at dd1600 where LE diverges** (plastic yielding regularises the state that made LE singular).
- LE and MCC use the **same Tuller retention** → the LE-vs-MCC comparison is clean (no retention confound).

## 2. Validation against the calibration anchor (Dixon 2023, EMDD≡ρ_d)

| ρ_d (g/cm³) | BGR-LE $P_s$ | Dixon (2023) experimental | dev |
|---|---|---|---|
| 1.40 | 5.03 MPa | 5.0 MPa | +0.6% |
| 1.60 | 14.44 MPa | 14.16 MPa | +2.0% |
| 1.80 | 40.86 MPa | 40.0 MPa | +2.2% |

The calibration **survives the port** to this branch (all within ~2%). Note BGR sits *above* the Villar/Lloret Eq.(7) line `exp(6.77ρ−9.07)` (1.5/5.8/22.5 MPa) because the calibration target is the EMDD≡ρ_d Dixon set, not the literal-axis Villar reading (see `[[feedback_emdd_dry_density_convention]]`). Both references are drawn on the figure.

## 3. Figures (in `figures/`, 220 dpi, branch+model in footnote)

| file | content | teams overlaid |
|---|---|---|
| `modelI_swelling_pressure_vs_density.png` | $P_s$ vs ρ_d; Villar Eq.(7) + Dixon squares + BGR LE/MCC | 15 |
| `modelI_suction_stress_path.png` | suction[log] vs mean stress, 3 ρ_d, LE+MCC | 15 |
| `modelI_permeability_vs_suction.png` | k[log] vs suction[log], 3 ρ_d | 15 |
| `modelIII_interteam.png` | mean stress vs time (centre) + gap closure vs time | 4 / 1 |
| `modelIV_interteam.png` | mean stress vs time + dry-density (clay vs pellet) | 3 |
| `modelVII_interteam.png` | void ratio vs time + axial stress vs time | 2 |

### Inter-team observations (descriptive)
- **Model I:** BGR LE/MCC land inside the team swelling-pressure cluster and on the Dixon anchors. BGR confined-swelling paths are near-vertical (stress builds only near saturation — the Π-path signature) vs more gradual team paths.
- **Model III:** BGR mean stress (~7–9 MPa) sits in the team band; MCC below LE (yielding).
- **Model IV:** **corrected to per-material K** (clay 85312.6 / pellet 13064 J/kg via `<medium id="1">`, pellet calibrated to the Dixon-consistent 0.350 MPa — see §6). Center mean stress dropped 9.12→**5.52 MPa**; dry-density evolution now shows the proper heterogeneity — clay (top, ρ_d 1.60→1.51) de-densifies while the pellet (bottom, 0.90→0.96) compacts under it (partial homogenisation). Still above the team band (~1–2 MPa): the DSM Π-path generates more swelling stress than the BBM-type team models (genuine model spread, not a setup error now).
- **Model VII:** BGR void ratio (0.75→1.0) and the loading-phase axial-stress spike (~5 MPa near 210 d) track the team cluster.

## 4. Caveats
- **Model III gap closure for BGR is a proxy** = |radial displacement u_r| at the gap interface (r=23 mm), not a contact-mechanics gap aperture (absent in OGS RM). Teams report true gap closure. Labeled as proxy on the figure.
- Team curves with no readable column for a quantity are silently absent (gap closure: only CTU-CU; VII: EPFL master-sheet not in the standard template) — see `figures/PROVENANCE.md`.
- BGR multi-element center probes are still largely *dry* (Sl≈0) at end-time — hydration from the boundary is slow under Tuller; the comparison is at the reported center location regardless.

## 5. Reproduce
```
cd /Users/vinaykumar/git/ogs/ms33_pdisj_aug_tuller_2026-06-08
# runs (grinder-watched): ./run_with_watchdog.sh <prj> <outdir> <max_wall_s> <stall_s> <omp>
HOME=/tmp python3 summarize_runs.py        # VTU -> reduced/*_history.csv (center probe)
HOME=/tmp python3 fig_modelI.py            # Model I figures
HOME=/tmp python3 fig_III_IV_VII.py        # III/IV/VII figures
```

## 6. Model IV — per-material pellet K (correction, 2026-06-08)

The pellet-clay model is a two-medium problem and the DSM swelling prefactor K
**must differ per material**; the branch PRJ applied a single global
`K=85312.6` (clay) to the pellet too, over-swelling it (center p=9.12 MPa).

**Fix.** Added a per-material `<medium id="1">` override inside
`<potential_exchange>` (the branch code parses these as partial overrides on the
global block — `CreateRichardsMechanicsProcess.cpp:623`, looked up by
`material_id`): pellet `K=13064 J/kg`, `n_s=0.32374`, `n_l0=6.59e-4` (the latter
two geometric, from the spec pellet ρ_d=0.9 → phi0=0.6763). Clay (MaterialID 0)
keeps the global K.

**Calibration target corrected (guardrail §12.1/§1.1).** The PRJ header listed
the pellet target as "Dixon ~0.051 MPa", but 0.051 is the *Villar Eq.7* value;
the whole suite (clay + dd1400/1600/1800) is calibrated to **Dixon EMDD≡ρ_d**
(σ=0.003·exp(5.2883·EMDD)), which at ρ_d=0.9 gives **0.350 MPa**. Vinay
confirmed the Dixon basis. K_pellet=13064 was secant-fit (single-element,
confined swelling, **swelling_stress** field — the σ0=−1.5e5 Pa seating dwarfs
the target at this density, so total mean stress is not the right metric here)
to 0.350 MPa (0.0% dev). K still extrapolated below Dixon's measured 1.4–1.8
range (edge case, approved). Calibration on a separate single-element case — the
Model IV deliverable is not derived from this K, so no fit-and-verify crossing (§2).

**Result.** Center mean stress 9.12→5.52 MPa; clay (top) ρ_d 1.60→1.51 and
pellet (bottom) 0.90→0.96 (partial homogenisation now visible). Scripts:
`calibrate_pellet_K.py`, run `calib_pellet/`.
