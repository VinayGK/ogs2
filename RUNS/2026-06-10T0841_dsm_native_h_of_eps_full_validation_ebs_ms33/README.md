<!-- RUN BUNDLE (CLAUDE.md §13) -->
**Bundle:** 2026-06-10T0841_dsm_native_h_of_eps_full_validation_ebs_ms33
**Run finished:** 2026-06-10 08:41 (+0200)
**Branch / commit:** `dsm_native_h_of_eps` @ `7ff8861847`
**Binary:** `~/git/build/h_of_eps_20260609/bin/ogs` — NOTE: its embedded version
string reads `...dsm_consolidate_on_film-...-gd98f5f83.dirty` because it was
built from the pre-commit working tree; the dirty content IS what was then
committed as `7ff8861847` (verified in-session: 28/28 unit tests on this
binary + dd1400 off-mode bit-for-bit vs a clean d98f5f8324 build).
**Mode:** `film_strain_coupling = off` everywhere (the validated default,
bit-for-bit the maxwell_conjugate baseline d98f5f8324).
**Full raw outputs (untracked, local):** `~/git/ogs/_full_outputs/2026-06-10T0841_full_validation/`
**Report / beamer:** TODO (this README is the run record).

---

# Full validation: MS33 suite + EBS measurement cases, one campaign

Requested by Vinay 2026-06-10 ("run the full EBS and MS33 tests, all in one
folder including results and comparison to measurements").

## MS33 (see `ms33/README_ms33.md` for the full record)

All six runs completed, zero rejected steps, total wall 165 s:

| model | result | anchor / teams |
|---|---|---|
| I dd1400 | σ_zz = −4.9218 MPa | Dixon (2023) target 4.922 (PRJ §12 header) |
| I dd1600 | σ_zz = −14.1602 MPa | target 14.161 |
| I dd1800 | σ_zz = −40.8448 MPa | target 40.86 |
| III gap2mm | mean stress T/C/B = 9.64/9.63/9.61 MPa | teams 2.6–17.1 (spread) |
| IV pellets | 6.49/6.49/1.21 MPa | teams 1.5–8.6 (spread) |
| VII free-swell | e_end = 1.4954 (e@200d 1.4995) | teams e_end ≈ 1.09 (over-swelling, known) |

Inter-team overlay figures regenerated with the canonical generator
(team-xlsx reading untouched): `ms33/figures/`.

## EBS measurement cases (see `ebs/README_ebs.md`)

First pass (as-found PRJs): mgr27 completed (binary parity with the
2026-06-09 record); mgr23 failed at its documented dt-floor point (partial
parity ≥4 decimals); ebs_task13 + 4/5 epfl_t33 PRJs were parse-BLOCKED by the
pre-rename schema (`vdw_augmentation_*`).

Schema-migrated re-run (`*_schema20260610.prj`, pure tag rename
`vdw_augmentation_prefactor→potential_augmentation_prefactor`,
`vdw_augmentation_decay_length→potential_augmentation_exponent`, values
unchanged; logs `ebs/run_*.log`, outputs `ebs/out_schema_rerun/`):

| PRJ | outcome |
|---|---|
| epfl beacon_t33_column_le | COMPLETED (12.6 s) |
| epfl path2_P2-1_dsm_mcc | COMPLETED (0.13 s — short restart-stage schedule, as authored) |
| epfl path2_P2-1_swellingpressure_dsm_mcc | COMPLETED (0.13 s, ditto) |
| ebs_task13 12a (0.5 MPa) | FAILED step 1 — MFront integration status −1 |
| ebs_task13 12b (6 MPa) | FAILED step 1 — MFront integration status −1 |
| epfl path1_P1-3_dsm_mcc | FAILED step ~62 (t≈9.8e5 s) — MFront status −1 |

The three failures are NOT schema artifacts: they reproduce the documented
DSM pathology analysed in
`~/tex/implementation_and_theory/implementation/film_pressure_imbibition_problem.tex`
— the Task-13 step-1 Π-magnitude blow-up (dry-IC facet, proposed fix = 1W
crystalline floor / Π cap) and the rate-triggered facet (proposed fix =
Route R, load → expel water). Both fixes are OPEN (Vinay's call); the
`micro_water_content_floor` mechanism (d98f5f8324) exists but the 1W value
is not set in these PRJs.

## Measurement comparisons (EBS, the cases' own metrics)

- **MGR27**: block−pellet density gap 0.208 g/cm³ vs measured 0.02 (the
  documented χ=0-family homogenisation deficit); axial p_sw 2.39 MPa vs
  measured ~1.2 (digitized ±0.3). Figures in `ebs/figures/`.
- **MGR23**: no final state (40 d of 210 d); partial at 34 d: gap 0.2389
  (experiment final 0.17), σ_zz,top 1.787 MPa (experiment final ~3.0).
- **EPFL T33 P2 homogeneous MCC**: e = 0.4695 at 20 MPa vs measured 0.57.
- **ebs_task13**: no model output to compare (failures above); measurement
  CSVs preserved in `ebs/measurements/`.

## Caveats

- All numbers above are stated side by side with their anchors; agreement
  judgements are Vinay's.
- MS33 III/IV/VII carry the Vinay-approved 20× spec permeability; Model I
  headers have an open TODO(Vinay) on E,ν provenance (see ms33 README).
- Model III completed in 86 s, far below the historical 1–2 h — reported as
  observed (full schedule reached, 438 steps, 0 rejected).
