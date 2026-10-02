# Superseded reference (moved 2026-09-30, IC-fix V0 re-cut; never deleted - CLAUDE.md 6.2/6.3)

`ms33_reference_dd1600_ts_855_t_17280000.000000.vtu` (md5 db714d4685daa1cb207af8a6807ac452) - the Reference reference of the gen5_conformant v2 campaign of 2026-09-08 (accepted steps 855), produced by the shipped binary
(`archive/dsm_native_Pi_fofnlev_branchtip_2026-08-11-49-gbed3e395`, md5 cea9d7d81972f732385b41a71e50f20e) with the gen-5 K table.
Superseded by `../ms33_reference_dd1600_ts_895_t_17280000.000000.vtu` (md5 6eb7fbda4358ac543e60f1d97edd1de1), the last frame of run `fx_Ref` of the IC-fix readiness suite (Shilpa's MBP
`~/ogs-models/scratch/2026-09-30_readiness_icfix/runs/fx_Ref/out/`, full copy on WD_elements `run_outputs_2026-09-30/shilpas-macbook-pro/...`, same md5): patched binary
(commit 686fcd6ef8, bin/ogs md5 23640465bff1ef6ce5c83d6c5dcd48d3) + the K_fix table (deck `ms33_reference_dd1600.prj`, md5 of the run deck 2354a64384c9f2f54f6634bcba09d19e; the committed deck has the same live content,
only header comments differ - `check_decks_vs_kfix.py` in the 2026-09-30_icfix_v0_recut record). Run: 895 accepted steps, 0 rejected, rc 0, 0 warning and 0 error lines; 1 thread, 20.3 s.

WHY: the IC fix (declared initial_micro_water_content realised in the previous micro porosity) changes the initial micro state and so the adaptive stepping and the end state;
the K table is refitted with the fixed binary (K_fix, ANCHORS_MS33_ModelI/icfix_2026-09-30/). vtkdiff matches output and reference by IDENTICAL file name (step index in the name: ts_855 -> ts_895),
and with the deck's own tolerances the old frame against the new one passes 4 of 11 checks, so the old reference cannot pass by construction. The refit calibrates Ps to the Dixon target, so
the agreement of the new frame with the target is calibration, not validation.

vtkdiff old -> new reference, MEASURED 2026-09-30 on the mac mini (vtkdiff from `~/git/build/ogs-dsm-active`, the only vtkdiff on the mini; repo tolerances of the deck's <test_definition>;
abs max norm and rel max norm per component / total where vtkdiff prints several; a field passes if either is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 8.031732832671547e-18, 7.173781189733744e-07 | 1.662916494343272e+02, 2.361027918781894e-03 | NO |
| saturation | 1e-8 / 1e-8 | 3.650148388020357e-09 | 4.292721726944075e-04 | yes |
| porosity | 1e-8 / 1e-7 | 3.821721960700630e-05 | 9.353840749427769e-05 | NO |
| transport_porosity | 1e-8 / 1e-7 | 1.309924592179448e-04 | 1.508981020311574e-02 | NO |
| micro_porosity | 1e-8 / 1e-7 | 1.194050113450529e-04 | 3.124124848976300e-04 | NO |
| micro_water_content | 1e-8 / 1e-7 | 7.019525346885214e-05 | 1.745956536789848e-04 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 1.062438705075692e-01 | 6.462277853223950e-05 | NO |
| sigma | 1e3 / 1e-2 | 3.753830077521503e+03, 1.626944940704852e+03, 3.753830077551305e+03, 3.866359493206586e-08 | 2.870920879189990e-04, 1.159246682856118e-04, 2.870920879212786e-04, 1.086499872463520e+02 | NO |
| swelling_stress | 1e3 / 1e-2 | 5.392668008822948e+03, 5.392668008822948e+03, 5.392668008822948e+03, 0.000000000000000e+00 | 4.484458850980795e-04, 4.484458850980795e-04, 4.484458850980795e-04, 0.000000000000000e+00 | yes |
| pressure | 1e3 / 1e-2 | 5.243116749107838e+03 | 2.120172782547364e-04 | yes |
| micro_pressure | 1e3 / 1e-2 | 3.413985016743839e+04 | 4.672445716914770e-04 | yes |

NOT run through ctest: the mini has no ctest and the patched binary was built with OGS_BUILD_TESTING=OFF. The new reference against the fx output it was copied from passes 11 of 11 (identity check, vtkdiff replay).
To run on the MBP: build this branch (gen5_icfix_conformant_2026-09-30) with OGS_BUILD_TESTING=ON and run `ctest -R "ANCHORS_MS33|ms33_" -j<cores>` in the build directory (Model IV: 8957.6 s on the MBP in the 2026-09-30 T4 run, registered RUNTIME 10736).
Record: ~/ogs-models/scratch/2026-09-30_icfix_v0_recut/OGS_GATE.md and ogs_step/vtkdiff_replay.md.
