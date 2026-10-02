# Superseded reference (moved 2026-09-30, IC-fix V0 re-cut; never deleted - CLAUDE.md 6.2/6.3)

`ms33_modelVII_freeswelling_ladder550_ts_1304_t_47520000.000000.vtu` (md5 255f08988a8b9fbbaf3d08be8c85c65d) - the VII ladder550 reference of the gen5_conformant v2 campaign of 2026-09-08 (accepted steps 1304 accepted / 0 rejected), produced by the shipped binary
(`archive/dsm_native_Pi_fofnlev_branchtip_2026-08-11-49-gbed3e395`, md5 cea9d7d81972f732385b41a71e50f20e) with the gen-5 K table.
Superseded by `../ms33_modelVII_freeswelling_ladder550_ts_1287_t_47520000.000000.vtu` (md5 d57025c5811e0df56b61f11977216f50), the last frame of run `fx_VII_550` of the IC-fix readiness suite (Shilpa's MBP
`~/ogs-models/scratch/2026-09-30_readiness_icfix/runs/fx_VII_550/out/`, full copy on WD_elements `run_outputs_2026-09-30/shilpas-macbook-pro/...`, same md5): patched binary
(commit 686fcd6ef8, bin/ogs md5 23640465bff1ef6ce5c83d6c5dcd48d3) + the K_fix table (deck `ms33_modelVII_freeswelling_ladder550.prj`, md5 of the run deck e138933daef81e5a610668c4838f6a6c; the committed deck has the same live content,
only header comments differ - `check_decks_vs_kfix.py` in the 2026-09-30_icfix_v0_recut record). Run: 1287 accepted steps, 0 rejected, rc 0, 0 warning and 0 error lines; OMP_NUM_THREADS 2, 466.5 s.

WHY: the IC fix (declared initial_micro_water_content realised in the previous micro porosity) changes the initial micro state and so the adaptive stepping and the end state;
the K table is refitted with the fixed binary (K_fix, ANCHORS_MS33_ModelI/icfix_2026-09-30/). vtkdiff matches output and reference by IDENTICAL file name (step index in the name: ts_1304 -> ts_1287),
and with the deck's own tolerances the old frame against the new one passes 4 of 11 checks, so the old reference cannot pass by construction. The refit calibrates Ps to the Dixon target, so
the agreement of the new frame with the target is calibration, not validation.

The two Model VII decks have separate reference files (different output prefixes); both are replaced.

vtkdiff old -> new reference, MEASURED 2026-09-30 on the mac mini (vtkdiff from `~/git/build/ogs-dsm-active`, the only vtkdiff on the mini; repo tolerances of the deck's <test_definition>;
abs max norm and rel max norm per component / total where vtkdiff prints several; a field passes if either is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 1.612462277424467e-07, 2.403717364921586e-07 | 9.115981262523045e-05, 1.676316803606925e-04 | NO |
| saturation | 1e-8 / 1e-8 | 2.038027578123158e-08 | 5.072021610582651e-05 | NO |
| porosity | 1e-8 / 1e-7 | 1.220959401537769e-05 | 2.302015403705564e-05 | NO |
| transport_porosity | 1e-8 / 1e-7 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| micro_porosity | 1e-8 / 1e-7 | 1.220959401537769e-05 | 2.302015403705564e-05 | NO |
| micro_water_content | 1e-8 / 1e-7 | 1.220959401537769e-05 | 2.302015403705564e-05 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 3.394267136263807e-02 | 2.599995209508103e-05 | NO |
| sigma | 1e3 / 1e-2 | 1.427128318661373e+03, 9.711825972880470e+02, 1.427128318662639e+03, 1.846468973983447e+02 | 7.511863777226863e+00, 2.265910014599000e-03, 1.846041722062379e+01, 5.441692804973800e+00 | NO |
| swelling_stress | 1e3 / 1e-2 | 1.385641934562474e+03, 1.385641934562474e+03, 1.385641934562474e+03, 0.000000000000000e+00 | 1.554989253346070e-04, 1.554989253346070e-04, 1.554989253346070e-04, 0.000000000000000e+00 | yes |
| pressure | 1e3 / 1e-2 | 3.468346760850400e+02 | 2.534496356832507e-05 | yes |
| micro_pressure | 1e3 / 1e-2 | 8.951079040877521e+02 | 4.311838156121477e-05 | yes |

NOT run through ctest: the mini has no ctest and the patched binary was built with OGS_BUILD_TESTING=OFF. The new reference against the fx output it was copied from passes 11 of 11 (identity check, vtkdiff replay).
To run on the MBP: build this branch (gen5_icfix_conformant_2026-09-30) with OGS_BUILD_TESTING=ON and run `ctest -R "ANCHORS_MS33|ms33_" -j<cores>` in the build directory (Model IV: 8957.6 s on the MBP in the 2026-09-30 T4 run, registered RUNTIME 10736).
Record: ~/ogs-models/scratch/2026-09-30_icfix_v0_recut/OGS_GATE.md and ogs_step/vtkdiff_replay.md.
