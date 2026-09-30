# Superseded reference (moved 2026-09-30, IC-fix V0 re-cut; never deleted - CLAUDE.md 6.2/6.3)

`ms33_modelIII_gapswitch_ts_1116_t_17280000.000000.vtu` (md5 1579f009e48a2b09ab11fd3c5e65f223) - the III gapswitch reference of the gen5_conformant v2 campaign of 2026-09-08 (accepted steps 1116 accepted / 5 rejected), produced by the shipped binary
(`archive/dsm_native_Pi_fofnlev_branchtip_2026-08-11-49-gbed3e395`, md5 cea9d7d81972f732385b41a71e50f20e) with the gen-5 K table.
Superseded by `../ms33_modelIII_gapswitch_ts_970_t_17280000.000000.vtu` (md5 020aa0821051c3ef993f93138f976b0a), the last frame of run `fx_III` of the IC-fix readiness suite (Shilpa's MBP
`~/ogs-models/scratch/2026-09-30_readiness_icfix/runs/fx_III/out/`, full copy on WD_elements `run_outputs_2026-09-30/shilpas-macbook-pro/...`, same md5): patched binary
(commit 686fcd6ef8, bin/ogs md5 23640465bff1ef6ce5c83d6c5dcd48d3) + the K_fix table (deck `ms33_modelIII_gapswitch.prj`, md5 of the run deck f35668d7116e62ecb4108b0636d7c0bc; the committed deck has the same live content,
only header comments differ - `check_decks_vs_kfix.py` in the 2026-09-30_icfix_v0_recut record). Run: 970 accepted steps, 0 rejected, rc 0, 0 warning and 0 error lines; OMP_NUM_THREADS 2 (runs single-threaded), 349.0 s.

WHY: the IC fix (declared initial_micro_water_content realised in the previous micro porosity) changes the initial micro state and so the adaptive stepping and the end state;
the K table is refitted with the fixed binary (K_fix, ANCHORS_MS33_ModelI/icfix_2026-09-30/). vtkdiff matches output and reference by IDENTICAL file name (step index in the name: ts_1116 -> ts_970),
and with the deck's own tolerances the old frame against the new one passes 4 of 11 checks, so the old reference cannot pass by construction. The refit calibrates Ps to the Dixon target, so
the agreement of the new frame with the target is calibration, not validation.

vtkdiff old -> new reference, MEASURED 2026-09-30 on the mac mini (vtkdiff from `~/git/build/ogs-dsm-active`, the only vtkdiff on the mini; repo tolerances of the deck's <test_definition>;
abs max norm and rel max norm per component / total where vtkdiff prints several; a field passes if either is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 1.097466763386873e-08, 9.593513325422473e-08 | 1.025172245965310e-05, 2.279205538925949e-01 | NO |
| saturation | 1e-8 / 1e-8 | 1.499176144158362e-08 | 4.642335246517624e-05 | NO |
| porosity | 1e-8 / 1e-7 | 6.330285866606999e-06 | 1.227385564542784e-05 | NO |
| transport_porosity | 1e-8 / 1e-7 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| micro_porosity | 1e-8 / 1e-7 | 6.330285866606999e-06 | 1.227385564542784e-05 | NO |
| micro_water_content | 1e-8 / 1e-7 | 6.330285866606999e-06 | 1.227385564542784e-05 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 1.759819470908042e-02 | 1.307262024601436e-05 | NO |
| sigma | 1e3 / 1e-2 | 9.408412666284712e+02, 5.830658113937825e+02, 1.054529591946164e+03, 3.874418822527514e+01 | 3.257802210035192e-03, 1.547117486557146e-04, 3.902825241221676e-03, 1.602037869868559e+01 | NO |
| swelling_stress | 1e3 / 1e-2 | 1.531345108553767e+03, 1.531345108553767e+03, 1.531345108553767e+03, 0.000000000000000e+00 | 1.739909114735753e-04, 1.739909114735753e-04, 1.739909114735753e-04, 0.000000000000000e+00 | yes |
| pressure | 1e3 / 1e-2 | 3.932153807133436e+02 | 2.321125597163290e-05 | yes |
| micro_pressure | 1e3 / 1e-2 | 6.828318453505635e+02 | 3.102545047902073e-05 | yes |

NOT run through ctest: the mini has no ctest and the patched binary was built with OGS_BUILD_TESTING=OFF. The new reference against the fx output it was copied from passes 11 of 11 (identity check, vtkdiff replay).
To run on the MBP: build this branch (gen5_icfix_conformant_2026-09-30) with OGS_BUILD_TESTING=ON and run `ctest -R "ANCHORS_MS33|ms33_" -j<cores>` in the build directory (Model IV: 8957.6 s on the MBP in the 2026-09-30 T4 run, registered RUNTIME 10736).
Record: ~/ogs-models/scratch/2026-09-30_icfix_v0_recut/OGS_GATE.md and ogs_step/vtkdiff_replay.md.
