# Superseded reference (moved 2026-09-30, IC-fix V0 re-cut; never deleted - CLAUDE.md 6.2/6.3)

`ms33_modelIV_pellets_ts_16527_t_17280000.000000.vtu` (md5 1bc0d9e9689920fc198cef2651746af0) - the IV pellets reference of the gen5_conformant v2 campaign of 2026-09-08 (accepted steps 16527 accepted / 0 rejected), produced by the shipped binary
(`archive/dsm_native_Pi_fofnlev_branchtip_2026-08-11-49-gbed3e395`, md5 cea9d7d81972f732385b41a71e50f20e) with the gen-5 K table.
Superseded by `../ms33_modelIV_pellets_ts_16494_t_17280000.000000.vtu` (md5 5da4f4e6a982bd192e78864180ce337b), the last frame of run `fx_IV` of the IC-fix readiness suite (Shilpa's MBP
`~/ogs-models/scratch/2026-09-30_readiness_icfix/runs/fx_IV/out/`, full copy on WD_elements `run_outputs_2026-09-30/shilpas-macbook-pro/...`, same md5): patched binary
(commit 686fcd6ef8, bin/ogs md5 23640465bff1ef6ce5c83d6c5dcd48d3) + the K_fix table (deck `ms33_modelIV_pellets.prj`, md5 of the run deck 2c8c493a53d2efb99e7c07458e616d16; the committed deck has the same live content,
only header comments differ - `check_decks_vs_kfix.py` in the 2026-09-30_icfix_v0_recut record). Run: 16494 accepted steps, 0 rejected, rc 0, 0 warning and 0 error lines; OMP_NUM_THREADS 6 (runs single-threaded), 16112.7 s on a loaded 10-core host (the record run took 7957.4 s on the MBP).

WHY: the IC fix (declared initial_micro_water_content realised in the previous micro porosity) changes the initial micro state and so the adaptive stepping and the end state;
the K table is refitted with the fixed binary (K_fix, ANCHORS_MS33_ModelI/icfix_2026-09-30/). vtkdiff matches output and reference by IDENTICAL file name (step index in the name: ts_16527 -> ts_16494),
and with the deck's own tolerances the old frame against the new one passes 3 of 11 checks, so the old reference cannot pass by construction. The refit calibrates Ps to the Dixon target, so
the agreement of the new frame with the target is calibration, not validation.

vtkdiff old -> new reference, MEASURED 2026-09-30 on the mac mini (vtkdiff from `~/git/build/ogs-dsm-active`, the only vtkdiff on the mini; repo tolerances of the deck's <test_definition>;
abs max norm and rel max norm per component / total where vtkdiff prints several; a field passes if either is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 1.979015005852940e-07, 8.742092987151784e-07 | 2.581637888250593e-04, 3.948839770456044e-03 | NO |
| saturation | 1e-8 / 1e-8 | 9.995805336920904e-07 | 1.429812499916363e-04 | NO |
| porosity | 1e-8 / 1e-7 | 5.930135851073270e-05 | 9.328078527622230e-05 | NO |
| transport_porosity | 1e-8 / 1e-7 | 5.983769670130923e-04 | 1.331861479520452e-01 | NO |
| micro_porosity | 1e-8 / 1e-7 | 5.390756085024706e-04 | 1.243982172944591e-03 | NO |
| micro_water_content | 1e-8 / 1e-7 | 2.138026408858851e-04 | 3.983791632202182e-04 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 1.648577766595736e-01 | 1.628214037706978e-04 | NO |
| sigma | 1e3 / 1e-2 | 8.462599714702228e+02, 1.185271838889457e+03, 8.541534632134717e+02, 1.674586669133569e+02 | 9.792575994329986e-04, 9.369843374424780e-04, 1.317573158581503e-03, 6.416666666666667e+00 | NO |
| swelling_stress | 1e3 / 1e-2 | 2.660784567527473e+03, 2.660784567527473e+03, 2.660784567527473e+03, 0.000000000000000e+00 | 3.215463465712948e-03, 3.215463465712948e-03, 3.215463465712948e-03, 0.000000000000000e+00 | yes |
| pressure | 1e3 / 1e-2 | 6.168508005160838e+02 | 7.144101165391435e-05 | yes |
| micro_pressure | 1e3 / 1e-2 | 2.245005470547825e+03 | 2.468269744620052e-04 | yes |

NOT run through ctest: the mini has no ctest and the patched binary was built with OGS_BUILD_TESTING=OFF. The new reference against the fx output it was copied from passes 11 of 11 (identity check, vtkdiff replay).
To run on the MBP: build this branch (gen5_icfix_conformant_2026-09-30) with OGS_BUILD_TESTING=ON and run `ctest -R "ANCHORS_MS33|ms33_" -j<cores>` in the build directory (Model IV: 8957.6 s on the MBP in the 2026-09-30 T4 run, registered RUNTIME 10736).
Record: ~/ogs-models/scratch/2026-09-30_icfix_v0_recut/OGS_GATE.md and ogs_step/vtkdiff_replay.md.
