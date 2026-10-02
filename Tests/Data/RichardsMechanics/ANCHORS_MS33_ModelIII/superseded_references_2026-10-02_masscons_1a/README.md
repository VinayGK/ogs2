# Superseded reference (moved 2026-10-02, mass-conserving candidate 2a; never deleted)

`ms33_modelIII_gapswitch_ts_970_t_17280000.000000.vtu` (md5 020aa0821051c3ef993f93138f976b0a): the III gapswitch reference of the IC-fix V0 campaign (candidate 1, branch
gen5_icfix_conformant_2026-09-30 at d64c47cb30; record binary bin/ogs md5 23640465bff1ef6ce5c83d6c5dcd48d3, commit 686fcd6ef8),
970 accepted steps (IC-FIX V0 paragraph of ProcessLib/RichardsMechanics/Tests.cmake).
Superseded on branch gen5_masscons_1a_conformant_2026-10-02 by `../ms33_modelIII_gapswitch_ts_822_t_17280000.000000.vtu` (md5 b472410c06a742b87586af3cfdba08a8), the t_end frame of the
candidate-2a record run `1a/runs/fx_III` (Shilpa's MBP `~/ogs-models/scratch/2026-10-02_candidate2_runs/1a/runs/fx_III/out/`;
copy on WD_elements `run_outputs_2026-10-02/shilpas-macbook-pro/2026-10-02_candidate2_runs/`, same md5 in the Shilpa and the
WD_elements md5 manifests). Binary: bin/ogs md5 fd2e4e76e51d490e45c15b3861c7471d, lib/libRichardsMechanics.dylib md5 31b02f0b7e28b25b9590c767d3f6b8aa, built from
e03c32a1a3 (gen5_masscons_1a_conformant_2026-10-02; the 1b branch has the same code). Run deck md5 e3f4f504348cbd941cbe2fca2ce08a1c;
the committed deck `ms33_modelIII_gapswitch.prj` has the same content apart from comments (c2/runs_prep/1a/CHECK_vs_c2branch.md).
Run: 822 accepted steps, 0 rejected, 6223 Newton iterations, rc 0, 0 error lines; wall 312 s
(OGS timer 311.602 s), OMP_NUM_THREADS 1, nice 10, Shilpa's MBP (M1 Max, 10 cores, one run per core).

WHY: the nine committed decks of this branch carry, inside <potential_exchange>, the KKT micro ceiling, F3, the kkt_active
tangents, the latched bishop_relperm gate, macro storage on phi_M, the T_m drop, darcy_relative_permeability_mobility = kirchhoff_element_mean (variant 1a) and
macro_storage_exact_time_levels = true. This changes the micro and macro water balances and so the adaptive stepping and
the end state. vtkdiff pairs output and reference by IDENTICAL file name (ts_970 -> ts_822), and with the deck's own
tolerances the old frame against the new one passes 3 of 11 checks, so the old reference cannot pass by construction.
The new reference registers this configuration's output for regression; it is not a validation of it.

vtkdiff old -> new reference, MEASURED 2026-10-02 on the mac mini (vtkdiff /Users/vinaykumar/git/build/kkt_v5_testing_20261002/bin/vtkdiff, md5 3633bcdd8fef6f2e86c5c27683396ea4;
tolerances of the deck's <test_definition>, not edited; a field passes if the abs or the rel max norm is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 1.484979456109047e-07, 6.882031398173013e-07 | 1.678729210468309e-04, 6.159507095784952e-01 | NO |
| saturation | 1e-8 / 1e-8 | 1.000041231424917e+00 | 2.268862536870462e+06 | NO |
| porosity | 1e-8 / 1e-7 | 3.688480696961705e-05 | 7.152056440488755e-05 | NO |
| transport_porosity | 1e-8 / 1e-7 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| micro_porosity | 1e-8 / 1e-7 | 3.688480696961705e-05 | 7.152056440488755e-05 | NO |
| micro_water_content | 1e-8 / 1e-7 | 3.688480696961705e-05 | 7.152056440488755e-05 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 1.025397633752618e-01 | 7.617050549107234e-05 | NO |
| sigma | 1e3 / 1e-2 | 4.175403599371610e+03, 4.678625637960155e+03, 5.724793248880072e+03, 4.253952299252041e+02 | 1.244093949104661e-02, 1.222626370843531e-03, 1.661800570584201e-02, 2.457950324971304e+03 | NO |
| swelling_stress | 1e3 / 1e-2 | 7.825069438699633e+03, 7.825069438699633e+03, 7.825069438699633e+03, 0.000000000000000e+00 | 8.798761417641812e-04, 8.798761417641812e-04, 8.798761417641812e-04, 0.000000000000000e+00 | yes |
| pressure | 1e3 / 1e-2 | 1.995661551343848e+07 | 4.264248842374473e+20 | NO |
| micro_pressure | 1e3 / 1e-2 | 3.582193163461983e+03 | 1.624568564723417e-04 | yes |

The new reference against the record-run frame it was copied from passes 11 of 11 (identity check).
ctest of this branch: run after this commit (P7c), recorded outside the repository; not stated here.
Record: ~/ogs-models/scratch/2026-10-02_candidate2/ctest/ (1a/vtkdiff_replay_full.md, 1a/p7a_reregister_1a.json).
