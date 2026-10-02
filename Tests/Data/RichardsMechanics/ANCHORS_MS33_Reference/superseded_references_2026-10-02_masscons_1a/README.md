# Superseded reference (moved 2026-10-02, mass-conserving candidate 2a; never deleted)

`ms33_reference_dd1600_ts_895_t_17280000.000000.vtu` (md5 6eb7fbda4358ac543e60f1d97edd1de1): the Reference reference of the IC-fix V0 campaign (candidate 1, branch
gen5_icfix_conformant_2026-09-30 at d64c47cb30; record binary bin/ogs md5 23640465bff1ef6ce5c83d6c5dcd48d3, commit 686fcd6ef8),
895 accepted steps (IC-FIX V0 paragraph of ProcessLib/RichardsMechanics/Tests.cmake).
Superseded on branch gen5_masscons_1a_conformant_2026-10-02 by `../ms33_reference_dd1600_ts_698_t_17280000.000000.vtu` (md5 fe7db2c52ab48d9bdae39a8cf4475bc9), the t_end frame of the
candidate-2a record run `1a/runs/fx_Ref` (Shilpa's MBP `~/ogs-models/scratch/2026-10-02_candidate2_runs/1a/runs/fx_Ref/out/`;
copy on WD_elements `run_outputs_2026-10-02/shilpas-macbook-pro/2026-10-02_candidate2_runs/`, same md5 in the Shilpa and the
WD_elements md5 manifests). Binary: bin/ogs md5 fd2e4e76e51d490e45c15b3861c7471d, lib/libRichardsMechanics.dylib md5 31b02f0b7e28b25b9590c767d3f6b8aa, built from
e03c32a1a3 (gen5_masscons_1a_conformant_2026-10-02; the 1b branch has the same code). Run deck md5 858cf89bec86028acda69b976a446e0d;
the committed deck `ms33_reference_dd1600.prj` has the same content apart from comments (c2/runs_prep/1a/CHECK_vs_c2branch.md).
Run: 698 accepted steps, 0 rejected, 4918 Newton iterations, rc 0, 0 error lines; wall 18 s
(OGS timer 16.9061 s), OMP_NUM_THREADS 1, nice 10, Shilpa's MBP (M1 Max, 10 cores, one run per core).

WHY: the nine committed decks of this branch carry, inside <potential_exchange>, the KKT micro ceiling, F3, the kkt_active
tangents, the latched bishop_relperm gate, macro storage on phi_M, the T_m drop, darcy_relative_permeability_mobility = kirchhoff_element_mean (variant 1a) and
macro_storage_exact_time_levels = true. This changes the micro and macro water balances and so the adaptive stepping and
the end state. vtkdiff pairs output and reference by IDENTICAL file name (ts_895 -> ts_698), and with the deck's own
tolerances the old frame against the new one passes 1 of 11 checks, so the old reference cannot pass by construction.
The new reference registers this configuration's output for regression; it is not a validation of it.

vtkdiff old -> new reference, MEASURED 2026-10-02 on the mac mini (vtkdiff /Users/vinaykumar/git/build/kkt_v5_testing_20261002/bin/vtkdiff, md5 3633bcdd8fef6f2e86c5c27683396ea4;
tolerances of the deck's <test_definition>, not edited; a field passes if the abs or the rel max norm is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 6.790393610491855e-18, 3.003211263457282e-05 | 9.775972660994873e+01, 1.318561073837982e-01 | NO |
| saturation | 1e-8 / 1e-8 | 1.000000729291232e+00 | 1.487885772765658e+07 | NO |
| porosity | 1e-8 / 1e-7 | 9.255741594501044e-04 | 2.144800553279663e-03 | NO |
| transport_porosity | 1e-8 / 1e-7 | 4.922075381138270e-02 | inf | NO |
| micro_porosity | 1e-8 / 1e-7 | 5.014632797082819e-02 | 1.311623852478780e-01 | NO |
| micro_water_content | 1e-8 / 1e-7 | 3.035393076084109e-02 | 7.548571996288783e-02 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 2.573096163270975e+00 | 1.630878032750014e-03 | NO |
| sigma | 1e3 / 1e-2 | 9.313763721899688e+04, 3.178000777287409e+04, 9.313763721900061e+04, 2.677396908235118e-08 | 6.407004191635675e-03, 2.263592578579953e-03, 6.407004191635933e-03, 1.242740812098612e+02 | NO |
| swelling_stress | 1e3 / 1e-2 | 1.421270699710250e+05, 1.421270699710250e+05, 1.421270699710250e+05, 0.000000000000000e+00 | 9.630259850915236e-03, 9.630259850915236e-03, 9.630259850915236e-03, 0.000000000000000e+00 | yes |
| pressure | 1e3 / 1e-2 | 5.110548958151349e+07 | 1.494569336228749e+20 | NO |
| micro_pressure | 1e3 / 1e-2 | 3.221635265995249e+06 | 6.079192533978318e-02 | NO |

The new reference against the record-run frame it was copied from passes 11 of 11 (identity check).
ctest of this branch: run after this commit (P7c), recorded outside the repository; not stated here.
Record: ~/ogs-models/scratch/2026-10-02_candidate2/ctest/ (1a/vtkdiff_replay_full.md, 1a/p7a_reregister_1a.json).
