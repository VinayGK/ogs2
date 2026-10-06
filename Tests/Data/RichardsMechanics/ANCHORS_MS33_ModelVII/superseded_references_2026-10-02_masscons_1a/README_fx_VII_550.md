# Superseded reference (moved 2026-10-02, mass-conserving candidate 2a; never deleted)

`ms33_modelVII_freeswelling_ladder550_ts_1287_t_47520000.000000.vtu` (md5 d57025c5811e0df56b61f11977216f50): the VII ladder550 reference of the IC-fix V0 campaign (candidate 1, branch
gen5_icfix_conformant_2026-09-30 at d64c47cb30; record binary bin/ogs md5 23640465bff1ef6ce5c83d6c5dcd48d3, commit 686fcd6ef8),
1287 accepted steps (IC-FIX V0 paragraph of ProcessLib/RichardsMechanics/Tests.cmake).
Superseded on branch gen5_masscons_1a_conformant_2026-10-02 by `../ms33_modelVII_freeswelling_ladder550_ts_1147_t_47520000.000000.vtu` (md5 144266cdbbb7e388a54d2cd0b79d1fa0), the t_end frame of the
candidate-2a record run `1a/runs/fx_VII_550` (Shilpa's MBP `~/ogs-models/scratch/2026-10-02_candidate2_runs/1a/runs/fx_VII_550/out/`;
copy on WD_elements `run_outputs_2026-10-02/shilpas-macbook-pro/2026-10-02_candidate2_runs/`, same md5 in the Shilpa and the
WD_elements md5 manifests). Binary: bin/ogs md5 fd2e4e76e51d490e45c15b3861c7471d, lib/libRichardsMechanics.dylib md5 31b02f0b7e28b25b9590c767d3f6b8aa, built from
e03c32a1a3 (gen5_masscons_1a_conformant_2026-10-02; the 1b branch has the same code). Run deck md5 6879ed405448707b827f2df329b39789;
the committed deck `ms33_modelVII_freeswelling_ladder550.prj` has the same content apart from comments (c2/runs_prep/1a/CHECK_vs_c2branch.md).
Run: 1147 accepted steps, 0 rejected, 6626 Newton iterations, rc 0, 0 error lines; wall 486 s
(OGS timer 482.737 s), OMP_NUM_THREADS 1, nice 10, Shilpa's MBP (M1 Max, 10 cores, one run per core).

WHY: the nine committed decks of this branch carry, inside <potential_exchange>, the KKT micro ceiling, F3, the kkt_active
tangents, the latched bishop_relperm gate, macro storage on phi_M, the T_m drop, darcy_relative_permeability_mobility = kirchhoff_element_mean (variant 1a) and
macro_storage_exact_time_levels = true. This changes the micro and macro water balances and so the adaptive stepping and
the end state. vtkdiff pairs output and reference by IDENTICAL file name (ts_1287 -> ts_1147), and with the deck's own
tolerances the old frame against the new one passes 3 of 11 checks, so the old reference cannot pass by construction.
The new reference registers this configuration's output for regression; it is not a validation of it.

The two Model VII decks have separate reference files (different output prefixes); both are replaced, one README each.

vtkdiff old -> new reference, MEASURED 2026-10-02 on the mac mini (vtkdiff /Users/vinaykumar/git/build/kkt_v5_testing_20261002/bin/vtkdiff, md5 3633bcdd8fef6f2e86c5c27683396ea4;
tolerances of the deck's <test_definition>, not edited; a field passes if the abs or the rel max norm is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 6.073157882882218e-07, 1.110149441872418e-06 | 3.468732837656468e-04, 6.537902698834110e-04 | NO |
| saturation | 1e-8 / 1e-8 | 1.000046397501458e+00 | 2.021527358037666e+06 | NO |
| porosity | 1e-8 / 1e-7 | 5.175155464187675e-05 | 9.753734509731914e-05 | NO |
| transport_porosity | 1e-8 / 1e-7 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| micro_porosity | 1e-8 / 1e-7 | 5.175155464187675e-05 | 9.753734509731914e-05 | NO |
| micro_water_content | 1e-8 / 1e-7 | 5.175155464187675e-05 | 9.753734509731914e-05 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 1.438693219045035e-01 | 1.102583519104161e-04 | NO |
| sigma | 1e3 / 1e-2 | 3.298351115424142e+03, 2.296969068806095e+03, 3.298351115423109e+03, 4.758026988379961e+02 | 9.350422777577249e+01, 5.554314632255485e-03, 6.660300604514288e+01, 2.216917796719618e+01 | NO |
| swelling_stress | 1e3 / 1e-2 | 7.038234491407871e+03, 7.038234491407871e+03, 7.038234491407871e+03, 0.000000000000000e+00 | 7.970194962799459e-04, 7.970194962799459e-04, 7.970194962799459e-04, 0.000000000000000e+00 | yes |
| pressure | 1e3 / 1e-2 | 1.883747181707468e+07 | 7.652785659880566e+18 | NO |
| micro_pressure | 1e3 / 1e-2 | 4.156741509597749e+03 | 2.004231337576580e-04 | yes |

The new reference against the record-run frame it was copied from passes 11 of 11 (identity check).
ctest of this branch: run after this commit (P7c), recorded outside the repository; not stated here.
Record: ~/ogs-models/scratch/2026-10-02_candidate2/ctest/ (1a/vtkdiff_replay_full.md, 1a/p7a_reregister_1a.json).
