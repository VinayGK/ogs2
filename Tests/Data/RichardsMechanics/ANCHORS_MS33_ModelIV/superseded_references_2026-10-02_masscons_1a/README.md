# Superseded reference (moved 2026-10-02, mass-conserving candidate 2a; never deleted)

`ms33_modelIV_pellets_ts_16494_t_17280000.000000.vtu` (md5 5da4f4e6a982bd192e78864180ce337b): the IV pellets reference of the IC-fix V0 campaign (candidate 1, branch
gen5_icfix_conformant_2026-09-30 at d64c47cb30; record binary bin/ogs md5 23640465bff1ef6ce5c83d6c5dcd48d3, commit 686fcd6ef8),
16494 accepted steps (IC-FIX V0 paragraph of ProcessLib/RichardsMechanics/Tests.cmake).
Superseded on branch gen5_masscons_1a_conformant_2026-10-02 by `../ms33_modelIV_pellets_ts_681_t_17280000.000000.vtu` (md5 273d003deb1d2ed69dc8b7f161bd0310), the t_end frame of the
candidate-2a record run `1a/runs/fx_IV` (Shilpa's MBP `~/ogs-models/scratch/2026-10-02_candidate2_runs/1a/runs/fx_IV/out/`;
copy on WD_elements `run_outputs_2026-10-02/shilpas-macbook-pro/2026-10-02_candidate2_runs/`, same md5 in the Shilpa and the
WD_elements md5 manifests). Binary: bin/ogs md5 fd2e4e76e51d490e45c15b3861c7471d, lib/libRichardsMechanics.dylib md5 31b02f0b7e28b25b9590c767d3f6b8aa, built from
e03c32a1a3 (gen5_masscons_1a_conformant_2026-10-02; the 1b branch has the same code). Run deck md5 5fc5129c44749ec2c15cba03debbf987;
the committed deck `ms33_modelIV_pellets.prj` has the same content apart from comments (c2/runs_prep/1a/CHECK_vs_c2branch.md).
Run: 681 accepted steps, 0 rejected, 4905 Newton iterations, rc 0, 0 error lines; wall 296 s
(OGS timer 294.499 s), OMP_NUM_THREADS 1, nice 10, Shilpa's MBP (M1 Max, 10 cores, one run per core).

WHY: the nine committed decks of this branch carry, inside <potential_exchange>, the KKT micro ceiling, F3, the kkt_active
tangents, the latched bishop_relperm gate, macro storage on phi_M, the T_m drop, darcy_relative_permeability_mobility = kirchhoff_element_mean (variant 1a) and
macro_storage_exact_time_levels = true. This changes the micro and macro water balances and so the adaptive stepping and
the end state. vtkdiff pairs output and reference by IDENTICAL file name (ts_16494 -> ts_681), and with the deck's own
tolerances the old frame against the new one passes 0 of 11 checks, so the old reference cannot pass by construction.
The new reference registers this configuration's output for regression; it is not a validation of it.

vtkdiff old -> new reference, MEASURED 2026-10-02 on the mac mini (vtkdiff /Users/vinaykumar/git/build/kkt_v5_testing_20261002/bin/vtkdiff, md5 3633bcdd8fef6f2e86c5c27683396ea4;
tolerances of the deck's <test_definition>, not edited; a field passes if the abs or the rel max norm is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 1.092264376733240e-04, 1.458052637870319e-03 | 2.310916175900777e-01, 4.451343835804923e+00 | NO |
| saturation | 1e-8 / 1e-8 | 1.001067283678774e+00 | 1.336685367931179e+06 | NO |
| porosity | 1e-8 / 1e-7 | 3.641414352332528e-02 | 6.020944413216326e-02 | NO |
| transport_porosity | 1e-8 / 1e-7 | 5.642505237972967e-01 | inf | NO |
| micro_porosity | 1e-8 / 1e-7 | 5.618766163113734e-01 | 6.623095627801960e+00 | NO |
| micro_water_content | 1e-8 / 1e-7 | 4.520304367719741e-01 | 2.321889715622367e+00 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 1.012313189948441e+02 | 1.016583373529480e-01 | NO |
| sigma | 1e3 / 1e-2 | 6.540375066200455e+05, 7.760729418593086e+05, 6.704432297874215e+05, 2.704121245795017e+05 | 8.769766035342486e-01, 1.074282937664130e+00, 1.359498229648874e+00, 2.826649913679531e+03 | NO |
| swelling_stress | 1e3 / 1e-2 | 7.124569536447899e+05, 7.124569536447899e+05, 7.124569536447899e+05, 0.000000000000000e+00 | 2.183434769126015e+00, 2.183434769126015e+00, 2.183434769126015e+00, 0.000000000000000e+00 | NO |
| pressure | 1e3 / 1e-2 | 1.531783102288760e+07 | 9.298584983855871e+21 | NO |
| micro_pressure | 1e3 / 1e-2 | 8.895539046325503e+06 | 2.194692156897514e+00 | NO |

The new reference against the record-run frame it was copied from passes 11 of 11 (identity check).
ctest of this branch: run after this commit (P7c), recorded outside the repository; not stated here.
Record: ~/ogs-models/scratch/2026-10-02_candidate2/ctest/ (1a/vtkdiff_replay_full.md, 1a/p7a_reregister_1a.json).
