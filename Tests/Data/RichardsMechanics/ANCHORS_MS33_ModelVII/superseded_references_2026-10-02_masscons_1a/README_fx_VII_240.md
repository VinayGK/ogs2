# Superseded reference (moved 2026-10-02, mass-conserving candidate 2a; never deleted)

`ms33_modelVII_freeswelling_ts_1006_t_20736000.000000.vtu` (md5 a30214d8c8505abe059bf89f9e3fa2ad): the VII 240 d reference of the IC-fix V0 campaign (candidate 1, branch
gen5_icfix_conformant_2026-09-30 at d64c47cb30; record binary bin/ogs md5 23640465bff1ef6ce5c83d6c5dcd48d3, commit 686fcd6ef8),
1006 accepted steps (IC-FIX V0 paragraph of ProcessLib/RichardsMechanics/Tests.cmake).
Superseded on branch gen5_masscons_1a_conformant_2026-10-02 by `../ms33_modelVII_freeswelling_ts_887_t_20736000.000000.vtu` (md5 4fca7b7e65eb509a91bae0a28cf07044), the t_end frame of the
candidate-2a record run `1a/runs/fx_VII_240` (Shilpa's MBP `~/ogs-models/scratch/2026-10-02_candidate2_runs/1a/runs/fx_VII_240/out/`;
copy on WD_elements `run_outputs_2026-10-02/shilpas-macbook-pro/2026-10-02_candidate2_runs/`, same md5 in the Shilpa and the
WD_elements md5 manifests). Binary: bin/ogs md5 fd2e4e76e51d490e45c15b3861c7471d, lib/libRichardsMechanics.dylib md5 31b02f0b7e28b25b9590c767d3f6b8aa, built from
e03c32a1a3 (gen5_masscons_1a_conformant_2026-10-02; the 1b branch has the same code). Run deck md5 ab779c4f5ba0ffc1e2545582841c8ebd;
the committed deck `ms33_modelVII_freeswelling.prj` has the same content apart from comments (c2/runs_prep/1a/CHECK_vs_c2branch.md).
Run: 887 accepted steps, 0 rejected, 6239 Newton iterations, rc 0, 0 error lines; wall 395 s
(OGS timer 390.192 s), OMP_NUM_THREADS 1, nice 10, Shilpa's MBP (M1 Max, 10 cores, one run per core).

WHY: the nine committed decks of this branch carry, inside <potential_exchange>, the KKT micro ceiling, F3, the kkt_active
tangents, the latched bishop_relperm gate, macro storage on phi_M, the T_m drop, darcy_relative_permeability_mobility = kirchhoff_element_mean (variant 1a) and
macro_storage_exact_time_levels = true. This changes the micro and macro water balances and so the adaptive stepping and
the end state. vtkdiff pairs output and reference by IDENTICAL file name (ts_1006 -> ts_887), and with the deck's own
tolerances the old frame against the new one passes 3 of 11 checks, so the old reference cannot pass by construction.
The new reference registers this configuration's output for regression; it is not a validation of it.

The two Model VII decks have separate reference files (different output prefixes); both are replaced, one README each.

vtkdiff old -> new reference, MEASURED 2026-10-02 on the mac mini (vtkdiff /Users/vinaykumar/git/build/kkt_v5_testing_20261002/bin/vtkdiff, md5 3633bcdd8fef6f2e86c5c27683396ea4;
tolerances of the deck's <test_definition>, not edited; a field passes if the abs or the rel max norm is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 9.947846858137602e-07, 2.497394219308849e-06 | 7.582419095622376e-04, 7.434695034185038e-04 | NO |
| saturation | 1e-8 / 1e-8 | 1.000046317702282e+00 | 2.025013106407723e+06 | NO |
| porosity | 1e-8 / 1e-7 | 1.668823107107142e-04 | 3.146115374658318e-04 | NO |
| transport_porosity | 1e-8 / 1e-7 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| micro_porosity | 1e-8 / 1e-7 | 1.668823107107142e-04 | 3.146115374658318e-04 | NO |
| micro_water_content | 1e-8 / 1e-7 | 1.668823107107142e-04 | 3.146115374658318e-04 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 4.639328237760765e-01 | 3.555272560848546e-04 | NO |
| sigma | 1e3 / 1e-2 | 4.118680492234431e+03, 2.957809165379615e+03, 4.118680492233907e+03, 7.760993118066472e+02 | 6.590411849136709e+02, 7.138716473217626e-03, 1.255795912008488e+02, 8.920238008993857e+01 | NO |
| swelling_stress | 1e3 / 1e-2 | 7.526304211862385e+03, 7.526304211862385e+03, 7.526304211862385e+03, 0.000000000000000e+00 | 8.559060580678680e-04, 8.559060580678680e-04, 8.559060580678680e-04, 0.000000000000000e+00 | yes |
| pressure | 1e3 / 1e-2 | 1.885307959786528e+07 | 1.275460483962684e+05 | NO |
| micro_pressure | 1e3 / 1e-2 | 1.877371650562063e+04 | 9.050622843776347e-04 | yes |

The new reference against the record-run frame it was copied from passes 11 of 11 (identity check).
ctest of this branch: run after this commit (P7c), recorded outside the repository; not stated here.
Record: ~/ogs-models/scratch/2026-10-02_candidate2/ctest/ (1a/vtkdiff_replay_full.md, 1a/p7a_reregister_1a.json).
