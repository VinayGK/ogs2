# Superseded reference (moved 2026-10-05, K-fix set with E = 52 MPa; never deleted)

`ms33_modelIII_gapswitch_ts_822_t_17280000.000000.vtu` (md5 b472410c06a742b87586af3cfdba08a8): the III gapswitch reference of Mass-conserving A (candidate 2a, branch
gen5_masscons_1a_conformant_2026-10-02 at fc14f19fb9), 822 accepted steps.
Superseded on branch gen5_masscons_1a_kfix_guard_E52_conformant_2026-10-05 by `../ms33_modelIII_gapswitch_ts_778_t_17280000.000000.vtu` (md5 70766c85cbc1ec3eb9bd71e2fd3dc22d), the t_end frame of the
record run `1a/suite/III` of the K-fix set (MacBook Pro, `~/ogs-models/scratch/2026-10-05_kfix_guard_recut/1a/runs/fx_III/out/`; the same md5 is printed in
CANONICAL_RESULTS_2026-10-05_gen5_masscons_1a_kfix_guard_v0.json, md5 ee36f00a77ffd89df0e1467fb4b1b047).
Binary of the record run: bin/ogs md5 ad1d8b16b52fba77874a0166391f109d, lib/libRichardsMechanics.dylib md5 7452094c2200330fb63724a238fc01f9 (built from 3739b7ec03 + guard_always_on.diff, md5 95e9c7c217325224690f6251b560e2a6;
its source tree equals the tree of the code commit of this branch, checked with git rev-parse). Run deck md5 8436ff20f8b36d2891e80233bbebbc93; the committed deck
`ms33_modelIII_gapswitch.prj` equals it with XML comments, test_definition and blank lines removed (E52 package ogs/apply_e52_decks.log.json).
Run: 778 accepted steps, 0 rejected, rc 0, 0 error lines; wall 177 s (OGS timer 172.594 s), OMP_NUM_THREADS 1, nice 10,
MacBook Pro (one run per core, other jobs on the machine).

WHY: the nine committed decks of this branch carry, in addition to the switch lines of Mass-conserving A, swelling_stress_k_level = true (the K-fix), and the
code carries the always-on corrected clamp guard of the live-K swelling tangent. This changes the swelling stress, the stepping and the end state.
vtkdiff pairs output and reference by IDENTICAL file name (ts_822 -> ts_778), and with the deck's own tolerances the old frame against the new one
passes 3 of 11 checks, so the old reference cannot pass by construction. The new reference registers this configuration's output for regression;
it is not a validation of it.

vtkdiff old -> new reference, MEASURED 2026-10-05 on the MacBook Pro (vtkdiff /Users/vinaykumar/git/build/dsm_kfix_guard_E52_ctest_20261005/bin/vtkdiff, md5 f8bd4c4a4f0eabe617bfa6a342b8b483;
tolerances of the deck's <test_definition>, not edited; a field passes if the abs or the rel max norm is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 4.159507889630817e-04, 1.891908717819359e-05 | 2.625870371102931e-01, 4.067907023369162e+03 | NO |
| saturation | 1e-8 / 1e-8 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| porosity | 1e-8 / 1e-7 | 1.828121898504947e-02 | 3.668594122075059e-02 | NO |
| transport_porosity | 1e-8 / 1e-7 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| micro_porosity | 1e-8 / 1e-7 | 1.828121898504947e-02 | 3.668594122075059e-02 | NO |
| micro_water_content | 1e-8 / 1e-7 | 1.828121898504947e-02 | 3.668594122075059e-02 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 5.082178877843694e+01 | 3.781783649407600e-02 | NO |
| sigma | 1e3 / 1e-2 | 4.777266188046940e+05, 1.185591674251766e+06, 4.777266188046933e+05, 1.094236135786433e+04 | 2.319498155692001e+06, 4.392582539045914e-01, 2.767906880903388e+06, 9.788546182511743e+11 | NO |
| swelling_stress | 1e3 / 1e-2 | 2.271237159617290e+06, 2.271237159617290e+06, 2.271237159617290e+06, 0.000000000000000e+00 | 3.397355879195015e-01, 3.397355879195015e-01, 3.397355879195015e-01, 0.000000000000000e+00 | NO |
| pressure | 1e3 / 1e-2 | 2.117197042756885e-11 | 3.030806552154300e+02 | yes |
| micro_pressure | 1e3 / 1e-2 | 1.796945340204798e+06 | 8.178219566981167e-02 | NO |

The new reference against the record-run frame it was copied from passes 11 of 11 (identity check).
ctest of this branch: run after this commit, recorded in the E52 package (ogs/CTEST.md); not stated here.
Record: ~/ogs-models/scratch/2026-10-05_kfix_guard_E52_package/ogs/ (vtkdiff_replay_full.md, e52_reregister.json).
