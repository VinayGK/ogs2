# Superseded reference (moved 2026-10-05, K-fix set with E = 52 MPa; never deleted)

`ms33_modelIV_pellets_ts_681_t_17280000.000000.vtu` (md5 273d003deb1d2ed69dc8b7f161bd0310): the IV pellets reference of Mass-conserving A (candidate 2a, branch
gen5_masscons_1a_conformant_2026-10-02 at fc14f19fb9), 681 accepted steps.
Superseded on branch gen5_masscons_1a_kfix_guard_E52_conformant_2026-10-05 by `../ms33_modelIV_pellets_ts_709_t_17280000.000000.vtu` (md5 cb081e83d5e15eb8a03b0584eaf5c7ac), the t_end frame of the
record run `1a/suite/IV` of the K-fix set (MacBook Pro, `~/ogs-models/scratch/2026-10-05_kfix_guard_recut/1a/runs/fx_IV/out/`; the same md5 is printed in
CANONICAL_RESULTS_2026-10-05_gen5_masscons_1a_kfix_guard_v0.json, md5 ee36f00a77ffd89df0e1467fb4b1b047).
Binary of the record run: bin/ogs md5 ad1d8b16b52fba77874a0166391f109d, lib/libRichardsMechanics.dylib md5 7452094c2200330fb63724a238fc01f9 (built from 3739b7ec03 + guard_always_on.diff, md5 95e9c7c217325224690f6251b560e2a6;
its source tree equals the tree of the code commit of this branch, checked with git rev-parse). Run deck md5 22d48f4fef385b72cb923cc49acce3b5; the committed deck
`ms33_modelIV_pellets.prj` equals it with XML comments, test_definition and blank lines removed (E52 package ogs/apply_e52_decks.log.json).
Run: 709 accepted steps, 0 rejected, rc 0, 0 error lines; wall 177 s (OGS timer 173.155 s), OMP_NUM_THREADS 1, nice 10,
MacBook Pro (one run per core, other jobs on the machine).

WHY: the nine committed decks of this branch carry, in addition to the switch lines of Mass-conserving A, swelling_stress_k_level = true (the K-fix), and the
code carries the always-on corrected clamp guard of the live-K swelling tangent. This changes the swelling stress, the stepping and the end state.
vtkdiff pairs output and reference by IDENTICAL file name (ts_681 -> ts_709), and with the deck's own tolerances the old frame against the new one
passes 3 of 11 checks, so the old reference cannot pass by construction. The new reference registers this configuration's output for regression;
it is not a validation of it.

vtkdiff old -> new reference, MEASURED 2026-10-05 on the MacBook Pro (vtkdiff /Users/vinaykumar/git/build/dsm_kfix_guard_E52_ctest_20261005/bin/vtkdiff, md5 f8bd4c4a4f0eabe617bfa6a342b8b483;
tolerances of the deck's <test_definition>, not edited; a field passes if the abs or the rel max norm is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 2.990532030167138e-04, 2.391245980778637e-04 | 3.316469976971954e-01, 4.580317544059191e+00 | NO |
| saturation | 1e-8 / 1e-8 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| porosity | 1e-8 / 1e-7 | 2.155830850672247e-02 | 3.768123632574086e-02 | NO |
| transport_porosity | 1e-8 / 1e-7 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| micro_porosity | 1e-8 / 1e-7 | 2.155830850672247e-02 | 3.768123632574086e-02 | NO |
| micro_water_content | 1e-8 / 1e-7 | 2.155830850672247e-02 | 3.768123632574086e-02 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 5.993209764868834e+01 | 6.134949860836706e-02 | NO |
| sigma | 1e3 / 1e-2 | 1.257812309147775e+05, 9.721483728941502e+05, 1.239156467400257e+05, 3.500008176159287e+04 | 1.605643317896983e-01, 2.891855665722005e-01, 2.778637218038300e-01, 9.489024987385609e+01 | NO |
| swelling_stress | 1e3 / 1e-2 | 1.968771948231543e+06, 1.968771948231543e+06, 1.968771948231543e+06, 0.000000000000000e+00 | 1.180352461698824e+00, 1.180352461698824e+00, 1.180352461698824e+00, 0.000000000000000e+00 | NO |
| pressure | 1e3 / 1e-2 | 1.235658907042448e-11 | 3.471269349096464e+03 | yes |
| micro_pressure | 1e3 / 1e-2 | 3.064492358225308e+06 | 2.333413196199667e-01 | NO |

The new reference against the record-run frame it was copied from passes 11 of 11 (identity check).
ctest of this branch: run after this commit, recorded in the E52 package (ogs/CTEST.md); not stated here.
Record: ~/ogs-models/scratch/2026-10-05_kfix_guard_E52_package/ogs/ (vtkdiff_replay_full.md, e52_reregister.json).
