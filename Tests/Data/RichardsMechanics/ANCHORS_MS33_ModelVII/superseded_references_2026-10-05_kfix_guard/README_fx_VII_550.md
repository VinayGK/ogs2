# Superseded reference (moved 2026-10-05, K-fix set with E = 52 MPa; never deleted)

`ms33_modelVII_freeswelling_ladder550_ts_1147_t_47520000.000000.vtu` (md5 144266cdbbb7e388a54d2cd0b79d1fa0): the VII ladder550 reference of Mass-conserving A (candidate 2a, branch
gen5_masscons_1a_conformant_2026-10-02 at fc14f19fb9), 1147 accepted steps.
Superseded on branch gen5_masscons_1a_kfix_guard_E52_conformant_2026-10-05 by `../ms33_modelVII_freeswelling_ladder550_ts_1150_t_47520000.000000.vtu` (md5 75a18a075e21ed04c366bafa95861759), the t_end frame of the
record run `1a/suite/VII_550` of the K-fix set (MacBook Pro, `~/ogs-models/scratch/2026-10-05_kfix_guard_recut/1a/runs/fx_VII_550/out/`; the same md5 is printed in
CANONICAL_RESULTS_2026-10-05_gen5_masscons_1a_kfix_guard_v0.json, md5 ee36f00a77ffd89df0e1467fb4b1b047).
Binary of the record run: bin/ogs md5 ad1d8b16b52fba77874a0166391f109d, lib/libRichardsMechanics.dylib md5 7452094c2200330fb63724a238fc01f9 (built from 3739b7ec03 + guard_always_on.diff, md5 95e9c7c217325224690f6251b560e2a6;
its source tree equals the tree of the code commit of this branch, checked with git rev-parse). Run deck md5 2b8f96577a3cda7bdab7ad0e522d2996; the committed deck
`ms33_modelVII_freeswelling_ladder550.prj` equals it with XML comments, test_definition and blank lines removed (E52 package ogs/apply_e52_decks.log.json).
Run: 1150 accepted steps, 0 rejected, rc 0, 0 error lines; wall 263 s (OGS timer 258.865 s), OMP_NUM_THREADS 1, nice 10,
MacBook Pro (one run per core, other jobs on the machine).

WHY: the nine committed decks of this branch carry, in addition to the switch lines of Mass-conserving A, swelling_stress_k_level = true (the K-fix), and the
code carries the always-on corrected clamp guard of the live-K swelling tangent. This changes the swelling stress, the stepping and the end state.
vtkdiff pairs output and reference by IDENTICAL file name (ts_1147 -> ts_1150), and with the deck's own tolerances the old frame against the new one
passes 3 of 11 checks, so the old reference cannot pass by construction. The new reference registers this configuration's output for regression;
it is not a validation of it.

The two Model VII decks have separate reference files (different output prefixes); both are replaced, one README each.

vtkdiff old -> new reference, MEASURED 2026-10-05 on the MacBook Pro (vtkdiff /Users/vinaykumar/git/build/dsm_kfix_guard_E52_ctest_20261005/bin/vtkdiff, md5 f8bd4c4a4f0eabe617bfa6a342b8b483;
tolerances of the deck's <test_definition>, not edited; a field passes if the abs or the rel max norm is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 3.927115117093260e-04, 1.091898543139138e-03 | 2.903760849444079e-01, 3.759321410989050e-01 | NO |
| saturation | 1e-8 / 1e-8 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| porosity | 1e-8 / 1e-7 | 2.352639922358968e-02 | 4.637523139838591e-02 | NO |
| transport_porosity | 1e-8 / 1e-7 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| micro_porosity | 1e-8 / 1e-7 | 2.352639922358968e-02 | 4.637523139838591e-02 | NO |
| micro_water_content | 1e-8 / 1e-7 | 2.352639922358968e-02 | 4.637523139838591e-02 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 6.540338984157984e+01 | 5.014490026919191e-02 | NO |
| sigma | 1e3 / 1e-2 | 1.670942466110903e+05, 1.348109841549186e+05, 1.670942466110888e+05, 1.538831585162885e+04 | 7.494220621196356e+02, 4.119143122387995e-01, 4.194317114155532e+02, 4.362470210263098e+02 | NO |
| swelling_stress | 1e3 / 1e-2 | 2.279059532984957e+06, 2.279059532984957e+06, 2.279059532984957e+06, 0.000000000000000e+00 | 3.434228723199930e-01, 3.434228723199930e-01, 3.434228723199930e-01, 0.000000000000000e+00 | NO |
| pressure | 1e3 / 1e-2 | 2.543644654915662e-10 | 3.461782028113853e+01 | yes |
| micro_pressure | 1e3 / 1e-2 | 2.135772454734053e+06 | 1.030613226581749e-01 | NO |

The new reference against the record-run frame it was copied from passes 11 of 11 (identity check).
ctest of this branch: run after this commit, recorded in the E52 package (ogs/CTEST.md); not stated here.
Record: ~/ogs-models/scratch/2026-10-05_kfix_guard_E52_package/ogs/ (vtkdiff_replay_full.md, e52_reregister.json).
