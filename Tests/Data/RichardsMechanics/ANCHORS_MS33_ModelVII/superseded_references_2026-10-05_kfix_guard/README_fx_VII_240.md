# Superseded reference (moved 2026-10-05, K-fix set with E = 52 MPa; never deleted)

`ms33_modelVII_freeswelling_ts_887_t_20736000.000000.vtu` (md5 4fca7b7e65eb509a91bae0a28cf07044): the VII 240 d reference of Mass-conserving A (candidate 2a, branch
gen5_masscons_1a_conformant_2026-10-02 at fc14f19fb9), 887 accepted steps.
Superseded on branch gen5_masscons_1a_kfix_guard_E52_conformant_2026-10-05 by `../ms33_modelVII_freeswelling_ts_881_t_20736000.000000.vtu` (md5 7a87de78bc24f62be0f6c55d7dd7c21c), the t_end frame of the
record run `1a/suite/VII_240` of the K-fix set (MacBook Pro, `~/ogs-models/scratch/2026-10-05_kfix_guard_recut/1a/runs/fx_VII_240/out/`; the same md5 is printed in
CANONICAL_RESULTS_2026-10-05_gen5_masscons_1a_kfix_guard_v0.json, md5 ee36f00a77ffd89df0e1467fb4b1b047).
Binary of the record run: bin/ogs md5 ad1d8b16b52fba77874a0166391f109d, lib/libRichardsMechanics.dylib md5 7452094c2200330fb63724a238fc01f9 (built from 3739b7ec03 + guard_always_on.diff, md5 95e9c7c217325224690f6251b560e2a6;
its source tree equals the tree of the code commit of this branch, checked with git rev-parse). Run deck md5 47c38261c244930382b5df028c6e25ef; the committed deck
`ms33_modelVII_freeswelling.prj` equals it with XML comments, test_definition and blank lines removed (E52 package ogs/apply_e52_decks.log.json).
Run: 881 accepted steps, 0 rejected, rc 0, 0 error lines; wall 217 s (OGS timer 213.695 s), OMP_NUM_THREADS 1, nice 10,
MacBook Pro (one run per core, other jobs on the machine).

WHY: the nine committed decks of this branch carry, in addition to the switch lines of Mass-conserving A, swelling_stress_k_level = true (the K-fix), and the
code carries the always-on corrected clamp guard of the live-K swelling tangent. This changes the swelling stress, the stepping and the end state.
vtkdiff pairs output and reference by IDENTICAL file name (ts_887 -> ts_881), and with the deck's own tolerances the old frame against the new one
passes 3 of 11 checks, so the old reference cannot pass by construction. The new reference registers this configuration's output for regression;
it is not a validation of it.

The two Model VII decks have separate reference files (different output prefixes); both are replaced, one README each.

vtkdiff old -> new reference, MEASURED 2026-10-05 on the MacBook Pro (vtkdiff /Users/vinaykumar/git/build/dsm_kfix_guard_E52_ctest_20261005/bin/vtkdiff, md5 f8bd4c4a4f0eabe617bfa6a342b8b483;
tolerances of the deck's <test_definition>, not edited; a field passes if the abs or the rel max norm is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 3.927018027319368e-04, 1.091388424212438e-03 | 2.900819772083431e-01, 3.754897954362432e-01 | NO |
| saturation | 1e-8 / 1e-8 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| porosity | 1e-8 / 1e-7 | 2.342925073455082e-02 | 4.618553465015297e-02 | NO |
| transport_porosity | 1e-8 / 1e-7 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| micro_porosity | 1e-8 / 1e-7 | 2.342925073455082e-02 | 4.618553465015297e-02 | NO |
| micro_water_content | 1e-8 / 1e-7 | 2.342925073455082e-02 | 4.618553465015297e-02 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 6.513331704205143e+01 | 4.992539050071700e-02 | NO |
| sigma | 1e3 / 1e-2 | 1.672585976397620e+05, 1.348826844205700e+05, 1.672585976397609e+05, 1.545508429881505e+04 | 1.370196862682042e+03, 4.115715828655541e-01, 1.600075566070819e+03, 1.950992603143551e+03 | NO |
| swelling_stress | 1e3 / 1e-2 | 2.278847966407308e+06, 2.278847966407308e+06, 2.278847966407308e+06, 0.000000000000000e+00 | 3.430500549145321e-01, 3.430500549145321e-01, 3.430500549145321e-01, 0.000000000000000e+00 | NO |
| pressure | 1e3 / 1e-2 | 1.535026346779902e+02 | 3.263439495151246e-01 | yes |
| micro_pressure | 1e3 / 1e-2 | 2.125440654127393e+06 | 1.024924037011253e-01 | NO |

The new reference against the record-run frame it was copied from passes 11 of 11 (identity check).
ctest of this branch: run after this commit, recorded in the E52 package (ogs/CTEST.md); not stated here.
Record: ~/ogs-models/scratch/2026-10-05_kfix_guard_E52_package/ogs/ (vtkdiff_replay_full.md, e52_reregister.json).
