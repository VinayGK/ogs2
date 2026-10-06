# Superseded reference (moved 2026-10-05, K-fix set with E = 52 MPa; never deleted)

`ms33_reference_dd1600_ts_698_t_17280000.000000.vtu` (md5 fe7db2c52ab48d9bdae39a8cf4475bc9): the Reference reference of Mass-conserving A (candidate 2a, branch
gen5_masscons_1a_conformant_2026-10-02 at fc14f19fb9), 698 accepted steps.
Superseded on branch gen5_masscons_1a_kfix_guard_E52_conformant_2026-10-05 by `../ms33_reference_dd1600_ts_591_t_17280000.000000.vtu` (md5 48f6f2034cd308e21bbda841908b3321), the t_end frame of the
record run `1a/suite/Ref` of the K-fix set (MacBook Pro, `~/ogs-models/scratch/2026-10-05_kfix_guard_recut/1a/runs/fx_Ref/out/`; the same md5 is printed in
CANONICAL_RESULTS_2026-10-05_gen5_masscons_1a_kfix_guard_v0.json, md5 ee36f00a77ffd89df0e1467fb4b1b047).
Binary of the record run: bin/ogs md5 ad1d8b16b52fba77874a0166391f109d, lib/libRichardsMechanics.dylib md5 7452094c2200330fb63724a238fc01f9 (built from 3739b7ec03 + guard_always_on.diff, md5 95e9c7c217325224690f6251b560e2a6;
its source tree equals the tree of the code commit of this branch, checked with git rev-parse). Run deck md5 985e9aa747f11fd43b1f0a7d3676f0e7; the committed deck
`ms33_reference_dd1600.prj` equals it with XML comments, test_definition and blank lines removed (E52 package ogs/apply_e52_decks.log.json).
Run: 591 accepted steps, 0 rejected, rc 0, 0 error lines; wall 10 s (OGS timer 7.24289 s), OMP_NUM_THREADS 1, nice 10,
MacBook Pro (one run per core, other jobs on the machine).

WHY: the nine committed decks of this branch carry, in addition to the switch lines of Mass-conserving A, swelling_stress_k_level = true (the K-fix), and the
code carries the always-on corrected clamp guard of the live-K swelling tangent. This changes the swelling stress, the stepping and the end state.
vtkdiff pairs output and reference by IDENTICAL file name (ts_698 -> ts_591), and with the deck's own tolerances the old frame against the new one
passes 3 of 11 checks, so the old reference cannot pass by construction. The new reference registers this configuration's output for regression;
it is not a validation of it.

vtkdiff old -> new reference, MEASURED 2026-10-05 on the MacBook Pro (vtkdiff /Users/vinaykumar/git/build/dsm_kfix_guard_E52_ctest_20261005/bin/vtkdiff, md5 f8bd4c4a4f0eabe617bfa6a342b8b483;
tolerances of the deck's <test_definition>, not edited; a field passes if the abs or the rel max norm is within tolerance):

| field | abs tol / rel tol | abs max | rel max | pass |
|---|---|---|---|---|
| displacement | 1e-9 / 1e-7 | 6.214704294161032e-18, 4.890712967372296e-04 | 4.877324753451560e+01, 5.516416410245739e+00 | NO |
| saturation | 1e-8 / 1e-8 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| porosity | 1e-8 / 1e-7 | 1.999284888131503e-02 | 4.896339630005455e-02 | NO |
| transport_porosity | 1e-8 / 1e-7 | 0.000000000000000e+00 | 0.000000000000000e+00 | yes |
| micro_porosity | 1e-8 / 1e-7 | 1.999284888131503e-02 | 4.896339630005455e-02 | NO |
| micro_water_content | 1e-8 / 1e-7 | 1.999284888131503e-02 | 4.896339630005455e-02 | NO |
| dry_density_solid | 1e-5 / 1e-8 | 5.558011989005604e+01 | 3.497180391709183e-02 | NO |
| sigma | 1e3 / 1e-2 | 1.689341262909701e+06, 3.222489190839939e+05, 1.689341262909705e+06, 4.161691726090547e-08 | 1.298754099211937e-01, 2.290443204761464e-02, 1.298754099211944e-01, 1.920606954055351e+02 | NO |
| swelling_stress | 1e3 / 1e-2 | 2.720177690451525e+06, 2.720177690451525e+06, 2.720177690451525e+06, 0.000000000000000e+00 | 2.261474172591371e-01, 2.261474172591371e-01, 2.261474172591371e-01, 0.000000000000000e+00 | NO |
| pressure | 1e3 / 1e-2 | 2.280850253512513e-10 | 1.445418807451887e+02 | yes |
| micro_pressure | 1e3 / 1e-2 | 1.754748937121565e+07 | 3.143966444521093e-01 | NO |

The new reference against the record-run frame it was copied from passes 11 of 11 (identity check).
ctest of this branch: run after this commit, recorded in the E52 package (ogs/CTEST.md); not stated here.
Record: ~/ogs-models/scratch/2026-10-05_kfix_guard_E52_package/ogs/ (vtkdiff_replay_full.md, e52_reregister.json).
