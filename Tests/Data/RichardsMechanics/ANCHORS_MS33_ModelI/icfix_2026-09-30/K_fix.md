# K_fix table (patched binary dsm-ic-fix-20260930), 2026-09-30

Produced by refit_icfix.py (bracket sweep of 13 K around the shipped knot, then log-linear secant rounds, one run per round) and kfix_table.py, stdlib only, on Shilpa's Mac. MEASURED = runs kept in calib/runs/. NOT ADOPTED: the values enter only the deck copies in decks_kfix/.

Ps = mean total stress -tr(sigma)/3 at t_end = 200 d [MPa], mean over the 4 nodes (equals node 0 and the swelling-stress mean to <1e-13 relative, see K_fix.json). Tolerance |Ps/target - 1| <= 1e-09 (the shipped gen-5 900 fit reached +2.2e-10; 1400/1600/1800 shipped at -6e-6, -7.4e-5, +7.7e-6).

| knot | target Ps [MPa] | K_ship [J/kg] | K_fix [J/kg] | K_fix/K_ship - 1 | Ps at K_fix [MPa] | dev. from target | patched deviation from target at K_ship [= fix-only shift + shipped-fit residual; CORRECTED 10:1x, part 1 mislabelled it 'the shift the fix causes'; split below] | iteration rounds (+13-point sweep) | failed runs |
|---|---|---|---|---|---|---|---|---|---|
| 900, floor 0.0 (Reference, III, VII) | 0.3500522009 | 24241.892204 | 24241.938995 | +1.9302e-06 | 0.350052201214 | +8.98e-10 | -1.92e-06 | 1 | 3 |
| 900, floor 0.08 (Model IV) | 0.3500522009 | 18469.144929 | 18469.1763367 | +1.7005e-06 | 0.350052200947 | +1.35e-10 | -1.69e-06 | 1 | 1 |
| 1400 | 4.9 | 46000 | 46000.3269694 | +7.1080e-06 | 4.9 | +9.99e-13 | -7.09e-06 | 2 | 0 |
| 1600 | 14.161 | 104689.9129 | 104698.192423 | +7.9086e-05 | 14.161 | -2.54e-12 | -7.83e-05 | 2 | 0 |
| 1800 | 40.6 | 265905.06 | 265909.813902 | +1.7878e-05 | 40.6 | +7.94e-13 | -1.78e-05 | 2 | 0 |

Controls (MEASURED, shipped binary, shipped K): dd900_f0 Ps 0.350052201 (dev +2.17e-10); dd900_f008 Ps 0.350052201 (dev +2.60e-10); dd1400 Ps 4.89997034 (dev -6.05e-06); dd1600 Ps 14.15994803 (dev -7.43e-05); dd1800 Ps 40.60031318 (dev +7.71e-06).


CORRECTION of part 1: the column 'Ps at K_ship, patched [dev.]' of K_fix.md and the README section-6 sentence 'the fix moves Ps at the shipped K by ...' gave the patched binary's DEVIATION FROM THE TARGET, not the shift the fix causes. The fix-only shift is the ratio of the two controls; the rest of the deviation is the miss of the SHIPPED fit, which the shipped procedure accepted (calibrate_maxwell_K.py default rel_tol 0.02).

| knot | K_fix/K_ship - 1 | fix-only shift of Ps at K_ship | shipped-fit residual (Ps_ship/target - 1) | patched deviation at K_ship | share of K change from the fix | share from the shipped miss |
|---|---|---|---|---|---|---|
| 900, floor 0.0 (Reference, III, VII, IX) | +1.930e-06 | -1.924e-06 | +2.172e-10 | -1.924e-06 | 100 % | -0 % |
| 900, floor 0.08 (Model IV) | +1.701e-06 | -1.694e-06 | +2.603e-10 | -1.694e-06 | 100 % | -0 % |
| 1400 | +7.108e-06 | -1.033e-06 | -6.053e-06 | -7.086e-06 | 15 % | 85 % |
| 1600 | +7.909e-05 | -3.998e-06 | -7.429e-05 | -7.828e-05 | 5 % | 95 % |
| 1800 | +1.788e-05 | -2.555e-05 | +7.714e-06 | -1.783e-05 | 143 % | -43 % |

Reading: at 900 (both cells) the K change is the fix (shipped fit was at 2e-10). At 1400 and 1600 about 85 % and 95 % of the K change closes the shipped fits' own misses, NOT the fix; at 1800 the two components have opposite sign (first-order shares outside 0..100 %). Every fx-vs-v2 difference near dd1600 (Reference, III, VII, IX) therefore mixes the fix with a recalibration whose K change at 1600 (+7.9e-5) is about 20 times the fix-only Ps shift (-4.0e-6). Shares are first order in log space (PREDICTED rule, arithmetic MEASURED).

Target convention: 900 knot target = the PRINTED 0.3500522009 (dd900 deck header; campaign_lib.TARGET_PRIMARY), as in part 1. The exact Dixon relation 0.003*exp(5.2883*0.9) = 0.3500522008569968 (= the record's model_I.cells.900.dixon_target_MPa). Deviations on both conventions:

| 900 cell | shipped dev vs printed | shipped dev vs exact | K_fix dev vs printed | K_fix dev vs exact |
|---|---|---|---|---|
| dd900_f0 | +2.17e-10 | +3.40e-10 | +8.98e-10 | +1.02e-09 |
| dd900_f008 | +2.60e-10 | +3.83e-10 | +1.35e-10 | +2.58e-10 |

Against the exact relation the floor-0.0 K_fix sits just outside the declared 1e-9 tolerance (declared against the printed value); the effect on K is about 1e-10 relative. NOT re-fitted to the exact value (Vinay's call which convention; 'ratified read-offs as printed' applies to 1400/1600/1800).

Control-deck wording (should_fix 6): the dd1400 control deck (md5 46fc278b, both ctl_ship_Kship and ctl_fix_Kship) is NOT byte-identical to the shipped deck (md5 a06ad902...): refit_icfix.py wrote the K string as '46000' instead of '46000.0' (ledger change line '46000.0 -> 46000'). Numerically the same K; the other four control decks are byte-identical to the shipped deck or differ only in the recorded live strings.

