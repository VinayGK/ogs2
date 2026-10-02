# IC-fix V0 calibration record (2026-09-30)

Copies (byte-identical) of the records of the readiness suite, `~/ogs-models/scratch/2026-09-30_readiness_icfix/` (Shilpa's MBP, mirror on the mini):
- `K_fix.json` (md5 e6f29bd492a5011c8f8a6f7688d6980e), `K_fix.md`, `K_fix_attribution.md`: the K table refitted with the patched binary (ogs commit 686fcd6ef8,
  bin/ogs md5 23640465bff1ef6ce5c83d6c5dcd48d3) by the shipped procedure. STATUS WORDING INSIDE THE COPIES ("NOT ADOPTED; values enter only the deck copies") is the
  state of the readiness lane at the time of writing; since 2026-09-30 ~22:10 the values are ADOPTED into the decks of branch gen5_icfix_conformant_2026-09-30
  (Vinay: K_fix accepted). The copies are not edited.
- `_calib_result_dd900_floor008_ModelIV_knot_icfix_2026-09-30.json`: the 900 knot of Model IV (floor-0.08 cell, gen-4 point). The four Model I cells have their
  records next to the decks: `../_calib_result_dd{900,1400,1600,1800}_icfix_2026-09-30.json` (same format as the shipped `_calib_result_dd*.json`, which stay as history).

Labels: MEASURED = K and Ps are read from the runs listed in K_fix.json (`calib/runs/<cell>__it0N`). Ps = target is calibration, not validation.
