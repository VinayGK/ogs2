# Superseded reference (moved 2026-09-08, campaign v2; never deleted — CLAUDE.md 6.2/6.3)

`ms33_modelVII_freeswelling_ts_973_t_20736000.000000.vtu` (md5 22c865d5d77d767296de28a0946437ce) — the Model VII free-swelling (240-d step-and-hold Q&A ladder, gen-5 900 knot) reference of the
gen5_conformant campaign of 2026-09-08 (morning; superseded_references_2026-09-08/ holds its predecessor).
Superseded by `../ms33_modelVII_freeswelling_ts_1023_t_20736000.000000.vtu` (md5 5509fa9e33df94a6b2456dfa63325f7e) from campaign v2 of
2026-09-08 (`/Users/vinaykumar/ogs-models/scratch/2026-09-08_0936_gen5_conformant_v2_inflight/runs/VII_240/out`, same md5), Vinay: "3x audit against specs ... results ground up".

WHY: the ONLY deck change is `<fixed_output_times>` set to the data-collection template's 5-day grid (audit
2026-09-08_spec_audit_3x LANE2 finding 3); physics, BCs, meshes and numerics are unchanged. Every listed output time
forces a step to land there, so the adaptive IterationNumberBased stepping takes a different path and the final
timestep index changes (ts_973 -> ts_1023); vtkdiff matches output and reference by IDENTICAL file name, so the
old reference cannot pass by construction. The two frames are the same physical time; their difference is the
time-stepping footprint only — MEASURED 2026-09-08 with vtkdiff (max abs / max componentwise rel, all 11 test fields)
in `/Users/vinaykumar/ogs-models/scratch/2026-09-08_0936_gen5_conformant_v2_inflight/OGS_TREE_LANE_2026-09-08.md`.
Note: the campaign-v2 grid for this 240-d deck keeps the 7 half-interval frames (207.5 ... 237.5 d). The ladder-550
variant added the same day (`../ms33_modelVII_freeswelling_ladder550.prj`, own reference ts_1304 at 550 d) is a
SUPPLEMENT and does not supersede this deck or this reference.
