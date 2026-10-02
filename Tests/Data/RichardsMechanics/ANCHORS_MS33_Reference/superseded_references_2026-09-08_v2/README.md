# Superseded reference (moved 2026-09-08, campaign v2; never deleted — CLAUDE.md 6.2/6.3)

`ms33_reference_dd1600_ts_771_t_17280000.000000.vtu` (md5 414d63ce70f9117a0311dcf1658142c0) — the Reference-configuration (confined dd1600, no gap) reference of the
gen5_conformant campaign of 2026-09-08 (morning; first tracked reference of this deck, tracked at commit 028c431fa0 — no earlier superseded folder exists here).
Superseded by `../ms33_reference_dd1600_ts_855_t_17280000.000000.vtu` (md5 db714d4685daa1cb207af8a6807ac452) from campaign v2 of
2026-09-08 (`/Users/vinaykumar/ogs-models/scratch/2026-09-08_0936_gen5_conformant_v2_inflight/runs/Ref/out`, same md5), Vinay: "3x audit against specs ... results ground up".

WHY: the ONLY deck change is `<fixed_output_times>` set to the data-collection template's 5-day grid (audit
2026-09-08_spec_audit_3x LANE2 finding 3); physics, BCs, meshes and numerics are unchanged. Every listed output time
forces a step to land there, so the adaptive IterationNumberBased stepping takes a different path and the final
timestep index changes (ts_771 -> ts_855); vtkdiff matches output and reference by IDENTICAL file name, so the
old reference cannot pass by construction. The two frames are the same physical time; their difference is the
time-stepping footprint only — MEASURED 2026-09-08 with vtkdiff (max abs / max componentwise rel, all 11 test fields)
in `/Users/vinaykumar/ogs-models/scratch/2026-09-08_0936_gen5_conformant_v2_inflight/OGS_TREE_LANE_2026-09-08.md`.
