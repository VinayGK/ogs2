# Superseded reference (moved 2026-09-08, campaign v2; never deleted — CLAUDE.md 6.2/6.3)

`ms33_modelIII_gapswitch_ts_943_t_17280000.000000.vtu` (md5 cbfd07f1813e9954935b327866382d17) — the Model III gap-switch (clay r = 23 mm, wall at r = 25 mm, gen-5 900 knot) reference of the
gen5_conformant campaign of 2026-09-08 (morning; superseded_references_2026-09-08/ holds its predecessor).
Superseded by `../ms33_modelIII_gapswitch_ts_1116_t_17280000.000000.vtu` (md5 1579f009e48a2b09ab11fd3c5e65f223) from campaign v2 of
2026-09-08 (`/Users/vinaykumar/ogs-models/scratch/2026-09-08_0936_gen5_conformant_v2_inflight/runs/III/out`, same md5), Vinay: "3x audit against specs ... results ground up".

WHY: the ONLY deck change is `<fixed_output_times>` set to the data-collection template's 5-day grid (audit
2026-09-08_spec_audit_3x LANE2 finding 3); physics, BCs, meshes and numerics are unchanged. Every listed output time
forces a step to land there, so the adaptive IterationNumberBased stepping takes a different path and the final
timestep index changes (ts_943 -> ts_1116); vtkdiff matches output and reference by IDENTICAL file name, so the
old reference cannot pass by construction. The two frames are the same physical time; their difference is the
time-stepping footprint only — MEASURED 2026-09-08 with vtkdiff (max abs / max componentwise rel, all 11 test fields)
in `/Users/vinaykumar/ogs-models/scratch/2026-09-08_0936_gen5_conformant_v2_inflight/OGS_TREE_LANE_2026-09-08.md`.
