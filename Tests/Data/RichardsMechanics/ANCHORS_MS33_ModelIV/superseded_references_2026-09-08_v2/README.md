# Superseded reference (moved 2026-09-08, campaign v2; never deleted — CLAUDE.md 6.2/6.3)

`ms33_modelIV_pellets_ts_16008_t_17280000.000000.vtu` (md5 27379f65629edabf3f1e8badf86b1458) — the Model IV concentric-pellets reference of the
gen5_conformant campaign of 2026-09-08 (morning; superseded_references_2026-09-08/ holds its predecessor, the stacked-layout
reference).
Superseded by `../ms33_modelIV_pellets_ts_16527_t_17280000.000000.vtu` (md5 1bc0d9e9689920fc198cef2651746af0) from campaign v2 of
2026-09-08 (`/Users/vinaykumar/ogs-models/scratch/2026-09-08_0936_gen5_conformant_v2_inflight/runs/IV/out`, same md5), Vinay: "3x audit against specs ... results ground up".

WHY: the ONLY deck change is `<fixed_output_times>` set to the data-collection template's 5-day grid (audit
2026-09-08_spec_audit_3x LANE2 finding 3); physics, BCs, mesh and numerics are unchanged. Every listed output time
forces a step to land there, so the adaptive IterationNumberBased stepping takes a different path and the final
timestep index changes (ts_16008 -> ts_16527); vtkdiff matches output and reference by IDENTICAL file name, so the
old reference cannot pass by construction. Run completed cleanly on the laptop after the poll in
OGS_TREE_LANE_2026-09-08.md's task 3 (16527 accepted steps, 0 rejected, 0 error lines, 7957.4 s at OMP 6). The two
frames are the same physical time; their difference is the time-stepping footprint only — MEASURED 2026-09-08 with
vtkdiff (max abs / max componentwise rel, all 11 test fields) in `/Users/vinaykumar/ogs-models/scratch/2026-09-08_0936_gen5_conformant_v2_inflight/OGS_TREE_LANE_2026-09-08.md`.

vtkdiff old (ts_16008) -> new (ts_16527), MEASURED 2026-09-08 on the mac mini (build maxwell_conjugate_20260824), abs max / rel max:
displacement 4.8e-8, 3.0e-7 / 7.3e-5, 5.8e-4; saturation 2.7e-7 / 2.2e-5; porosity 3.1e-5 / 5.1e-5; transport_porosity 1.7e-4 / 5.5e-2
(at its floor); micro_porosity 1.5e-4 / 3.6e-4; micro_water_content 6.0e-5 / 1.1e-4; dry_density_solid 8.6e-2 / 7.8e-5; sigma 1.8e2 Pa
/ 1.4e1 (near-zero shear); swelling_stress 2.4e2 Pa / 9.5e-4; pressure 11.6 Pa / 1.1e-5; micro_pressure 1.7e3 Pa / 2.7e-4. Spec-frame
probes 2.397472 / 2.391977 / 2.373192 MPa (04:37 run: 2.397459 / 2.391964 / 2.373176), pellet centre 1.290556 (1.290552).
