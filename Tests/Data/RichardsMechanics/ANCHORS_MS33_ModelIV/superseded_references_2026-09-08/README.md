# Superseded reference (moved 2026-09-08, never deleted — CLAUDE.md 6.2/6.3)

`ms33_modelIV_pellets_ts_21500_t_17280000.000000.vtu` — the gen-4 Model IV reference (2026-08-31, commit 7ec39ecf4c) of the
STACKED layout (pellets z 0–35 mm under the block z 35–70 mm), which is NOT the spec's geometry (spec Fig4pellets.jpg: concentric).
Superseded by `../ms33_modelIV_pellets_ts_16008_t_17280000.000000.vtu` from the gen5_conformant campaign of 2026-09-08 (Vinay:
"approve all five, run them as one campaign"): CONCENTRIC mesh `ms33_pellets_concentric_r25_h70.vtu` (clay core r <= 15 mm, pellet
annulus 15–25 mm), gen-4 900 knot 18469.144929 retained (the gen-5 knot stalls on this mesh at 27.28 d — see the deck). Old
values (clay-local 70/40.25/10.5, stacked probes): 4.393549 / 4.361817 / 0.892442 MPa; new (spec frame 64.75/35/5.25, all in the
clay core): 2.397459 / 2.391964 / 2.373176 MPa, pellet centre (20,35) 1.290552 MPa.
Run card: ~/git/eurad-anchors/runs/2026-09-08_*_gen5_conformant_2026-09-08_*/README.md.
