#!/usr/bin/env python3
"""
BEACON WP5.3 MGR — cross-team comparison data (compiled 2026-06-01).

PROVENANCE DISCIPLINE (CLAUDE.md §1.1):
  * EXACT      = tabulated in the deliverable text (cited table) or verified from our own VTU.
  * DIGITIZED  = read by eye off a published figure (cited figure + PDF page). APPROXIMATE,
                 ~+-0.02 g/cm3 on density, ~+-0.3 MPa on stress. NOT exact. Never to be used
                 as a material parameter or assertion value — comparison/plotting only.

Geometry convention (all): distance-to-hydration z in [0,10] cm; hydrated layer at z=0.
  MGR22/23: pellets 0-5 cm (bottom, hydrated), block 5-10 cm (top).
  MGR27   : block   0-5 cm (bottom, hydrated), pellets 5-10 cm (top)  [inverted].
Initial dry densities (EXACT, D5.6 text L231/L348): pellet 1.30, block 1.60 g/cm3 -> gap 0.30.
"""

# ---- EXPERIMENT (EXACT; D5.6 final-state tables) -----------------------------------------
# D5.6 L422-426 (MGR22/23), L456-460 (MGR27).
EXPERIMENT = {
  # test : dict(pellet, block, gap, src)   densities in g/cm3
  "MGR22": dict(pellet=1.35, block=1.51, gap=0.16, src="D5.6 Tab(L422-426) EXACT"),
  "MGR23": dict(pellet=1.34, block=1.51, gap=0.17, src="D5.6 Tab(L422-426) EXACT"),
  "MGR27": dict(pellet=1.434, block=1.454, gap=0.02, src="D5.6 Tab(L456-460) EXACT"),
}
# Axial swelling pressure (measured, friction-AFFECTED real sample):
EXPERIMENT_PSW = {  # MPa
  "MGR23": dict(psw=3.0, src="D5.6 Fig 3-140 'Reference Data' p.137 DIGITIZED ~+-0.3"),
  "MGR27": dict(psw=1.2, src="D5.7 Fig 5.3-13a 'experimental data' p.251 DIGITIZED ~+-0.3"),
}

# ---- UPC / CODE_BRIGHT  (BExM double-structure, live macro WRC + f_c/f_s; the reference) ---
# "captures very satisfactorily the final state ... in terms of dry density in all cases" (D5.7 5.3.7).
UPC = {
  "MGR23": dict(pellet=1.39, block=1.45, gap=0.06, src="D5.7 Fig 5.3-11a p.249 DIGITIZED"),
  "MGR27": dict(pellet=1.44, block=1.44, gap=0.02, src="D5.7 Fig 5.3-16a p.253 DIGITIZED"),
}
UPC_PSW = {"MGR27": dict(psw=2.9, src="D5.7 Fig 5.3-13a 'model' p.251 DIGITIZED (OVER-predicts; no friction)")}

# ---- BGR-old (Vinay & Steffen, BEACON submission; LE skeleton + chi(S)={1 if S=1 else 0}) --
# Same Bishop-cutoff weak HM coupling our native DSM uses. "density jump not homogenised"
# (D5.6 L4776); "high jump still present ... far more homogenised in the experiment" (p.136).
BGR_OLD = {
  "MGR22": dict(pellet=1.30, block=1.575, gap=0.275, src="D5.6 Fig 3-138 p.137 DIGITIZED (sim t=155d)"),
  "MGR23": dict(pellet=1.315, block=1.560, gap=0.245, src="D5.6 Fig 3-142 p.138 DIGITIZED (sim t=270d)"),
  "MGR27": dict(pellet=1.325, block=1.565, gap=0.24, src="D5.6 Fig 3-132/133 p.132-133 DIGITIZED (sim t=600d)"),
}
BGR_OLD_PSW = {  # under-predicts: weak chi=0 coupling, swelling law too weak
  "MGR22": dict(psw=1.3, src="D5.6 Fig 3-137 p.136 DIGITIZED (sim, vs exp ~3.0)"),
  "MGR23": dict(psw=1.4, src="D5.6 Fig 3-140 p.137 DIGITIZED (sim, vs exp ~3.0)"),
}

# ---- NATIVE DSM (OURS, 2026-06-01; Pi-path vdW disjoining, chi=0 BishopsSaturationCutoff) --
# VERIFIED from VTU (results/cal_base, mgr27_predict, mgr27_predict_calk0). EXACT (model output).
# block lambda=0; pellet k0 6.5e-16 lambda9; K pellet 3979 / block 25268 J/kg (Villar-frozen).
NATIVE_DSM = {
  # 'base'   = block k0 5e-19 ; 'maxk' = block k0 5e-15 (x10000, diagnostic upper bound)
  "MGR23": dict(base=dict(pellet=1.350, block=1.541, gap=0.192),
                maxk=dict(pellet=1.352, block=1.538, gap=0.186),
                src="VTU results/cal_base, cal_x10k EXACT(model)"),
  "MGR27": dict(base=dict(pellet=1.353, block=1.538, gap=0.185, psw=2.79),
                maxk=dict(pellet=1.345, block=1.547, gap=0.202, psw=3.21),
                src="VTU results/mgr27_predict, mgr27_predict_calk0 EXACT(model)"),
}

if __name__ == "__main__":
    import json
    print(json.dumps(dict(EXPERIMENT=EXPERIMENT, UPC=UPC, BGR_OLD=BGR_OLD,
                          NATIVE_DSM=NATIVE_DSM), indent=1, default=str))
