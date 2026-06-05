#!/usr/bin/env python3
"""EPFL BEACON Task 3.3 — e vs sigma_v, model (current_porosity_split, K re-fit to
Dixon MX-80) vs measured (Ferrari et al., Acta Geotechnica 2022, Fig 9 / D3.3 Tab 4-1).
Model from VTU: e = phi/(1-phi), sigma_ax = -sigma_yy/1e6 (compression-positive).
Mode: micro_solid_volume_fraction_mode=current_porosity_split (nS=1-n_l, micro-only).
K = 56183 J/kg @ rho_d=1500 (log-interp of recalibrated 36796@1400 / 85776@1600)."""
import glob, os, re
import numpy as np, meshio
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def traj(od, smin=1e-3):
    fs = sorted(glob.glob(os.path.join(ROOT, "results", od, "*.vtu")),
                key=lambda f: float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)) if re.search(r"_t_([0-9.]+)\.vtu", f) else -1)
    S, E = [], []
    for f in fs:
        m = meshio.read(f)
        phi = float(np.asarray(m.point_data["porosity"]).mean())
        sax = -np.asarray(m.point_data["sigma"]).mean(axis=0)[1] / 1e6
        if sax > smin:
            S.append(sax); E.append(phi / (1 - phi))
    return np.array(S), np.array(E)

S_p2, E_p2 = traj("run_beacon_t33_path2_P2-1_dsm_mcc")
S_le, E_le = traj("run_p1_LE", smin=-1); ipk = int(np.argmax(E_le))
S_p1sw, E_p1sw = S_le[:ipk+1], E_le[:ipk+1]
S_p1c, E_p1c = traj("run_p1_MCCrestart", smin=-1)

mP1_s = [0.021, 0.021, 0.6, 1.0, 2.0, 3.24, 5.0, 10.0, 20.0]
mP1_e = [0.84,  2.34,  2.30, 2.05, 1.45, 1.10, 0.92, 0.72, 0.56]
mP2_s = [0.021, 0.5, 3.5, 5.0, 10.0, 20.0]
mP2_e = [0.84,  0.85, 0.85, 0.78, 0.67, 0.57]

fig, ax = plt.subplots(figsize=(7.6, 5.8))
ax.plot(mP1_s, mP1_e, "-o", color="#c0392b", lw=1.6, ms=5, mec="k", mew=0.4, zorder=5, label="measured P1 (free swell A-B-C-D)")
ax.plot(mP2_s, mP2_e, "--s", color="#e08e0b", lw=1.6, ms=5, mec="k", mew=0.4, zorder=5, label="measured P2 (const. vol A-B'-C')")
ax.plot(S_p2, E_p2, "-o", color="#1f5fa8", ms=4, lw=1.6, zorder=3, label="model P2 (MCC)")
ax.plot(S_p1sw, E_p1sw, "-^", color="#5fa8d3", ms=4, lw=1.4, zorder=3, label="model P1 free-swell (LE, A->B)")
ax.plot(S_p1c, E_p1c, "-D", color="#7b3fa0", ms=4, lw=1.5, zorder=3, label="model P1 compression (MCC, B->D)")
for n,(s,e) in {"A":(0.021,0.84),"B":(0.021,2.34),"B'":(3.5,0.86),"C":(3.24,1.10),"C'":(20.0,0.57),"D":(20.0,0.56)}.items():
    ax.annotate(n,(s,e),textcoords="offset points",xytext=(6,5),fontsize=11,fontweight="bold",zorder=6)
ax.text(0.012, 0.34,
        "micro-only disjoining  nS = 1 - n_l  (current-porosity-split)\n"
        "K re-fit to Dixon MX-80 under the new mode:\n"
        "K = 56183 J/kg @ rho_d=1500 (log-interp 36796@1400 / 85776@1600)",
        fontsize=7.6, color="#1f4e79",
        bbox=dict(boxstyle="round,pad=0.3", fc="#eef5ff", ec="#1f5fa8", lw=0.8))
ax.set_xscale("log"); ax.set_xlim(0.01, 60); ax.set_ylim(0.3, 2.6)
ax.set_xlabel("total vertical stress  $\\sigma_v$  [MPa]", fontsize=11)
ax.set_ylabel("void ratio  $e$  [-]", fontsize=11)
ax.set_title("BEACON Task 3.3 granular MX-80 — $e$ vs $\\sigma_v$\n"
             "measured (Ferrari et al. 2022) vs model (current-porosity-split, K re-fit to Dixon)", fontsize=9.4)
ax.grid(True, which="both", ls=":", alpha=0.4); ax.legend(fontsize=8, loc="upper right", framealpha=0.95)
fig.tight_layout()
out = os.path.join(ROOT, "figures", "epfl_cps_both_paths.png")
fig.savefig(out, dpi=150)
print("wrote", out)
print(f"P1 free-swell peak e = {E_le[ipk]:.3f} (meas B 2.34)")
print(f"P1 MCC compression endpoint e = {E_p1c[-1]:.3f} @ {S_p1c[-1]:.1f} MPa (meas D 0.56)")
print(f"P2 endpoint e = {E_p2[-1]:.3f} @ {S_p2[-1]:.1f} MPa (meas C' 0.57)")
