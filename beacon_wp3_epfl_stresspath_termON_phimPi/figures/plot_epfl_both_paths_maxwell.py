#!/usr/bin/env python3
"""
EPFL BEACON Task 3.3 — BOTH stress paths in e-vs-sigma_v space, on the
maxwell-conjugate binary (full-Pi closure + Maxwell-conjugate term).

Adapted from beacon_wp3_epfl_repro_2026-05-31/figures/plot_epfl_both_paths.py.
MEASURED loci: Ferrari et al., Acta Geotechnica 2022, Fig 9 / D3.3 Tab 4-1
  (digitised; cited in caption). UNCHANGED from the original.
MODEL (measured-from-VTU): this run's results/ dirs. e = phi/(1-phi),
  sigma_ax = -sigma_yy/1e6 (axial, soil-mechanics compression-positive).

INTEGRITY NOTE drawn on the figure: the Maxwell-conjugate term is INERT on every
EPFL path here (coded gate p' >= Pi never opens: max p'~10 MPa < Pi~22-28 MPa),
so these curves are the full-Pi baseline, NOT a demonstration of the term.
"""
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

# --- model trajectories (this run) ---
S_p2, E_p2 = traj("run_beacon_t33_path2_P2-1_dsm_mcc")            # P2 full A-B'-C'
S_le, E_le = traj("run_p1_LE", smin=-1)                           # P1 LE: swell then compress
ipk = int(np.argmax(E_le))                                       # free-swell peak = B
S_p1sw, E_p1sw = S_le[:ipk+1], E_le[:ipk+1]                      # P1 free-swell branch A->B
S_p1c, E_p1c = traj("run_p1_MCCrestart", smin=-1)                # P1 MCC compression B->C->D

# --- MEASURED loci (Ferrari et al. 2022 Fig 9 / D3.3 Tab 4-1) — unchanged ---
mP1_s = [0.021, 0.021, 0.6, 1.0, 2.0, 3.24, 5.0, 10.0, 20.0]
mP1_e = [0.84,  2.34,  2.30, 2.05, 1.45, 1.10, 0.92, 0.72, 0.56]
mP2_s = [0.021, 0.5, 3.5, 5.0, 10.0, 20.0]
mP2_e = [0.84,  0.85, 0.85, 0.78, 0.67, 0.57]

fig, ax = plt.subplots(figsize=(7.6, 5.8))

# measured
ax.plot(mP1_s, mP1_e, "-",  color="#c0392b", lw=1.6, zorder=4, label="measured P1 (free swell $A\\!-\\!B\\!-\\!C\\!-\\!D$)")
ax.plot(mP1_s, mP1_e, "o",  color="#c0392b", ms=5, zorder=5, mec="k", mew=0.4)
ax.plot(mP2_s, mP2_e, "--", color="#e08e0b", lw=1.6, zorder=4, label="measured P2 (const. volume $A\\!-\\!B'\\!-\\!C'$)")
ax.plot(mP2_s, mP2_e, "s",  color="#e08e0b", ms=5, zorder=5, mec="k", mew=0.4)

# model
ax.plot(S_p2, E_p2, "-o", color="#1f5fa8", ms=4, lw=1.6, zorder=3, label="model P2 (full-$\\Pi$, MCC)")
ax.plot(S_p1sw, E_p1sw, "-^", color="#5fa8d3", ms=4, lw=1.4, zorder=3, label="model P1 free-swell (LE, $A\\!\\to\\!B$)")
ax.plot(S_p1c, E_p1c, "-D", color="#7b3fa0", ms=4, lw=1.5, zorder=3, label="model P1 compression (MCC, $B\\!\\to\\!D$)")

# stage-point labels (measured)
pts = {"A": (0.021, 0.84), "B": (0.021, 2.34), "B'": (3.5, 0.86),
       "C": (3.24, 1.10), "C'": (20.0, 0.57), "D": (20.0, 0.56)}
for n, (s, e) in pts.items():
    ax.annotate(n, (s, e), textcoords="offset points", xytext=(6, 5), fontsize=11, fontweight="bold", zorder=6)

# honest annotations
ax.annotate("model free swell caps at\n$e\\approx%.2f$ (measured 2.34)\n— small-strain limit" % E_le[ipk],
            xy=(0.021, E_le[ipk]), xytext=(0.03, 1.5), fontsize=8, color="#1f4e79",
            arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.0))
ax.text(0.012, 0.34,
        "Maxwell-conjugate term INERT on all paths\n(coded gate $p'\\!\\geq\\!\\Pi$ never opens: $p'\\!\\sim\\!10$ MPa $<\\Pi\\!\\sim\\!25$ MPa)\n"
        "$\\Rightarrow$ these curves = full-$\\Pi$ baseline; term needs gate rescale ($p'\\!\\geq\\!\\phi_m\\Pi$)",
        fontsize=7.4, color="#7a2d00",
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff3e6", ec="#e08e0b", lw=0.8))

ax.set_xscale("log"); ax.set_xlim(0.01, 60); ax.set_ylim(0.3, 2.6)
ax.set_xlabel("total vertical stress  $\\sigma_v$  [MPa]", fontsize=11)
ax.set_ylabel("void ratio  $e$  [-]", fontsize=11)
ax.set_title("BEACON Task 3.3 granular MX-80 — both stress paths, $e$ vs $\\sigma_v$\n"
             "measured (Ferrari et al. 2022) vs model (maxwell-conjugate binary, full-$\\Pi$; term inert)", fontsize=9.2)
ax.grid(True, which="both", ls=":", alpha=0.4)
ax.legend(fontsize=8, loc="upper right", framealpha=0.95)
fig.tight_layout()
out = os.path.join(ROOT, "figures", "epfl_both_paths_maxwell.png")
fig.savefig(out, dpi=150); fig.savefig(out.replace(".png", ".pdf"))
print("wrote", out)
print(f"model P1 free-swell peak e = {E_le[ipk]:.3f} (meas B 2.34)")
print(f"model P1 MCC compression endpoint e = {E_p1c[-1]:.3f} @ {S_p1c[-1]:.1f} MPa (meas D 0.56)")
print(f"model P2 endpoint e = {E_p2[-1]:.3f} @ {S_p2[-1]:.1f} MPa (meas C' 0.57)")
