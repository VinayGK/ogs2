#!/usr/bin/env python3
"""
EPFL BEACON Task 3.3 — BOTH stress paths in e-vs-sigma_v space, published style
(D3.3 Fig 4-6 / Acta Geotech 2022 Fig 9). P1 (A-B-C-D free swell then load) AND
P2 (A-B'-C' const-volume swell then load), measured vs native-DSM model.

MEASURED: D3.3 Table 4-1 Experiment row + Acta Geotech 2022 (Ferrari et al.) Fig 9
  published e-log(sigma) loci; cited in caption.
MODEL (measured-from-VTU, results/rerun_*): P2 full A-B'-C' (MCC); P1 LE free-swell
  variant (MCC fails at TS1). e=phi/(1-phi), sigma_ax = -sigma_yy/1e6.
Honest result shown: model P1 free swell reaches e~1.07 only (measured 2.34) — the
small-strain / no-large-macro-opening limit, drawn not hidden.
"""
import glob, os, re
import numpy as np, meshio
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def traj(od, smin=1e-3):
    fs = sorted(glob.glob(os.path.join(ROOT, od, "*.vtu")),
                key=lambda f: float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)) if re.search(r"_t_([0-9.]+)\.vtu", f) else -1)
    S, E = [], []
    for f in fs:
        m = meshio.read(f)
        phi = float(np.asarray(m.point_data["porosity"]).mean())
        sax = -np.asarray(m.point_data["sigma"]).mean(axis=0)[1] / 1e6
        if sax > smin:
            S.append(sax); E.append(phi / (1 - phi))
    return np.array(S), np.array(E)

S_p2, E_p2 = traj("results/rerun_beacon_t33_path2_P2-1_dsm_mcc")
S_p1, E_p1 = traj("results/rerun_p1_LE")
# P1 compression: MCC restart from the LE-achieved free-swell state (e=1.069), not LE.
S_p1c, E_p1c = traj("results/rerun_p1_MCCrestart")

# MEASURED loci (D3.3 Tab 4-1 / Acta Fig 9)
mP1_s = [0.021, 0.021, 0.6, 1.0, 2.0, 3.24, 5.0, 10.0, 20.0]
mP1_e = [0.84,  2.34,  2.30, 2.05, 1.45, 1.10, 0.92, 0.72, 0.56]
mP2_s = [0.021, 0.5, 3.5, 5.0, 10.0, 20.0]
mP2_e = [0.84,  0.85, 0.85, 0.78, 0.67, 0.57]

fig, ax = plt.subplots(figsize=(7.4, 5.6))

ax.plot(mP1_s, mP1_e, "-",  color="#c0392b", lw=1.6, zorder=4,
        label="measured P1 ($A\\!-\\!B\\!-\\!C\\!-\\!D$, free swell)")
ax.plot(mP1_s, mP1_e, "o",  color="#c0392b", ms=5, zorder=5, mec="k", mew=0.4)
ax.plot(mP2_s, mP2_e, "--", color="#e08e0b", lw=1.6, zorder=4,
        label="measured P2 ($A\\!-\\!B'\\!-\\!C'$, const. volume)")
ax.plot(mP2_s, mP2_e, "s",  color="#e08e0b", ms=5, zorder=5, mec="k", mew=0.4)

ax.plot(S_p2, E_p2, "-o", color="#1f5fa8", ms=4, lw=1.6, zorder=3,
        label="DSM model P2 ($A\\!-\\!B'\\!-\\!C'$, MCC)")
# P1 model: LE free-swell A->B (vertical at ~21 kPa, e 0.83 -> 1.07), then MCC compression B->C->D
ax.plot([0.021, 0.021], [0.83, E_p1.max()], "-^", color="#5fa8d3", ms=4, lw=1.4, zorder=3,
        label="DSM model P1 free-swell (LE, $A\\!\\to\\!B$)")
ax.plot(S_p1c, E_p1c, "-D", color="#7b3fa0", ms=4, lw=1.5, zorder=3,
        label="DSM model P1 compression (MCC, $B\\!\\to\\!C\\!\\to\\!D$)")

pts = {"A": (0.021, 0.84), "B": (0.021, 2.34), "B'": (3.5, 0.85),
       "C": (3.24, 1.10), "C'": (20.0, 0.57), "D": (20.0, 0.56)}
for n, (s, e) in pts.items():
    ax.annotate(n, (s, e), textcoords="offset points", xytext=(6, 5),
                fontsize=11, fontweight="bold", zorder=6)

ax.annotate("model $A\\to B$ free swell:\n$e$ 0.83$\\to$1.07 only\n(measured 2.34; small-strain caps it)\n"
            "then MCC compresses $B\\to D$",
            xy=(0.021, 1.07), xytext=(0.028, 1.62), fontsize=7.8, color="#1f4e79",
            arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.0))
ax.annotate("", xy=(3.5, 1.10), xytext=(3.5, 0.85),
            arrowprops=dict(arrowstyle="<->", color="k", lw=1.0))
ax.text(3.9, 0.97, "$\\Delta e\\approx0.26$\n(path dep.)", fontsize=8)

ax.set_xscale("log"); ax.set_xlim(0.01, 60); ax.set_ylim(0.3, 2.6)
ax.set_xlabel("total vertical stress  $\\sigma_v$  [MPa]", fontsize=11)
ax.set_ylabel("void ratio  $e$  [-]", fontsize=11)
ax.set_title("BEACON Task 3.3 granular MX-80 — both stress paths, $e$ vs $\\sigma_v$\n"
             "measured (Ferrari et al. 2022; D3.3 Tab 4-1) vs native DSM (Pi-path)", fontsize=9.5)
ax.grid(True, which="both", ls=":", alpha=0.4)
ax.legend(fontsize=8, loc="upper right", framealpha=0.95)
fig.tight_layout()
out = os.path.join(ROOT, "figures", "epfl_both_paths.png")
fig.savefig(out, dpi=150); fig.savefig(out.replace(".png", ".pdf"))
print("wrote", out)
print("model P1 free-swell e: 0.83 -> %.3f (meas A->B 0.84->2.34)" % E_p1.max())
print("model P1 MCC compression endpoint e=%.3f at %.1f MPa (meas D 0.56)" % (E_p1c[-1], S_p1c[-1]))
print("model P2 endpoint e=%.3f at %.1f MPa (meas C' 0.57)" % (E_p2[-1], S_p2[-1]))
