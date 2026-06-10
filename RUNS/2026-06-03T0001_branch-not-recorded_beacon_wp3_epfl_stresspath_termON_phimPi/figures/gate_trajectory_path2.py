#!/usr/bin/env python3
"""
EPFL path2 (P2-1) on the SCALED (p'>=phi_m*Pi) Maxwell binary: the gate trajectory.
Shows why the rescaled term still does not bite -- p' hugs the phi_m*Pi operating
point (open only at the strain-free start and a marginal endpoint), far below the
intrinsic Pi. term-ON == term-OFF bit-for-bit (max|Δ|=0 over all fields/steps).
"""
import glob, os, re
import numpy as np, meshio
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
run = os.path.join(ROOT, "results", "run_beacon_t33_path2_P2-1_dsm_mcc")
def tkey(f):
    m = re.search(r"_t_([0-9.]+)\.vtu", f); return float(m.group(1)) if m else -1
fs = sorted(glob.glob(run + "/*.vtu"), key=tkey)
t, pconf, gate, Pi = [], [], [], []
for f in fs:
    m = meshio.read(f); pd = m.point_data; sig = np.asarray(pd["sigma"])
    t.append(tkey(f) / 86400.0)
    pconf.append(float((-(sig[:, 0] + sig[:, 1] + sig[:, 2]) / 3).mean() / 1e6))
    mp = float(np.asarray(pd["micro_pressure"]).mean() / 1e6)
    phiM = float(np.asarray(pd["micro_porosity"]).mean())
    gate.append(phiM * mp); Pi.append(mp)
t = np.array(t); pconf = np.array(pconf); gate = np.array(gate); Pi = np.array(Pi)

fig, ax = plt.subplots(figsize=(7.6, 5.0))
ax.plot(t, Pi, "-", color="#999999", lw=1.3, label=r"$\Pi$ (intrinsic micro disjoining) = coded gate")
ax.plot(t, gate, "-s", color="#e08e0b", ms=3, lw=1.6, label=r"$\phi_m\Pi$ (REV gate, option 1)")
ax.plot(t, pconf, "-o", color="#1f5fa8", ms=3, lw=1.8, label=r"$p'=-\mathrm{tr}(\sigma')/3$ (confining)")
# shade where the REV gate is open
op = pconf >= gate
ax.fill_between(t, 0, 1, where=op, transform=ax.get_xaxis_transform(),
                color="#5fa8d3", alpha=0.15, label="REV gate OPEN")
ax.annotate("open but $\\varepsilon_v\\approx0$\n(term $\\propto\\varepsilon_v$ = 0)",
            xy=(8, 1.5), fontsize=8, color="#1f4e79", ha="center")
ax.annotate("open, but marginal &\nlate (lagged gate misses)",
            xy=(237, 11.5), xytext=(180, 17), fontsize=8, color="#1f4e79",
            arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.0))
ax.text(70, 19.5,
        "term-ON $\\equiv$ term-OFF bit-for-bit (max$|\\Delta|=0$):\n"
        "$p'$ hugs the $\\phi_m\\Pi$ operating point through the whole path",
        fontsize=8.5, color="#7a2d00",
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff3e6", ec="#e08e0b", lw=0.8))
ax.set_xlabel("time  [d]  (const-volume swell $\\to$ compression)", fontsize=11)
ax.set_ylabel("stress  [MPa]", fontsize=11)
ax.set_ylim(0, 30)
ax.set_title("EPFL path2 on the scaled binary: $p'$ grazes the $\\phi_m\\Pi$ gate, never reaches $\\Pi$",
             fontsize=10)
ax.grid(True, ls=":", alpha=0.4); ax.legend(fontsize=8, loc="center right")
fig.tight_layout()
out = os.path.join(ROOT, "figures", "gate_trajectory_path2.png")
fig.savefig(out, dpi=150); fig.savefig(out.replace(".png", ".pdf"))
print("wrote", out)
print(f"p' range [{pconf.min():.2f},{pconf.max():.2f}] | phi_m*Pi range [{gate.min():.2f},{gate.max():.2f}] | Pi range [{Pi.min():.1f},{Pi.max():.1f}]")
print(f"REV-gate-open steps: {int(op.sum())}/{len(op)} (start ε_v≈0; end marginal)")
