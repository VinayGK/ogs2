#!/usr/bin/env python3
"""
EPFL BEACON Task 3.3 — void ratio vs total vertical stress, in the published style
(D3.3 Fig 4-6 / Acta Geotech 2022 Fig 9: e on linear y, total vertical stress on log x,
points A/B/B'/C/C'/D marked). Native-DSM (Pi-path) + MFront MCC, lambda=9 perm rerun.

Model curves (measured-from-VTU, results/rerun_*):
  - P2-1 full A-B'-C': isochoric swell to B' then oedometric compression.
  - homogeneous-compression: B'-C' alone from the saturated state (overlay).
Measured anchors (Acta Geotech 2022 Table/Figs; D3.3 Table 4-1) plotted as points.
NO fitted/invented values: model points are extracted; measured points are the
published Experiment-row values, cited in the caption.
"""
import glob, re, os
import numpy as np, meshio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def traj(od):
    fs = glob.glob(os.path.join(ROOT, od, "*.vtu"))
    gt = lambda f: float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)) if re.search(r"_t_([0-9.]+)\.vtu", f) else -1
    fs = sorted(fs, key=gt)
    S, E = [], []
    for f in fs:
        m = meshio.read(f)
        phi = float(np.asarray(m.point_data["porosity"]).mean())
        sig = np.asarray(m.point_data["sigma"]).mean(axis=0)
        sax = -sig[1] / 1e6
        if sax > 1e-3:
            S.append(sax); E.append(phi / (1 - phi))
    return np.array(S), np.array(E)

S_p2, E_p2 = traj("results/rerun_beacon_t33_path2_P2-1_dsm_mcc")
S_hc, E_hc = traj("results/rerun_beacon_t33_path2_compression_homogeneous_mcc")

# MEASURED Task-3.3 anchors (D3.3 Table 4-1 Experiment row; Acta Geotech 2022).
# A as-poured; B free-swell; B' swelling pressure; C path-1 reload pt; C'/D 20 MPa endpoints.
meas = {
    "A":  (0.021, 0.84),    # as-poured, sigma_a=21 kPa
    "B":  (0.021, 2.34),    # free swell (P1)
    "B'": (3.5,   0.85),    # swelling pressure (P2), Ps 3.2-3.6 MPa
    "C":  (3.24,  1.10),    # P1 reload point (D3.3:2293)
    "C'": (20.0,  0.57),    # P2 endpoint
    "D":  (20.0,  0.56),    # P1 endpoint
}

fig, ax = plt.subplots(figsize=(7.0, 5.2))

# model P2-1 full path
ax.plot(S_p2, E_p2, "-o", color="#1f5fa8", ms=4, lw=1.6,
        label=r"DSM model, P2-1 ($A\to B'\to C'$)")
# model homogeneous compression overlay (dashed)
ax.plot(S_hc, E_hc, "--s", color="#3a8f4f", ms=3.5, lw=1.3,
        label=r"DSM model, homogenised $B'\!-\!C'$")

# measured points
for name, (s, e) in meas.items():
    red = name in ("B'", "C'")
    ax.plot(s, e, "^" if red else "o", color="#c0392b", ms=8, zorder=5,
            markeredgecolor="k", markeredgewidth=0.5)
    dx = 1.15 if name in ("A",) else 1.0
    ax.annotate(name, (s, e), textcoords="offset points",
                xytext=(6, 6), fontsize=11, fontweight="bold")
ax.plot([], [], "^", color="#c0392b", ms=8, markeredgecolor="k",
        markeredgewidth=0.5, label="measured (Ferrari et al. 2022; D3.3 Tab 4-1)")

# the path-dependency gap annotation at the swelling stress (~3.5 MPa)
ax.annotate("", xy=(3.5, 1.10), xytext=(3.5, 0.85),
            arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.2))
ax.text(3.9, 0.97, r"$\Delta e\approx0.26$", color="#c0392b", fontsize=9)

ax.set_xscale("log")
ax.set_xlim(0.01, 60)
ax.set_ylim(0.3, 2.6)
ax.set_xlabel("total vertical stress  $\\sigma_v$  [MPa]", fontsize=11)
ax.set_ylabel("void ratio  $e$  [-]", fontsize=11)
ax.set_title("BEACON Task 3.3 — granular MX-80, $e$ vs $\\sigma_v$\n"
             "native DSM (Pi-path + MCC, $k(\\phi_M)^9$) vs measured", fontsize=10)
ax.grid(True, which="both", ls=":", alpha=0.4)
ax.legend(fontsize=8.5, loc="upper right", framealpha=0.95)
fig.tight_layout()
out = os.path.join(ROOT, "figures", "epfl_e_sigma.png")
fig.savefig(out, dpi=150)
fig.savefig(out.replace(".png", ".pdf"))
print("wrote", out)
print("\nmodel P2-1 endpoint: sigma=%.1f MPa, e=%.3f" % (S_p2[-1], E_p2[-1]))
print("measured C'/D at 20 MPa: e=0.56-0.57")
print("model B' (swelling pressure): %.2f MPa (measured 3.2-3.6)" %
      (S_p2[(np.abs(E_p2-0.85)).argmin()] if len(S_p2) else float('nan')))
