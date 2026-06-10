#!/usr/bin/env python3
"""
MGR23 team-style comparison: final dry-density profile vs distance from hydration face.
Layout follows D5.7 Fig 5.3-11a / D4.1 Fig 5-10: rho_d (g/cm3) on x or y vs sample height;
here height (mm, 0 = hydration/pellet end) on x, rho_d on y, with the INITIAL step profile
(pellet 1.30 / block 1.60) as reference and the MEASURED final (D4.1 Table 5-3) overlaid.

Model curves: native DSM (Pi-path) measured-from-VTU. rho_d = 2700*(1-phi)/1000.
Honest result: model retains a near-initial STEP at the 50 mm interface (little
homogenisation) vs the measured flattened profile -> the missing macro-redistribution.
"""
import glob, os, re
import numpy as np, meshio
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RHO_S = 2700.0

def profile(od, prefix):
    fs = sorted(glob.glob(os.path.join(ROOT, "results", od, prefix + "_*.vtu")),
                key=lambda f: float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)))
    m = meshio.read(fs[-1]); y = m.points[:, 1]; phi = np.asarray(m.point_data["porosity"])
    zs = np.array(sorted(set(np.round(y, 6))))
    rd = np.array([RHO_S * (1 - phi[np.abs(y - z) < 1e-7].mean()) / 1000 for z in zs])
    return zs * 1000, rd  # mm, g/cm3

z_c, rd_c = profile("mgr23c", "mgr23c_powerlaw")     # block 5e-21
z_d, rd_d = profile("mgr23d", "mgr23d_block100x")    # block 5e-19 (100x)

# INITIAL step (as-built): pellet 0-50mm = 1.30, block 50-100mm = 1.60 (D4.1 Tab 5-3 MGR23)
z_init = [0, 50, 50, 100]; rd_init = [1.30, 1.30, 1.60, 1.60]
# MEASURED final (D4.1 Tab 5-3 MGR23 layer averages): pellet 1.34 (mid 25mm), block 1.51 (mid 75mm)
# D4.1 Fig 5-10 shows a near-linear flattened profile; plot the two layer-average points + trend.
z_meas = [25, 75]; rd_meas = [1.34, 1.51]

fig, ax = plt.subplots(figsize=(7.2, 5.2))

ax.plot(rd_init, z_init, ":", color="gray", lw=1.8, label="initial (as-built step)")
ax.plot(rd_c, z_c, "-o", color="#1f5fa8", ms=4, lw=1.5, label="DSM model (block $k$=5e-21)")
ax.plot(rd_d, z_d, "-s", color="#3a8f4f", ms=4, lw=1.5, label="DSM model (block $k$=5e-19, 100$\\times$)")
ax.plot(rd_meas, z_meas, "^", color="#c0392b", ms=11, mec="k", mew=0.6,
        label="measured final (D4.1 Tab 5-3)")
# measured layer-average bands (pellet 0-50, block 50-100)
ax.plot([1.34, 1.34], [0, 50], "-", color="#c0392b", lw=2.5, alpha=0.5)
ax.plot([1.51, 1.51], [50, 100], "-", color="#c0392b", lw=2.5, alpha=0.5)

ax.axhline(50, color="k", ls="--", lw=0.7, alpha=0.5)
ax.text(1.62, 51, "pellet | block interface", fontsize=7.5, va="bottom")
ax.text(1.62, 5,  "PELLET (hydration end)", fontsize=8, color="#555")
ax.text(1.62, 95, "BLOCK (top)", fontsize=8, color="#555", va="top")

ax.set_xlabel("dry density  $\\rho_d$  [g/cm$^3$]", fontsize=11)
ax.set_ylabel("distance from hydration face  [mm]", fontsize=11)
ax.set_xlim(1.25, 1.68); ax.set_ylim(0, 100)
ax.set_title("BEACON WP5.3 MGR23 — final dry-density profile\n"
             "native DSM (FEBEX) vs measured (CIEMAT)", fontsize=10)
ax.grid(True, ls=":", alpha=0.4)
ax.legend(fontsize=8, loc="lower left", framealpha=0.95)
fig.tight_layout()
out = os.path.join(ROOT, "figures", "mgr_drydensity_profile.png")
fig.savefig(out, dpi=150); fig.savefig(out.replace(".png", ".pdf"))
print("wrote", out)
print("model step at interface (block5e-19): pellet~%.2f block~%.2f (gap %.2f)" %
      (rd_d[z_d < 50].mean(), rd_d[z_d > 50].mean(), rd_d[z_d > 50].mean() - rd_d[z_d < 50].mean()))
print("measured: pellet 1.34 block 1.51 (gap 0.17); initial gap 0.30")
