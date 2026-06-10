#!/usr/bin/env python3
"""EBS validation bundle 2026-06-10 — measurement-comparison figures, re-generated
on the h_of_eps_20260609 binary outputs.

Adapted ONLY in input paths from the cases' own scripts:
  * MGR profile panel + team data: validation_2026-06-09/mgr27/results/team_comparison_data.py
    (data module, EXACT vs DIGITIZED provenance as documented there) and
    beacon_wp5_mgr_repro_2026-06-01/make_team_comparison.py (profile() and panel layout).
  * EPFL measured loci: beacon_wp3_epfl_stresspath_2026-06-02/figures/
    plot_epfl_both_paths_maxwell.py (Ferrari et al. 2022 Fig 9 / D3.3 Tab 4-1, digitized,
    copied UNCHANGED from that script).
No measurement value is introduced here that is not in those case-own files.
"""
import sys, glob, re, json, os
import numpy as np, meshio
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

EBS = "/Users/vinaykumar/git/ogs/RUNS/_INPROGRESS_full_validation/ebs"
sys.path.insert(0, "/Users/vinaykumar/git/ogs/validation_2026-06-09/mgr27/results")
from team_comparison_data import EXPERIMENT, UPC, BGR_OLD, NATIVE_DSM, EXPERIMENT_PSW

RHO_S = 2700.0  # kg/m3, FEBEX solid density as used by the case's own profile script

def vtus(d, prefix):
    return sorted(glob.glob(f"{d}/{prefix}_*.vtu"),
                  key=lambda f: float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)))

def profile(vtu, nbin=25):
    m = meshio.read(vtu); y = m.points[:, 1]
    phi = np.asarray(m.point_data["porosity"])
    zc = np.linspace(0, 0.10, nbin + 1); zmid = 0.5 * (zc[:-1] + zc[1:]); rd = []
    for i in range(nbin):
        s = (y >= zc[i]) & (y <= zc[i + 1])
        rd.append(RHO_S * (1 - phi[s].mean()) / 1000 if s.any() else np.nan)
    return zmid * 100, np.array(rd)  # cm, g/cm3

def layer_means(vtu):
    m = meshio.read(vtu); y = m.points[:, 1]
    phi = np.asarray(m.point_data["porosity"])
    rd = RHO_S * (1 - phi) / 1000
    return rd[y < 0.05].mean(), rd[y > 0.05].mean()

def szz_top(vtu):
    m = meshio.read(vtu); y = m.points[:, 1]
    sig = np.asarray(m.point_data["sigma"])
    top = np.abs(y - y.max()) < 1e-9
    return -sig[top, 1].mean() / 1e6  # MPa, compression positive

# ---------------- MGR27 (completed run) ----------------
d27 = f"{EBS}/out_mgr27/run_mgr27_le_calk0"
f27 = vtus(d27, "mgr27_predict_calk0")
z27, rd27 = profile(f27[-1])
blk27, pel27 = layer_means(f27[-1])   # MGR27: block at bottom (z<5), pellets top
gap27 = blk27 - pel27
psw27 = szz_top(f27[-1])
e = EXPERIMENT["MGR27"]

fig, ax = plt.subplots(figsize=(5.6, 5.3))
ax.plot([1.60, 1.60], [0, 5], "--", c="firebrick", lw=1.2, label="initial (1.30/1.60)")
ax.plot([1.30, 1.30], [5, 10], "--", c="firebrick", lw=1.2)
ax.plot(rd27, z27, "-o", c="#1f77b4", ms=3, lw=1.8,
        label="Native DSM final (h_of_eps 2026-06-10, VTU)", zorder=5)
ax.plot([e["pellet"], e["block"]], [7.5, 2.5], "s", c="k", ms=9,
        label="Experiment final (exact, layer avg)", zorder=6)
u = UPC["MGR27"]; b = BGR_OLD["MGR27"]
ax.plot([u["pellet"], u["block"]], [7.5, 2.5], "^", c="#2ca02c", ms=8, mfc="none",
        mew=1.6, label="UPC BExM (digitized)")
ax.plot([b["pellet"], b["block"]], [7.5, 2.5], "v", c="#7f7f7f", ms=8, mfc="none",
        mew=1.6, label="BGR-old (digitized)")
ax.axhline(5, c="gray", lw=0.8, ls=":")
ax.set_title("MGR27 (hydrate z=0; block / pellets)\nh_of_eps_20260609 re-run", fontsize=10)
ax.set_xlabel("dry density  [g/cm$^3$]"); ax.set_ylabel("distance to hydration surface  [cm]")
ax.set_xlim(1.2, 1.7); ax.set_ylim(0, 10); ax.grid(alpha=0.3)
ax.legend(fontsize=7.4, loc="upper center", framealpha=0.95)
fig.tight_layout(); fig.savefig(f"{EBS}/figures/mgr27_profile_team_20260610.png", dpi=155)
plt.close(fig)

# ---------------- MGR23 (failed run — partial-state consistency vs case's own partial json) ----------------
d23 = f"{EBS}/out_mgr23/run_mgr23_column_calibrated"
f23 = vtus(d23, "sw_A_blk0")
rows = []
for f in f23:
    t_d = float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)) / 86400.0
    pel, blk = layer_means(f)   # MGR23: pellets at bottom (z<5)
    rows.append(dict(t_d=t_d, pel_rd=pel, blk_rd=blk, gap=blk - pel, szz_top=szz_top(f)))
ref = json.load(open("/Users/vinaykumar/git/ogs/validation_2026-06-09/mgr23/results/mgr23_reduced_partial.json"))

fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.4))
axs[0].plot([r["t_d"] for r in rows], [r["gap"] for r in rows], "-o", c="#1f77b4",
            label="this re-run (calibrated, partial)")
axs[0].plot([r["t_d"] for r in ref], [r["gap"] for r in ref], "s", c="#7f7f7f", mfc="none",
            label="case record 2026-06-09 (partial json)")
axs[0].axhline(EXPERIMENT["MGR23"]["gap"], ls="--", c="k", lw=1,
               label=f'experiment FINAL gap {EXPERIMENT["MGR23"]["gap"]} (D5.6, exact)')
axs[0].set_xlabel("time [d]"); axs[0].set_ylabel("dry-density gap block−pellet [g/cm$^3$]")
axs[0].set_title("MGR23 — gap evolution (run FAILED at t≈40 d)", fontsize=9.5)
axs[0].legend(fontsize=7.5); axs[0].grid(alpha=0.3)
axs[1].plot([r["t_d"] for r in rows], [r["szz_top"] for r in rows], "-o", c="#1f77b4")
axs[1].plot([r["t_d"] for r in ref], [r["szz_top"] for r in ref], "s", c="#7f7f7f", mfc="none")
axs[1].axhline(EXPERIMENT_PSW["MGR23"]["psw"], ls="--", c="k", lw=1,
               label=f'experiment axial p_sw {EXPERIMENT_PSW["MGR23"]["psw"]} MPa (digitized ±0.3)')
axs[1].set_xlabel("time [d]"); axs[1].set_ylabel("axial stress at top [MPa]")
axs[1].set_title("MGR23 — axial swelling pressure (partial)", fontsize=9.5)
axs[1].legend(fontsize=7.5); axs[1].grid(alpha=0.3)
fig.suptitle("MGR23 h_of_eps_20260609 re-run — PARTIAL (dt-floor failure at t≈3.47e6 s); "
             "no final-state comparison possible", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(f"{EBS}/figures/mgr23_partial_consistency_20260610.png", dpi=155)
plt.close(fig)

# ---------------- EPFL T33 — only the homogeneous P2 compression case parsed/ran ----------------
# measured P2 locus copied UNCHANGED from the case-family script plot_epfl_both_paths_maxwell.py
# (Ferrari et al., Acta Geotechnica 2022, Fig 9 / D3.3 Tab 4-1, digitized):
mP2_s = [0.021, 0.5, 3.5, 5.0, 10.0, 20.0]
mP2_e = [0.84, 0.85, 0.85, 0.78, 0.67, 0.57]

dep = f"{EBS}/out_epfl_t33/run_beacon_t33_path2_compression_homogeneous_mcc"
S, E = [], []
for f in vtus(dep, "beacon_t33_path2_compression_homog"):
    m = meshio.read(f)
    phi = float(np.asarray(m.point_data["porosity"]).mean())
    sax = -np.asarray(m.point_data["sigma"]).mean(axis=0)[1] / 1e6
    if sax > 1e-3:
        S.append(sax); E.append(phi / (1 - phi))
fig, ax = plt.subplots(figsize=(6.4, 5.0))
ax.semilogx(mP2_s, mP2_e, "s--", c="k", ms=7, label="measured P2 locus (Ferrari 2022 Fig 9, digitized)")
ax.semilogx(S, E, "-o", c="#1f77b4", ms=4, label="model: P2 homogeneous MCC compression (this run)")
ax.set_xlabel("axial stress  $\\sigma_v$  [MPa]"); ax.set_ylabel("void ratio  $e$")
ax.set_title("EPFL T3.3 — Path-2 compression (only runnable sub-case)\n"
             "4/5 PRJs BLOCKED at parse on h_of_eps_20260609", fontsize=9.5)
ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{EBS}/figures/epfl_t33_p2_compression_20260610.png", dpi=155)
plt.close(fig)

# ---------------- numbers ----------------
out = dict(
    MGR27=dict(model=dict(pellet=round(pel27, 3), block=round(blk27, 3),
                          gap=round(gap27, 3), psw_MPa=round(psw27, 2)),
               experiment=dict(pellet=e["pellet"], block=e["block"], gap=e["gap"],
                               psw_MPa=EXPERIMENT_PSW["MGR27"]["psw"]),
               prior_native_dsm_calk0=NATIVE_DSM["MGR27"]["maxk"]),
    MGR23=dict(partial_this_run=rows, partial_case_record=ref,
               experiment_final=EXPERIMENT["MGR23"]),
    EPFL_P2_homog=dict(model_final_sigma_MPa=round(S[-1], 3), model_final_e=round(E[-1], 4),
                       measured_locus_sigma=mP2_s, measured_locus_e=mP2_e),
)
with open(f"{EBS}/figures/comparison_numbers.json", "w") as f:
    json.dump(out, f, indent=1)
print(json.dumps(out["MGR27"], indent=1))
print("MGR23 partial rows (this run):")
for r in rows: print(" ", {k: round(v, 4) for k, v in r.items()})
print("EPFL P2 homog final:", out["EPFL_P2_homog"]["model_final_sigma_MPa"], "MPa, e =",
      out["EPFL_P2_homog"]["model_final_e"])
