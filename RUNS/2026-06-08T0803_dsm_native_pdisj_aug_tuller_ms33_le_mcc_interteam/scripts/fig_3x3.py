#!/usr/bin/env python3
"""MS33 BGR DSM (pdisj_aug_tuller) vs EURAD-2 teams — single 3x3 panel.

Reuses the authoritative loaders in fig_modelI.py + fig_III_IV_VII.py
(team time-series from the Family-A master sheets; BGR from reduced/*.csv of
THIS submission branch). No new data path, no fabrication.

Layout (user spec):
  row1: (1) swelling pressure vs dry density   (2) suction vs mean stress   (3) Model III mean stress vs t
  row2: (4) Model IV mean stress vs t          (5) Model III gap closure    (6) Model IV dry density
  row3: (7) Model VII loading cycle e-sigma_a  (8) VII void ratio vs t       (9) VII axial stress vs t
  ( (4): user listed 'model 3 mean stress' twice; rendered as Model IV mean stress to
    avoid a duplicate of (3) and to pair with the Model IV DD panel. Easy to flip. )
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUN = "/Users/vinaykumar/git/ogs/ms33_pdisj_aug_tuller_2026-06-08"
sys.path.insert(0, RUN)
import fig_modelI as M1
import fig_III_IV_VII as M3

GREY = "#cfcfcf"; RED = "#cc2222"; GRN = "#2a8a2a"
DD = M1.DD; DDCOL = M1.DDCOL

fig, AX = plt.subplots(3, 3, figsize=(17, 14.5))
fig.suptitle("MS33 — BGR DSM (native $\\Pi$-path + vdW augmentation + Tuller macro-WRC) vs EURAD-2 teams",
             fontsize=15, fontweight="bold", y=0.992)

# ===== (1) swelling pressure vs dry density =====
ax = AX[0, 0]
teamsI = M1.team_files_modelI()
rr = np.linspace(1.30, 1.90, 100)
ax.plot(rr, np.exp(6.77*rr - 9.07), "-", color="#888", lw=1.3, label="Villar/Lloret Eq.(7)")
first = True; dixon = {}
for team, path in teamsI.items():
    try: t = M1.team_modelI(path)
    except Exception: t = None
    if not t: continue
    xs, ys = [], []
    for k, b in t["blocks"].items():
        if len(b["mean"]): xs.append(DD[k]); ys.append(b["mean"][-1])
    if xs:
        ax.plot(xs, ys, "-o", color="#bbb", lw=0.8, ms=3, label="other teams" if first else None); first = False
    if t["dixon"]:
        for dd_, sp_ in zip(t["dixon"]["dd"], t["dixon"]["sp"]): dixon[round(float(dd_),3)] = float(sp_)
if dixon:
    dx = sorted(dixon); ax.plot(dx, [dixon[x] for x in dx], "s", color="#444", ms=8,
                                mfc="#ccc", mec="#333", label="Dixon (2023) exp.")
le = {k: M1.bgr_modelI(f"LE_I_dd{k}") for k in DD}
xs = [DD[k] for k in DD if le[k] is not None]
ys = [M1.saturated_endpoint(le[k])[0] for k in DD if le[k] is not None]
ax.plot(xs, ys, "-o", color=RED, lw=2.2, ms=9, label="BGR DSM-LE")
for x, y in zip(xs, ys): ax.annotate(f"{y:.1f}", (x, y), textcoords="offset points", xytext=(6,6), fontsize=8, color=RED)
mcc = {k: M1.bgr_modelI(f"MCC_I_dd{k}") for k in ["1400","1600"]}
xm = [DD[k] for k in ["1400","1600"] if mcc[k] is not None]
ym = [M1.saturated_endpoint(mcc[k])[0] for k in ["1400","1600"] if mcc[k] is not None]
ax.plot(xm, ym, "--D", color=GRN, lw=2.0, ms=8, label="BGR DSM-MCC")
ax.set_yscale("log"); ax.set_xlabel("dry density $\\rho_d$ [g/cm³]")
ax.set_ylabel("swelling pressure $P_s$ [MPa]")
ax.set_title("(1) Model I — swelling pressure vs $\\rho_d$", fontsize=11, fontweight="bold")
ax.grid(True, which="both", ls=":", color="#ddd"); ax.legend(fontsize=7, loc="upper left")

# ===== (2) suction vs mean stress path =====
ax = AX[0, 1]
first = True
for team, path in teamsI.items():
    try: t = M1.team_modelI(path)
    except Exception: t = None
    if not t: continue
    for k, b in t["blocks"].items():
        if len(b["mean"]):
            ax.plot(b["mean"], b["suction"], "-", color=GREY, lw=0.8, label="other teams" if first else None); first = False
for tag, style, lab in [("LE", dict(ls="-", lw=2.0), "LE"), ("MCC", dict(ls="--", lw=1.8), "MCC")]:
    for k in (["1400","1600","1800"] if tag=="LE" else ["1400","1600"]):
        d = M1.bgr_modelI(f"{tag}_I_dd{k}")
        if d is None: continue
        ax.plot(d["mean_stress_MPa"], d["suction_MPa"], color=DDCOL[k], label=f"BGR {lab} $\\rho_d$={DD[k]}", **style)
ax.set_yscale("log"); ax.set_ylim(0.1, 110); ax.set_xlim(0, 60)
ax.set_xlabel("mean stress $p$ [MPa]"); ax.set_ylabel("suction $s$ [MPa]")
ax.set_title("(2) Model I — suction–mean-stress path", fontsize=11, fontweight="bold")
ax.grid(True, which="both", ls=":", color="#ddd"); ax.legend(fontsize=6.5, ncol=2, loc="upper right")

# ---- helper for vs-time team overlays ----
def overlay_team_time(ax, model, needle, loc):
    first = True
    for team, path in M3.team_files(model).items():
        df = M3.load_master(path, model); r = M3.central_series(df, needle, loc=loc)
        if r:
            ax.plot(r[0], r[1], "-", color=GREY, lw=0.9, label="other teams" if first else None); first = False

# ===== (3) Model III mean stress vs time =====
ax = AX[0, 2]
overlay_team_time(ax, "III", "mean stress", "central")
d = M3.bgr_csv("LE_III");  ax.plot(d["time_days"], d["mean_stress_MPa"], "-",  color=RED, lw=2.2, label="BGR LE")
d = M3.bgr_csv("MCC_III"); ax.plot(d["time_days"], d["mean_stress_MPa"], "--", color=GRN, lw=2.0, label="BGR MCC")
ax.set_xlim(0, 200); ax.set_xlabel("time [d]"); ax.set_ylabel("mean stress $p$ [MPa]")
ax.set_title("(3) Model III — mean stress vs time (centre)", fontsize=11, fontweight="bold")
ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)

# ===== (4) Model IV mean stress vs time =====
ax = AX[1, 0]
overlay_team_time(ax, "IV", "mean stress", "central")
d = M3.bgr_csv("LE_IV"); ax.plot(d["time_days"], d["mean_stress_MPa"], "-", color=RED, lw=2.2, label="BGR LE (clay centre)")
ax.set_xlim(0, 200); ax.set_xlabel("time [d]"); ax.set_ylabel("mean stress $p$ [MPa]")
ax.set_title("(4) Model IV — mean stress vs time (centre)", fontsize=11, fontweight="bold")
ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)

# ===== (5) Model III gap closure vs time =====
ax = AX[1, 1]
overlay_team_time(ax, "III", "gap closure", "")
for lab, outdir, pre, col, ls in [("BGR LE", f"{RUN}/LE/ModelIII/out","ms33_modelIII_gap2mm",RED,"-"),
                                   ("BGR MCC", f"{RUN}/MCC/ModelIII/out","ms33_modelIII_gap2mm_mcc_native",GRN,"--")]:
    if os.path.isdir(outdir):
        T, Y = M3.bgr_probe_series(outdir, pre, "displacement", (0.023, 0.040), comp=0)
        if len(T): ax.plot(T, np.abs(Y)*1e3, ls, color=col, lw=2.0, label=f"{lab} $|u_r|$ @gap")
ax.set_xlim(0, 200); ax.set_xlabel("time [d]"); ax.set_ylabel("gap closure / $|u_r|$ [mm]")
ax.set_title("(5) Model III — gap closure vs time (r=23 mm)", fontsize=11, fontweight="bold")
ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)

# ===== (6) Model IV dry density vs time =====
ax = AX[1, 2]
overlay_team_time(ax, "IV", "dry density", "")
for zlab, rz, col, ls in [("clay (top, z=52mm)",(0.0,0.0525),RED,"-"), ("pellet (bottom, z=18mm)",(0.0,0.0175),"#cc8822","--")]:
    outdir = f"{RUN}/LE/ModelIV/out_perK"
    if os.path.isdir(outdir):
        T, Y = M3.bgr_probe_series(outdir, "ms33_modelIV_pellets", "dry_density_solid", rz)
        if len(T): ax.plot(T, Y/1000.0, ls, color=col, lw=2.0, label=f"BGR LE {zlab}")
ax.set_xlim(0, 200); ax.set_xlabel("time [d]"); ax.set_ylabel("dry density $\\rho_d$ [g/cm³]")
ax.set_title("(6) Model IV — dry-density evolution (clay vs pellet)", fontsize=11, fontweight="bold")
ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)

# ===== (7) Model VII loading cycle: e vs axial stress =====
ax = AX[2, 0]
first = True
for team, path in M3.team_files("VII").items():
    df = M3.load_master(path, "VII")
    re_ = M3.central_series(df, "void ratio", loc="")
    ra_ = M3.central_series(df, "axial stress", loc="central")
    if re_ and ra_:
        # align axial onto void-ratio time base
        a_on_e = np.interp(re_[0], ra_[0], ra_[1])
        ax.plot(a_on_e, re_[1], "-", color=GREY, lw=0.9, label="other teams" if first else None); first = False
d = M3.bgr_csv("LE_VII")
ax.plot(d["axial_stress_MPa"], d["void_ratio"], "-", color=RED, lw=2.2, label="BGR LE (centre)")
ax.set_xlabel("axial stress $\\sigma_a$ [MPa]"); ax.set_ylabel("void ratio $e$ [–]")
ax.set_title("(7) Model VII — loading cycle ($e$–$\\sigma_a$)", fontsize=11, fontweight="bold")
ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)

# ===== (8) Model VII void ratio vs time =====
ax = AX[2, 1]
overlay_team_time(ax, "VII", "void ratio", "")
d = M3.bgr_csv("LE_VII"); ax.plot(d["time_days"], d["void_ratio"], "-", color=RED, lw=2.2, label="BGR LE (centre)")
ax.set_xlim(0, 240); ax.set_xlabel("time [d]"); ax.set_ylabel("void ratio $e$ [–]")
ax.set_title("(8) Model VII — void ratio vs time", fontsize=11, fontweight="bold")
ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)

# ===== (9) Model VII axial stress vs time =====
ax = AX[2, 2]
overlay_team_time(ax, "VII", "axial stress", "central")
d = M3.bgr_csv("LE_VII"); ax.plot(d["time_days"], d["axial_stress_MPa"], "-", color=RED, lw=2.2, label="BGR LE (centre)")
ax.set_xlim(0, 240); ax.set_xlabel("time [d]"); ax.set_ylabel("axial stress $\\sigma_a$ [MPa]")
ax.set_title("(9) Model VII — axial stress vs time", fontsize=11, fontweight="bold")
ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)

fig.text(0.5, 0.004, f"Simulation branch: {M3.BRANCH}  ·  BGR reduced/*.csv (this submission) + EURAD-2 Family-A team sheets  ·  2026-06-08",
         ha="center", fontsize=7, color="#555")
fig.tight_layout(rect=[0, 0.012, 1, 0.978])
out = f"{RUN}/figures/ms33_3x3_interteam.png"
fig.savefig(out, dpi=200); print("wrote", out)
