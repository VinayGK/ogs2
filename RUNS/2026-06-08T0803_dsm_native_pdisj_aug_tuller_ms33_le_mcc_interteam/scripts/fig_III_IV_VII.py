#!/usr/bin/env python3
"""Models III/IV/VII inter-team time-series comparison — BGR pdisj_aug_tuller vs EURAD-2 teams.
Team time-series read from Family-A master sheets by header-string matching (per TEAM_DATA_MAP.md).
BGR series from reduced center-point CSVs (verified) + extra VTU probes (pellet zone, gap location).
NO invented numbers.
"""
import os, glob, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meshio

RUN = "/Users/vinaykumar/git/ogs/ms33_pdisj_aug_tuller_2026-06-08"
ROOTA = "/Users/vinaykumar/tex/eurad2_MS34/MSXX/g_Support_Section_Data Collection"
FIG = f"{RUN}/figures"; os.makedirs(FIG, exist_ok=True)
BRANCH = "dsm_native_pdisj_aug_tuller (OGS 6.5.8)"

def footnote(fig, model, out="ms33_pdisj_aug_tuller_2026-06-08/{LE,MCC}"):
    fig.text(0.5, 0.012,
             f"Simulation branch: {BRANCH}  ·  Model {model}  ·  BGR DSM (native $\\Pi$-path + vdW augmentation + Tuller macro-WRC)  ·  2026-06-08  ·  out: {out}",
             ha="center", va="bottom", fontsize=6.2, color="#555555")
def norm(s): return re.sub(r"\s+", " ", str(s).strip().lower())

# ---------- Family-A time-series reader ----------
def load_master(path, model):
    xl = pd.ExcelFile(path)
    key = {"III":"model_iii","IV":"model_iv","VII":"model_vii"}[model]
    master = next((s for s in xl.sheet_names if norm(s).startswith(key)), None)
    if master is None: return None
    return pd.read_excel(path, sheet_name=master, header=None)

def _hdr(df, c):
    """concatenated header text (rows 1..6) for column c, normalized."""
    return " | ".join(norm(df.iat[r, c]) for r in range(1, 7) if r < len(df))

def central_series(df, needle, hr=None, loc="central"):
    """Return (time_days, value) for a quantity group.
    Searches header rows 1..6 for `needle`; for main groups prefers the
    row-3 'central' column; for single-location appended groups takes the match.
    Time column = nearest left column whose header says time/days. Data = dropna."""
    if df is None or len(df) < 8: return None
    n = df.shape[1]
    r3 = [norm(df.iat[3, c]) for c in range(n)]
    def is_time(c):
        h = _hdr(df, c); return ("time" in h) or ("(days)" in h)
    cands = [c for c in range(n) if needle in _hdr(df, c) and not is_time(c)]
    if not cands: return None
    pick = None
    if loc == "central":
        cc = [c for c in cands if r3[c] == "central"]
        if cc: pick = cc[0]
    if pick is None: pick = cands[0]
    tcol = max([c for c in range(pick) if is_time(c)], default=None)
    if tcol is None: return None
    t = pd.to_numeric(df[tcol], errors="coerce")
    y = pd.to_numeric(df[pick], errors="coerce")
    ok = t.notna() & y.notna() & (t >= 0) & (t < 1e4)
    if ok.sum() < 2: return None
    return t[ok].to_numpy(), y[ok].to_numpy()

def team_files(model):
    sub = {"III":"Model_III","IV":"Model_IV","VII":"Model_VII"}[model]
    files = {}
    for d in sorted(glob.glob(f"{ROOTA}/*_DATA/{sub}")):
        team = os.path.basename(os.path.dirname(d)).replace("_DATA","")
        cands = [f for f in glob.glob(f"{d}/*.xlsx") if "teamname" not in os.path.basename(f).lower()]
        if cands: files[team] = cands[0]
    return files

# ---------- BGR ----------
def bgr_csv(label):
    f = f"{RUN}/reduced/{label}_history.csv"
    return pd.read_csv(f) if os.path.exists(f) else None

def ts_of(fn):
    m = re.search(r"_ts_(\d+)_t_([0-9.]+)\.vtu$", fn); return (int(m.group(1)), float(m.group(2)))

def bgr_probe_series(outdir, prefix, field, rz, comp=None, reduce_mean=False):
    """Time series of a field at the mesh node nearest (r,z)=rz, from a VTU series."""
    vtus = sorted(glob.glob(f"{outdir}/{prefix}_ts_*.vtu"), key=lambda f: ts_of(f)[0])
    T, Y = [], []
    idx = None
    for v in vtus:
        _, t = ts_of(v); m = meshio.read(v); pts = m.points
        if idx is None:
            idx = int(np.argmin((pts[:,0]-rz[0])**2 + (pts[:,1]-rz[1])**2))
        if field not in m.point_data: continue
        a = m.point_data[field]
        val = a[idx] if a.ndim==1 else a[idx, comp if comp is not None else 0]
        T.append(t/86400.0); Y.append(float(val))
    return np.array(T), np.array(Y)

GREY = "#cfcfcf"

# ============== MODEL III: stresses + gap closure ==============
def fig_III():
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.0))
    # (a) mean stress vs time at center
    ax = axs[0]; first=True
    for team, path in team_files("III").items():
        df = load_master(path, "III"); r = central_series(df, "mean stress")
        if r: ax.plot(r[0], r[1], "-", color=GREY, lw=0.9, label="other teams" if first else None); first=False
    d = bgr_csv("LE_III");  ax.plot(d["time_days"], d["mean_stress_MPa"], "-", color="#cc2222", lw=2.2, label="BGR LE (center)")
    d = bgr_csv("MCC_III"); ax.plot(d["time_days"], d["mean_stress_MPa"], "--", color="#2a8a2a", lw=2.0, label="BGR MCC (center)")
    ax.set_xlabel("time [d]"); ax.set_ylabel("mean stress $p$ [MPa]"); ax.set_xlim(0,200)
    ax.set_title("Model III — mean stress vs time (centre)"); ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)
    # (b) gap closure vs time: teams report 'gap closure'; BGR proxy = radial displacement u_r at gap (r=0.023,z=0.040)
    ax = axs[1]; first=True
    for team, path in team_files("III").items():
        df = load_master(path, "III"); r = central_series(df, "gap closure", loc="")
        if r: ax.plot(r[0], r[1], "-", color=GREY, lw=0.9, label="other teams" if first else None); first=False
    for lab, outdir, pre, col in [("BGR LE", f"{RUN}/LE/ModelIII/out","ms33_modelIII_gap2mm","#cc2222"),
                                   ("BGR MCC", f"{RUN}/MCC/ModelIII/out","ms33_modelIII_gap2mm_mcc_native","#2a8a2a")]:
        T,Y = bgr_probe_series(outdir, pre, "displacement", (0.023,0.040), comp=0)
        if len(T): ax.plot(T, np.abs(Y)*1e3, "-" if "LE" in lab else "--", color=col, lw=2.0, label=f"{lab} |u_r| @gap")
    ax.set_xlabel("time [d]"); ax.set_ylabel("gap closure / $|u_r|$ [mm]"); ax.set_xlim(0,200)
    ax.set_title("Model III — gap closure vs time (gap interface r=23 mm)"); ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)
    footnote(fig, "III"); fig.subplots_adjust(bottom=0.13, wspace=0.25)
    out=f"{FIG}/modelIII_interteam.png"; fig.savefig(out, dpi=220); plt.close(fig); print("wrote", out)

# ============== MODEL IV: mean stress + dry density (two zones) ==============
def fig_IV():
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.0))
    ax = axs[0]; first=True
    for team, path in team_files("IV").items():
        df = load_master(path, "IV"); r = central_series(df, "mean stress")
        if r: ax.plot(r[0], r[1], "-", color=GREY, lw=0.9, label="other teams" if first else None); first=False
    d = bgr_csv("LE_IV"); ax.plot(d["time_days"], d["mean_stress_MPa"], "-", color="#cc2222", lw=2.2, label="BGR LE (clay centre)")
    ax.set_xlabel("time [d]"); ax.set_ylabel("mean stress $p$ [MPa]"); ax.set_xlim(0,200)
    ax.set_title("Model IV — mean stress vs time (centre)"); ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)
    # dry density two zones: clay (r=0,z=0.035) + pellet (r=0.020,z=0.035)
    ax = axs[1]; first=True
    for team, path in team_files("IV").items():
        df = load_master(path, "IV"); r = central_series(df, "dry density", hr=4, loc="")
        if r: ax.plot(r[0], r[1], "-", color=GREY, lw=0.9, label="other teams" if first else None); first=False
    # vertical stack: clay (MaterialID 0) on top z[0.036,0.069], pellet (id 1) bottom z[0.001,0.034]
    for zlab, rz, col, ls in [("clay (top, z=52mm)",(0.0,0.0525),"#cc2222","-"), ("pellet (bottom, z=18mm)",(0.0,0.0175),"#cc8822","--")]:
        T,Y = bgr_probe_series(f"{RUN}/LE/ModelIV/out_perK","ms33_modelIV_pellets","dry_density_solid", rz)
        if len(T): ax.plot(T, Y/1000.0, ls, color=col, lw=2.0, label=f"BGR LE {zlab}")
    ax.set_xlabel("time [d]"); ax.set_ylabel("dry density $\\rho_d$ [g/cm³]"); ax.set_xlim(0,200)
    ax.set_title("Model IV — dry-density evolution (clay vs pellet)"); ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)
    footnote(fig, "IV"); fig.subplots_adjust(bottom=0.13, wspace=0.25)
    out=f"{FIG}/modelIV_interteam.png"; fig.savefig(out, dpi=220); plt.close(fig); print("wrote", out)

# ============== MODEL VII: void ratio + axial stress ==============
def fig_VII():
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.0))
    ax = axs[0]; first=True
    for team, path in team_files("VII").items():
        df = load_master(path, "VII"); r = central_series(df, "void ratio", loc="")
        if r: ax.plot(r[0], r[1], "-", color=GREY, lw=0.9, label="other teams" if first else None); first=False
    d = bgr_csv("LE_VII"); ax.plot(d["time_days"], d["void_ratio"], "-", color="#cc2222", lw=2.2, label="BGR LE (centre)")
    ax.set_xlabel("time [d]"); ax.set_ylabel("void ratio $e$ [–]"); ax.set_xlim(0,240)
    ax.set_title("Model VII — void ratio vs time (free swelling)"); ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)
    ax = axs[1]; first=True
    for team, path in team_files("VII").items():
        df = load_master(path, "VII"); r = central_series(df, "axial stress")
        if r: ax.plot(r[0], r[1], "-", color=GREY, lw=0.9, label="other teams" if first else None); first=False
    d = bgr_csv("LE_VII"); ax.plot(d["time_days"], d["axial_stress_MPa"], "-", color="#cc2222", lw=2.2, label="BGR LE (centre)")
    ax.set_xlabel("time [d]"); ax.set_ylabel("axial stress $\\sigma_a$ [MPa]"); ax.set_xlim(0,240)
    ax.set_title("Model VII — axial stress vs time"); ax.grid(True, ls=":", color="#ddd"); ax.legend(fontsize=7)
    footnote(fig, "VII"); fig.subplots_adjust(bottom=0.13, wspace=0.25)
    out=f"{FIG}/modelVII_interteam.png"; fig.savefig(out, dpi=220); plt.close(fig); print("wrote", out)

if __name__ == "__main__":
    # quick reader self-test
    import sys
    if "--test" in sys.argv:
        for model in ["III","IV","VII"]:
            for team, path in team_files(model).items():
                df = load_master(path, model)
                for q in (["mean stress","gap closure"] if model=="III" else
                          ["mean stress","dry density"] if model=="IV" else ["void ratio","axial stress"]):
                    hr = 4 if q=="dry density" else 6
                    r = central_series(df, q, hr=hr, loc="" if q in ("gap closure","dry density","void ratio") else "central")
                    print(f"{model:4} {team:18} {q:14} {'OK n=%d'%len(r[0]) if r else 'MISS'}")
        sys.exit(0)
    fig_III(); fig_IV(); fig_VII()
