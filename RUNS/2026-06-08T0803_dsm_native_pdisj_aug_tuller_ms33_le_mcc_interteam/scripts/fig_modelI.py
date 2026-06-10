#!/usr/bin/env python3
"""Model I inter-team comparison figures — BGR pdisj_aug_tuller runs vs other EURAD-2 teams.
Reuses the column conventions from TEAM_DATA_MAP.md (header-string matching, not fixed offsets).
NO invented numbers: BGR curves come from the reduced solver-output CSVs; team curves from xlsx.
"""
import os, glob, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUN = "/Users/vinaykumar/git/ogs/ms33_pdisj_aug_tuller_2026-06-08"
ROOTA = "/Users/vinaykumar/tex/eurad2_MS34/MSXX/g_Support_Section_Data Collection"
FIG = f"{RUN}/figures"; os.makedirs(FIG, exist_ok=True)
BRANCH = "dsm_native_pdisj_aug_tuller (OGS 6.5.8)"
DD = {"1400": 1.40, "1600": 1.60, "1800": 1.80}
DDCOL = {"1400": "#d4a017", "1600": "#000000", "1800": "#1f4fd4"}  # gold / black / blue

def footnote(fig, model, out="ms33_pdisj_aug_tuller_2026-06-08/{LE,MCC}"):
    fig.text(0.5, 0.012,
             f"Simulation branch: {BRANCH}  ·  Model {model}  ·  BGR DSM (native $\\Pi$-path + vdW augmentation + Tuller macro-WRC)  ·  2026-06-08  ·  out: {out}",
             ha="center", va="bottom", fontsize=6.2, color="#555555")

# ---------- BGR data (reduced solver output) ----------
def bgr_modelI(label):
    f = f"{RUN}/reduced/{label}_history.csv"
    if not os.path.exists(f): return None
    d = pd.read_csv(f)
    return d

def saturated_endpoint(d):
    """swelling pressure = mean stress at the most-saturated row (max saturation, min |suction|)."""
    i = d["saturation"].idxmax()
    # prefer the last row at/after full saturation
    sat = d[d["saturation"] >= 0.999]
    if len(sat): i = sat.index[-1]
    return float(d.loc[i, "mean_stress_MPa"]), float(d.loc[i, "saturation"]), float(d.loc[i, "time_days"])

# ---------- Family-A team reader (Model I) ----------
def norm(s): return re.sub(r"\s+", " ", str(s).strip().lower())

def load_master_modelI(path):
    xl = pd.ExcelFile(path)
    master = next((s for s in xl.sheet_names if s.lower().startswith("model_i")), None)
    if master is None: return None
    return pd.read_excel(path, sheet_name=master, header=None)

def hdr_cols(df, row, needle):
    r = df.iloc[row].astype(str).map(norm)
    return [i for i, v in r.items() if needle in v]

def team_modelI(path):
    """Return {dd_key: {'suction':[], 'mean':[], 'perm':[], 'kind':'eff'/'tot'}} + dixon."""
    df = load_master_modelI(path)
    if df is None or len(df) < 6: return None
    HR = 4
    sucts = sorted(hdr_cols(df, HR, "suction"))
    perms = sorted(hdr_cols(df, HR, "permeability"))
    # the 3 density mean-stress columns: header has BOTH 'mean' and 'stress'
    row = df.iloc[HR].astype(str).map(norm)
    means = sorted([i for i, v in row.items() if "mean" in v and "stress" in v])
    out = {"blocks": {}, "dixon": None}
    # iterate the mean-stress columns directly (one per density) -> robust to extra
    # suction/permeability columns (BGE-CU 7-col, ULIEGE 8-col strides).
    nextmean = means[1:] + [10**9]
    for k, mc, mc_next in zip(["1400", "1600", "1800"], means[:3], nextmean[:3]):
        sc = max([s for s in sucts if s < mc], default=None)          # suction just left of mean
        pc = min([p for p in perms if mc < p < mc_next], default=None)  # perm in this block
        if sc is None: continue
        sub = df.iloc[HR+1:].copy()
        s = pd.to_numeric(sub[sc], errors="coerce")
        m = pd.to_numeric(sub[mc], errors="coerce")
        p = pd.to_numeric(sub[pc], errors="coerce") if pc is not None else pd.Series(np.nan, index=sub.index)
        ok = s.notna() & m.notna()
        kind = "tot" if "total" in norm(df.iloc[HR, mc]) else "eff"
        out["blocks"][k] = dict(suction=s[ok].to_numpy(), mean=m[ok].to_numpy(),
                                perm=p[ok].to_numpy() if pc is not None else None, kind=kind)
    # Dixon experimental block
    spc = hdr_cols(df, HR, "sweling pressure") or hdr_cols(df, HR, "swelling pressure")
    if spc:
        spc = spc[0]
        # suction col = nearest 'suction' header to the left; dd col = nearest to the right
        sc = max([c for c in hdr_cols(df, HR, "suction") if c < spc], default=None)
        ddc = min([c for c in hdr_cols(df, HR, "dry density") if c > spc], default=spc+1)
        sub = df.iloc[HR+1:]
        sp = pd.to_numeric(sub[spc], errors="coerce")
        ddv = pd.to_numeric(sub[ddc], errors="coerce")
        ok = sp.notna() & ddv.notna()
        if ok.any():
            out["dixon"] = dict(sp=sp[ok].to_numpy(), dd=ddv[ok].to_numpy())
    return out

def team_files_modelI():
    files = {}
    for d in sorted(glob.glob(f"{ROOTA}/*_DATA/Model_I")):
        team = os.path.basename(os.path.dirname(d)).replace("_DATA", "")
        cands = [f for f in glob.glob(f"{d}/*.xlsx") if "teamname" not in os.path.basename(f).lower()]
        if cands: files[team] = cands[0]
    return files

# ===================== FIGURE I-1: swelling pressure vs dry density =====================
def fig_I1():
    teams = team_files_modelI()
    fig, ax = plt.subplots(figsize=(7.8, 5.2))
    # Villar/Lloret Eq.(7)
    rr = np.linspace(1.30, 1.90, 100)
    ax.plot(rr, np.exp(6.77*rr - 9.07), "-", color="#888", lw=1.4, zorder=2,
            label="Villar/Lloret Eq.(7)")
    # team saturated swelling pressures (last row of mean-stress per block)
    first = True
    dixon_pts = {}
    for team, path in teams.items():
        try: t = team_modelI(path)
        except Exception: t = None
        if not t: continue
        xs, ys = [], []
        for k, b in t["blocks"].items():
            if len(b["mean"]):
                xs.append(DD[k]); ys.append(b["mean"][-1])  # saturated endpoint (suction->0)
        if xs:
            ax.plot(xs, ys, "-", color="#bbbbbb", lw=0.9, marker="o", ms=3,
                    zorder=3, label="other teams" if first else None)
            first = False
        if t["dixon"]:
            for dd_, sp_ in zip(t["dixon"]["dd"], t["dixon"]["sp"]):
                dixon_pts[round(float(dd_), 3)] = float(sp_)
    # Dixon experimental medians (grey squares)
    if dixon_pts:
        dx = sorted(dixon_pts); dy = [dixon_pts[x] for x in dx]
        ax.plot(dx, dy, "s", color="#444", ms=8, mfc="#cccccc", mec="#333",
                zorder=5, label="Dixon (2023) experimental")
    # BGR LE
    le = {k: bgr_modelI(f"LE_I_dd{k}") for k in DD}
    xs, ys = [], []
    for k in DD:
        if le[k] is not None:
            p, sat, td = saturated_endpoint(le[k]); xs.append(DD[k]); ys.append(p)
    ax.plot(xs, ys, "-o", color="#cc2222", lw=2.2, ms=9, zorder=6, label="BGR DSM-LE (this branch)")
    for x, y in zip(xs, ys):
        ax.annotate(f"{y:.1f}", (x, y), textcoords="offset points", xytext=(6, 6), fontsize=8, color="#cc2222")
    # BGR MCC
    mcc = {k: bgr_modelI(f"MCC_I_dd{k}") for k in ["1400", "1600"]}
    xm, ym = [], []
    for k in ["1400", "1600"]:
        if mcc[k] is not None:
            p, sat, td = saturated_endpoint(mcc[k]); xm.append(DD[k]); ym.append(p)
    ax.plot(xm, ym, "--D", color="#2a8a2a", lw=2.0, ms=8, zorder=6, label="BGR DSM-MCC (this branch)")
    ax.set_yscale("log")
    ax.set_xlabel("dry density  $\\rho_d$  [g/cm³]")
    ax.set_ylabel("swelling pressure  $P_s$  [MPa]  (mean stress at saturation)")
    ax.set_title("MS33 Model I — swelling pressure vs dry density: BGR vs EURAD-2 teams")
    ax.grid(True, which="both", ls=":", color="#ddd")
    ax.legend(fontsize=8, framealpha=0.95, loc="upper left")
    footnote(fig, "I")
    fig.subplots_adjust(bottom=0.12)
    out = f"{FIG}/modelI_swelling_pressure_vs_density.png"
    fig.savefig(out, dpi=220); plt.close(fig)
    print("wrote", out, "| BGR LE", list(zip(xs, [round(v,2) for v in ys])),
          "| BGR MCC", list(zip(xm, [round(v,2) for v in ym])), "| dixon", dixon_pts, "| teams", len(teams))

# ===================== FIGURE I-2: suction vs mean stress paths =====================
def fig_I2():
    teams = team_files_modelI()
    fig, ax = plt.subplots(figsize=(7.8, 5.2))
    first = True
    for team, path in teams.items():
        try: t = team_modelI(path)
        except Exception: t = None
        if not t: continue
        for k, b in t["blocks"].items():
            if len(b["mean"]):
                ax.plot(b["mean"], b["suction"], "-", color="#cfcfcf", lw=0.8, zorder=2,
                        label="other teams" if first else None); first = False
    # BGR LE + MCC per density
    for tag, style, lab in [("LE", dict(ls="-", lw=2.0), "LE"), ("MCC", dict(ls="--", lw=1.8), "MCC")]:
        for k in (["1400","1600","1800"] if tag=="LE" else ["1400","1600"]):
            d = bgr_modelI(f"{tag}_I_dd{k}")
            if d is None: continue
            ax.plot(d["mean_stress_MPa"], d["suction_MPa"], color=DDCOL[k], zorder=6,
                    label=f"BGR {lab} ρ_d={DD[k]}", **style)
    ax.set_yscale("log"); ax.set_ylim(0.1, 110); ax.set_xlim(0, 60)
    ax.set_xlabel("mean stress  $p$  [MPa]  (compression +)")
    ax.set_ylabel("suction  $s$  [MPa]")
    ax.set_title("MS33 Model I — suction–mean-stress path: BGR vs EURAD-2 teams")
    ax.grid(True, which="both", ls=":", color="#ddd")
    ax.legend(fontsize=7, framealpha=0.95, ncol=2, loc="upper right")
    footnote(fig, "I"); fig.subplots_adjust(bottom=0.12)
    out = f"{FIG}/modelI_suction_stress_path.png"; fig.savefig(out, dpi=220); plt.close(fig); print("wrote", out)

# ===================== FIGURE I-3: permeability vs suction =====================
def fig_I3():
    teams = team_files_modelI()
    fig, ax = plt.subplots(figsize=(7.8, 5.2))
    first = True
    for team, path in teams.items():
        try: t = team_modelI(path)
        except Exception: t = None
        if not t: continue
        for k, b in t["blocks"].items():
            if b["perm"] is not None and np.isfinite(b["perm"]).any():
                ax.plot(b["suction"], b["perm"], "-", color="#cfcfcf", lw=0.8, zorder=2,
                        label="other teams" if first else None); first = False
    for k in ["1400","1600","1800"]:
        d = bgr_modelI(f"LE_I_dd{k}")
        if d is None or not np.isfinite(d["permeability_m2"]).any(): continue
        ax.plot(d["suction_MPa"], d["permeability_m2"], "-", color=DDCOL[k], lw=2.0, zorder=6,
                label=f"BGR LE ρ_d={DD[k]}")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 110)
    ax.set_xlabel("suction  $s$  [MPa]"); ax.set_ylabel("intrinsic permeability  $k$  [m²]")
    ax.set_title("MS33 Model I — permeability vs suction: BGR vs EURAD-2 teams")
    ax.grid(True, which="both", ls=":", color="#ddd"); ax.legend(fontsize=7, framealpha=0.95)
    footnote(fig, "I"); fig.subplots_adjust(bottom=0.12)
    out = f"{FIG}/modelI_permeability_vs_suction.png"; fig.savefig(out, dpi=220); plt.close(fig); print("wrote", out)

if __name__ == "__main__":
    fig_I1(); fig_I2(); fig_I3()
