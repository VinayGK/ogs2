#!/usr/bin/env python3
"""Reduce MS33 pdisj_aug_tuller run VTU series to center-point history CSVs.

Conventions (verbatim from Vinay's plot_ms33_standard / plot_interteam_panels):
  suction      s = -pressure / 1e6                         [MPa]
  mean stress  p = -(sxx+syy+szz)/3 / 1e6   (compression+) [MPa]
  radial       = -sigma_xx/1e6,  axial = -sigma_yy/1e6     [MPa]
  swelling     = -(ss_xx+ss_yy+ss_zz)/3 / 1e6              [MPa]
  void ratio   e = phi/(1-phi)
Probe = mesh node nearest specimen centre (r=0, z=H/2). Model I (4 nodes) = mean.
NO calibration, NO invented values — pure reduction of solver output.
"""
import glob, re, os, sys, csv
import numpy as np
import meshio

RUN = "/Users/vinaykumar/git/ogs/ms33_pdisj_aug_tuller_2026-06-08"

# (label, out_dir, prefix)
RUNS = [
    ("LE_I_dd1400",  f"{RUN}/LE/ModelI/dd1400_st_out",  "ms33_modelI_dd1400"),
    ("LE_I_dd1600",  f"{RUN}/LE/ModelI/dd1600_st_out",  "ms33_modelI_dd1600"),   # partial (fails ~22d)
    ("LE_I_dd1800",  f"{RUN}/LE/ModelI/dd1800_st_out",  "ms33_modelI_dd1800"),
    ("LE_III",       f"{RUN}/LE/ModelIII/out",          "ms33_modelIII_gap2mm"),
    ("LE_IV",        f"{RUN}/LE/ModelIV/out_perK",       "ms33_modelIV_pellets"),  # per-material pellet K (13064)
    ("LE_VII",       f"{RUN}/LE/ModelVII/out",          "ms33_modelVII_freeswelling"),
    ("MCC_I_dd1400", f"{RUN}/MCC/ModelI/dd1400_out",    "ms33_modelI_dd1400_mcc_native"),
    ("MCC_I_dd1600", f"{RUN}/MCC/ModelI/dd1600_out",    "ms33_modelI_dd1600_mcc_native"),
    ("MCC_III",      f"{RUN}/MCC/ModelIII/out",         "ms33_modelIII_gap2mm_mcc_native"),
]

def ts_of(fn):
    m = re.search(r"_ts_(\d+)_t_([0-9.]+)\.vtu$", fn)
    return (int(m.group(1)), float(m.group(2))) if m else (None, None)

def probe_idx(pts):
    zc = 0.5 * (pts[:,1].min() + pts[:,1].max())
    d = (pts[:,0]-0.0)**2 + (pts[:,1]-zc)**2
    return int(np.argmin(d))

def scal(pd, name, idx, single_elem):
    if name not in pd: return np.nan
    a = pd[name]
    if a.ndim == 1:
        return float(a.mean() if single_elem else a[idx])
    return float(a[:,0].mean() if single_elem else a[idx,0])

def vec(pd, name, idx, single_elem):
    if name not in pd: return None
    a = pd[name]
    return a.mean(axis=0) if single_elem else a[idx]

os.makedirs(f"{RUN}/reduced", exist_ok=True)
manifest = []
for label, d, pre in RUNS:
    vtus = sorted(glob.glob(f"{d}/{pre}_ts_*.vtu"), key=lambda f: ts_of(f)[0])
    if not vtus:
        manifest.append((label, 0, np.nan, "NO_OUTPUT")); continue
    rows = []
    for v in vtus:
        _, t = ts_of(v)
        m = meshio.read(v); pd = m.point_data; pts = m.points
        se = (len(pts) <= 4)
        idx = probe_idx(pts)
        sig = vec(pd, "sigma", idx, se)
        ss  = vec(pd, "swelling_stress", idx, se)
        phi = scal(pd, "porosity", idx, se)
        rows.append(dict(
            time_days = t/86400.0,
            suction_MPa = -scal(pd,"pressure",idx,se)/1e6,
            mean_stress_MPa = (-(sig[0]+sig[1]+sig[2])/3.0/1e6) if sig is not None else np.nan,
            radial_stress_MPa = (-sig[0]/1e6) if sig is not None else np.nan,
            axial_stress_MPa  = (-sig[1]/1e6) if sig is not None else np.nan,
            swelling_stress_MPa = (-(ss[0]+ss[1]+ss[2])/3.0/1e6) if ss is not None else np.nan,
            saturation = scal(pd,"saturation",idx,se),
            porosity = phi,
            void_ratio = (phi/(1.0-phi)) if phi==phi else np.nan,
            permeability_m2 = scal(pd,"intrinsic_permeability",idx,se),
            dry_density = scal(pd,"dry_density_solid",idx,se),
        ))
    out = f"{RUN}/reduced/{label}_history.csv"
    with open(out,"w",newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    last = rows[-1]
    manifest.append((label, len(rows), last["time_days"], f"p={last['mean_stress_MPa']:.3f} s={last['suction_MPa']:.3f} Sl={last['saturation']:.3f}"))

print(f"{'run':14} {'nframes':>7} {'t_end_d':>9}  final")
for label, n, te, info in manifest:
    print(f"{label:14} {n:7d} {te:9.2f}  {info}")
