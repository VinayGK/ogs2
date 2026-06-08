#!/usr/bin/env python3
"""Extract MS33 Model II (Friction_Top) LE run in the team-comparable quantities.

Headline comparison (team xlsx Model_II chart-sheet "Axial stress (comparison)"):
  AXIAL stress sigma_zz at Top(0,0,0.070), Central(0,0,0.040), Bottom(0,0,0.010), MPa vs day.
Also reported: mean, radial stress; suction; porosity; permeability at the 3 points.
Sign: OGS sigma compression negative; teams tabulate magnitudes (positive). Output both.
OGS axisym sigma Kelvin order: [rr, zz, thth, sqrt2*rz]. mean=(rr+zz+thth)/3.
"""
import glob, os, re, sys, csv
import numpy as np
import meshio

RUN_DIR = os.path.dirname(os.path.abspath(__file__))
PREFIX = "wf2_modelII_LE"
OBS = {"Top": 0.070, "Central": 0.040, "Bottom": 0.010}
DAY = 86400.0

def ts_time(fn):
    m = re.search(r"_ts_(\d+)_t_([0-9.]+)\.vtu$", fn)
    return (int(m.group(1)), float(m.group(2)))

def nearest_node(pts, zt):
    d = np.hypot(pts[:, 0] - 0.0, pts[:, 1] - zt)
    return int(d.argmin())

files = sorted(glob.glob(os.path.join(RUN_DIR, PREFIX + "_ts_*.vtu")), key=ts_time)
if not files:
    print("NO VTU FILES", file=sys.stderr); sys.exit(2)
m0 = meshio.read(files[0]); pts = m0.points
nidx = {k: nearest_node(pts, z) for k, z in OBS.items()}
for k, z in OBS.items():
    i = nidx[k]
    print(f"# {k} z={z:.3f} -> node {i} ({pts[i,0]:.4f},{pts[i,1]:.4f})", file=sys.stderr)

rows = []
for fn in files:
    ts, t = ts_time(fn)
    pd = meshio.read(fn).point_data
    row = {"timestep": ts, "time_d": round(t / DAY, 4)}
    for k in OBS:
        i = nidx[k]
        sig = np.asarray(pd["sigma"][i]).ravel()
        s_rr, s_zz, s_tt = sig[0], sig[1], sig[2]
        mean = (s_rr + s_zz + s_tt) / 3.0
        p = float(np.asarray(pd["pressure"][i]).ravel()[0])
        por = float(np.asarray(pd["porosity"][i]).ravel()[0])
        sat = float(np.asarray(pd["saturation"][i]).ravel()[0])
        perm = float(np.asarray(pd["intrinsic_permeability"][i]).ravel()[0])
        row[f"{k}_axial_mag_MPa"] = -s_zz / 1e6
        row[f"{k}_radial_mag_MPa"] = -s_rr / 1e6
        row[f"{k}_mean_mag_MPa"] = -mean / 1e6
        row[f"{k}_suction_MPa"] = -p / 1e6
        row[f"{k}_porosity"] = por
        row[f"{k}_saturation"] = sat
        row[f"{k}_perm_m2"] = perm
    rows.append(row)

out_csv = os.path.join(RUN_DIR, "wf2_modelII_LE_extract.csv")
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

last = rows[-1]
print(f"\nLE final state at day {last['time_d']:.1f}  (axial/mean/radial = team magnitude, MPa):")
print(f"{'point':8} {'axial':>8} {'mean':>8} {'radial':>8} {'suction':>9} {'porosity':>9}")
for k in OBS:
    print(f"{k:8} {last[k+'_axial_mag_MPa']:8.3f} {last[k+'_mean_mag_MPa']:8.3f} "
          f"{last[k+'_radial_mag_MPa']:8.3f} {last[k+'_suction_MPa']:9.4f} {last[k+'_porosity']:9.4f}")
print(f"\nFull CSV: {out_csv}")
