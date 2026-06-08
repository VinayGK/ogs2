#!/usr/bin/env python3
"""Extract MS33 Model II (Friction_Top) run in the TEAM-COMPARABLE quantities.

Headline comparison quantity (EURAD-2 MS33 spec, Model II):
  AXIAL (vertical) stress sigma_zz on TOP of the domain (0,0,0.070) in MPa vs time (days).
Also: porosity at top + bottom; and the full set the team xlsx tabulate
  (mean, radial, axial stress; suction; porosity; permeability) at
  Top(0,0,0.070) Central(0,0,0.040) Bottom(0,0,0.010).

Sign: OGS sigma compression negative. Teams report stress magnitudes (positive in
the xlsx). We output sigma_zz BOTH raw (OGS, negative) and as magnitude (-sigma_zz)
so the overlay can match the team convention. Suction = -pressure (MPa).

OGS axisym sigma Kelvin 4-comp order: [sigma_rr, sigma_zz, sigma_thetatheta, sqrt2*sigma_rz].
Radial stress reported by teams = sigma_rr; axial = sigma_zz; mean = (rr+zz+thth)/3.
"""
import glob, os, re, sys
import numpy as np
import meshio
import csv

RUN_DIR = os.path.dirname(os.path.abspath(__file__))
PREFIX = "wf_ms33_modelII_friction_top"
OBS = {"Top": 0.070, "Central": 0.040, "Bottom": 0.010}
DAY = 86400.0

def ts_time(fn):
    m = re.search(r"_ts_(\d+)_t_([0-9.]+)\.vtu$", fn)
    return (int(m.group(1)), float(m.group(2)))

def nearest_node(pts, zt):
    d = np.hypot(pts[:, 0] - 0.0, pts[:, 1] - zt)
    return int(d.argmin())

def main():
    files = sorted(glob.glob(os.path.join(RUN_DIR, PREFIX + "_ts_*.vtu")), key=ts_time)
    if not files:
        print("NO VTU FILES", file=sys.stderr); sys.exit(2)
    # node indices from first file
    m0 = meshio.read(files[0]); pts = m0.points
    nidx = {k: nearest_node(pts, z) for k, z in OBS.items()}
    for k, z in OBS.items():
        i = nidx[k]
        print(f"# {k} target z={z:.3f} -> node {i} at ({pts[i,0]:.4f},{pts[i,1]:.4f})", file=sys.stderr)

    rows = []
    for fn in files:
        ts, t = ts_time(fn)
        m = meshio.read(fn); pd = m.point_data
        row = {"timestep": ts, "time_s": t, "time_d": t / DAY}
        for k in OBS:
            i = nidx[k]
            sig = np.asarray(pd["sigma"][i]).ravel()
            s_rr, s_zz, s_tt = sig[0], sig[1], sig[2]
            mean = (s_rr + s_zz + s_tt) / 3.0
            p = float(np.asarray(pd["pressure"][i]).ravel()[0])
            por = float(np.asarray(pd["porosity"][i]).ravel()[0])
            sat = float(np.asarray(pd["saturation"][i]).ravel()[0])
            perm = float(np.asarray(pd["intrinsic_permeability"][i]).ravel()[0])
            nl = float(np.asarray(pd["micro_water_content"][i]).ravel()[0]) if "micro_water_content" in pd else float("nan")
            row[f"{k}_sigma_zz_MPa"] = s_zz / 1e6          # raw OGS (negative=compression)
            row[f"{k}_axial_mag_MPa"] = -s_zz / 1e6        # team magnitude convention
            row[f"{k}_radial_mag_MPa"] = -s_rr / 1e6
            row[f"{k}_mean_mag_MPa"] = -mean / 1e6
            row[f"{k}_suction_MPa"] = -p / 1e6
            row[f"{k}_porosity"] = por
            row[f"{k}_saturation"] = sat
            row[f"{k}_perm_m2"] = perm
            row[f"{k}_micro_n_l"] = nl
        rows.append(row)

    # write full CSV
    out_csv = os.path.join(RUN_DIR, "wf_modelII_extract.csv")
    cols = list(rows[0].keys())
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

    # headline CSV: time_d, axial stress on top (magnitude, MPa) + porosity top/bottom
    head_csv = os.path.join(RUN_DIR, "wf_modelII_axial_top_vs_time.csv")
    with open(head_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_d", "axial_stress_top_MPa", "porosity_top", "porosity_bottom"])
        for r in rows:
            w.writerow([f"{r['time_d']:.4f}", f"{r['Top_axial_mag_MPa']:.6f}",
                        f"{r['Top_porosity']:.6f}", f"{r['Bottom_porosity']:.6f}"])

    last = rows[-1]
    print("\n=== SUMMARY ===")
    print(f"VTUs: {len(files)}  last ts={last['timestep']} t={last['time_d']:.3f} d")
    print(f"full CSV : {out_csv}")
    print(f"headline : {head_csv}")
    print("\n--- final values @ t={:.2f} d ---".format(last["time_d"]))
    for k in OBS:
        print(f"  {k:7s}: axial(mag)={last[k+'_axial_mag_MPa']:.4f} MPa  radial(mag)={last[k+'_radial_mag_MPa']:.4f}  mean(mag)={last[k+'_mean_mag_MPa']:.4f}  suction={last[k+'_suction_MPa']:.4f} MPa  poro={last[k+'_porosity']:.4f}  sat={last[k+'_saturation']:.4f}  n_l={last[k+'_micro_n_l']:.4e}")
    print("\n--- HEADLINE: axial stress on TOP (magnitude) vs time ---")
    for r in rows:
        if abs(r['time_d'] % 20.0) < 0.6 or r is rows[-1]:
            print(f"  t={r['time_d']:7.2f} d  sigma_zz_top(mag)={r['Top_axial_mag_MPa']:8.4f} MPa  (raw {r['Top_sigma_zz_MPa']:8.4f})")

if __name__ == "__main__":
    main()
