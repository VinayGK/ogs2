#!/usr/bin/env python3
"""
BEACON Task 3.3 — extract void ratio e and axial total stress sigma_zz vs time
from the native-DSM oedometer runs, for comparison against the measured
e-vs-sigma_axial curves (Acta Geotech 2022, Ferrari et al., Figs 9/11).

e = phi/(1-phi)            (phi = macro+micro total porosity field 'porosity')
sigma_zz (total)           OGS 'sigma' component 1 (Pa, compression negative)
sigma_axial_total = -sigma_zz/1e6   (MPa, soil-mechanics positive-compression)
back-pressure 20 kPa is constant; effective ~ total here for plotting.

Single-element mesh -> all integration points identical; take cell mean.
Writes one CSV per path to results/.
NO expected values asserted here (pure extraction); comparison/plot is separate.
"""
import glob, os, re, csv, sys
import numpy as np
import meshio

ROOT = os.path.dirname(os.path.abspath(__file__))

def cell_mean(m, name):
    if name in m.cell_data:
        return float(np.concatenate([np.asarray(b).ravel() for b in m.cell_data[name]]).mean())
    if name in m.point_data:
        return float(np.asarray(m.point_data[name]).mean())
    return None

def sigma_zz(m):
    # 'sigma' is a tensor: components [xx, yy, zz, xy] (2D axisym -> yy is axial z)
    # OGS RM 2D axisymmetric: comp order (rr, zz, thetatheta, rz). Axial = component 1 (zz).
    for key in ("sigma",):
        if key in m.cell_data:
            arr = np.asarray(m.cell_data[key][0])
            arr = arr.reshape(arr.shape[0], -1).mean(axis=0)
            return arr
        if key in m.point_data:
            arr = np.asarray(m.point_data[key]).mean(axis=0)
            return arr
    return None

def extract(path_dir, out_csv):
    fs = glob.glob(os.path.join(ROOT, "results", path_dir, "*.vtu"))
    rows = []
    for f in fs:
        mobj = re.search(r"_t_([0-9.]+)\.vtu$", os.path.basename(f))
        if not mobj:
            continue
        t = float(mobj.group(1))
        m = meshio.read(f)
        phi = cell_mean(m, "porosity")
        e = phi/(1.0-phi) if phi is not None else None
        s = sigma_zz(m)
        szz = float(s[1]) if s is not None and len(s) > 1 else None
        sat = cell_mean(m, "saturation")
        micro_p = cell_mean(m, "micro_pressure")
        swell = None
        ss = m.cell_data.get("swelling_stress") or m.point_data.get("swelling_stress")
        if ss is not None:
            ssa = np.asarray(ss[0] if isinstance(ss, list) else ss).reshape(-1)
            if ssa.size > 1:
                swell = float(ssa.reshape(-1, 4 if ssa.size % 4 == 0 else ssa.size)[...,1].mean()) if ssa.size>=4 else float(ssa.mean())
        rows.append(dict(t=t, e=e, phi=phi,
                         sigma_zz_Pa=szz,
                         sigma_axial_MPa=(-szz/1e6 if szz is not None else None),
                         saturation=sat, micro_pressure_Pa=micro_p))
    rows.sort(key=lambda r: r["t"])
    with open(os.path.join(ROOT, "results", out_csv), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"{path_dir}: {len(rows)} rows -> results/{out_csv}")
    for r in rows:
        print(f"  t={r['t']:>12.0f}  e={r['e']:.4f}  sig_ax={r['sigma_axial_MPa']:>8.3f} MPa  "
              f"S={r['saturation']:.3f}  micro_p={r['micro_pressure_Pa']:.3e}")
    return rows

if __name__ == "__main__":
    for d, out in (("p1", "path1_e_sigma.csv"), ("p2", "path2_e_sigma.csv")):
        if glob.glob(os.path.join(ROOT, "results", d, "*.vtu")):
            extract(d, out)
        else:
            print(f"{d}: no VTUs yet")
