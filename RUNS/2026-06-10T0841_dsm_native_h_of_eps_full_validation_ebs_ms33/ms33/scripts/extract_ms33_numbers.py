#!/usr/bin/env python3
"""Extract MS33 validation-bundle numbers from the 2026-06-10 h_of_eps runs.
Probe map per task spec: Top (r=0,z=0.07), Central (r=0,z=0.04), Bottom (r=0,z=0.01).
Mean stress = -tr(sigma)/3 [MPa], compression positive.
sigma components (2D axisym Kelvin order): [0]=rr, [1]=zz(axial), [2]=hoop, [3]=rz.
"""
import glob, re, json
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

O = "/Users/vinaykumar/git/ogs/RUNS/_INPROGRESS_full_validation/ms33"

def load(p):
    r = vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(p); r.Update(); return r.GetOutput()

def tkey(f): return float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1))

def final_vtu(slug):
    fs = sorted(glob.glob(f"{O}/out_{slug}/*.vtu"), key=tkey)
    return fs[-1]

def vtu_at(slug, t):
    for f in sorted(glob.glob(f"{O}/out_{slug}/*.vtu"), key=tkey):
        if abs(tkey(f) - t) < 1.0:
            return f
    raise SystemExit(f"no frame at t={t} for {slug}")

res = {}

# ---- Model I: final sigma_zz (axial, comp [1]) mean over 4 nodes, MPa
for dd in (1400, 1600, 1800):
    s = f"modelI_dd{dd}"
    f = final_vtu(s)
    g = load(f); sig = vtk_to_numpy(g.GetPointData().GetArray("sigma"))
    res[s] = {
        "vtu": f, "t": tkey(f),
        "sigma_zz_mean_MPa": float(np.mean(sig[:, 1])) / 1e6,
        "sigma_rr_mean_MPa": float(np.mean(sig[:, 0])) / 1e6,
        "mean_stress_comp_pos_MPa": float(np.mean(-(sig[:,0]+sig[:,1]+sig[:,2])/3))/1e6,
    }

# ---- Models III, IV: mean stress at probes
def probe(g, r_m, z_m):
    pts = vtk_to_numpy(g.GetPoints().GetData())
    return int(np.argmin(np.linalg.norm(pts - np.array([r_m, z_m, 0.0]), axis=1)))

for s in ("modelIII_gap2mm", "modelIV_pellets"):
    f = final_vtu(s)
    g = load(f); sig = vtk_to_numpy(g.GetPointData().GetArray("sigma"))
    pts = vtk_to_numpy(g.GetPoints().GetData())
    d = {"vtu": f, "t": tkey(f)}
    for name, z in (("Top", 0.07), ("Central", 0.04), ("Bottom", 0.01)):
        i = probe(g, 0.0, z)
        sgm = sig[i]
        d[name] = {"node_xyz": pts[i].tolist(),
                   "mean_stress_MPa": float(-(sgm[0]+sgm[1]+sgm[2])/3)/1e6}
    res[s] = d

# ---- Model VII: void ratio from domain-mean porosity at t=17280000 and final t=20736000
s = "modelVII_freeswelling"
res[s] = {}
for tag, t in (("end_free_swelling_t17280000", 17280000.0), ("end_t20736000", 20736000.0)):
    f = vtu_at(s, t)
    g = load(f)
    phi = float(np.mean(vtk_to_numpy(g.GetPointData().GetArray("porosity"))))
    res[s][tag] = {"vtu": f, "phi_mean": phi, "void_ratio_e": phi/(1-phi)}

print(json.dumps(res, indent=2))
with open(f"{O}/extracted_numbers.json", "w") as fh:
    json.dump(res, fh, indent=2)
