#!/usr/bin/env python3
"""Extract central-point (r=0, y=0.040 m) time series for Model IX from a PVD.
Mean stress reported COMPRESSION-POSITIVE (= -(sxx+syy+szz)/3) to match team xlsx sign.
Also reports suction (=-pressure) and intrinsic permeability kxx.
Usage: python3 extract_central.py <pvd> <out.json>
"""
import sys, json, xml.etree.ElementTree as ET, os
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

pvd, out = sys.argv[1], sys.argv[2]
base = os.path.dirname(os.path.abspath(pvd))
tree = ET.parse(pvd)
ds = [(float(d.get("timestep")), d.get("file")) for d in tree.iter("DataSet")]

def read_vtu(fn):
    r = vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(os.path.join(base, fn)); r.Update()
    return r.GetOutput()

# locate central node from t=0 mesh: r=x≈0, y≈0.040
g0 = read_vtu(ds[0][1])
pts = vtk_to_numpy(g0.GetPoints().GetData())
target = np.array([0.0, 0.040])
d = np.hypot(pts[:,0]-target[0], pts[:,1]-target[1])
ni = int(np.argmin(d))
print("central node idx", ni, "coord", pts[ni], "dist", d[ni])

res = {"t_day": [], "mean_MPa": [], "suction_MPa": [], "kxx_m2": [], "porosity": []}
for t, fn in ds:
    g = read_vtu(fn)
    pd = g.GetPointData()
    sig = vtk_to_numpy(pd.GetArray("sigma"))[ni]          # [xx,yy,zz,xy] tension-positive
    p   = vtk_to_numpy(pd.GetArray("pressure"))[ni]
    mean_comp = -(sig[0]+sig[1]+sig[2])/3.0 / 1e6         # compression-positive, MPa
    res["t_day"].append(t/86400.0)
    res["mean_MPa"].append(round(float(mean_comp),5))
    res["suction_MPa"].append(round(float(-p)/1e6,5))
    try:
        k = vtk_to_numpy(pd.GetArray("intrinsic_permeability"))[ni]
        kxx = float(k[0]) if hasattr(k,"__len__") else float(k)
        res["kxx_m2"].append(kxx)
    except Exception:
        res["kxx_m2"].append(None)
    try:
        res["porosity"].append(round(float(vtk_to_numpy(pd.GetArray("porosity"))[ni]),5))
    except Exception:
        res["porosity"].append(None)

json.dump(res, open(out,"w"), indent=1)
print("CENTRAL mean stress (compression+) MPa:")
for t,m,s in zip(res["t_day"],res["mean_MPa"],res["suction_MPa"]):
    print("  t=%7.2f d  mean=%8.4f MPa  suction=%9.4f MPa"%(t,m,s))
