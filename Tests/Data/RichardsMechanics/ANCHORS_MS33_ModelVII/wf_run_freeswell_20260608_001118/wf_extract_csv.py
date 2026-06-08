#!/usr/bin/env python3
"""Final extraction: write CSV of void ratio vs mean/axial stress at 3 axis probe
points for every output frame. Comparison quantity = void ratio (e) vs mean stress (MPa),
compression-positive (team convention, EURAD-2 MS33 Model VII 'Free Swelling' sheet)."""
import sys, glob, os, re, csv
import numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy
scr=sys.argv[1]
probes=[('Top',0.0,0.070),('Central',0.0,0.040),('Bottom',0.0,0.010)]
def outputs(scr):
    o=[]
    for f in glob.glob(os.path.join(scr,"ms33_modelVII_freeswelling_ts_*_t_*.vtu")):
        m=re.search(r"_ts_(\d+)_t_([0-9.]+)\.vtu",f)
        if m: o.append((int(m.group(1)),float(m.group(2)),f))
    return sorted(o)
def probe(f):
    r=vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(f); r.Update(); g=r.GetOutput()
    pts=vtk.vtkPoints()
    for _,x,z in probes: pts.InsertNextPoint(x,z,0.0)
    poly=vtk.vtkPolyData(); poly.SetPoints(pts)
    pf=vtk.vtkProbeFilter(); pf.SetInputData(poly); pf.SetSourceData(g); pf.Update()
    pd=pf.GetOutput().GetPointData()
    A=lambda n:(vtk_to_numpy(pd.GetArray(n)) if pd.GetArray(n) is not None else None)
    poro,sig,pres,sat=A('porosity'),A('sigma'),A('pressure'),A('saturation')
    rows=[]
    for i,(name,_,_) in enumerate(probes):
        phi=float(poro[i]); e=phi/(1-phi)
        s_rr,s_zz,s_th=float(sig[i][0]),float(sig[i][1]),float(sig[i][2])
        p=float(pres[i]); Sw=float(sat[i]); chi=min(Sw,1.0); biot=1.0
        srr_t=s_rr+biot*chi*p; szz_t=s_zz+biot*chi*p; sth_t=s_th+biot*chi*p
        mean_eff=(s_rr+s_zz+s_th)/3.0; mean_tot=(srr_t+szz_t+sth_t)/3.0
        rows.append(dict(loc=name,phi=phi,e=e,Sw=Sw,p_MPa=p/1e6,
            axial_eff_MPa=-s_zz/1e6, axial_tot_MPa=-szz_t/1e6,
            radial_eff_MPa=-s_rr/1e6, radial_tot_MPa=-srr_t/1e6,
            mean_eff_MPa=-mean_eff/1e6, mean_tot_MPa=-mean_tot/1e6))
    return rows
o=outputs(scr)
out_csv=os.path.join(scr,"wf_modelVII_void_vs_stress.csv")
with open(out_csv,'w',newline='') as fh:
    w=csv.writer(fh)
    w.writerow(['timestep','time_s','time_d','location','porosity','void_ratio',
                'saturation','pressure_MPa','axial_eff_MPa','axial_tot_MPa',
                'radial_eff_MPa','radial_tot_MPa','mean_eff_MPa','mean_tot_MPa'])
    for ts,t,f in o:
        for r in probe(f):
            w.writerow([ts,t,t/86400.0,r['loc'],f"{r['phi']:.6f}",f"{r['e']:.6f}",
                        f"{r['Sw']:.6f}",f"{r['p_MPa']:.4f}",f"{r['axial_eff_MPa']:.5f}",
                        f"{r['axial_tot_MPa']:.5f}",f"{r['radial_eff_MPa']:.5f}",
                        f"{r['radial_tot_MPa']:.5f}",f"{r['mean_eff_MPa']:.5f}",f"{r['mean_tot_MPa']:.5f}"])
print("WROTE",out_csv,"frames=",len(o),"rows=",len(o)*3)
# print summary at key locations (Central) for each frame
print("\n-- Central probe summary (void e vs mean/axial stress, comp+ MPa) --")
for ts,t,f in o:
    for r in probe(f):
        if r['loc']=='Central':
            print(f"t={t/86400.0:7.2f}d e={r['e']:.4f} Sw={r['Sw']:.3f} axial_tot={r['axial_tot_MPa']:8.4f} mean_eff={r['mean_eff_MPa']:8.4f} mean_tot={r['mean_tot_MPa']:8.4f}")
