#!/usr/bin/env python3
"""Extract MS33 Model VII free-swell run in the team comparison quantity:
void ratio (e) vs mean/axial stress (MPa), at axis probe points Top/Central/Bottom.
Team convention: stress positive in COMPRESSION, MPa. OGS sigma is tension-positive,
components (axisym) = [s_rr, s_zz, s_thth, s_rz]; effective stress.
Mean stress (total) and axial (total) reported compression-positive in MPa.
Probe points (m): Top (0,0.070), Central (0,0.040), Bottom (0,0.010)."""
import sys, glob, os, re
import numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy

scr=sys.argv[1]
probes={'Top':(0.0,0.070),'Central':(0.0,0.040),'Bottom':(0.0,0.010)}

def list_outputs(scr):
    fs=glob.glob(os.path.join(scr,"ms33_modelVII_freeswelling_ts_*_t_*.vtu"))
    out=[]
    for f in fs:
        m=re.search(r"_ts_(\d+)_t_([0-9.]+)\.vtu",f)
        if m: out.append((int(m.group(1)), float(m.group(2)), f))
    return sorted(out)

def probe_file(f):
    r=vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(f); r.Update()
    g=r.GetOutput()
    # build probe point set
    pts=vtk.vtkPoints()
    for name in probes:
        x,z=probes[name]; pts.InsertNextPoint(x,z,0.0)
    poly=vtk.vtkPolyData(); poly.SetPoints(pts)
    pf=vtk.vtkProbeFilter(); pf.SetInputData(poly); pf.SetSourceData(g); pf.Update()
    res=pf.GetOutput(); pdout=res.GetPointData()
    def arr(n): 
        a=pdout.GetArray(n); return vtk_to_numpy(a) if a is not None else None
    poro=arr('porosity'); sig=arr('sigma'); pres=arr('pressure'); sat=arr('saturation')
    sw=arr('swelling_stress'); disp=arr('displacement')
    rows={}
    for i,name in enumerate(probes):
        phi=float(poro[i])
        e=phi/(1.0-phi)
        s=sig[i]  # [s_rr, s_zz, s_thth, s_rz] tension-positive, effective
        s_rr,s_zz,s_th=float(s[0]),float(s[1]),float(s[2])
        p=float(pres[i]); Sw=float(sat[i])
        # Bishop chi = Sw (BishopsSaturationCutoff cutoff=1 -> chi=Sw until Sw=1)
        chi=min(Sw,1.0)
        # total stress = effective + biot*chi*p*I (p<0 suction); tension-positive
        biot=1.0
        s_rr_tot=s_rr+biot*chi*p
        s_zz_tot=s_zz+biot*chi*p
        s_th_tot=s_th+biot*chi*p
        mean_eff=(s_rr+s_zz+s_th)/3.0
        mean_tot=(s_rr_tot+s_zz_tot+s_th_tot)/3.0
        # compression-positive MPa
        rows[name]=dict(
            phi=phi, e=e, Sw=Sw, p_MPa=p/1e6,
            axial_eff_comp_MPa=-s_zz/1e6, axial_tot_comp_MPa=-s_zz_tot/1e6,
            radial_eff_comp_MPa=-s_rr/1e6, radial_tot_comp_MPa=-s_rr_tot/1e6,
            mean_eff_comp_MPa=-mean_eff/1e6, mean_tot_comp_MPa=-mean_tot/1e6,
        )
    return rows

if __name__=='__main__':
    outs=list_outputs(scr)
    sel=sys.argv[2:] if len(sys.argv)>2 else None
    for ts,t,f in outs:
        if sel and str(ts) not in sel: continue
        rows=probe_file(f)
        td=t/86400.0
        for name in probes:
            r=rows[name]
            print(f"ts={ts:5d} t={td:7.2f}d {name:8s} e={r['e']:.4f} phi={r['phi']:.4f} Sw={r['Sw']:.4f} "
                  f"p={r['p_MPa']:8.3f}MPa | axial_tot={r['axial_tot_comp_MPa']:8.4f} mean_tot={r['mean_tot_comp_MPa']:8.4f} "
                  f"axial_eff={r['axial_eff_comp_MPa']:8.4f} mean_eff={r['mean_eff_comp_MPa']:8.4f} (comp+,MPa)")
