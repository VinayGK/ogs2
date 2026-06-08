#!/usr/bin/env python3
"""wf_extract.py  (scratch; MS33 Model III convergence drive, 2026-06-08)

Extract the EURAD-2 MS33 Model III team-comparable quantities from an OGS run,
in the exact quantity + units the teams report in
  .../g_Support_Section_Data Collection/<TEAM>_DATA/Model_III/Model_III_*.xlsx
  sheet "Model_III":  Final mean / radial / axial stress [MPa], final porosity,
  final dry density [g/cm3], time(d) to reach s=0, initial/final permeability [m2]
at four probe locations (mm, axisym x=r y=z):
  Top(0,70) Central(0,40) Bottom(0,10) Clay-gap intersection(23,40).

Sign convention (anchors-ms33-workflow skill): OGS sigma is signed, tension +.
Team "stress" columns are COMPRESSION-POSITIVE, so we negate:
  radial  stress = -sigma_xx   (xx == rr in axisym)
  axial   stress = -sigma_yy   (yy == zz)
  hoop    stress = -sigma_zz   (zz == theta-theta, the 3rd comp in OGS 2D axisym)
  mean    stress = -(sigma_xx+sigma_yy+sigma_zz)/3
suction[MPa] = -pressure/1e6.  dry_density: g/cm3 = dry_density_solid/1000.

Usage: wf_extract.py <prefix.pvd | final.vtu> [--pvd]
"""
import sys, os, glob, re
import numpy as np
import pyvista as pv

MESH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'ms33_gap2mm_r27_h70.vtu')

PROBES = {  # name: (r[m], z[m])
    'Top':            (0.000, 0.070),
    'Central':        (0.000, 0.040),
    'Bottom':         (0.000, 0.010),
    'ClayGap_r23z40': (0.023, 0.040),
}

def probe_nodes():
    m = pv.read(MESH); pts = m.points
    out = {}
    for name, (r, z) in PROBES.items():
        d = np.hypot(pts[:, 0]-r, pts[:, 1]-z)
        out[name] = int(np.argmin(d))
    return out, pts

def find_final_vtu(arg):
    if arg.endswith('.vtu'):
        return arg
    # pvd -> last DataSet timestep
    if arg.endswith('.pvd'):
        txt = open(arg).read()
        files = re.findall(r'file="([^"]+)"', txt)
        times = [float(t) for t in re.findall(r'timestep="([^"]+)"', txt)]
        if files:
            base = os.path.dirname(os.path.abspath(arg))
            # pick max time
            i = int(np.argmax(times)) if times else len(files)-1
            return os.path.join(base, files[i]), (times[i] if times else None)
    # else: treat as prefix, glob
    cands = sorted(glob.glob(arg+'*_ts_*_t_*.vtu'))
    if not cands:
        raise SystemExit(f'no vtu found for {arg}')
    # sort by timestep number
    def tsnum(p):
        m = re.search(r'_ts_(\d+)_', p); return int(m.group(1)) if m else -1
    cands.sort(key=tsnum)
    return cands[-1]

def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    arg = sys.argv[1]
    res = find_final_vtu(arg)
    tfinal = None
    if isinstance(res, tuple):
        vtu, tfinal = res
    else:
        vtu = res
    nodes, pts = probe_nodes()
    m = pv.read(vtu)
    # time from filename if not from pvd
    if tfinal is None:
        mm = re.search(r'_t_([0-9.]+)\.vtu', vtu)
        if mm:
            tfinal = float(mm.group(1))
    sig = m.point_data['sigma']            # (n,4): xx,yy,zz,xy
    p   = m.point_data['pressure']
    sat = m.point_data['saturation']
    por = m.point_data['porosity']
    dds = m.point_data['dry_density_solid']
    perm = m.point_data['intrinsic_permeability']  # (n,4)
    sw  = m.point_data.get('swelling_stress')
    print(f'# VTU: {vtu}')
    print(f'# t_final = {tfinal} s = {tfinal/86400:.3f} d' if tfinal else '# t_final unknown')
    print(f'# sign: stresses COMPRESSION-POSITIVE (team convention); suction=-p/1e6 MPa')
    print()
    hdr = f'{"probe":16s} {"mean[MPa]":>10s} {"radial[MPa]":>11s} {"axial[MPa]":>10s} {"hoop[MPa]":>9s} {"swell_mean[MPa]":>15s} {"suction[MPa]":>12s} {"S_L":>7s} {"poros":>7s} {"dd[g/cm3]":>9s} {"k[m2]":>10s}'
    print(hdr)
    rows = {}
    for name, nid in nodes.items():
        sxx, syy, szz = sig[nid, 0], sig[nid, 1], sig[nid, 2]
        mean = -(sxx+syy+szz)/3.0/1e6
        rad  = -sxx/1e6
        ax   = -syy/1e6
        hoop = -szz/1e6
        swm = (-(sw[nid,0]+sw[nid,1]+sw[nid,2])/3.0/1e6) if sw is not None else float('nan')
        suc = -p[nid]/1e6
        kk = perm[nid, 0]
        print(f'{name:16s} {mean:10.4f} {rad:11.4f} {ax:10.4f} {hoop:9.4f} {swm:15.4f} {suc:12.4f} {sat[nid]:7.4f} {por[nid]:7.4f} {dds[nid]/1000:9.4f} {kk:10.3e}')
        rows[name] = dict(mean=mean, radial=rad, axial=ax, hoop=hoop, swell_mean=swm,
                          suction=suc, S_L=float(sat[nid]), porosity=float(por[nid]),
                          dry_density_g_cm3=float(dds[nid]/1000), perm=float(kk))
    print()
    print('# Team comparison headline = FINAL MEAN STRESS [MPa] at Central (z=40 mm).')
    print(f'# our Central mean stress = {rows["Central"]["mean"]:.4f} MPa')
    return rows

if __name__ == '__main__':
    main()
