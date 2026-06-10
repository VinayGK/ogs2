#!/usr/bin/env python3
"""
MGR23 MCC alpha_M (exchange-rate) sweep (Vinay 2026-06-01): does slowing the micro<->macro exchange
gentle the Pi-path swelling-stress onset enough for the semi-explicit MCC to integrate THROUGH it
(and reach pellet yield/collapse)? Lower alpha_M -> slower interlayer fill -> gentler dsigma/dt.
Hold everything else (constE, FEBEX pc_pel=0.4 so the pellet genuinely yields, gentle ramp 10d,
dt_max 0.1d). DIAGNOSTIC sweep (alpha_M is phenomenological/unidentified); a low alpha_M that only
"works" by freezing the kinetics is reported as such (gap won't close for the right reason).
LE reference: gap 0.192 (pellet 1.350 / block 1.541); elastic pellet cap ~1.353 (>1.36 => yielded).
"""
import os, re, subprocess, glob
import numpy as np, meshio

ROOT="/Users/vinaykumar/git/ogs/beacon_wp5_mgr_repro_2026-06-01"; MODEL=ROOT+"/model"
OGS="/Users/vinaykumar/git/build/dsm-native-omp-release/bin/ogs"; RHO_S=2700.0
BASE=open(f"{MODEL}/mgr23_mcc.prj").read()
# reset pc_pel to the FEBEX-sourced 0.4 MPa (current file has 1.0); sweep alpha_M
BASE=re.sub(r'(<name>pc0_pel</name><type>Constant</type><value>)[^<]+(</value>)', r'\g<1>0.4e6\g<2>', BASE)
assert "0.4e6" in BASE

ALPHAS=[("1e-13","1e-13"),("1e-14","1e-14"),("1e-15","1e-15"),("1e-16","1e-16"),("1e-17","1e-17")]

def run(tag, aM):
    s=re.sub(r'<mass_exchange_coefficient>[^<]+</mass_exchange_coefficient>',
             f'<mass_exchange_coefficient>{aM}</mass_exchange_coefficient>', BASE)
    s=re.sub(r"<prefix>[^<]+</prefix>", f"<prefix>{tag}</prefix>", s)
    pf=f"{MODEL}/{tag}.prj"; open(pf,"w").write(s)
    import xml.etree.ElementTree as ET; ET.parse(pf)
    od=f"{ROOT}/results/{tag}"; os.makedirs(od, exist_ok=True)
    r=subprocess.run([OGS,f"{tag}.prj","-o",od], cwd=MODEL, capture_output=True, text=True)
    log=r.stdout+r.stderr
    done="Simulation completed" in log
    tt=re.findall(r"Time: ([0-9.eE+]+)", log)
    lastt=float(tt[-1].rstrip('.'))/86400 if tt else -1
    rej=log.count("rejected")
    return done,lastt,rej,od

def metrics(od,tag):
    fs=sorted([f for f in glob.glob(f"{od}/{tag}_*.vtu") if re.search(r'_t_([0-9.]+)\.vtu',f)],
              key=lambda f:float(re.search(r'_t_([0-9.]+)\.vtu',f).group(1)))
    if not fs: return None
    m=meshio.read(fs[-1]); y=m.points[:,1]
    def L(fld,lo,hi):
        a=np.asarray(m.point_data[fld]);s=(y>=lo)&(y<=hi);return float(a[s].mean())
    pe=RHO_S*(1-L("porosity",0,0.0499))/1000; bl=RHO_S*(1-L("porosity",0.0501,0.10))/1000
    return dict(pe=pe,bl=bl,gap=bl-pe)

out=["=== MGR23 MCC alpha_M sweep (constE, FEBEX pc_pel=0.4, ramp10d, dt0.1d) ===",
     "LE ref: gap 0.192 (pel 1.350/blk 1.541); pellet elastic cap ~1.353 (rho_d>1.36 => YIELDED)",
     f"{'alpha_M':>8}{'done':>6}{'lastT[d]':>9}{'rej':>5}{'pel_rd':>8}{'blk_rd':>8}{'gap':>7}  verdict"]
for tag,aM in ALPHAS:
    cname=f"mcc_aM_{tag}"
    done,lt,rej,od=run(cname,aM); mm=metrics(od,cname)
    if mm:
        yld = "YIELDED/collapsed" if mm['pe']>1.36 else "elastic-only"
        v = f"{yld}; gap {mm['gap']:.3f}" if done else f"WALLED @{lt:.1f}d"
        out.append(f"{aM:>8}{str(done):>6}{lt:>9.1f}{rej:>5}{mm['pe']:>8.3f}{mm['bl']:>8.3f}{mm['gap']:>7.3f}  {v}")
    else:
        out.append(f"{aM:>8}{str(done):>6}{lt:>9.1f}{rej:>5}   NO VTU")
txt="\n".join(out); open(f"{ROOT}/results/ALPHAM_MCC_SWEEP.txt","w").write(txt+"\n"); print(txt)
