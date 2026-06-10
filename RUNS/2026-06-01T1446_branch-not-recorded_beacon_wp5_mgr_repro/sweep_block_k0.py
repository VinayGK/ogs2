#!/usr/bin/env python3
"""
MGR23 block-permeability CALIBRATION sweep (Vinay directive 2026-06-01):
"make the block permeability in the calibrative exercise as high as necessary to agree to
experiments, then use that calibration for the predictive exercise."

Calibrate block k0 (lambda=0 held) on MGR23 (pellet bottom = wet end) to match measured
homogenisation; freeze the matched k0 for the MGR27 blind prediction (separate script).

GUARDRAIL NOTE (CLAUDE.md §12.5/§1.1): block k0 here is a DIAGNOSTIC CALIBRATION value tuned to
the BEACON D5.6 MGR23 homogenisation target (cited), per explicit user directive. NOT a literature
constant. Block swelling K stays Villar-frozen (25268 J/kg); only the transport knob moves.
This is calibrate-on-MGR23 / predict-on-MGR27 (different geometry) => NOT fit-and-verify in one
test (§2 satisfied). The MGR23 match itself is 'calibrated, not validated'; MGR27 is the validation.

MEASURED MGR23 (D5.6 Tab 2-5 / EXTRACTION_MGR.md:49-51): pellet 1.30->1.34, block 1.60->1.51,
gap 0.30->0.17 (~43% closure). Geometry: pellet BOTTOM (y 0-0.05), block TOP (y 0.05-0.10).
"""
import os, re, subprocess, glob
import numpy as np, meshio

ROOT  = "/Users/vinaykumar/git/ogs/beacon_wp5_mgr_repro_2026-06-01"
MODEL = ROOT + "/model"
OGS   = "/Users/vinaykumar/git/build/dsm-native-omp-release/bin/ogs"
RHO_S = 2700.0
BASE  = open(f"{MODEL}/sw_A_blk0.prj").read()
assert BASE.count("5.0e-19 5.0e-19") == 1, "block k0 string not unique"   # pellet is 6.5e-16

# block k0 sweep at lambda=0: base, x10, x100, x1000, x10000 (spans above pellet 6.5e-16)
CASES = [("base", "5.0e-19"), ("x10", "5.0e-18"), ("x100", "5.0e-17"),
         ("x1k", "5.0e-16"), ("x10k", "5.0e-15")]

def make(tag, k0):
    s = BASE.replace("5.0e-19 5.0e-19", f"{k0} {k0}")
    s = re.sub(r"<prefix>[^<]+</prefix>", f"<prefix>{tag}</prefix>", s)
    pf = f"{MODEL}/{tag}.prj"; open(pf, "w").write(s)
    return pf

def run(tag):
    od = f"{ROOT}/results/{tag}"; os.makedirs(od, exist_ok=True)
    r = subprocess.run([OGS, f"{tag}.prj", "-o", od], cwd=MODEL, capture_output=True, text=True)
    log = r.stdout + r.stderr
    lastt = re.findall(r"Time: ([0-9.]+)", log)
    return ("Simulation completed" in log), (float(lastt[-1])/86400 if lastt else -1)

def metrics(tag):
    fs = sorted(glob.glob(f"{ROOT}/results/{tag}/{tag}_*.vtu"),
                key=lambda f: float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)))
    if not fs: return None
    m = meshio.read(fs[-1]); y = m.points[:, 1]
    def Lr(fld, lo, hi):
        a = np.asarray(m.point_data[fld]); s = (y >= lo) & (y <= hi); return float(a[s].mean())
    pe = RHO_S*(1-Lr("porosity", 0, 0.0499))/1000      # pellet = BOTTOM (MGR23)
    bl = RHO_S*(1-Lr("porosity", 0.0501, 0.10))/1000   # block  = TOP
    bnl = Lr("micro_water_content", 0.0501, 0.10)
    pnl = Lr("micro_water_content", 0, 0.0499)
    top = np.abs(y-y.max()) < 1e-7
    szz = -np.asarray(m.point_data["sigma"])[top, 1].mean()/1e6
    return dict(pe=pe, bl=bl, gap=bl-pe, bnl=bnl, pnl=pnl, szz=szz)

out = ["=== MGR23 block-k0 CALIBRATION sweep (lambda=0 held; pellet k0 6.5e-16 lambda=9 frozen) ===",
       "MEASURED MGR23: pellet 1.34, block 1.51, gap 0.17 (~43% closure); init gap 0.30",
       f"{'case':<6}{'block_k0':>10}{'done':>6}{'lastT':>7}{'pel_rd':>8}{'blk_rd':>8}{'gap':>7}{'blk_nl':>8}{'pel_nl':>8}{'szz':>7}"]
for tag, k0 in CASES:
    cname = f"cal_{tag}"
    make(cname, k0); done, lt = run(cname); mm = metrics(cname)
    if mm:
        out.append(f"{tag:<6}{k0:>10}{str(done):>6}{lt:>7.0f}{mm['pe']:>8.3f}{mm['bl']:>8.3f}"
                   f"{mm['gap']:>7.3f}{mm['bnl']:>8.3f}{mm['pnl']:>8.3f}{mm['szz']:>7.2f}")
    else:
        out.append(f"{tag:<6}{k0:>10}{str(done):>6}{lt:>7.0f}   NO VTU")
txt = "\n".join(out)
open(f"{ROOT}/results/BLOCK_K0_CAL.txt", "w").write(txt+"\n")
print(txt)
