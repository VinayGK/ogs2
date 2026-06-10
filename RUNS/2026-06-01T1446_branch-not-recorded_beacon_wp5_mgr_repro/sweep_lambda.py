#!/usr/bin/env python3
"""
MGR permeability-EXPONENT (lambda) study (Vinay 2026-06-01).
k = k0*(phi_M/phi_M0)^lambda per region. Base = mgr23e (both k0 x100, alpha_M 1e-13,
lambda 9/9). Diagnostic sweep to find WHICH region's macro-porosity throttle gates the
homogenisation; lambda values are a sweep, NOT asserted material constants (CLAUDE.md
§1.1 sweep-to-locate-threshold). A matched pair would be 'calibrated, not validated' (§2).

Study A: pellet lambda = 9 (held), block lambda in {9,3,1,0}
Study B: block  lambda = 9 (held), pellet lambda in {9,3,1,0}
(case 9/9 is the base mgr23e; run once.)
Each case: replace the two positional <exponents> blocks (pellet first, block second),
run, extract @last-time: density gap, block n_l, szz_top. Verified VTU numbers only.
"""
import os, re, subprocess, glob
import numpy as np, meshio

ROOT = "/Users/vinaykumar/git/ogs/beacon_wp5_mgr_repro_2026-06-01"
MODEL = ROOT + "/model"
OGS = "/Users/vinaykumar/git/build/dsm-native-omp-release/bin/ogs"
RHO_S = 2700.0
BASE = open(f"{MODEL}/mgr23e_both100x.prj").read()
assert BASE.count("<exponents>9 9</exponents>") == 2

def make(pel_lam, blk_lam, tag):
    # replace the TWO exponents blocks POSITIONALLY via a counting regex sub
    # (sequential str.replace fails when pel_lam stays '9' — it re-hits the 1st block).
    repl = [f"<exponents>{pel_lam} {pel_lam}</exponents>",   # 1st = pellet
            f"<exponents>{blk_lam} {blk_lam}</exponents>"]   # 2nd = block
    it = iter(repl)
    s = re.sub(r"<exponents>[^<]+</exponents>", lambda m: next(it), BASE, count=2)
    s = s.replace("<prefix>mgr23e_both100x</prefix>", f"<prefix>{tag}</prefix>")
    # verify positional correctness via XML
    import xml.etree.ElementTree as ET
    pf = f"{MODEL}/{tag}.prj"
    open(pf, "w").write(s)
    r = ET.parse(pf).getroot()
    exps = [p.findtext("exponents") for med in r.find("media").findall("medium")
            for p in med.find("properties").findall("property") if p.findtext("name") == "permeability"]
    assert exps == [f"{pel_lam} {pel_lam}", f"{blk_lam} {blk_lam}"], f"{tag}: {exps}"
    return pf

def run(tag):
    od = f"{ROOT}/results/{tag}"
    os.makedirs(od, exist_ok=True)
    r = subprocess.run([OGS, f"{tag}.prj", "-o", od], cwd=MODEL, capture_output=True, text=True)
    log = r.stdout + r.stderr
    done = "Simulation completed" in log
    lastt = re.findall(r"Time: ([0-9.]+)", log)
    return done, (float(lastt[-1])/86400 if lastt else -1)

def metrics(tag):
    fs = sorted(glob.glob(f"{ROOT}/results/{tag}/{tag}_*.vtu"),
                key=lambda f: float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)))
    if not fs: return None
    m = meshio.read(fs[-1]); y = m.points[:, 1]
    def Lr(fld, lo, hi):
        a = np.asarray(m.point_data[fld]); s = (y >= lo) & (y <= hi); return float(a[s].mean())
    pe = RHO_S*(1-Lr("porosity", 0, 0.0499))/1000
    bl = RHO_S*(1-Lr("porosity", 0.0501, 0.10))/1000
    bnl = Lr("micro_water_content", 0.0501, 0.10)
    top = np.abs(y-y.max()) < 1e-7
    szz = -np.asarray(m.point_data["sigma"])[top, 1].mean()/1e6
    t = float(re.search(r"_t_([0-9.]+)\.vtu", fs[-1]).group(1))/86400
    return dict(t=t, pe=pe, bl=bl, gap=bl-pe, bnl=bnl, szz=szz)

CASES = [("A_blk9","9","9"),  # = base 9/9
         ("A_blk3","9","3"), ("A_blk1","9","1"), ("A_blk0","9","0"),
         ("B_pel3","3","9"), ("B_pel1","1","9"), ("B_pel0","0","9")]
out = ["=== MGR lambda (perm-exponent) study — base mgr23e (k0 x100 both, alpha_M 1e-13) ===",
       "MEASURED: gap 0.17 ; szz ~3 MPa",
       f"{'case':<10}{'pel_lam':>8}{'blk_lam':>8}{'done':>6}{'lastT[d]':>9}{'gap':>7}{'blk_nl':>8}{'szz[MPa]':>9}"]
for tag, pl, bl in CASES:
    pf = make(pl, bl, f"sw_{tag}")
    done, lt = run(f"sw_{tag}")
    mm = metrics(f"sw_{tag}")
    if mm:
        out.append(f"{tag:<10}{pl:>8}{bl:>8}{str(done):>6}{lt:>9.1f}{mm['gap']:>7.3f}{mm['bnl']:>8.3f}{mm['szz']:>9.3f}")
    else:
        out.append(f"{tag:<10}{pl:>8}{bl:>8}{str(done):>6}{lt:>9.1f}   NO VTU")
txt = "\n".join(out)
open(f"{ROOT}/results/LAMBDA_SWEEP.txt", "w").write(txt+"\n")
print(txt)
