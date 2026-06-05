#!/usr/bin/env python3
"""Recalibrate K to the Dixon MX-80 Ps(rho_d) anchors under
micro_solid_volume_fraction_mode = current_porosity_split.

Per-density secant on K (Ps ~ linear in K). On a saturation-transition linear-solver
failure (singular at S_L->1, from the not-yet-consistent dnS/dn_l chain in split mode)
the dt is reduced and K is KEPT (the failure is numerical, not K-too-high).
Calibration target = Dixon (2023) MX-80, sigma_swell = 0.003*exp(5.2883*rho_d[Mg/m3]) MPa
(prj-header anchors 4.922/14.161/40.86 at 1.40/1.60/1.80). Run per density, in parallel.
"""
import sys, os, re, glob, subprocess
import numpy as np, meshio

OGS = "/Users/vinaykumar/git/build/pdisj_maxwell_revref_20260605/bin/ogs"
HERE = os.path.dirname(os.path.realpath(__file__))
TARGET = {1400: 4.922, 1600: 14.161, 1800: 40.86}        # MPa, Dixon (prj headers)
K0     = {1400: 35625.4, 1600: 85312.6, 1800: 224610.0}  # reference-mode K (start)

def set_param(prj, tag, val):
    s = open(prj, encoding="iso-8859-1").read()
    s = re.sub(fr"<{tag}>[^<]+</{tag}>", f"<{tag}>{val}</{tag}>", s)
    open(prj, "w", encoding="iso-8859-1").write(s)

def run(prj, outdir):
    os.makedirs(outdir, exist_ok=True)
    with open(f"{outdir}/run.log", "w") as f:
        r = subprocess.run([OGS, "-o", outdir, "-l", "warn", prj],
                           stdout=f, stderr=subprocess.STDOUT, timeout=2400)
    return r.returncode == 0

def psw(outdir, prefix):
    fs = glob.glob(f"{outdir}/{prefix}_ts_*_t_17280000.000000.vtu")
    if not fs:
        return None
    m = meshio.read(sorted(fs)[-1]); sig = m.point_data["sigma"]
    return -float(np.mean((sig[:, 0] + sig[:, 1] + sig[:, 2]) / 3.0)) / 1e6

def fit(rho):
    prj = os.path.join(HERE, f"ms33_modelI_dd{rho}.prj"); prefix = f"ms33_modelI_dd{rho}"
    tgt = TARGET[rho]; K = K0[rho]; dt = 0.5; p = None; err = None
    for step in range(8):
        set_param(prj, "vdw_augmentation_prefactor", f"{K:.8g}")
        set_param(prj, "initial_dt", f"{dt:g}")
        o = f"/tmp/cal_dd{rho}_s{step}"
        ok = run(prj, o); p = psw(o, prefix) if ok else None
        if p is None:
            dt *= 0.3
            print(f"dd{rho} step{step} K={K:.1f} dt->{dt:g} FAILED(numerical), keep K", flush=True)
            continue
        err = (p - tgt) / tgt
        print(f"dd{rho} step{step} K={K:.1f} dt={dt:g} p_sw={p:.4f} tgt={tgt} err={err*100:+.2f}%", flush=True)
        if abs(err) < 0.005:
            break
        K = K * tgt / p
    set_param(prj, "vdw_augmentation_prefactor", f"{K:.8g}")  # leave calibrated K in prj
    print(f"dd{rho} FINAL K={K:.6g} dt={dt:g} p_sw={p} tgt={tgt} err={None if err is None else round(err*100,3)}%", flush=True)

if __name__ == "__main__":
    fit(int(sys.argv[1]))
