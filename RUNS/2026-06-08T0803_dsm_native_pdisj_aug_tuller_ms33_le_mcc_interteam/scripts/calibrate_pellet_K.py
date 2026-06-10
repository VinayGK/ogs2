#!/usr/bin/env python3
"""Secant-fit the pellet-zone DSM swelling prefactor K so the single-element
confined swelling pressure at rho_d=900 (phi0=0.6763) matches the Dixon (2023)
EMDD==rho_d extrapolated target Ps=0.051 MPa. Same single-element procedure as
Model I per-density calibration. NO assertion is made on a Model-IV deliverable
derived from this K (calibration on a separate single-element case; not circular)."""
import os, re, glob, subprocess, sys
import numpy as np, meshio

RUN = "/Users/vinaykumar/git/ogs/ms33_pdisj_aug_tuller_2026-06-08"
OGS = "/Users/vinaykumar/git/build/ufz_integration_20260602/bin/ogs"
CAL = f"{RUN}/calib_pellet"
TMPL = open(f"{CAL}/ms33_pellet_dd900.tmpl.prj").read()
TARGET = 0.3501   # MPa, Dixon(2023) median 0.003*exp(5.2883*0.9), EMDD==rho_d (suite basis)
RELTOL = 0.02
env = dict(os.environ, HOME="/tmp", MPLCONFIGDIR="/tmp/mplcache", XDG_CACHE_HOME="/tmp", OMP_NUM_THREADS="1")

def Ps_of_K(K):
    od = f"{CAL}/out_K{K:.4g}"
    os.makedirs(od, exist_ok=True)
    prj = f"{CAL}/run.prj"
    open(prj, "w").write(TMPL.replace("__KPELLET__", repr(float(K))))
    subprocess.run([OGS, prj, "-o", od], env=env, cwd=CAL,
                   stdout=open(f"{od}/run.log","w"), stderr=subprocess.STDOUT)
    vtus = sorted(glob.glob(f"{od}/pellet_dd900_ts_*.vtu"),
                  key=lambda f: int(re.search(r"_ts_(\d+)_", f).group(1)))
    if not vtus: return None
    # most-saturated frame
    # Calibrate on the SWELLING_STRESS field (the DSM swelling-induced stress =
    # the seating-independent swelling pressure). Total mean stress is floored by
    # the sigma0 = -1.5e5 Pa (0.15 MPa) seating, which dwarfs the 0.051 MPa target
    # at rho_d=900; swelling_stress isolates the swelling contribution (Dixon analog).
    best = None
    for v in vtus:
        m = meshio.read(v); s = float(m.point_data["saturation"].mean())
        ss = m.point_data["swelling_stress"].mean(axis=0)
        psw = -(ss[0]+ss[1]+ss[2])/3.0/1e6
        if best is None or s >= best[0]:
            best = (s, psw)
    return best[1]

# secant in log K (swelling pressure ~ linear in K; seed near header estimate scaled up)
K = [8000.0, 14000.0]
P = [Ps_of_K(K[0]), Ps_of_K(K[1])]
print(f"seed: K={K[0]:.4g} -> Ps={P[0]:.4f} | K={K[1]:.4g} -> Ps={P[1]:.4f} (target {TARGET})")
hist = list(zip(K, P))
for it in range(10):
    lo = abs(P[-1]-TARGET)/TARGET
    if lo <= RELTOL:
        break
    # secant on log Ps vs log K
    lk0, lk1 = np.log(K[-2]), np.log(K[-1])
    lp0, lp1 = np.log(P[-2]), np.log(P[-1])
    lt = np.log(TARGET)
    lk = lk1 + (lt-lp1)*(lk1-lk0)/(lp1-lp0)
    Kn = float(np.exp(lk)); Pn = Ps_of_K(Kn)
    K.append(Kn); P.append(Pn); hist.append((Kn, Pn))
    print(f"  it{it}: K={Kn:.5g} -> Ps={Pn:.5f} MPa  (dev {abs(Pn-TARGET)/TARGET*100:.1f}%)")

print("\nCALIBRATED K_pellet = %.5g J/kg  -> Ps=%.5f MPa (target %.3f, dev %.1f%%)"
      % (K[-1], P[-1], TARGET, abs(P[-1]-TARGET)/TARGET*100))
open(f"{CAL}/K_pellet.txt","w").write("%.6g\n"%K[-1])
