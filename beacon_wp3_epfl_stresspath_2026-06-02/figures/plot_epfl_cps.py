#!/usr/bin/env python3
"""EPFL Task 3.3 e-sigma_v: model (current_porosity_split, K re-fit to Dixon) vs
measured (Ferrari et al. 2022). P1: LE free-swell (capped by small-strain at e~1.0,
cannot reach measured B=2.34) + MCC compression run from the IMPOSED measured B=2.34.
P2 (confined): split, K=56183. e=phi/(1-phi), sigma_ax=-sigma_yy/1e6."""
import glob, os, re
import numpy as np, meshio
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def traj(od, smin=1e-3):
    fs=sorted(glob.glob(os.path.join(ROOT,"results",od,"*.vtu")),
              key=lambda f:(float(re.search(r"_t_([0-9.]+)\.vtu",f).group(1)) if re.search(r"_t_([0-9.]+)\.vtu",f) else -1))
    S,E=[],[]
    for f in fs:
        m=meshio.read(f); phi=float(np.asarray(m.point_data["porosity"]).mean())
        sax=-np.asarray(m.point_data["sigma"]).mean(axis=0)[1]/1e6
        if sax>smin: S.append(sax); E.append(phi/(1-phi))
    return np.array(S),np.array(E)
S_p2,E_p2=traj("run_beacon_t33_path2_P2-1_dsm_mcc")
S_le,E_le=traj("run_p1_LE",smin=-1); ipk=int(np.argmax(E_le)); S_sw,E_sw=S_le[:ipk+1],E_le[:ipk+1]
S_c,E_c=traj("run_p1_MCCrestart",smin=-1)
mP1_s=[0.021,0.021,0.6,1.0,2.0,3.24,5.0,10.0,20.0]; mP1_e=[0.84,2.34,2.30,2.05,1.45,1.10,0.92,0.72,0.56]
mP2_s=[0.021,0.5,3.5,5.0,10.0,20.0]; mP2_e=[0.84,0.85,0.85,0.78,0.67,0.57]
fig,ax=plt.subplots(figsize=(8,6))
ax.plot(mP1_s,mP1_e,"-o",color="#c0392b",lw=1.6,ms=5,mec="k",mew=0.4,zorder=5,label="measured P1 (free swell A-B-C-D)")
ax.plot(mP2_s,mP2_e,"--s",color="#e08e0b",lw=1.6,ms=5,mec="k",mew=0.4,zorder=5,label="measured P2 (const. vol A-B'-C')")
ax.plot(S_sw,E_sw,"-^",color="#5fa8d3",ms=4,lw=1.4,zorder=3,label="model P1 free-swell (LE; capped e~1.0, small-strain)")
ax.plot(S_c,E_c,"-D",color="#7b3fa0",ms=4,lw=1.5,zorder=3,label="model P1 compression (MCC from B=2.34, $\lambda$=0.31 fit to D)")
ax.plot(S_p2,E_p2,"-o",color="#1f5fa8",ms=4,lw=1.7,zorder=4,label="model P2 (split, K=56183)")
for n,(s,e) in {"A":(0.021,0.84),"B":(0.021,2.34),"C":(3.24,1.10),"C'":(20.0,0.57),"D":(20.0,0.56)}.items():
    ax.annotate(n,(s,e),textcoords="offset points",xytext=(6,5),fontsize=11,fontweight="bold",zorder=6)
ax.text(0.012,0.34,
        "micro-only disjoining nS=1-n_l (current-porosity-split); K re-fit to Dixon\n"
        "P1: LE free-swell capped at e~1.0 (82% swell strain > small-strain limit,\n"
        "both E-soften and K-crank diverge); MCC from IMPOSED B=2.34,\n"
        "lambda recalibrated 0.077->0.31 to match measured D (Cc~0.59); endpoint e=0.56",
        fontsize=7.4,color="#1f4e79",bbox=dict(boxstyle="round,pad=0.3",fc="#eef5ff",ec="#1f5fa8",lw=0.8))
ax.set_xscale("log"); ax.set_xlim(0.01,60); ax.set_ylim(0.3,2.6)
ax.set_xlabel("$\\sigma_v$ [MPa]"); ax.set_ylabel("void ratio $e$ [-]")
ax.set_title("EPFL Task 3.3 MX-80 — model (current-porosity-split, K re-fit) vs Ferrari (2022)\nP1 = LE free-swell + MCC compression from imposed B; P2 confined",fontsize=9.2)
ax.grid(True,which="both",ls=":",alpha=0.4); ax.legend(fontsize=7.6,loc="upper right")
fig.tight_layout(); out=os.path.join(ROOT,"figures","epfl_cps_both_paths.png"); fig.savefig(out,dpi=150)
print("wrote",out)
print(f"P1 MCC compression from B=2.34: end e={E_c[-1]:.3f} @ {S_c[-1]:.1f} MPa (meas D 0.56)")
print(f"P2 end e={E_p2[-1]:.3f} @ {S_p2[-1]:.1f} MPa (meas C' 0.57)")
