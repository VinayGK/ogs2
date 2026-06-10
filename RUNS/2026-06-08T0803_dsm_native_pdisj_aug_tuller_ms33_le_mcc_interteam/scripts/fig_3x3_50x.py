#!/usr/bin/env python3
"""MS33 3x3 inter-team panel — IV & VII from the 50x high-perm runs (/tmp/ms33_LE_50x),
I & III from the existing submission reduced CSVs. Authoritative style (vs-time, grey
team curves, BGR red/green). Reuses fig_modelI.py / fig_III_IV_VII.py loaders + the
summarize_runs.py centre-probe conventions. No invented numbers.
"""
import os, sys, glob, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meshio

RUN = "/Users/vinaykumar/git/ogs/ms33_pdisj_aug_tuller_2026-06-08"
SX  = "/tmp/ms33_LE_50x"
RED50 = "/tmp/ms33_50x_reduced"; os.makedirs(RED50, exist_ok=True)
sys.path.insert(0, RUN)
import fig_modelI as M1
import fig_III_IV_VII as M3
GREY="#cfcfcf"; RED="#cc2222"; GRN="#2a8a2a"; DD=M1.DD; DDCOL=M1.DDCOL

# ---- centre-probe reducer (verbatim conventions from summarize_runs.py) ----
def ts_of(fn):
    m=re.search(r"_ts_(\d+)_t_([0-9.]+)\.vtu$",fn); return (int(m.group(1)),float(m.group(2))) if m else (None,None)
def probe_idx(pts):
    zc=0.5*(pts[:,1].min()+pts[:,1].max()); return int(np.argmin((pts[:,0])**2+(pts[:,1]-zc)**2))
def _s(pd,n,i):
    if n not in pd: return np.nan
    a=pd[n]; return float(a[i] if a.ndim==1 else a[i,0])
def _v(pd,n,i):
    if n not in pd: return None
    a=pd[n]; return a[i]
def reduce_run(outdir, prefix, label):
    vtus=sorted(glob.glob(f"{outdir}/{prefix}_ts_*.vtu"), key=lambda f: ts_of(f)[0])
    rows=[]
    for v in vtus:
        _,t=ts_of(v); m=meshio.read(v); pd_=m.point_data; pts=m.points; i=probe_idx(pts)
        sig=_v(pd_,"sigma",i); phi=_s(pd_,"porosity",i)
        rows.append(dict(time_days=t/86400.0,
            suction_MPa=-_s(pd_,"pressure",i)/1e6,
            mean_stress_MPa=(-(sig[0]+sig[1]+sig[2])/3/1e6) if sig is not None else np.nan,
            axial_stress_MPa=(-sig[1]/1e6) if sig is not None else np.nan,
            saturation=_s(pd_,"saturation",i), porosity=phi,
            void_ratio=(phi/(1-phi)) if phi==phi else np.nan,
            dry_density=_s(pd_,"dry_density_solid",i)))
    df=pd.DataFrame(rows); df.to_csv(f"{RED50}/{label}_history.csv",index=False)
    print(f"reduced {label}: n={len(df)}  e_end={df['void_ratio'].iloc[-1]:.3f}  p_end={df['mean_stress_MPa'].iloc[-1]:.3f}")
    return df

IV50  = reduce_run(f"{SX}/ModelIV/out",  "ms33_modelIV_pellets",      "LE_IV_50x")
VII50 = reduce_run(f"{SX}/ModelVII/out", "ms33_modelVII_freeswelling","LE_VII_50x")

def overlay_team_time(ax, model, needle, loc):
    first=True
    for team,path in M3.team_files(model).items():
        df=M3.load_master(path,model); r=M3.central_series(df,needle,loc=loc)
        if r: ax.plot(r[0],r[1],"-",color=GREY,lw=0.9,label="other teams" if first else None); first=False

fig,AX=plt.subplots(3,3,figsize=(17,14.5))
fig.suptitle("MS33 — BGR DSM vs EURAD-2 teams  (Models III, IV & VII = 50$\\times$ high-permeability runs, equilibrated)",
             fontsize=14, fontweight="bold", y=0.992)

# (1) swelling pressure vs dry density (Model I, existing)
ax=AX[0,0]; teamsI=M1.team_files_modelI()
rr=np.linspace(1.30,1.90,100); ax.plot(rr,np.exp(6.77*rr-9.07),"-",color="#888",lw=1.3,label="Villar/Lloret Eq.(7)")
first=True; dixon={}
for team,path in teamsI.items():
    try: t=M1.team_modelI(path)
    except Exception: t=None
    if not t: continue
    xs=[DD[k] for k,b in t["blocks"].items() if len(b["mean"])]; ys=[b["mean"][-1] for k,b in t["blocks"].items() if len(b["mean"])]
    if xs: ax.plot(xs,ys,"-o",color="#bbb",lw=0.8,ms=3,label="other teams" if first else None); first=False
    if t["dixon"]:
        for dd_,sp_ in zip(t["dixon"]["dd"],t["dixon"]["sp"]): dixon[round(float(dd_),3)]=float(sp_)
if dixon:
    dx=sorted(dixon); ax.plot(dx,[dixon[x] for x in dx],"s",color="#444",ms=8,mfc="#ccc",mec="#333",label="Dixon (2023) exp.")
le={k:M1.bgr_modelI(f"LE_I_dd{k}") for k in DD}
xs=[DD[k] for k in DD if le[k] is not None]; ys=[M1.saturated_endpoint(le[k])[0] for k in DD if le[k] is not None]
ax.plot(xs,ys,"-o",color=RED,lw=2.2,ms=9,label="BGR DSM-LE")
for x,y in zip(xs,ys): ax.annotate(f"{y:.1f}",(x,y),textcoords="offset points",xytext=(6,6),fontsize=8,color=RED)
ax.set_yscale("log"); ax.set_xlabel("dry density $\\rho_d$ [g/cm³]"); ax.set_ylabel("swelling pressure $P_s$ [MPa]")
ax.set_title("(1) Model I — swelling pressure vs $\\rho_d$",fontsize=11,fontweight="bold")
ax.grid(True,which="both",ls=":",color="#ddd"); ax.legend(fontsize=7,loc="upper left")

# (2) suction vs mean stress (Model I, existing)
ax=AX[0,1]; first=True
for team,path in teamsI.items():
    try: t=M1.team_modelI(path)
    except Exception: t=None
    if not t: continue
    for k,b in t["blocks"].items():
        if len(b["mean"]): ax.plot(b["mean"],b["suction"],"-",color=GREY,lw=0.8,label="other teams" if first else None); first=False
for k in ["1400","1600","1800"]:
    d=M1.bgr_modelI(f"LE_I_dd{k}")
    if d is not None: ax.plot(d["mean_stress_MPa"],d["suction_MPa"],"-",color=DDCOL[k],lw=2.0,label=f"BGR LE $\\rho_d$={DD[k]}")
ax.set_yscale("log"); ax.set_ylim(0.1,110); ax.set_xlim(0,60)
ax.set_xlabel("mean stress $p$ [MPa]"); ax.set_ylabel("suction $s$ [MPa]")
ax.set_title("(2) Model I — suction–mean-stress path",fontsize=11,fontweight="bold")
ax.grid(True,which="both",ls=":",color="#ddd"); ax.legend(fontsize=6.5,ncol=2,loc="upper right")

# (3) Model III mean stress vs time — 50x equilibrated
ax=AX[0,2]; overlay_team_time(ax,"III","mean stress","central")
III50=pd.read_csv("/tmp/ms33_III_50x/LE_III_50x_history.csv")
ax.plot(III50["t"],III50["p"],"-",color=RED,lw=2.2,label="BGR LE 50$\\times$ (centre)")
ax.set_xlim(0,200); ax.set_xlabel("time [d]"); ax.set_ylabel("mean stress $p$ [MPa]")
ax.set_title("(3) Model III — mean stress vs time (50$\\times$)",fontsize=11,fontweight="bold")
ax.grid(True,ls=":",color="#ddd"); ax.legend(fontsize=7)

# (4) Model IV mean stress vs time — 50x
ax=AX[1,0]; overlay_team_time(ax,"IV","mean stress","central")
ax.plot(IV50["time_days"],IV50["mean_stress_MPa"],"-",color=RED,lw=2.2,label="BGR LE 50$\\times$ (clay centre)")
ax.set_xlim(0,200); ax.set_xlabel("time [d]"); ax.set_ylabel("mean stress $p$ [MPa]")
ax.set_title("(4) Model IV — mean stress vs time (50$\\times$)",fontsize=11,fontweight="bold")
ax.grid(True,ls=":",color="#ddd"); ax.legend(fontsize=7)

# (5) Model III gap closure vs time — 50x equilibrated
ax=AX[1,1]; overlay_team_time(ax,"III","gap closure","")
outdir="/tmp/ms33_III_50x/out"
if os.path.isdir(outdir):
    T,Y=M3.bgr_probe_series(outdir,"ms33_modelIII_gap2mm","displacement",(0.023,0.040),comp=0)
    if len(T): ax.plot(T,np.abs(Y)*1e3,"-",color=RED,lw=2.0,label="BGR LE 50$\\times$ $|u_r|$ @gap")
ax.set_xlim(0,200); ax.set_xlabel("time [d]"); ax.set_ylabel("gap closure / $|u_r|$ [mm]")
ax.set_title("(5) Model III — gap closure vs time (50$\\times$, r=23 mm)",fontsize=11,fontweight="bold")
ax.grid(True,ls=":",color="#ddd"); ax.legend(fontsize=7)

# (6) Model IV dry density vs time — 50x (clay + pellet probes)
ax=AX[1,2]; overlay_team_time(ax,"IV","dry density","")
for zlab,rz,col,ls in [("clay (top, z=52mm)",(0.0,0.0525),RED,"-"),("pellet (bottom, z=18mm)",(0.0,0.0175),"#cc8822","--")]:
    T,Y=M3.bgr_probe_series(f"{SX}/ModelIV/out","ms33_modelIV_pellets","dry_density_solid",rz)
    if len(T): ax.plot(T,Y/1000.0,ls,color=col,lw=2.0,label=f"BGR LE 50$\\times$ {zlab}")
ax.set_xlim(0,200); ax.set_xlabel("time [d]"); ax.set_ylabel("dry density $\\rho_d$ [g/cm³]")
ax.set_title("(6) Model IV — dry-density evolution (50$\\times$)",fontsize=11,fontweight="bold")
ax.grid(True,ls=":",color="#ddd"); ax.legend(fontsize=7)

# (7) Model VII loading cycle e vs axial — 50x
ax=AX[2,0]; first=True
for team,path in M3.team_files("VII").items():
    df=M3.load_master(path,"VII"); re_=M3.central_series(df,"void ratio",loc=""); ra_=M3.central_series(df,"axial stress",loc="central")
    if re_ and ra_:
        a_on_e=np.interp(re_[0],ra_[0],ra_[1]); ax.plot(a_on_e,re_[1],"-",color=GREY,lw=0.9,label="other teams" if first else None); first=False
ax.plot(VII50["axial_stress_MPa"],VII50["void_ratio"],"-",color=RED,lw=2.2,label="BGR LE 50$\\times$ (centre)")
ax.set_xlabel("axial stress $\\sigma_a$ [MPa]"); ax.set_ylabel("void ratio $e$ [–]")
ax.set_title("(7) Model VII — loading cycle ($e$–$\\sigma_a$, 50$\\times$)",fontsize=11,fontweight="bold")
ax.grid(True,ls=":",color="#ddd"); ax.legend(fontsize=7)

# (8) Model VII void ratio vs time — 50x
ax=AX[2,1]; overlay_team_time(ax,"VII","void ratio","")
ax.plot(VII50["time_days"],VII50["void_ratio"],"-",color=RED,lw=2.2,label="BGR LE 50$\\times$ (centre)")
ax.set_xlim(0,240); ax.set_xlabel("time [d]"); ax.set_ylabel("void ratio $e$ [–]")
ax.set_title("(8) Model VII — void ratio vs time (50$\\times$)",fontsize=11,fontweight="bold")
ax.grid(True,ls=":",color="#ddd"); ax.legend(fontsize=7)

# (9) Model VII axial stress vs time — 50x
ax=AX[2,2]; overlay_team_time(ax,"VII","axial stress","central")
ax.plot(VII50["time_days"],VII50["axial_stress_MPa"],"-",color=RED,lw=2.2,label="BGR LE 50$\\times$ (centre)")
ax.set_xlim(0,240); ax.set_xlabel("time [d]"); ax.set_ylabel("axial stress $\\sigma_a$ [MPa]")
ax.set_title("(9) Model VII — axial stress vs time (50$\\times$)",fontsize=11,fontweight="bold")
ax.grid(True,ls=":",color="#ddd"); ax.legend(fontsize=7)

fig.text(0.5,0.004,"BGR: Model I = submission reduced CSVs (g9bd1c32b);  Model III = 50× high-perm (g9bd1c32b, /tmp/ms33_III_50x);  Models IV/VII = 50× high-perm (g295dc649, /tmp/ms33_LE_50x).  Teams: EURAD-2 Family-A sheets.  2026-06-08",
         ha="center",fontsize=6.2,color="#555")
fig.tight_layout(rect=[0,0.012,1,0.978])
out=f"{RUN}/figures/ms33_3x3_50x_equil.png"; fig.savefig(out,dpi=200); print("wrote",out)
