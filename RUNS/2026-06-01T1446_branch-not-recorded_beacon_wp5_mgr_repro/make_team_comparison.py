#!/usr/bin/env python3
"""BEACON WP5.3 MGR — cross-team 1:1 comparison figures (Vinay 2026-06-01).
Native DSM (ours) vs Experiment (exact) vs UPC BExM (digitized) vs BGR-old (digitized).
Solid/filled = EXACT (tabulated or our VTU). Hatched/open = DIGITIZED from D5.6/D5.7 figures (approx)."""
import sys, glob, re
sys.path.insert(0, "results")
import numpy as np, meshio
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from team_comparison_data import EXPERIMENT, UPC, BGR_OLD, NATIVE_DSM, EXPERIMENT_PSW, UPC_PSW, BGR_OLD_PSW

ROOT="/Users/vinaykumar/git/ogs/beacon_wp5_mgr_repro_2026-06-01"; RHO_S=2700.0

def profile(vtudir, prefix, nbin=25):
    fs=sorted(glob.glob(f"{ROOT}/results/{vtudir}/{prefix}_*.vtu"),
              key=lambda f:float(re.search(r"_t_([0-9.]+)\.vtu",f).group(1)))
    m=meshio.read(fs[-1]); y=m.points[:,1]
    zc=np.linspace(0,0.10,nbin+1); zmid=0.5*(zc[:-1]+zc[1:]); rd=[]
    phi=np.asarray(m.point_data["porosity"])
    for i in range(nbin):
        s=(y>=zc[i])&(y<=zc[i+1])
        rd.append(RHO_S*(1-phi[s].mean())/1000 if s.any() else np.nan)
    return zmid*100, np.array(rd)   # cm, g/cm3

# our verified final profiles
z23, rd23 = profile("cal_base","cal_base")              # MGR23 pellet@bottom
z27, rd27 = profile("mgr27_predict","mgr27_predict_blk0")# MGR27 block@bottom

# ============================ FIG 1 — HERO: density-gap 1:1 ============================
fig,ax=plt.subplots(figsize=(8.4,5.2))
tests=["MGR23","MGR27"]; x=np.arange(len(tests)); w=0.2
teams=[("Experiment", EXPERIMENT, "k", False),
       ("UPC BExM (ref.)", UPC, "#2ca02c", True),
       ("BGR-old (χ=0)", BGR_OLD, "#7f7f7f", True),
       ("Native DSM (ours, χ=0)", None, "#1f77b4", False)]
for j,(name,data,col,hatch) in enumerate(teams):
    xs=x+(j-1.5)*w
    if name.startswith("Native"):
        gaps=[NATIVE_DSM[t]["base"]["gap"] for t in tests]
        maxk=[NATIVE_DSM[t]["maxk"]["gap"] for t in tests]
        ax.bar(xs,gaps,w,color=col,label=name,zorder=3)
        # perm-sensitivity range base->maxk (marker)
        for xi,g,mk in zip(xs,gaps,maxk):
            ax.plot([xi,xi],[g,mk],color="darkblue",lw=1.4,zorder=4)
            ax.plot(xi,mk,marker="_",color="darkblue",ms=11,mew=2,zorder=4)
    else:
        gaps=[data[t]["gap"] for t in tests]
        ax.bar(xs,gaps,w,color=col,label=name,zorder=3,
               hatch=("///" if hatch else None),
               edgecolor=("white" if hatch else col),alpha=(0.9 if not hatch else 1))
ax.axhline(0.30,ls="--",c="firebrick",lw=1.3,zorder=1)
ax.text(1.46,0.305,"initial gap 0.30 (1.60/1.30)",color="firebrick",fontsize=8,va="bottom",ha="right")
ax.set_xticks(x); ax.set_xticklabels(tests,fontsize=11)
ax.set_ylabel("final dry-density gap  block − pellet  [g/cm$^3$]",fontsize=10)
ax.set_title("BEACON MGR homogenisation — final density gap (lower = more homogenised)",fontsize=10.5)
ax.set_ylim(0,0.33)
ax.annotate("experiment\nhomogenises",xy=(1-0.30,0.02),xytext=(0.35,0.13),fontsize=8.5,
            ha="center",arrowprops=dict(arrowstyle="->",color="k",lw=0.8))
ax.annotate("χ=0 family fails\n(BGR-old & ours)",xy=(1+0.15,0.20),xytext=(1.42,0.255),fontsize=8.5,
            ha="center",color="#444",arrowprops=dict(arrowstyle="->",color="#444",lw=0.8))
h,l=ax.get_legend_handles_labels()
h.append(Patch(facecolor="white",edgecolor="gray",hatch="///",label="hatched = digitized (approx.)"))
ax.legend(handles=h,fontsize=8.2,loc="upper left",ncol=1,framealpha=0.95)
ax.grid(axis="y",alpha=0.3)
fig.text(0.5,0.005,"Exact: Experiment (D5.6 tables), Native DSM (our VTU). Digitized ±~0.02: UPC (D5.7 Figs 5.3-11/16), "
         "BGR-old (D5.6 Figs 3-138/142/132). Native DSM bar=base k$_0$; cap=block k$_0$×10$^4$.",
         fontsize=6.6,ha="center",color="#555")
fig.tight_layout(rect=[0,0.02,1,1]); fig.savefig(f"{ROOT}/figures/mgr_gap_1to1.png",dpi=155); plt.close(fig)

# ===================== FIG 2 — participating-team dry-density profiles =====================
fig,axs=plt.subplots(1,2,figsize=(10.2,5.3),sharey=True)
def panel(ax,test,z,rd,pel_bottom):
    # initial step (exact)
    if pel_bottom:  # MGR23: pellet 0-5 (1.30), block 5-10 (1.60)
        ax.plot([1.30,1.30],[0,5],"--",c="firebrick",lw=1.2)
        ax.plot([1.60,1.60],[5,10],"--",c="firebrick",lw=1.2,label="initial (1.30/1.60)")
    else:           # MGR27: block 0-5 (1.60), pellet 5-10 (1.30)
        ax.plot([1.60,1.60],[0,5],"--",c="firebrick",lw=1.2,label="initial (1.30/1.60)")
        ax.plot([1.30,1.30],[5,10],"--",c="firebrick",lw=1.2)
    ax.plot(rd,z,"-o",c="#1f77b4",ms=3,lw=1.8,label="Native DSM final (ours, VTU)",zorder=5)
    # experiment final layer averages (EXACT) as points at layer mid-depth
    e=EXPERIMENT[test]
    if pel_bottom: pz,bz=2.5,7.5
    else:          bz,pz=2.5,7.5
    ax.plot([e["pellet"],e["block"]],[pz,bz],"s",c="k",ms=9,label="Experiment final (exact, layer avg)",zorder=6)
    # UPC + BGR-old final (digitized) layer markers
    u=UPC[test]; b=BGR_OLD[test]
    ax.plot([u["pellet"],u["block"]],[pz,bz],"^",c="#2ca02c",ms=8,mfc="none",mew=1.6,label="UPC BExM (digitized)")
    ax.plot([b["pellet"],b["block"]],[pz,bz],"v",c="#7f7f7f",ms=8,mfc="none",mew=1.6,label="BGR-old (digitized)")
    ax.axhline(5,c="gray",lw=0.8,ls=":")
    lay="pellets / block" if pel_bottom else "block / pellets"
    ax.set_title(f"{test}   (hydrate z=0; {lay})",fontsize=10)
    ax.set_xlabel("dry density  [g/cm$^3$]",fontsize=10); ax.set_xlim(1.2,1.7); ax.grid(alpha=0.3)
panel(axs[0],"MGR23",z23,rd23,pel_bottom=True)
panel(axs[1],"MGR27",z27,rd27,pel_bottom=False)
axs[0].set_ylabel("distance to hydration surface  [cm]",fontsize=10); axs[0].set_ylim(0,10)
axs[0].legend(fontsize=7.6,loc="upper center",framealpha=0.95)
fig.suptitle("Native DSM as a participating team — final dry-density profile vs experiment, UPC, BGR-old",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.96]); fig.savefig(f"{ROOT}/figures/mgr_profile_team.png",dpi=155); plt.close(fig)

print("WROTE figures/mgr_gap_1to1.png  figures/mgr_profile_team.png")
print("native DSM MGR23 final profile pellet(z<5) mean=%.3f block(z>5) mean=%.3f"%(np.nanmean(rd23[z23<5]),np.nanmean(rd23[z23>5])))
print("native DSM MGR27 final profile block(z<5) mean=%.3f pellet(z>5) mean=%.3f"%(np.nanmean(rd27[z27<5]),np.nanmean(rd27[z27>5])))
