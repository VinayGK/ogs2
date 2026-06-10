#!/usr/bin/env python3
"""Emit figures/PROVENANCE.md: every curve -> exact source + sample values."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import fig_modelI as M1
import fig_III_IV_VII as M3

RUN = "/Users/vinaykumar/git/ogs/ms33_pdisj_aug_tuller_2026-06-08"
L = []
def w(s=""): L.append(s)

w("# Figure provenance — every plotted curve and its source")
w()
w("BGR curves are reductions of solver output (no fitting). Team curves are read")
w("verbatim from the EURAD-2 data-collection `.xlsx` (Family-A master sheet, columns")
w("located by header-string match per `TEAM_DATA_MAP.md`). Sample = (first, last) (x,y).")
w()

# ---- Model I ----
w("## Model I  (`modelI_*` figures)")
w("BGR source: `reduced/{LE,MCC}_I_dd*_history.csv` (center probe; suction=-p/1e6, p=-(σxx+σyy+σzz)/3/1e6).")
for k in ["1400","1600","1800"]:
    d = M1.bgr_modelI(f"LE_I_dd{k}")
    if d is not None:
        p,sat,td = M1.saturated_endpoint(d)
        w(f"- BGR LE ρ_d={M1.DD[k]}: n={len(d)}, saturated P_s={p:.3f} MPa @ t={td:.0f} d (Sl={sat:.3f})")
for k in ["1400","1600"]:
    d = M1.bgr_modelI(f"MCC_I_dd{k}")
    if d is not None:
        p,sat,td = M1.saturated_endpoint(d)
        w(f"- BGR MCC ρ_d={M1.DD[k]}: n={len(d)}, saturated P_s={p:.3f} MPa @ t={td:.0f} d")
w()
w("Team Model I files read (saturated mean-stress endpoint + path + perm):")
for team, path in M1.team_files_modelI().items():
    try: t = M1.team_modelI(path)
    except Exception as e: w(f"- {team}: READ ERROR {e}"); continue
    if not t: w(f"- {team}: no master sheet"); continue
    eps = {k: (round(float(b['mean'][-1]),2) if len(b['mean']) else None) for k,b in t['blocks'].items()}
    dix = "yes" if t['dixon'] else "no"
    w(f"- **{team}** (`{os.path.basename(path)}`): sat P_s by ρ_d {eps}; Dixon block: {dix}")
w()

# ---- III / IV / VII ----
def dump(model, quants):
    w(f"## Model {model}  (`model{model}_interteam.png`)")
    for team, path in M3.team_files(model).items():
        df = M3.load_master(path, model)
        line = [f"**{team}** (`{os.path.basename(path)}`):"]
        for q, loc in quants:
            r = M3.central_series(df, q, loc=loc)
            if r:
                line.append(f"{q}: n={len(r[0])} ({r[1][0]:.3g}→{r[1][-1]:.3g})")
            else:
                line.append(f"{q}: —")
        w("- " + "  |  ".join(line))
    w()

dump("III", [("mean stress","central"),("gap closure","")])
w("BGR III: mean stress from `reduced/{LE,MCC}_III_history.csv` (center); gap closure = "
  "**proxy** |u_r| at node nearest (r=0.023, z=0.040) m from VTU `displacement[0]`.")
w()
dump("IV", [("mean stress","central"),("dry density","")])
w("BGR IV: per-material pellet K (clay 85312.6 / pellet 13064 J/kg via <medium id=1>, Dixon 0.350 MPa). "
  "Mean stress center CSV (out_perK); dry density two zones from VTU `dry_density_solid` at the on-axis "
  "zone centroids (0,0.0525) clay [MaterialID 0, top] and (0,0.0175) pellet [MaterialID 1, bottom], ÷1000 kg/m³→g/cm³.")
w()
dump("VII", [("void ratio",""),("axial stress","central")])
w("BGR VII: void ratio e=φ/(1−φ) and axial stress −σyy from `reduced/LE_VII_history.csv` (center).")
w()
w("### Data gaps (team quantity absent / non-standard sheet → silently omitted from plot)")
w("- Model III gap closure: only CTU-CU exposes a `Gap closure (mm)` column; ICL/UCLM/ULIEGE omit it.")
w("- Model VII: EPFL master sheet not in the standard time-series template (read returns none).")

open(f"{RUN}/figures/PROVENANCE.md","w").write("\n".join(L)+"\n")
print("wrote figures/PROVENANCE.md")
