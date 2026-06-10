#!/bin/bash
# =====================================================================
# BEACON Task 3.3 — COLUMN BVP test (Vinay 2026-06-01)
# Question: does a lab-dimension wetting-front column lift the phi_M->0 floor
# that the single element showed? (My claim: it distributes but does not lift it,
# because a sharp Pi over-homogenises every material point regardless of mesh.)
#
# Isolation: LE skeleton (not MCC) -- the question is HYDRAULIC (does macro
# survive the front); LE removes the MCC integrator wall so hydration runs to
# saturation and phi_M(z,t) is readable. Factorial isolation per Vinay's method.
#
# Geometry: axisymmetric column, lab oedometer R=17.5mm x H=12.5mm, 8x40 quad.
# Hydration: suction ramp 100->0 MPa on the BOTTOM face only -> vertical front.
# Confinement: u_r=0 both walls, u_z=0 base, top free (drained oedometer).
# Then: freeze state, read phi_M profile along the column axis.
#
# Writes a single verdict file: results/column/VERDICT.txt
# =====================================================================
set +e
MODEL=/Users/vinaykumar/git/ogs/beacon_wp3_epfl_repro_2026-05-31/model
RES=/Users/vinaykumar/git/ogs/beacon_wp3_epfl_repro_2026-05-31/results/column
OGS=/Users/vinaykumar/git/build/dsm-native-omp-release/bin/ogs
BIN=/Users/vinaykumar/git/build/dsm-native-omp-release/bin
V=$RES/VERDICT.txt
mkdir -p "$RES"
echo "=== COLUMN BVP run $(date) ===" > "$V"

cd "$MODEL" || { echo "no model dir" >> "$V"; exit 1; }

# 1. mesh (idempotent)
if [ ! -f col_r175_h125.vtu ]; then
  "$BIN/generateStructuredMesh" -e quad --lx 0.0175 --ly 0.0125 --nx 8 --ny 40 -o col_r175_h125.vtu >> "$V" 2>&1
fi
python3 make_column_boundaries.py >> "$V" 2>&1
echo "-- meshes --" >> "$V"
ls -1 col_r175_h125*.vtu >> "$V" 2>&1

# 2. run
echo "-- running ogs --" >> "$V"
"$OGS" beacon_t33_column_le.prj -o "$RES" >> "$RES/run.log" 2>&1
echo "ogs EXIT=$?" >> "$V"
grep -E "simulation succeeded|terminated with error|integration failed" "$RES/run.log" 2>/dev/null | tail -2 >> "$V"
grep -o "Time: [0-9.]*" "$RES/run.log" 2>/dev/null | tail -1 >> "$V"

# 3. phi_M profile along axis at the last available step
echo "-- phi_M(z) profile --" >> "$V"
python3 - >> "$V" 2>&1 <<'PY'
import glob, os, re
import numpy as np, meshio
RES="/Users/vinaykumar/git/ogs/beacon_wp3_epfl_repro_2026-05-31/results/column"
fs=glob.glob(RES+"/*.vtu")
if not fs:
    print("NO OUTPUT VTU"); raise SystemExit
def t(f):
    m=re.search(r"_t_([0-9.]+)\.vtu",os.path.basename(f)); return float(m.group(1)) if m else -1
fs=sorted(fs,key=t)
last=fs[-1]
m=meshio.read(last)
pts=m.points
def pd(n): return np.asarray(m.point_data[n]) if n in m.point_data else None
phi=pd("porosity"); nl=pd("micro_water_content"); sat=pd("saturation"); pr=pd("pressure")
y=pts[:,1]
# sample along the axis (x~0): bin by height
order=np.argsort(y)
print("file:",os.path.basename(last)," t=",t(last))
print(f"{'z[mm]':>7} {'S':>6} {'phi':>7} {'n_l':>7} {'phi_M':>8} {'suction[MPa]':>12}")
ys=np.unique(np.round(y,6))
for yy in ys[::4]:
    sel=np.abs(y-yy)<1e-7
    ph=float(np.mean(phi[sel])); n=float(np.mean(nl[sel]))
    phiM=(ph-n)/(1-n) if (1-n)!=0 else float('nan')
    s=float(np.mean(sat[sel])); p=float(np.mean(pr[sel]))
    suc=-p/1e6 if p<0 else 0.0
    print(f"{yy*1000:7.2f} {s:6.3f} {ph:7.4f} {n:7.4f} {phiM:8.4f} {suc:12.3f}")
# verdict numbers
phiM_all=(phi-nl)/(1-nl)
print("\nphi_M: min=%.4f max=%.4f mean=%.4f  (phi_M0 ~ 0.42)"%(phiM_all.min(),phiM_all.max(),phiM_all.mean()))
frac_collapsed=float(np.mean(phiM_all<0.02))
print("fraction of column with phi_M<0.02 (macro gone): %.1f%%"%(100*frac_collapsed))
print("\nVERDICT: if max phi_M ~0 everywhere -> floor NOT lifted (sharp-Pi over-homogenises, mesh-independent).")
print("         if phi_M survives behind the front -> floor IS lifted (gap recoverable by BVP).")
PY
echo "=== END ===" >> "$V"
