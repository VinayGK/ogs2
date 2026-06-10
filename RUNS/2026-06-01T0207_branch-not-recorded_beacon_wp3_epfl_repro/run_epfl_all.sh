#!/bin/bash
# Re-run all EPFL Task-3.3 PRJs with the lambda=9 porosity power-law permeability
# (Vinay "run the epfl" 2026-06-01). One result subdir per PRJ; summary -> results/EPFL_RERUN_SUMMARY.txt
set +e
M=/Users/vinaykumar/git/ogs/beacon_wp3_epfl_repro_2026-05-31/model
R=/Users/vinaykumar/git/ogs/beacon_wp3_epfl_repro_2026-05-31/results
OGS=/Users/vinaykumar/git/build/dsm-native-omp-release/bin/ogs
S=$R/EPFL_RERUN_SUMMARY.txt
cd "$M" || exit 1
echo "=== EPFL Task-3.3 re-run (lambda=9 power-law perm) $(date) ===" > "$S"

for prj in beacon_t33_path1_P1-3_dsm_mcc \
           beacon_t33_path2_P2-1_dsm_mcc \
           beacon_t33_path2_P2-1_swellingpressure_dsm_mcc \
           beacon_t33_path2_compression_homogeneous_mcc \
           beacon_t33_column_le ; do
  od="$R/rerun_$prj"
  rm -rf "$od"; mkdir -p "$od"
  "$OGS" "$prj.prj" -o "$od" > "$od/run.log" 2>&1
  rc=$?
  status=$(grep -oE "Simulation completed|terminated with error|integration failed with status -1" "$od/run.log" 2>/dev/null | tail -1)
  lastt=$(grep -o "Time: [0-9.]*" "$od/run.log" 2>/dev/null | tail -1)
  nv=$(ls "$od"/*.vtu 2>/dev/null | wc -l | tr -d ' ')
  echo "" >> "$S"
  echo "### $prj : rc=$rc | $status | $lastt | vtus=$nv" >> "$S"
done

echo "" >> "$S"
echo "=== COLUMN phi_M(z) at final step (the headline: does throttled k preserve macro?) ===" >> "$S"
python3 - >> "$S" 2>&1 <<'PY'
import glob,re,numpy as np,meshio
od="/Users/vinaykumar/git/ogs/beacon_wp3_epfl_repro_2026-05-31/results/rerun_beacon_t33_column_le"
fs=sorted(glob.glob(od+"/*.vtu"),key=lambda f:float(re.search(r'_t_([0-9.]+)\.vtu',f).group(1)) if re.search(r'_t_([0-9.]+)\.vtu',f) else -1)
if not fs:
    print("  no column vtus"); raise SystemExit
m=meshio.read(fs[-1]); y=m.points[:,1]
phi=np.asarray(m.point_data["porosity"]); nl=np.asarray(m.point_data["micro_water_content"])
sat=np.asarray(m.point_data["saturation"])
phiM=(phi-nl)/(1-nl)
print(f"  final file: {fs[-1].split('/')[-1]}")
print(f"  {'z[mm]':>6} {'S':>6} {'phi':>7} {'n_l':>7} {'phi_M':>8}")
for zz in sorted(set(np.round(y,5)))[::4]:
    s=np.abs(y-zz)<1e-7
    print(f"  {zz*1000:6.2f} {sat[s].mean():6.3f} {phi[s].mean():7.4f} {nl[s].mean():7.4f} {phiM[s].mean():8.4f}")
print(f"\n  phi_M overall: min={phiM.min():.4f} max={phiM.max():.4f} mean={phiM.mean():.4f} (phi_M0~0.42)")
print(f"  PREVIOUS (KozenyCarman, unlimited supply): phi_M -> 0.008 mean (floor held)")
print(f"  => if phi_M now SURVIVES, the k(phi_M)^9 throttle gates the floor (same as MGR finding)")
PY
echo "=== END ===" >> "$S"
