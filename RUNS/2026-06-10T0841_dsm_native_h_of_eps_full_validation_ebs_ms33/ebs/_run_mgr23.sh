#!/bin/bash
# mgr23: no runner script; prior records show all three PRJs attempted (calibrated/robust/bicg)
EBS=/Users/vinaykumar/git/ogs/RUNS/_INPROGRESS_full_validation/ebs
OGS=/Users/vinaykumar/git/build/h_of_eps_20260609/bin/ogs
M=/Users/vinaykumar/git/ogs/validation_2026-06-09/mgr23/model
export OMP_NUM_THREADS=3
cd "$M" || exit 1
for prj in mgr23_column_calibrated mgr23_column_robust mgr23_column_bicg; do
  od=$EBS/out_mgr23/run_$prj; mkdir -p "$od"
  echo "=== START $prj $(date)"; t0=$SECONDS
  "$OGS" "$prj.prj" -o "$od" > "$od/run.log" 2>&1
  echo "=== END $prj rc=$? elapsed=$((SECONDS-t0))s $(date)"
  tail -3 "$od/run.log"
done
