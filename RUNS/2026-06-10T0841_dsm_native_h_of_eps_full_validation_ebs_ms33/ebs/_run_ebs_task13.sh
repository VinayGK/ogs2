#!/bin/bash
# ebs_task13: Stage-1 GAP tests, 0.5 MPa (12a) and 6 MPa (12b) suction — both PRJs in stg1/
EBS=/Users/vinaykumar/git/ogs/RUNS/_INPROGRESS_full_validation/ebs
OGS=/Users/vinaykumar/git/build/h_of_eps_20260609/bin/ogs
D=/Users/vinaykumar/git/ogs/validation_2026-06-09/ebs_task13/stg1
export OMP_NUM_THREADS=3
cd "$D" || exit 1
for prj in 12a_t-13_MCC_0.5MPa 12b_t-13_MCC_6MPa; do
  od=$EBS/out_ebs_task13/run_$prj; mkdir -p "$od"
  echo "=== START $prj $(date)"; t0=$SECONDS
  "$OGS" "$prj.prj" -o "$od" > "$od/run.log" 2>&1
  echo "=== END $prj rc=$? elapsed=$((SECONDS-t0))s $(date)"
  tail -3 "$od/run.log"
done
