#!/bin/bash
# mgr27: single PRJ mgr27_le_calk0.prj (the case's only model; previously completed)
EBS=/Users/vinaykumar/git/ogs/RUNS/_INPROGRESS_full_validation/ebs
OGS=/Users/vinaykumar/git/build/h_of_eps_20260609/bin/ogs
M=/Users/vinaykumar/git/ogs/validation_2026-06-09/mgr27/model
export OMP_NUM_THREADS=3
cd "$M" || exit 1
od=$EBS/out_mgr27/run_mgr27_le_calk0; mkdir -p "$od"
echo "=== START mgr27_le_calk0 $(date)"; t0=$SECONDS
"$OGS" mgr27_le_calk0.prj -o "$od" > "$od/run.log" 2>&1
echo "=== END mgr27_le_calk0 rc=$? elapsed=$((SECONDS-t0))s $(date)"
tail -3 "$od/run.log"
