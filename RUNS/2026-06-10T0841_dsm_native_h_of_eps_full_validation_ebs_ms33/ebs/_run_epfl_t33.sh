#!/bin/bash
# epfl_t33: per its own run_epfl_t33.sh — 5 PRJs sequential; only binary+outdir+OMP adapted
EBS=/Users/vinaykumar/git/ogs/RUNS/_INPROGRESS_full_validation/ebs
OGS=/Users/vinaykumar/git/build/h_of_eps_20260609/bin/ogs
M=/Users/vinaykumar/git/ogs/validation_2026-06-09/epfl_t33/model
export OMP_NUM_THREADS=3
cd "$M" || exit 1
for prj in beacon_t33_column_le \
           beacon_t33_path1_P1-3_dsm_mcc \
           beacon_t33_path2_P2-1_dsm_mcc \
           beacon_t33_path2_P2-1_swellingpressure_dsm_mcc \
           beacon_t33_path2_compression_homogeneous_mcc ; do
  od=$EBS/out_epfl_t33/run_$prj; mkdir -p "$od"
  echo "=== START $prj $(date)"; t0=$SECONDS
  "$OGS" "$prj.prj" -o "$od" > "$od/run.log" 2>&1
  echo "=== END $prj rc=$? elapsed=$((SECONDS-t0))s $(date)"
  tail -3 "$od/run.log"
done
