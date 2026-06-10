#!/bin/bash
# EPFL Task-3.3 stress-path suite re-run with the MAXWELL-CONJUGATE binary.
# (Vinay "run ... EPFL test (the stress path dependency)" 2026-06-02; chose Maxwell-conjugate only.)
#
# Physics generation note: inputs ported from beacon_wp3_epfl_repro_2026-05-31/model
# by stripping the retired <accumulate_swelling_contributions> tag (switch retired
# 2026-06-01; the full-Pi disjoining closure accumulates the swelling eigenstress by
# default). NO parameter changed — K=44200 (Dixon 2023 MX-80, rho_d=1500) and all
# vdW/micro params carry over. The binary additionally carries the Maxwell-conjugate
# micro-potential term (residual only; analytic Jacobian block not yet wired).
set +e
NEW=/Users/vinaykumar/git/ogs/beacon_wp3_epfl_stresspath_2026-06-02
M=$NEW/model
R=$NEW/results
OGS=/Users/vinaykumar/git/build/maxwell-conjugate-20260602/bin/ogs
S=$R/EPFL_MAXWELL_SUMMARY.txt
export OMP_NUM_THREADS=4
cd "$M" || exit 1
echo "=== EPFL Task-3.3 on maxwell-conjugate-20260602 (full-Pi closure + Maxwell term) ===" > "$S"
echo "binary: $OGS" >> "$S"
echo "started: $(date)" >> "$S"

for prj in beacon_t33_column_le \
           beacon_t33_path1_P1-3_dsm_mcc \
           beacon_t33_path2_P2-1_dsm_mcc \
           beacon_t33_path2_P2-1_swellingpressure_dsm_mcc \
           beacon_t33_path2_compression_homogeneous_mcc ; do
  od="$R/run_$prj"; rm -rf "$od"; mkdir -p "$od"
  t0=$SECONDS
  "$OGS" "$prj.prj" -o "$od" > "$od/run.log" 2>&1
  rc=$?
  el=$((SECONDS-t0))
  status=$(grep -oE "Simulation completed|terminated with error|integration failed with status|read 1 time\(s\) less" "$od/run.log" | tail -1)
  lastt=$(grep -oE "Time: [0-9.eE+-]+" "$od/run.log" | tail -1)
  nv=$(ls "$od"/*.vtu 2>/dev/null | wc -l | tr -d ' ')
  echo "" >> "$S"
  echo "### $prj : rc=$rc ${el}s | ${status:-NO-STATUS} | ${lastt:-no-time} | vtus=$nv" >> "$S"
done
echo "" >> "$S"
echo "=== DONE $(date) ===" >> "$S"
