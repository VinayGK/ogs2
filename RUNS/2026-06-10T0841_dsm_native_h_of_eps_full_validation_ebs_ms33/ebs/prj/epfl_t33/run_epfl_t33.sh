#!/bin/bash
# EPFL/Beacon Task-3.3 stress-path-dependency suite on the mc_20260608 binary.
# PRJs reused verbatim from beacon_wp3_epfl_stresspath_2026-06-02 (maxwell-conjugate port);
# they parse unchanged on mc_20260608 (dsm_consolidate_on_film) — NO config adapted.
# K=44200 J/kg (Dixon 2023 MX-80 block locus, rho_d=1500 log-interp; per-PRJ §12 header).
set +e
NEW=/Users/vinaykumar/git/ogs/validation_2026-06-09/epfl_t33
M=$NEW/model; R=$NEW/results
OGS=/Users/vinaykumar/git/build/mc_20260608/bin/ogs
S=$R/EPFL_T33_SUMMARY.txt
export OMP_NUM_THREADS=1
cd "$M" || exit 1
echo "=== EPFL Task-3.3 on mc_20260608 (dsm_consolidate_on_film) ===" > "$S"
echo "binary: $OGS" >> "$S"; echo "started: $(date)" >> "$S"
for prj in beacon_t33_column_le \
           beacon_t33_path1_P1-3_dsm_mcc \
           beacon_t33_path2_P2-1_dsm_mcc \
           beacon_t33_path2_P2-1_swellingpressure_dsm_mcc \
           beacon_t33_path2_compression_homogeneous_mcc ; do
  od="$R/run_$prj"; rm -rf "$od"; mkdir -p "$od"
  t0=$SECONDS
  "$OGS" "$prj.prj" -o "$od" > "$od/run.log" 2>&1
  rc=$?; el=$((SECONDS-t0))
  status=$(grep -oE "Simulation completed|terminated with error|integration failed" "$od/run.log" | tail -1)
  lastt=$(grep -oE "Time: [0-9.eE+-]+" "$od/run.log" | tail -1)
  nv=$(ls "$od"/*.vtu 2>/dev/null | wc -l | tr -d ' ')
  echo "" >> "$S"
  echo "### $prj : rc=$rc ${el}s | ${status:-NO-STATUS} | ${lastt:-no-time} | vtus=$nv" >> "$S"
done
echo "" >> "$S"; echo "=== DONE $(date) ===" >> "$S"
