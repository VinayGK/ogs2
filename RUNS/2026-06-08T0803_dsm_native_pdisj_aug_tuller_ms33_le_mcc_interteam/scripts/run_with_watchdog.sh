#!/bin/bash
# Grinder-aware OGS runner.
# Usage: run_with_watchdog.sh PRJ OUTDIR [MAX_WALL_S] [STALL_S] [OMP]
# Kills the sim if wall-clock exceeds MAX_WALL_S (grinder) or the log file
# stops growing for STALL_S (hung). Otherwise lets it finish naturally.
OGS=/Users/vinaykumar/git/build/ufz_integration_20260602/bin/ogs
PRJ="$1"; OUT="$2"; MAXW="${3:-900}"; STALL="${4:-180}"; OMP="${5:-2}"
export HOME=/tmp MPLCONFIGDIR=/tmp/mplcache XDG_CACHE_HOME=/tmp OMP_NUM_THREADS="$OMP"
mkdir -p "$OUT"
echo "START $(date +%T) prj=$(basename "$PRJ") maxwall=${MAXW}s stall=${STALL}s omp=$OMP"
"$OGS" "$PRJ" -o "$OUT" > "$OUT/run.log" 2>&1 &
PID=$!
start=$(date +%s); last_mtime=$start
while kill -0 $PID 2>/dev/null; do
  sleep 15
  now=$(date +%s)
  lm=$(stat -f %m "$OUT/run.log" 2>/dev/null || echo $now)
  wall=$((now-start)); idle=$((now-lm))
  cnt=$(ls "$OUT"/*.vtu 2>/dev/null | wc -l | tr -d ' ')
  simt=$(grep -oE 't = ?[0-9.eE+]+ s' "$OUT/run.log" 2>/dev/null | tail -1)
  echo "  ...wall=${wall}s idle=${idle}s vtus=$cnt last=$simt"
  if [ $wall -ge $MAXW ]; then echo "WATCHDOG_KILL wall ${wall}s>=${MAXW}s vtus=$cnt"; kill -9 $PID 2>/dev/null; wait $PID 2>/dev/null; echo "KILLED_WALL"; exit 0; fi
  if [ $idle -ge $STALL ]; then echo "WATCHDOG_KILL hung idle ${idle}s>=${STALL}s vtus=$cnt"; kill -9 $PID 2>/dev/null; wait $PID 2>/dev/null; echo "KILLED_STALL"; exit 0; fi
done
wait $PID; rc=$?
echo "OGS_EXIT $rc vtus=$(ls "$OUT"/*.vtu 2>/dev/null|wc -l) $(grep -c 'Simulation completed' "$OUT/run.log" 2>/dev/null|tr -d ' ')_completed"
