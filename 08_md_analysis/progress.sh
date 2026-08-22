#!/bin/bash
# =====================================================================
# progress.sh - live status of the dual-site docking and apo MD
# Run:  bash ~/Research/progress.sh          (single snapshot)
#       watch -n 30 bash ~/Research/progress.sh   (refresh every 30 s)
# =====================================================================
R=~/Research
D=$R/dual_site
APO=$R/04_Molecular_Dynamics/Apo/md_run

echo "════════════════════════════════════════════════════════════"
echo "  $(date '+%H:%M:%S')"
echo "════════════════════════════════════════════════════════════"

# ---------------------------------------------------------- docking
echo
echo "── DUAL-SITE DOCKING ──"
NPY=$(pgrep -fc "dual_site_dock" 2>/dev/null | head -1 || echo 0)
NVINA=$(pgrep -fc "vina --receptor" 2>/dev/null | head -1 || echo 0)
echo "  python processes : $NPY      vina processes : $NVINA"

if [ "$NPY" -eq 0 ]; then
  echo "  status : NOT RUNNING"
else
  MODE=$(ps -o args= -p "$(pgrep -f dual_site_dock | head -1)" 2>/dev/null \
         | grep -q -- "--test" && echo "TEST" || echo "FULL")
  START=$(ps -o lstart= -p "$(pgrep -f dual_site_dock | head -1)" 2>/dev/null)
  ELAPSED=$(ps -o etime= -p "$(pgrep -f dual_site_dock | head -1)" 2>/dev/null | tr -d ' ')
  CPU=$(ps -o %cpu= -p $(pgrep -f dual_site_dock | tr '\n' ',' | sed 's/,$//') \
        2>/dev/null | awk '{s+=$1} END{printf "%.0f", s}')
  echo "  mode    : $MODE"
  echo "  elapsed : $ELAPSED"
  echo "  total CPU across workers : ${CPU}%"
  if [ "${CPU%.*}" -lt 100 ] 2>/dev/null; then
    echo "  ⚠  low CPU - docking is being starved (likely by the MD job)"
  fi
fi

for f in "$D/dual_site_results_test.csv" "$D/dual_site_results.csv"; do
  if [ -f "$f" ]; then
    N=$(( $(wc -l < "$f") - 1 ))
    [ "$N" -lt 0 ] && N=0
    SZ=$(du -h "$f" | cut -f1)
    echo "  $(basename "$f") : $N rows ($SZ)"
    if [ "$N" -gt 0 ] && [ "$(basename "$f")" = "dual_site_results.csv" ]; then
      SECS=$(ps -o etimes= -p "$(pgrep -f dual_site_dock | head -1)" 2>/dev/null | tr -d ' ')
      if [ -n "$SECS" ] && [ "$SECS" -gt 0 ]; then
        awk -v n="$N" -v s="$SECS" 'BEGIN{
          r=n/(s/60);
          printf "     rate %.1f/min   ", r;
          if (r>0) printf "ETA %.1f h for 6319\n", (6319-n)/r/60; else print ""
        }'
      fi
    fi
  fi
done
echo "  note: results flush every 100 molecules, so a --test 50 run"
echo "        writes nothing until it finishes"

# ---------------------------------------------------------- apo MD
echo
echo "── APO MD ──"
if pgrep -f "gmx mdrun" >/dev/null; then
  if [ -f "$APO/md_prod.log" ]; then
    STEP=$(grep -E "^ *Step +Time" -A1 "$APO/md_prod.log" | tail -1 | awk '{print $1}')
    if [ -n "$STEP" ]; then
      awk -v s="$STEP" 'BEGIN{
        printf "  step %d / 25000000  (%.1f%%)  = %.1f ns of 50\n",
               s, 100*s/25000000, s*0.002/1000
      }'
    fi
  fi
  PERF=$(grep -E "Performance:" "$APO/md_prod.log" 2>/dev/null | tail -1 | awk '{print $2}')
  [ -n "$PERF" ] && echo "  performance : $PERF ns/day"
  echo "  status : RUNNING"
else
  echo "  status : NOT RUNNING"
  [ -f "$APO/md_prod.gro" ] && echo "  md_prod.gro exists -> COMPLETE" \
                            || echo "  no md_prod.gro -> stopped early, resume needed"
fi

# ---------------------------------------------------------- load
echo
echo "── SYSTEM ──"
echo "  load average :$(cut -d' ' -f1-3 /proc/loadavg | sed 's/^/ /')   (12 cores)"
free -h | awk 'NR==2{printf "  memory       : %s used of %s\n", $3, $2}'
echo "════════════════════════════════════════════════════════════"
