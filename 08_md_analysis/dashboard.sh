#!/bin/bash
# =====================================================================
# dashboard.sh - combined progress for everything currently running
#
# Usage:
#   bash ~/Research/dashboard.sh              single snapshot
#   watch -n 60 bash ~/Research/dashboard.sh  refresh every 60s (Ctrl+C to exit)
# =====================================================================
R=~/Research

bar() {
  # bar PERCENT WIDTH
  local pct=$1 width=${2:-30}
  local filled=$(awk -v p="$pct" -v w="$width" 'BEGIN{printf "%d", p/100*w}')
  local empty=$((width - filled))
  printf "["
  [ "$filled" -gt 0 ] && printf "%${filled}s" | tr ' ' '#'
  [ "$empty" -gt 0 ] && printf "%${empty}s" | tr ' ' '-'
  printf "] %5.1f%%" "$pct"
}

echo "════════════════════════════════════════════════════════════════"
echo "  DASHBOARD   $(date '+%Y-%m-%d %H:%M:%S')"
echo "════════════════════════════════════════════════════════════════"

# ---------------------------------------------------------- MD Stage 2
echo
echo "── MD STAGE 2 (replicates 2 and 3) ──────────────────────────────"
LEADS=("LOA22-B1" "LOA22-B2" "LOA22-B3")
REPS=(2 3)
TOTAL_STAGES=$(( ${#LEADS[@]} * ${#REPS[@]} ))   # 6 lead/replicate pairs
DONE_STAGES=0
CURRENT_LABEL=""; CURRENT_PCT=0

for REP in "${REPS[@]}"; do
  for LEAD in "${LEADS[@]}"; do
    RUN=$R/05_MD_v2/$LEAD/rep$REP
    LABEL="$LEAD rep$REP"
    if [ -f "$RUN/md_prod.gro" ]; then
      echo "  $LABEL   production  COMPLETE   100.0%"
      DONE_STAGES=$((DONE_STAGES + 1))
    elif [ -f "$RUN/md_prod.log" ]; then
      STEP=$(grep -A1 "^ *Step *Time" "$RUN/md_prod.log" 2>/dev/null | tail -1 | awk '{print $1}')
      if [ -n "$STEP" ] && [ "$STEP" -gt 0 ] 2>/dev/null; then
        PCT=$(awk -v s="$STEP" 'BEGIN{printf "%.1f", 100*s/50000000}')
        NS=$(awk -v s="$STEP" 'BEGIN{printf "%.1f", s*0.002/1000}')
        PERF=$(grep "Performance:" "$RUN/md_prod.log" 2>/dev/null | tail -1 | awk '{print $2}')
        echo "  $LABEL   production  $(bar "$PCT" 20)  ${NS} ns/100${PERF:+  (${PERF} ns/day)}"
        CURRENT_LABEL="$LABEL production"; CURRENT_PCT="$PCT"
        DONE_STAGES=$(awk -v d="$DONE_STAGES" -v p="$PCT" 'BEGIN{print d+p/100}')
      else
        echo "  $LABEL   production  starting..."
        CURRENT_LABEL="$LABEL production (starting)"; CURRENT_PCT=0
      fi
    elif [ -f "$RUN/npt.log" ] && ! [ -f "$RUN/npt.gro" ]; then
      echo "  $LABEL   NPT equilibrating..."
      CURRENT_LABEL="$LABEL NPT"; CURRENT_PCT=0
      DONE_STAGES=$(awk -v d="$DONE_STAGES" 'BEGIN{print d+0.15}')
    elif [ -f "$RUN/nvt.log" ] && ! [ -f "$RUN/nvt.gro" ]; then
      echo "  $LABEL   NVT equilibrating..."
      CURRENT_LABEL="$LABEL NVT"; CURRENT_PCT=0
      DONE_STAGES=$(awk -v d="$DONE_STAGES" 'BEGIN{print d+0.03}')
    else
      echo "  $LABEL   not started"
    fi
  done
done

OVERALL=$(awk -v d="$DONE_STAGES" -v t="$TOTAL_STAGES" 'BEGIN{printf "%.1f", 100*d/t}')
echo
echo "  OVERALL STAGE 2:  $(bar "$OVERALL" 30)"
if [ -n "$CURRENT_LABEL" ]; then
  echo "  currently on: $CURRENT_LABEL"
fi

# ---------------------------------------------------------- Arm C eval
echo
echo "── ARM C DUAL-SITE EVALUATION ────────────────────────────────────"
LOG=$R/armC_eval.log
RESULTS=$R/arm_c/armC_dualsite_results.csv
if [ -f "$LOG" ]; then
  N=$(grep -oE '^[0-9]+/[0-9]+' "$LOG" | tail -1 | cut -d/ -f1)
  TOT=$(grep -oE '^[0-9]+/[0-9]+' "$LOG" | tail -1 | cut -d/ -f2)
  if [ -f "$RESULTS" ]; then
    N=$(( $(wc -l < "$RESULTS") - 1 ))
  fi
  if [ -n "$N" ] && [ -n "$TOT" ] && [ "$TOT" -gt 0 ] 2>/dev/null; then
    PCT=$(awk -v n="$N" -v t="$TOT" 'BEGIN{printf "%.1f", 100*n/t}')
    echo "  $(bar "$PCT" 30)   $N / $TOT molecules"

    START=$(stat -c %Y "$LOG" 2>/dev/null)
    NOW=$(date +%s)
    if [ -n "$START" ] && [ "$N" -gt 0 ]; then
      ELAPSED_MIN=$(( (NOW - START) / 60 ))
      if [ "$ELAPSED_MIN" -gt 0 ]; then
        RATE=$(awk -v n="$N" -v m="$ELAPSED_MIN" 'BEGIN{printf "%.1f", n/m}')
        REMAIN=$(awk -v n="$N" -v t="$TOT" -v r="$RATE" 'BEGIN{
          if (r>0) printf "%.0f", (t-n)/r; else print "?"}')
        echo "  rate: ${RATE} molecules/min   ETA: ~${REMAIN} min remaining"
      fi
    fi
  else
    echo "  waiting for first progress line..."
  fi
  if grep -q "THREE-ARM COMPARISON" "$LOG" 2>/dev/null; then
    echo "  >>> COMPLETE - comparison table is in the log <<<"
  fi
else
  echo "  not started"
fi

# ---------------------------------------------------------- system
echo
echo "── SYSTEM ─────────────────────────────────────────────────────────"
echo "  load average:$(cut -d' ' -f1-3 /proc/loadavg | sed 's/^/ /')   (12 cores)"
free -h | awk 'NR==2{printf "  memory: %s used of %s\n", $3, $2}'
NGMX=$(pgrep -fc "gmx mdrun" 2>/dev/null)
NVINA=$(pgrep -fc "vina --receptor" 2>/dev/null)
echo "  processes: gmx mdrun=${NGMX:-0}   vina=${NVINA:-0}"

echo "════════════════════════════════════════════════════════════════"
