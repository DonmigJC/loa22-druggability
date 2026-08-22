#!/bin/bash
for L in LOA22-B1 LOA22-B2 LOA22-B3; do
  for R in 1 2 3; do
    D=~/Research/05_MD_v2/$L/rep$R/mmpbsa
    if [ -f "$D/RESULTS_v2.dat" ] || [ -f "$D/FINAL_RESULTS.dat" ]; then
      F=$([ -f "$D/RESULTS_v2.dat" ] && echo "$D/RESULTS_v2.dat" || echo "$D/FINAL_RESULTS.dat")
      PB=$(grep -A1 "POISSON BOLTZMANN" "$F" | grep "ΔTOTAL" | awk '{print $2}')
      GB=$(grep -A1 "GENERALIZED BORN" "$F" | grep "ΔTOTAL" | awk '{print $2}')
      printf "  %-12s rep%s  DONE   PB %8s   GB %8s kcal/mol\n" "$L" "$R" "${PB:-?}" "${GB:-?}"
    elif [ -d "$D" ]; then
      LOG=$(ls -t ~/Research/mmpbsa_${L#LOA22-}_rep$R.log ~/Research/mmpbsa_*rep$R.log 2>/dev/null | head -1)
      S=$(grep -E "Beginning|calculating" "$LOG" 2>/dev/null | tail -1 | sed 's/\[INFO *\] *//')
      printf "  %-12s rep%s  running  %s\n" "$L" "$R" "${S:-starting}"
    else
      printf "  %-12s rep%s  not started\n" "$L" "$R"
    fi
  done
done
echo "  sander processes: $(pgrep -c sander)"
