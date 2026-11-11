#!/bin/bash
# =====================================================================
# Resumable wrapper around md_docked_replicates.sh and apo_replicates.sh.
# Written 5 Nov 2026 after a Windows power-off killed the queue at 4.37 ns.
#
# Why it exists: both driver scripts skip a replicate only if md_prod.gro
# exists. A replicate interrupted mid-production has none, so re-running a
# driver repeats NVT, NPT and production from zero and abandons the
# checkpoint. This wrapper first resumes any partial production from its
# checkpoint, then calls the unchanged drivers for what remains, in order.
#
# Safe to re-run after any interruption. Refuses to start if a queue or an
# mdrun is already live, so it can never create a second mdrun.
#   setsid nohup /home/donmigcage/Research/05_MD_v2/run_queue_resumable.sh \
#       >> /home/donmigcage/Research/05_MD_v2/queue.log 2>&1 < /dev/null &
# =====================================================================
set -uo pipefail
R=/home/donmigcage/Research
exec 9>"$R/05_MD_v2/.queue.lock"
flock -n 9 || { echo "$(date '+%F %T') queue already running, exiting"; exit 0; }
if pgrep -x gmx >/dev/null; then echo "$(date '+%F %T') a gmx process is live, exiting"; exit 1; fi

set +u; source /usr/local/gromacs/bin/GMXRC >/dev/null 2>&1; set -u
export GMX=/usr/local/gromacs/bin/gmx
export GMXRUN="-ntmpi 1 -ntomp 6 -nb gpu -pme gpu -bonded gpu -pin on"
echo "$(date '+%F %T') queue start, $($GMX --version | grep -m1 'GROMACS version')"

resume_partial () {
  local d=$1
  [ -f "$d/md_prod.tpr" ] && [ -f "$d/md_prod.cpt" ] && [ ! -f "$d/md_prod.gro" ] || return 0
  echo "$(date '+%F %T') resuming $d from checkpoint"
  ( cd "$d" && "$GMX" mdrun $GMXRUN -deffnm md_prod -cpi md_prod.cpt -append >> md_prod_resume.log 2>&1 )
  local rc=$?
  if [ $rc -ne 0 ] || [ ! -f "$d/md_prod.gro" ]; then
    echo "$(date '+%F %T') !! resume failed in $d, rc=$rc. Queue stopped."; exit 1
  fi
  echo "$(date '+%F %T') finished $d"
}

for C in LOA22-B1 LOA22-B2 LOA22-B3; do
  for P in rep2 rep3; do resume_partial "$R/05_MD_v2/$C/$P"; done
done
for P in rep1 rep2 rep3; do resume_partial "$R/06_MD_apo_v2/$P"; done

cd "$R/05_MD_v2"
./md_docked_replicates.sh "$R/05_MD_v2" >> replicates.log 2>&1 \
  || { echo "$(date '+%F %T') !! docked driver failed, see replicates.log. Queue stopped."; exit 1; }
echo "$(date '+%F %T') docked driver finished"

cd "$R/06_MD_apo_v2"
./apo_replicates.sh >> apo_replicates.log 2>&1 \
  || { echo "$(date '+%F %T') !! apo driver failed, see apo_replicates.log. Queue stopped."; exit 1; }
echo "$(date '+%F %T') queue complete"
