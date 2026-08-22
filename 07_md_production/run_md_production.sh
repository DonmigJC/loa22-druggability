#!/bin/bash
# =====================================================================
# PRODUCTION MD — STAGED
#
# Stage 1 (default): replicate 1 for all three systems.   ~307 ns, ~4 days
# Stage 2:           replicates 2 and 3.                  ~615 ns, ~8 days
#
# Each replicate is INDEPENDENTLY equilibrated with its own velocity
# seed, so the three runs do not share initial conditions.  Replicates
# branched from a single equilibrated structure would test only
# trajectory divergence; independent equilibration tests the protocol.
#
# Every stage checkpoints, so an interruption costs only the work since
# the last checkpoint.  Re-running this script resumes automatically.
#
# Usage
#   bash ~/Research/run_md_production.sh            # stage 1 (rep 1)
#   bash ~/Research/run_md_production.sh 2 3        # replicates 2 and 3
#   bash ~/Research/run_md_production.sh 1 2 3      # everything
# =====================================================================
set -u
R=~/Research
MD=$R/05_MD_v2
LEADS=("LOA22-B1" "LOA22-B2" "LOA22-B3")

declare -A SEED=( [1]=42 [2]=1337 [3]=24601 )

REPS=("$@")
[ ${#REPS[@]} -eq 0 ] && REPS=(1)

NT=6             # OpenMP threads
GPUFLAGS="-nb gpu -pme gpu -bonded gpu -pin on"

echo "############################################################"
echo "  PRODUCTION MD   $(date)"
echo "  replicates: ${REPS[*]}"
echo "  systems   : ${LEADS[*]}"
echo "############################################################"

# helper: run a stage only if its output is absent; resume from .cpt
run_stage() {
  local dir=$1 name=$2 label=$3
  cd "$dir" || return 1
  if [ -f "${name}.gro" ]; then
    echo "      $label already complete"
    return 0
  fi
  local resume=""
  if [ -f "${name}.cpt" ]; then
    resume="-cpi ${name}.cpt"
    echo "      $label resuming from checkpoint"
  else
    echo "      $label starting"
  fi
  gmx mdrun -deffnm "$name" $resume -ntmpi 1 -ntomp $NT $GPUFLAGS \
      >> "${name}_run.log" 2>&1
  if [ ! -f "${name}.gro" ]; then
    echo "      !! $label did not finish - see ${dir}/${name}_run.log"
    tail -12 "${name}_run.log" | sed 's/^/         /'
    return 1
  fi
  local perf
  perf=$(grep -E "Performance:" "${name}.log" 2>/dev/null | tail -1 | awk '{print $2}')
  echo "      $label complete${perf:+  (${perf} ns/day)}"
  return 0
}

for REP in "${REPS[@]}"; do
  S=${SEED[$REP]}
  echo
  echo "############################################################"
  echo "  REPLICATE $REP   (velocity seed $S)"
  echo "############################################################"

  for LEAD in "${LEADS[@]}"; do
    BUILD=$MD/$LEAD/build
    RUN=$MD/$LEAD/rep$REP
    echo
    echo "  ---- $LEAD  replicate $REP ----"

    if [ ! -f "$BUILD/em.gro" ]; then
      echo "      !! build incomplete for $LEAD - skipping"
      continue
    fi

    mkdir -p "$RUN"; cd "$RUN" || continue

    # link the shared build inputs
    for f in topol.top index.ndx posre.itp "posre_${LEAD}.itp" \
             "${LEAD}_atomtypes.itp" "${LEAD}_moleculetype.itp" \
             "${LEAD}_GMX.itp" em.gro; do
      [ -f "$BUILD/$f" ] && ln -sfn "$BUILD/$f" "$f"
    done
    for d in "$BUILD"/*.ff; do
      [ -d "$d" ] && ln -sfn "$d" "$(basename "$d")"
    done

    # per-replicate mdp files: only gen_seed differs between replicates
    for m in nvt npt md_prod; do
      sed "s/^gen_seed .*/gen_seed                 = $S/" \
          "$BUILD/$m.mdp" > "$m.mdp"
    done

    # ---------------- NVT ------------------------------------------
    if [ ! -f nvt.tpr ]; then
      gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top \
          -n index.ndx -o nvt.tpr > grompp_nvt.log 2>&1
      if [ ! -f nvt.tpr ]; then
        echo "      !! grompp nvt failed"
        grep -B2 -A6 -i "error" grompp_nvt.log | head -20 | sed 's/^/         /'
        continue
      fi
    fi
    run_stage "$RUN" nvt "NVT 500 ps (restrained)" || continue

    # ---------------- NPT ------------------------------------------
    if [ ! -f npt.tpr ]; then
      gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top \
          -n index.ndx -o npt.tpr > grompp_npt.log 2>&1
      if [ ! -f npt.tpr ]; then
        echo "      !! grompp npt failed"
        grep -B2 -A6 -i "error" grompp_npt.log | head -20 | sed 's/^/         /'
        continue
      fi
    fi
    run_stage "$RUN" npt "NPT 2 ns (restrained)" || continue

    # check equilibration before committing to production
    if [ -f npt.edr ]; then
      printf 'Temperature\nPressure\nDensity\n\n' | \
        gmx energy -f npt.edr -o npt_stats.xvg > energy.log 2>&1
      grep -E "^Temperature|^Pressure|^Density" energy.log \
        | awk '{printf "      %-14s %10.2f +/- %8.2f\n", $1, $2, $3}'
    fi

    # ---------------- production -----------------------------------
    if [ ! -f md_prod.tpr ]; then
      gmx grompp -f md_prod.mdp -c npt.gro -t npt.cpt -p topol.top \
          -n index.ndx -o md_prod.tpr > grompp_prod.log 2>&1
      if [ ! -f md_prod.tpr ]; then
        echo "      !! grompp production failed"
        grep -B2 -A6 -i "error" grompp_prod.log | head -20 | sed 's/^/         /'
        continue
      fi
    fi
    echo "      production 100 ns starting $(date '+%H:%M')"
    run_stage "$RUN" md_prod "production 100 ns" || continue

    echo "      >>> $LEAD replicate $REP COMPLETE"
  done
done

echo
echo "############################################################"
echo "  ALL REQUESTED RUNS FINISHED   $(date)"
echo "############################################################"
echo
echo "Completed trajectories:"
for REP in "${REPS[@]}"; do
  for LEAD in "${LEADS[@]}"; do
    f=$MD/$LEAD/rep$REP/md_prod.xtc
    if [ -f "$f" ]; then
      printf "  %-12s rep %s   %s\n" "$LEAD" "$REP" "$(du -h "$f" | cut -f1)"
    else
      printf "  %-12s rep %s   incomplete\n" "$LEAD" "$REP"
    fi
  done
done
echo
echo "Next: bash ~/Research/analyse_md_v2.sh"
