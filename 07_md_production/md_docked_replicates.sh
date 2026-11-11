#!/usr/bin/env bash
# =============================================================================
# Two additional replicates of the DOCKED-POSE molecular dynamics, so that each
# compound has three independent trajectories, as Protein Science requires.
#
# Design: identical to how the Arg142-pose replicates were produced. Every input
# is reused unchanged from rep1 and only the velocity-generation seed in nvt.mdp
# differs, so each replicate re-equilibrates independently from the same
# energy-minimized structure. Seeds 1337 and 24601 are the same values used for
# rep2_arg142pose and rep3_arg142pose, which keeps the two pose series matched.
#
# Run this on the Linux machine that holds the MD tree (the one whose paths
# appear in the existing .xvg headers as /home/donmigcage/Research/05_MD_v2),
# not on Windows. GROMACS 2021.6 must be on PATH.
#
# Usage:
#   ./md_docked_replicates.sh /home/donmigcage/Research/05_MD_v2
#   ./md_docked_replicates.sh /path/to/05_MD_v2 LOA22-B2        # one compound
# =============================================================================
set -euo pipefail

MDROOT="${1:?usage: $0 <path to 05_MD_v2> [compound ...]}"
shift || true
COMPOUNDS=("$@")
if [ ${#COMPOUNDS[@]} -eq 0 ]; then
    COMPOUNDS=(LOA22-B1 LOA22-B2 LOA22-B3)
fi

# replicate label -> gen_seed
declare -A SEEDS=( [rep2]=1337 [rep3]=24601 )

# Files copied from rep1 into each new replicate directory. These are the
# build products; nothing here depends on the trajectory.
INPUTS=(em.gro topol.top posre.itp index.ndx npt.mdp md_prod.mdp)

# mdrun options. Edit GMXRUN for your hardware. The existing rep1 run reached
# 120.3 ns/day on 117,185 atoms, so 100 ns takes about 20 hours per replicate.
GMX="${GMX:-gmx}"
GMXRUN="${GMXRUN:--ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu}"

for CMPD in "${COMPOUNDS[@]}"; do
    SRC="$MDROOT/$CMPD/rep1"
    [ -d "$SRC" ] || { echo "SKIP $CMPD: no rep1 at $SRC" >&2; continue; }

    for REP in rep2 rep3; do
        DST="$MDROOT/$CMPD/${REP}"
        SEED="${SEEDS[$REP]}"

        if [ -f "$DST/md_prod.gro" ]; then
            echo "== $CMPD/$REP already finished, skipping"
            continue
        fi

        echo "=============================================================="
        echo "== $CMPD $REP   gen_seed = $SEED"
        echo "=============================================================="
        mkdir -p "$DST"
        cd "$DST"

        # --- inputs, copied verbatim from rep1 -------------------------------
        for f in "${INPUTS[@]}"; do
            [ -f "$SRC/$f" ] && cp -n "$SRC/$f" .
        done
        # ligand topology include and the force field directory
        cp -n "$SRC/${CMPD}_GMX.itp" . 2>/dev/null || true
        cp -n "$SRC/posre_${CMPD}.itp" . 2>/dev/null || true
        [ -d charmm36-jul2022.ff ] || cp -r "$SRC/charmm36-jul2022.ff" .

        # --- the only thing that differs between replicates ------------------
        sed -e "s/^\( *gen_seed *= *\).*/\1${SEED}/" \
            -e "s/^\( *gen_vel *= *\).*/\1yes/" \
            "$SRC/nvt.mdp" > nvt.mdp
        grep -E "gen_vel|gen_temp|gen_seed" nvt.mdp

        # --- NVT, 500 ps with position restraints ----------------------------
        "$GMX" grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top \
                      -n index.ndx -o nvt.tpr -maxwarn 1 2>&1 | tee grompp_nvt.log
        "$GMX" mdrun $GMXRUN -deffnm nvt 2>&1 | tee nvt_run.log

        # --- NPT, 2 ns with position restraints ------------------------------
        "$GMX" grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top \
                      -n index.ndx -o npt.tpr -maxwarn 1 2>&1 | tee grompp_npt.log
        "$GMX" mdrun $GMXRUN -deffnm npt 2>&1 | tee npt_run.log

        # --- production, 100 ns unrestrained ---------------------------------
        "$GMX" grompp -f md_prod.mdp -c npt.gro -t npt.cpt -p topol.top \
                      -n index.ndx -o md_prod.tpr -maxwarn 1 2>&1 | tee grompp_prod.log
        "$GMX" mdrun $GMXRUN -deffnm md_prod 2>&1 | tee md_prod_run.log

        echo "== $CMPD $REP done"
    done
done

echo
echo "All requested replicates finished. Next: ./md_replicate_analysis.sh $MDROOT"
