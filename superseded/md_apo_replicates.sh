#!/usr/bin/env bash
# =============================================================================
# DO NOT RUN THIS AGAINST THE v1 APO BUILD. Superseded 5 November 2026.
# =============================================================================
#
# This script was written to add two replicates to the existing ligand-free run
# in 04_Molecular_Dynamics/Apo/md_run, by reusing that run's build products and
# changing only the velocity seed. Checking the build afterwards showed that
# approach to be wrong, so the script now refuses to run against it.
#
# The v1 apo build does not match the protocol that Materials and Methods 4.6
# states, and the project's own build_md_systems.sh already records every one of
# these as a defect it was written to fix:
#
#   topol.top includes charmm36-jul2022.ff/spce.itp, so the water is SPC/E.
#   CHARMM36 is parameterised against CHARMM-modified TIP3P.      FIX 3
#
#   [ molecules ] carries SOL 38227 and CL 2, that is two chloride ions for
#   charge neutrality and no added salt, against a stated 150 mM.  FIX 4
#   The v2 ligand systems carry NA 109 and CL 112.
#
#   Equilibration was 200 ps NVT and 1 ns NPT at 300 K, against 500 ps NVT
#   and 2 ns NPT at 310 K for every ligand trajectory.             FIX 6
#
# Copying those inputs into rep2 and rep3 would reproduce all three defects in
# two new trajectories and make the apo series internally consistent but still
# unreproducible from the stated methods. Seeding is not the problem here.
#
# The correct route is to build three matched v2 apo replicates with
# build_md_systems.sh, which applies FIX 3, FIX 4, FIX 6 and FIX 8, and to drop
# the v1 run from the paper while leaving it on disk and in the deposit.
#
# Kept in the repository as a record of what was tried and why it was abandoned.
#
# Original usage, retained for reference:
#   ./md_apo_replicates.sh /home/donmigcage/Research/04_Molecular_Dynamics/Apo
# =============================================================================
set -eo pipefail

# --- guard: refuse any apo build whose water model or salt does not match 4.6
_guard_src="${1:-}/md_run"
if [ -f "$_guard_src/topol.top" ]; then
    if grep -q 'spce\.itp' "$_guard_src/topol.top"; then
        echo "REFUSING TO RUN: $_guard_src/topol.top includes spce.itp." >&2
        echo "That is SPC/E water. CHARMM36 requires CHARMM-modified TIP3P," >&2
        echo "and Materials and Methods 4.6 states TIP3P. See the header of" >&2
        echo "this script and build_md_systems.sh FIX 3." >&2
        echo "Build matched v2 apo replicates with build_md_systems.sh instead." >&2
        exit 2
    fi
    if ! grep -qE '^\s*NA\s+[0-9]+' "$_guard_src/topol.top"; then
        echo "REFUSING TO RUN: $_guard_src/topol.top has no sodium." >&2
        echo "Methods 4.6 states 150 mM ionic strength. See FIX 4." >&2
        exit 2
    fi
fi

APOROOT="${1:?usage: $0 <path to the Apo directory>}"
SRC="$APOROOT/md_run"
[ -d "$SRC" ] || { echo "no existing apo run at $SRC" >&2; exit 1; }

declare -A SEEDS=( [rep2]=1337 [rep3]=24601 )

# build products reused unchanged. The apo tree names the restraint file
# loa22_posre.itp rather than posre.itp, and may or may not carry an index.ndx.
INPUTS=(em.gro topol.top loa22_posre.itp posre.itp index.ndx npt.mdp md_prod.mdp)

GMX="${GMX:-gmx}"
GMXRUN="${GMXRUN:--ntmpi 1 -ntomp 6 -nb gpu -pme gpu -bonded gpu -pin on}"

echo "apo source      : $SRC"
echo "gmx             : $("$GMX" --version 2>/dev/null | grep -m1 'GROMACS version')"
echo "mdrun options   : $GMXRUN"
grep -iE "^\s*(nsteps|ref_t)" "$SRC/md_prod.mdp" | sed 's/^/production      : /'

for REP in rep2 rep3; do
    DST="$APOROOT/$REP"
    SEED="${SEEDS[$REP]}"

    if [ -f "$DST/md_prod.gro" ]; then
        echo "== apo/$REP already finished, skipping"
        continue
    fi

    echo "=============================================================="
    echo "== apo $REP   gen_seed = $SEED"
    echo "=============================================================="
    mkdir -p "$DST"
    cd "$DST"

    for f in "${INPUTS[@]}"; do
        if [ -e "$SRC/$f" ]; then cp -n "$SRC/$f" . 2>/dev/null || true; fi
    done
    [ -e charmm36-jul2022.ff ] || cp -r "$SRC/charmm36-jul2022.ff" . 2>/dev/null || true

    # the only thing that differs between replicates
    sed -e "s/^\( *gen_seed *= *\).*/\1${SEED}/" \
        -e "s/^\( *gen_vel *= *\).*/\1yes/" \
        "$SRC/nvt.mdp" > nvt.mdp
    grep -iE "gen_vel|gen_temp|gen_seed|ref_t" nvt.mdp | sed 's/^/   /'

    # restraint file name as referenced by topol.top
    RESTR=em.gro
    NDX=()
    [ -e index.ndx ] && NDX=(-n index.ndx)

    "$GMX" grompp -f nvt.mdp -c em.gro -r "$RESTR" -p topol.top "${NDX[@]}" \
                  -o nvt.tpr -maxwarn 2 2>&1 | tee grompp_nvt.log | tail -3
    "$GMX" mdrun $GMXRUN -deffnm nvt 2>&1 | tee nvt_run.log | tail -3

    "$GMX" grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top "${NDX[@]}" \
                  -o npt.tpr -maxwarn 2 2>&1 | tee grompp_npt.log | tail -3
    "$GMX" mdrun $GMXRUN -deffnm npt 2>&1 | tee npt_run.log | tail -3

    "$GMX" grompp -f md_prod.mdp -c npt.gro -t npt.cpt -p topol.top "${NDX[@]}" \
                  -o md_prod.tpr -maxwarn 2 2>&1 | tee grompp_prod.log | tail -3
    "$GMX" mdrun $GMXRUN -deffnm md_prod 2>&1 | tee md_prod_run.log | tail -3

    echo "== apo $REP done"
done

echo
echo "All apo replicates finished. The apo series is now rep1 (the existing"
echo "md_run directory), rep2 and rep3, seeds 42, 1337 and 24601."
