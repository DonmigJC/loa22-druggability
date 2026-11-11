#!/usr/bin/env bash
# =============================================================================
# build_apo_v2.sh - a ligand-free (apo) system matched to the v2 protocol,
# plus three replicates.
#
# WHY THIS EXISTS
#
# The apo run in 04_Molecular_Dynamics/Apo/md_run is a v1 build and does not
# match what Materials and Methods 4.6 states. Its topol.top includes
# charmm36-jul2022.ff/spce.itp, so the water is SPC/E rather than the
# CHARMM-modified TIP3P that CHARMM36 is parameterised against, and its
# [ molecules ] block carries SOL 38227 with CL 2 and no sodium, so there is no
# added salt against a stated 150 mM. It equilibrated for 200 ps NVT and 1 ns
# NPT at 300 K, against 500 ps and 2 ns at 310 K for every ligand trajectory,
# and it is a single replicate. build_md_systems.sh records all four as FIX 3,
# FIX 4, FIX 6 and FIX 8.
#
# Adding replicates to that build would reproduce every one of those defects,
# so md_apo_replicates.sh now refuses to do it. This script builds the apo
# system again under the same protocol as the ligand systems instead.
#
# WHAT IS THE SAME AS build_md_systems.sh
#   force field and water     CHARMM36-jul2022 with -water tip3p      FIX 3
#   ionic strength            genion -neutral -conc 0.15              FIX 4
#   non-bonded settings       force-switched LJ, no DispCorr          FIX 2
#   equilibration             500 ps NVT and 2 ns NPT, restrained     FIX 1, FIX 6
#   seeds                     42, 1337, 24601, recorded               FIX 5, FIX 8
#   grompp tolerance          zero, except the one pre-genion call    FIX 7, FIX 11
#   box                       dodecahedron, 1.2 nm solute to boundary
#   production                100 ns unrestrained at 310 K and 1 bar
#
# WHAT DIFFERS, AND ONLY BECAUSE THERE IS NO LIGAND
#   no ligand merge, so the solute is protein.gro rather than complex.gro
#   no ligand topology include and no POSRES_LIG block
#   tc-grps is "Protein Non-Protein" rather than "Protein_LIG Water_and_ions",
#     which are default groups, so no custom index.ndx is needed and grompp
#     runs without -n. This is what the v1 apo run also did.
#
# Usage
#   ./build_apo_v2.sh                 build, then run all three replicates
#   ./build_apo_v2.sh --build-only    build and minimise, run nothing
#
# Analysis afterwards is analyse_apo.sh, unchanged. It builds its own index
# from md_prod.tpr and computes the residue 115 to 150 binding-site quantities,
# so it works on these replicates without modification.
# =============================================================================
set -eo pipefail

R=${RESEARCH_ROOT:-$HOME/Research}
# 06_MD_apo_v2 is the path run_queue_resumable.sh already expects: it resumes
# $R/06_MD_apo_v2/rep{1,2,3} from a checkpoint and then calls
# $R/06_MD_apo_v2/apo_replicates.sh. Do not move this without changing that
# wrapper, or an interrupted apo replicate will be restarted from NVT instead
# of resumed.
APO=${APO_ROOT:-$R/06_MD_apo_v2}
FF=$R/04_Molecular_Dynamics/charmm36-jul2022.ff
PROD_NS=${PROD_NS:-100}
REPLICATES=(1 2 3)
SEEDS=(42 1337 24601)
BUILD_ONLY=0
[ "${1:-}" = "--build-only" ] && BUILD_ONLY=1

GMX="${GMX:-gmx}"
GMXRUN="${GMXRUN:--ntmpi 1 -ntomp 6 -nb gpu -pme gpu -bonded gpu -pin on}"
# Minimisation is CPU-only in build_md_systems.sh, which used -ntomp 12 -pin on.
# Kept separate so the build can be tested while another mdrun holds the GPU and
# the pinned cores. Steepest descent converges to the same structure either way.
EM_RUN="${EM_RUN:--ntmpi 1 -ntomp 12 -pin on}"

# the same receptor the ligand systems were built from
RECPDB=""
for c in "$R/loa22_alphafold2_raw.pdb" \
         "/mnt/c/Users/donmi/Desktop/Loa22_Research/04_structure/loa22_alphafold2_raw.pdb"; do
    [ -f "$c" ] && { RECPDB="$c"; break; }
done
[ -n "$RECPDB" ] || { echo "receptor pdb not found" >&2; exit 1; }
[ -d "$FF" ] || { echo "force field not found at $FF" >&2; exit 1; }

echo "receptor : $RECPDB"
echo "gmx      : $("$GMX" --version 2>/dev/null | grep -m1 'GROMACS version')"
echo "output   : $APO"
echo "mdrun    : $GMXRUN"

# ---------------------------------------------------------------- mdp
write_mdp() {
  local dir=$1 seed=$2

  local NB="; FIX 2 - CHARMM36 requires force-switched LJ, no DispCorr
cutoff-scheme            = Verlet
coulombtype              = PME
pme_order                = 4
fourierspacing           = 0.16
rcoulomb                 = 1.2
vdwtype                  = cutoff
vdw-modifier             = force-switch
rvdw-switch              = 1.0
rvdw                     = 1.2
DispCorr                 = no"

  local CONS="constraint_algorithm     = lincs
constraints              = h-bonds
lincs_iter               = 1
lincs_order              = 4"

  # the only substantive difference from the ligand systems: default
  # coupling groups, because there is no ligand to merge into the solute group
  local TC="tcoupl                   = V-rescale
tc-grps                  = Protein Non-Protein
tau_t                    = 0.1     0.1
ref_t                    = 310     310"

  local PC="pcoupl                   = C-rescale
pcoupltype               = isotropic
tau_p                    = 2.0
ref_p                    = 1.0
compressibility          = 4.5e-5"

  cat > "$dir/em.mdp" <<EOF
; Energy minimisation
integrator               = steep
emtol                    = 500.0
emstep                   = 0.01
nsteps                   = 100000
nstlist                  = 10
$NB
pbc                      = xyz
EOF

  cat > "$dir/ions.mdp" <<EOF
integrator               = steep
emtol                    = 1000.0
nsteps                   = 5000
nstlist                  = 10
$NB
pbc                      = xyz
EOF

  cat > "$dir/nvt.mdp" <<EOF
; NVT equilibration - 500 ps, RESTRAINED   (FIX 1, FIX 6)
define                   = -DPOSRES
integrator               = md
nsteps                   = 250000
dt                       = 0.002
nstxout-compressed       = 5000
nstenergy                = 1000
nstlog                   = 1000
continuation             = no
$CONS
nstlist                  = 10
$NB
$TC
pcoupl                   = no
pbc                      = xyz
gen_vel                  = yes
gen_temp                 = 310
gen_seed                 = $seed
EOF

  cat > "$dir/npt.mdp" <<EOF
; NPT equilibration - 2 ns, RESTRAINED   (FIX 1, FIX 6)
define                   = -DPOSRES
integrator               = md
nsteps                   = 1000000
dt                       = 0.002
nstxout-compressed       = 5000
nstenergy                = 1000
nstlog                   = 1000
continuation             = yes
$CONS
nstlist                  = 10
$NB
$TC
$PC
refcoord_scaling         = com
pbc                      = xyz
gen_vel                  = no
EOF

  cat > "$dir/md_prod.mdp" <<EOF
; Production - ${PROD_NS} ns, unrestrained
integrator               = md
nsteps                   = $((PROD_NS * 500000))
dt                       = 0.002
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstxout-compressed       = 5000
compressed-x-grps        = System
nstenergy                = 5000
nstlog                   = 5000
continuation             = yes
$CONS
nstlist                  = 20
$NB
$TC
$PC
pbc                      = xyz
gen_vel                  = no
EOF
}

# ---------------------------------------------------------------- build
BUILD=$APO/build
if [ -f "$BUILD/em.gro" ]; then
    echo
    echo "################ build already present at $BUILD, skipping the build"
    # The mdp files are written by write_mdp during a fresh build, so a build
    # directory produced by any other route may not have them. They are a pure
    # function of the settings above, so regenerate them unconditionally rather
    # than letting the replicate loop copy files that are not there.
    write_mdp "$BUILD" 42
    echo "  mdp files regenerated in $BUILD"
    for f in em.gro topol.top posre.itp nvt.mdp npt.mdp md_prod.mdp; do
        [ -f "$BUILD/$f" ] && echo "     OK   $f" || { echo "     MISS $f"; exit 1; }
    done
    echo "  verifying the existing build against Methods 4.6"
    grep -q 'tip3p\.itp' "$BUILD/topol.top" \
        && echo "     OK   water is TIP3P" || { echo "     FAIL water model"; exit 1; }
    grep -q 'spce\.itp' "$BUILD/topol.top" \
        && { echo "     FAIL spce.itp present"; exit 1; }
    _na=$(awk '/^\[ molecules \]/{f=1} f&&/^[[:space:]]*NA[[:space:]]/{print $2}' "$BUILD/topol.top")
    [ -n "$_na" ] && [ "${_na:-0}" -gt 10 ] \
        && echo "     OK   salt present, NA $_na" || { echo "     FAIL no added sodium"; exit 1; }
else
    echo
    echo "################ building the apo system"
    mkdir -p "$BUILD"; cd "$BUILD"
    ln -sfn "$FF" charmm36-jul2022.ff

    echo "  [1/5] protein topology (CHARMM36 + TIP3P)        FIX 3"
    printf '1\n1\n' | "$GMX" pdb2gmx -f "$RECPDB" \
        -o protein.gro -p topol.top -i posre.itp \
        -ff charmm36-jul2022 -water tip3p -ignh -ter \
        > pdb2gmx.log 2>&1
    [ -f protein.gro ] || { echo "     !! pdb2gmx failed, see $BUILD/pdb2gmx.log"; exit 1; }
    grep -E "Total charge" pdb2gmx.log | tail -1 | sed 's/^/     /'
    # no ligand merge: the solute is the protein alone

    echo "  [2/5] box and solvation"
    "$GMX" editconf -f protein.gro -o box.gro -c -d 1.2 -bt dodecahedron \
        > editconf.log 2>&1
    "$GMX" solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top \
        > solvate.log 2>&1
    grep -E "^Number of solvent" solvate.log | sed 's/^/     /'

    write_mdp "$BUILD" 42

    echo "  [3/5] adding 150 mM NaCl and neutralising        FIX 4"
    # FIX 11: -maxwarn 1 ONLY here, for the single expected net-charge warning
    "$GMX" grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1 \
        > grompp_ions.log 2>&1
    [ -f ions.tpr ] || {
        echo "     !! grompp failed beyond the single tolerated warning (FIX 11)"
        grep -A5 -iE "error|warning" grompp_ions.log | head -30 | sed 's/^/     /'
        exit 1; }
    if grep -qi "^WARNING" grompp_ions.log && ! grep -qi "net charge" grompp_ions.log; then
        echo "     !! unexpected warning under -maxwarn 1 (FIX 11)"; exit 1
    fi
    printf 'SOL\n' | "$GMX" genion -s ions.tpr -o ions.gro -p topol.top \
        -pname NA -nname CL -neutral -conc 0.15 > genion.log 2>&1
    grep -E "Replacing|Number of" genion.log | tail -3 | sed 's/^/     /'

    echo "  [4/5] energy minimisation"
    "$GMX" grompp -f em.mdp -c ions.gro -p topol.top -o em.tpr \
        > grompp_em.log 2>&1
    [ -f em.tpr ] || { echo "     !! grompp failed"; \
        grep -iE "error" grompp_em.log | head -10 | sed 's/^/     /'; exit 1; }
    "$GMX" mdrun -deffnm em $EM_RUN > em_run.log 2>&1
    tail -25 em.log | grep -E "Steepest|Potential|Maximum force" | sed 's/^/     /'

    echo "  [5/5] verifying the build against Methods 4.6"
    ok=1
    grep -q 'tip3p\.itp' topol.top && echo "     OK   water is TIP3P" || { echo "     FAIL water model"; ok=0; }
    grep -q 'spce\.itp'  topol.top && { echo "     FAIL spce.itp present"; ok=0; }
    na=$(awk '/^\[ molecules \]/{f=1} f&&/^[[:space:]]*NA[[:space:]]/{print $2}' topol.top)
    cl=$(awk '/^\[ molecules \]/{f=1} f&&/^[[:space:]]*CL[[:space:]]/{print $2}' topol.top)
    if [ -n "$na" ] && [ "${na:-0}" -gt 10 ]; then
        echo "     OK   salt present, NA $na and CL ${cl:-0}"
    else
        echo "     FAIL no added sodium, NA ${na:-0} CL ${cl:-0}"; ok=0
    fi
    sed -n '/\[ molecules \]/,$p' topol.top | sed 's/^/       /'
    [ "$ok" -eq 1 ] || { echo "     build does not match Methods 4.6, stopping"; exit 1; }
    echo "  BUILD COMPLETE -> $BUILD"
fi

[ "$BUILD_ONLY" -eq 1 ] && { echo; echo "--build-only, stopping here."; exit 0; }

# ---------------------------------------------------------------- replicates
for i in "${!REPLICATES[@]}"; do
    REP="rep${REPLICATES[$i]}"
    SEED="${SEEDS[$i]}"
    DST="$APO/$REP"

    if [ -f "$DST/md_prod.gro" ]; then
        echo "== apo $REP already finished, skipping"; continue
    fi

    echo
    echo "=============================================================="
    echo "== apo $REP   gen_seed = $SEED   ${PROD_NS} ns at 310 K"
    echo "=============================================================="
    mkdir -p "$DST"; cd "$DST"
    for f in em.gro topol.top posre.itp npt.mdp md_prod.mdp; do
        cp -n "$BUILD/$f" . 2>/dev/null || true
    done
    [ -e charmm36-jul2022.ff ] || ln -sfn "$FF" charmm36-jul2022.ff
    sed -e "s/^\( *gen_seed *= *\).*/\1${SEED}/" "$BUILD/nvt.mdp" > nvt.mdp
    grep -iE "gen_vel|gen_temp|gen_seed|ref_t|tc-grps" nvt.mdp | sed 's/^/   /'

    "$GMX" grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr \
        2>&1 | tee grompp_nvt.log | tail -2
    "$GMX" mdrun $GMXRUN -deffnm nvt 2>&1 | tee nvt_run.log | tail -2

    "$GMX" grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr \
        2>&1 | tee grompp_npt.log | tail -2
    "$GMX" mdrun $GMXRUN -deffnm npt 2>&1 | tee npt_run.log | tail -2

    "$GMX" grompp -f md_prod.mdp -c npt.gro -t npt.cpt -p topol.top -o md_prod.tpr \
        2>&1 | tee grompp_prod.log | tail -2
    "$GMX" mdrun $GMXRUN -deffnm md_prod 2>&1 | tee md_prod_run.log | tail -2

    echo "== apo $REP done"
done

echo
echo "All apo replicates finished: rep1, rep2, rep3 at seeds 42, 1337, 24601."
echo "Next: analyse_apo.sh in each replicate directory, then recompute the apo"
echo "row of Supplementary Table S12 and the section 2.6 sentence."
