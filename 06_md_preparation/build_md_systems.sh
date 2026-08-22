#!/bin/bash
# =====================================================================
# MD SYSTEM BUILD AND SIMULATION v2 — CORRECTED PROTOCOL
#
# Every parameter defect from the original runs is fixed here, and each
# fix is annotated with the problem it addresses.
#
#   FIX 1  Position restraints ARE activated during equilibration.
#          Original: `define = -DPOSRES` was absent from nvt.mdp and
#          npt.mdp, leaving the complex exposed to solvent shock from the
#          first femtosecond.  The presence of `refcoord_scaling = com`
#          showed restraints were intended.
#
#   FIX 2  CHARMM36 non-bonded settings are correct.
#          Original: plain cutoff with `DispCorr = EnerPres`.  CHARMM36
#          is parameterised with force-switched Lennard-Jones and NO
#          dispersion correction.  A typo (`dispersion-correction`
#          instead of `DispCorr`) also meant the setting was silently
#          ignored during NVT but active during NPT and production.
#
#   FIX 3  Water model matches the force field.
#          Original: SPC/E.  CHARMM36 is parameterised against
#          CHARMM-modified TIP3P.
#
#   FIX 4  Physiological ionic strength.
#          Original: 2 chloride ions only, i.e. zero added salt, while
#          the manuscript claimed 150 mM NaCl.  Electrostatic screening
#          differs substantially between the two.
#
#   FIX 5  All random seeds are fixed and recorded.
#          Original: `gen_seed = -1`, so velocity generation was not
#          reproducible even in principle.
#
#   FIX 6  Longer equilibration: 500 ps NVT + 2 ns NPT.
#          Original: 100 ps each, short for a ~120,000-atom system.
#
#   FIX 7  No `-maxwarn`, with one narrow, deliberate exception (see
#          FIX 11).  Every OTHER grompp warning must be read, not
#          suppressed.  Original used `-maxwarn 2` throughout, which
#          concealed the FIX 2 typo.
#
#   FIX 8  Three replicates per system with different velocity seeds,
#          addressing the n = 1 limitation.
#
#   FIX 9  Starting coordinates come from the exhaustiveness-32 DOCKED
#          POSE, not from an independently embedded ligand.
#
#   FIX 10 Ligand [ atomtypes ] must be read before ANY [ moleculetype ]
#          opens, including the protein's own - already defined by this
#          point in pdb2gmx's topol.top. The old insertion point (just
#          before [ system ]) put the ACPYPE itp's atomtypes block
#          after that moleculetype, which grompp rejects as "Invalid
#          order for directive atomtypes". The ligand #include now
#          goes right after the forcefield #include instead.
#
#   FIX 11 The ions.mdp grompp call, and ONLY that call, uses
#          `-maxwarn 1` for exactly one expected warning: PME with a
#          net charge. The system is deliberately un-neutralised at
#          this point - genion runs right after, using this same tpr,
#          to place the counter-ions. The warning is guaranteed to
#          appear every time and is not a sign of a defect, unlike the
#          FIX 2 typo that blanket -maxwarn was originally hiding. If
#          grompp reports more than 1 warning here, the build still
#          stops - only the single anticipated one is tolerated. Every
#          other grompp call (em, nvt, npt, production) remains at
#          zero tolerance per FIX 7.
#
# Run:  bash ~/Research/build_md_systems.sh
# =====================================================================
set -u
R=~/Research
MD=$R/05_MD_v2
FF=$R/04_Molecular_Dynamics/charmm36-jul2022.ff
LEADS=("LOA22-B1" "LOA22-B2" "LOA22-B3")

# locate the receptor PDB (it may live on the Windows side)
RECPDB=""
for c in "$R/loa22_alphafold2_raw.pdb" \
         "/mnt/c/Users/donmi/Desktop/Research/loa22_alphafold2_raw.pdb" \
         "$R/04_Molecular_Dynamics/loa22_alphafold2_raw.pdb" \
         "$R/04_Molecular_Dynamics/loa22_clean.pdb"; do
  if [ -f "$c" ]; then RECPDB="$c"; break; fi
done
if [ -z "$RECPDB" ]; then echo "!! receptor PDB not found"; exit 1; fi
echo "receptor: $RECPDB"

# locate the CHARMM36 force field directory
if [ ! -d "$FF" ]; then
  ALT=$(find "$R" -maxdepth 3 -name "charmm36*.ff" -type d 2>/dev/null | head -1)
  if [ -n "$ALT" ]; then FF="$ALT"; else
    echo "!! charmm36 force field directory not found"; exit 1; fi
fi
echo "force field: $FF"
REPLICATES=(1 2 3)
SEEDS=(42 1337 24601)          # FIX 5 / FIX 8
PROD_NS=100

echo "########## MD BUILD v2  $(date) ##########"

# ---------------------------------------------------------------- mdp
write_mdp() {
  local dir=$1 seed=$2
  # --- shared blocks -------------------------------------------------
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

  local TC="tcoupl                   = V-rescale
tc-grps                  = Protein_LIG Water_and_ions
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
define                   = -DPOSRES -DPOSRES_LIG
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
define                   = -DPOSRES -DPOSRES_LIG
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
for LEAD in "${LEADS[@]}"; do
  echo
  echo "################ $LEAD ################"
  WD=$MD/$LEAD
  ACP=$WD/${LEAD}.acpype
  if [ ! -f "$ACP/${LEAD}_GMX.itp" ]; then
    echo "  !! parameters missing - run prep_md_leads.py first"; continue
  fi

  BUILD=$WD/build
  mkdir -p "$BUILD"; cd "$BUILD" || continue
  ln -sfn "$FF" charmm36-jul2022.ff

  # --- protein topology  (FIX 3: -water tip3p) ----------------------
  echo "  [1/7] protein topology (CHARMM36 + TIP3P)"
  printf '1\n1\n' | gmx pdb2gmx -f "$RECPDB" \
      -o protein.gro -p topol.top -i posre.itp \
      -ff charmm36-jul2022 -water tip3p -ignh -ter \
      > pdb2gmx.log 2>&1
  if [ ! -f protein.gro ]; then
    echo "     !! pdb2gmx failed - see $BUILD/pdb2gmx.log"; continue
  fi
  grep -E "Total charge" pdb2gmx.log | tail -1 | sed 's/^/     /'

  # --- merge the docked ligand into the complex  (FIX 9) ------------
  echo "  [2/7] complex from the DOCKED pose"
  python3 - "$BUILD" "$WD/${LEAD}_posed.gro" <<'PY'
import sys
build, liggro = sys.argv[1], sys.argv[2]
def rd(p):
    L = open(p).read().splitlines()
    n = int(L[1]); return L[0], n, L[2:2+n], L[2+n]
_, np_, pa, pbox = rd(f"{build}/protein.gro")
_, nl, la, _     = rd(liggro)
with open(f"{build}/complex.gro","w") as f:
    f.write("Loa22 + ligand (docked pose)\n")
    f.write(f"{np_+nl}\n")
    for l in pa: f.write(l+"\n")
    for l in la: f.write(l+"\n")
    f.write(pbox+"\n")
print(f"     protein {np_} + ligand {nl} = {np_+nl} atoms")
PY

  # --- topology: include ligand + its restraints  (FIX 1) -----------
  echo "  [3/7] topology assembly"
  cp "$ACP/${LEAD}_GMX.itp" .
  cp "$ACP/posre_${LEAD}.itp" . 2>/dev/null || \
     cp "$ACP/posre.itp" "posre_${LEAD}.itp" 2>/dev/null
  python3 - "$BUILD" "$LEAD" <<'PY'
import sys, re
build, lead = sys.argv[1], sys.argv[2]
p = f"{build}/topol.top"
t = open(p).read()
inc = (f'\n; ligand topology\n#include "{lead}_GMX.itp"\n'
       f'#ifdef POSRES_LIG\n#include "posre_{lead}.itp"\n#endif\n\n')
# FIX 10: insert right after the forcefield #include (before the
# protein's own [ moleculetype ], already defined later in this file)
# instead of before [ system ], so the ligand itp's [ atomtypes ]
# block - which ACPYPE puts at its top - is read before ANY
# moleculetype opens. grompp enforces this globally, not per-molecule.
m = re.search(r'#include\s+"[^"]+forcefield\.itp"\s*\n', t)
if not m:
    sys.exit(f"FIX 10: forcefield #include not found in {p}")
t = t[:m.end()] + inc + t[m.end():]
t = t.rstrip() + f"\n{lead}                 1\n"
open(p, "w").write(t)
print("     ligand topology (atomtypes before protein moleculetype) and POSRES_LIG block inserted")
PY
  if ! grep -q "${LEAD}_GMX.itp" topol.top; then
    echo "     !! topology insertion failed (FIX 10) - see traceback above"; continue
  fi

  # --- box, solvate, ionise  (FIX 3, FIX 4) -------------------------
  echo "  [4/7] box and solvation"
  gmx editconf -f complex.gro -o box.gro -c -d 1.2 -bt dodecahedron \
      > editconf.log 2>&1
  gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top \
      > solvate.log 2>&1
  grep -E "^Number of solvent" solvate.log | sed 's/^/     /'

  write_mdp "$BUILD" 42
  echo "  [5/7] adding 150 mM NaCl and neutralising"
  # FIX 11: -maxwarn 1 ONLY here, for the single expected pre-genion
  # net-charge warning. Every other grompp call in this script stays
  # at zero tolerance (FIX 7).
  gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1 \
      > grompp_ions.log 2>&1
  if [ ! -f ions.tpr ]; then
    echo "     !! grompp failed - exceeded the single tolerated warning (FIX 11)"
    grep -A5 -iE "error|warning" grompp_ions.log | head -30 | sed 's/^/     /'
    continue
  fi
  # confirm the one tolerated warning is actually the expected one,
  # not something unrelated slipping through under the same allowance
  if grep -qi "^WARNING" grompp_ions.log && \
     ! grep -qi "net charge" grompp_ions.log; then
    echo "     !! unexpected warning under -maxwarn 1 (FIX 11) - check grompp_ions.log"
    continue
  fi
  printf 'SOL\n' | gmx genion -s ions.tpr -o ions.gro -p topol.top \
      -pname NA -nname CL -neutral -conc 0.15 > genion.log 2>&1
  grep -E "Replacing|Number of" genion.log | tail -3 | sed 's/^/     /'

  # --- minimisation -------------------------------------------------
  echo "  [6/7] energy minimisation"
  gmx grompp -f em.mdp -c ions.gro -p topol.top -o em.tpr \
      > grompp_em.log 2>&1
  if [ ! -f em.tpr ]; then
    echo "     !! grompp failed"; grep -iE "error" grompp_em.log | head -10 \
      | sed 's/^/     /'; continue
  fi
  gmx mdrun -deffnm em -ntmpi 1 -ntomp 12 -pin on > em_run.log 2>&1
  tail -25 em.log | grep -E "Steepest|Potential|Maximum force" | sed 's/^/     /'

  # --- index group for temperature coupling -------------------------
  echo "  [7/7] index groups"
  NG=$(printf 'q\n' | gmx make_ndx -f em.gro -o tmp.ndx 2>/dev/null \
       | grep -cE "^ *[0-9]+ ")
  printf '1 | 13\nname %s Protein_LIG\n"Water" | "Ion"\nname %s Water_and_ions\nq\n' \
      "$NG" "$((NG+1))" | gmx make_ndx -f em.gro -o index.ndx \
      > make_ndx.log 2>&1
  grep -E "Protein_LIG|Water_and_ions" make_ndx.log | tail -2 | sed 's/^/     /'
  echo "     (verify these groups before equilibration)"

  echo "  BUILD COMPLETE for $LEAD -> $BUILD"
done

echo
echo "########## BUILD PHASE DONE $(date) ##########"
echo
echo "Check each build/em.log for convergence, and index.ndx for the"
echo "Protein_LIG and Water_and_ions groups, THEN run:"
echo "    bash ~/Research/run_md_production.sh"
