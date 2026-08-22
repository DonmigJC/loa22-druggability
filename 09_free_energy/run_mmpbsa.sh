#!/bin/bash
# =====================================================================
# MM-PBSA BINDING FREE ENERGY
#
# What this calculates
# --------------------
# Docking scores are empirical predictions fitted to other proteins.
# MM-PBSA computes binding energy from the simulation itself:
#
#   dG_bind = G(complex) - G(protein) - G(ligand)
#
# averaged over many frames, where each G combines molecular-mechanics
# energy (electrostatics, van der Waals) with a continuum solvation term
# (Poisson-Boltzmann for polar, surface-area for non-polar).
#
# It is still approximate - entropy is usually omitted, and absolute
# values are unreliable - but RELATIVE comparisons between ligands on
# the same target are meaningful, and it comes with an error estimate.
#
# Why per-residue decomposition matters here
# ------------------------------------------
# It reports how much each individual amino acid contributes.  That
# converts the bridging finding from geometric ("the ligand sits 3.5 A
# from Asp121") into energetic ("Asp121 contributes X kcal/mol"), which
# is the difference between showing a molecule is NEAR the pharmacophore
# and showing it is BOUND TO it.
#
# Frame sampling
# --------------
# Each frame needs a Poisson-Boltzmann solve taking seconds.  All 5,001
# frames would take weeks.  This uses 100 evenly spaced frames from the
# last 50 ns (the equilibrated portion), which is standard practice and
# converges well.
#
# Usage:
#   bash ~/Research/run_mmpbsa.sh LOA22-B2 1
#   bash ~/Research/run_mmpbsa.sh LOA22-B2 2
#   bash ~/Research/run_mmpbsa.sh LOA22-B2 3
# =====================================================================
set -u

LEAD=${1:-LOA22-B2}
REP=${2:-1}

R=~/Research
RUN=$R/05_MD_v2/$LEAD/rep$REP
OUT=$RUN/mmpbsa
BUILD=$R/05_MD_v2/$LEAD/build

echo "############################################################"
echo "  MM-PBSA   $LEAD  replicate $REP   $(date)"
echo "############################################################"

for f in "$RUN/md_prod.tpr" "$RUN/md_prod.xtc" "$BUILD/topol.top"; do
  [ -f "$f" ] || { echo "  missing: $f"; exit 1; }
done

mkdir -p "$OUT"; cd "$OUT" || exit 1

# link the topology and its includes
ln -sfn "$BUILD/topol.top" topol.top
for f in "$BUILD"/*.itp; do [ -f "$f" ] && ln -sfn "$f" "$(basename "$f")"; done
for d in "$BUILD"/*.ff; do [ -d "$d" ] && ln -sfn "$d" "$(basename "$d")"; done
ln -sfn "$RUN/md_prod.tpr" md_prod.tpr
ln -sfn "$RUN/md_prod.xtc" md_prod.xtc

# ---------------------------------------------------------------- 1
echo
echo "=== 1. Index groups (protein and ligand separately) ==="
gi() { awk -v g="$2" 'BEGIN{c=0} /^\[/{n=$0; gsub(/[][]/,"",n);
       gsub(/^ +| +$/,"",n); if(n==g){print c; exit} c++}' "$1"; }

printf 'q\n' | gmx make_ndx -f md_prod.tpr -o mmpbsa.ndx > mkndx.log 2>&1
P=$(gi mmpbsa.ndx Protein)
L=$(gi mmpbsa.ndx "$LEAD")
[ -z "$L" ] && L=$(gi mmpbsa.ndx Other)
echo "  Protein group = $P    Ligand group = $L"
if [ -z "$P" ] || [ -z "$L" ]; then
  echo "  !! group lookup failed"; grep '^\[' mmpbsa.ndx | nl -v0 -w3; exit 1
fi

# ---------------------------------------------------------------- 2
echo
echo "=== 2. Preparing a PBC-corrected, protein-centred trajectory ==="
if [ ! -f mmpbsa_traj.xtc ]; then
  printf '0\n' | gmx trjconv -s md_prod.tpr -f md_prod.xtc -n mmpbsa.ndx \
      -o whole.xtc -pbc whole > /dev/null 2>&1
  printf '0\n' | gmx trjconv -s md_prod.tpr -f whole.xtc -n mmpbsa.ndx \
      -o nojump.xtc -pbc nojump > /dev/null 2>&1
  rm -f whole.xtc
  printf '%s\n0\n' "$P" | gmx trjconv -s md_prod.tpr -f nojump.xtc \
      -n mmpbsa.ndx -o mmpbsa_traj.xtc -center -ur compact > /dev/null 2>&1
  rm -f nojump.xtc
fi
echo "  done"

# ---------------------------------------------------------------- 3
echo
echo "=== 3. MM-PBSA input file ==="
cat > mmpbsa.in <<'EOF'
Binding free energy with per-residue decomposition

&general
  sys_name             = "Loa22_dual_site"
  startframe           = 2501
  endframe             = 5001
  interval             = 25
  verbose              = 2
  forcefields          = "leaprc.protein.ff14SB,leaprc.gaff2"
  PBRadii              = 3
/

&gb
  igb                  = 5
  saltcon              = 0.150
/

&pb
  istrng               = 0.150
  inp                  = 1
  radiopt              = 1
/

&decomp
  idecomp              = 2
  dec_verbose          = 0
  print_res            = "within 6"
/
EOF
echo "  frames 2501-5001 step 25 = 101 frames from the last 50 ns"
echo "  GB and PB both computed; per-residue decomposition within 6 A"

# ---------------------------------------------------------------- 4
echo
echo "=== 4. Running MM-PBSA (expect several hours) ==="
echo "    started $(date)"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate mmpbsa

gmx_MMPBSA -O -i mmpbsa.in \
    -cs md_prod.tpr -ci mmpbsa.ndx -cg "$P" "$L" \
    -ct mmpbsa_traj.xtc -cp topol.top \
    -o FINAL_RESULTS.dat -eo FINAL_RESULTS.csv \
    -do FINAL_DECOMP.dat -deo FINAL_DECOMP.csv \
    -nogui 2>&1 | tee gmx_mmpbsa_run.log

conda activate base

# ---------------------------------------------------------------- 5
echo
echo "=== 5. Summary ==="
if [ -f FINAL_RESULTS.dat ]; then
  echo
  echo "  --- binding free energy ---"
  grep -A 6 "DELTA TOTAL" FINAL_RESULTS.dat | head -20 | sed 's/^/  /'
  echo
  echo "  --- contribution of the pharmacophoric residues ---"
  if [ -f FINAL_DECOMP.dat ]; then
    grep -E "^\s*(ASP|ARG|HIS)\s*(121|142|119)" FINAL_DECOMP.dat \
      | head -10 | sed 's/^/  /'
    echo
    echo "  --- five largest contributors ---"
    awk '/Total Energy Decomposition/,/^$/' FINAL_DECOMP.dat \
      | awk 'NF>5 && $1!~/^[A-Z][a-z]/' | sort -k18 -n | head -5 | sed 's/^/  /'
  fi
  echo
  echo "  full results : $OUT/FINAL_RESULTS.dat"
  echo "  decomposition: $OUT/FINAL_DECOMP.dat"
else
  echo "  !! no results produced - check gmx_mmpbsa_run.log"
  tail -25 gmx_mmpbsa_run.log | sed 's/^/  /'
fi

echo
echo "########## COMPLETE $(date) ##########"
