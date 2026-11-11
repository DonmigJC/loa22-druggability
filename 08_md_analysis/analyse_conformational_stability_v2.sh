#!/bin/bash
# =====================================================================
# CONFORMATIONAL STABILITY, v2 TRAJECTORIES ONLY
#
# Why this script exists
# ---------------------
# Supplementary Table S12 as submitted reports the v1 simulations for
# every row, including the three labelled LOA22-B1, LOA22-B2 and
# LOA22-B3. Checked 5 November 2026, its four value pairs match these
# files to two decimals:
#
#   Apo        17.62 +/- 4.50 | 1.21 +/- 0.26
#              04_Molecular_Dynamics/Apo/md_run/analysis
#   LOA22-B1   20.67 +/- 6.52 | 1.09 +/- 0.20
#              04_Molecular_Dynamics/Lead_1/md_run/analysis_v2
#   LOA22-B2   14.60 +/- 2.43 | 1.32 +/- 0.27
#              04_Molecular_Dynamics/Lead_4/md_run/analysis_v2
#   LOA22-B3   16.22 +/- 5.16 | 2.43 +/- 0.59
#              04_Molecular_Dynamics/Lead_5/md_run/analysis_v2
#
# All four have 5001 frames, that is 50 ns at 10 ps sampling. The v2
# runs give 10001 frames. The "analysis_v2" directory name refers to a
# re-analysis of v1 trajectories, not to the v2 simulations.
#
# The v1 systems used SPC/E water, two chloride ions and no added salt,
# loa22_clean.pdb, 300 K and 50 ns, and the v1 leads embedded the ligand
# independently rather than using the exhaustiveness-32 docked pose.
# build_md_systems.sh documents each of those as a defect, FIX 3, FIX 4,
# FIX 6 and FIX 9. Methods section 4.6 describes only the v2 protocol,
# so nothing in the stated methods reproduces Table S12 as submitted.
#
# This script recomputes the same two statistics on the v2 trajectories
# alone, three replicates per system, so Table S12 becomes reproducible
# from the Methods as written.
#
# Method is copied verbatim from analyse_apo.sh: BindingSite is the
# backbone atoms N, CA, C and O of residues 115 to 150, the trajectory
# gets the same three PBC passes, and the same gmx rms calls run against
# md_prod.tpr. Index groups are resolved by name, never by position.
#
# The backbone column
# -------------------
# Table S12's backbone column uses the GROMACS default Backbone group, which
# is all 195 residues. That group is dominated by the two disordered termini.
# Residues 1 to 47 have a mean pLDDT of 46.8 and residues 191 to 195 fall to
# 31.8, while residues 78 to 186 are all at or above pLDDT 90 with a minimum
# of 90.69.
#
# Measured on the same trajectories Table S12 reports, the difference is not
# a detail. Whole protein against OmpA domain 78 to 186, in Angstrom:
#
#   Apo        17.62 +/- 4.50   against   1.18 +/- 0.22
#   LOA22-B1   20.67 +/- 6.52   against   1.70 +/- 0.22
#   LOA22-B2   14.60 +/- 2.43   against   2.80 +/- 0.58
#   LOA22-B3   16.22 +/- 5.16   against   2.68 +/- 0.57
#
# The whole-protein column also inverts the ranking. LOA22-B1 has the largest
# whole-protein value of the three and the smallest domain value, so a reader
# would take the most stable system for the least stable one. That is why this
# script computes both and reports the domain-restricted figure first.
#
# Usage:
#   bash analyse_conformational_stability_v2.sh            # all systems
#   bash analyse_conformational_stability_v2.sh apo        # apo only
# Re-running is safe: a replicate that already has rmsd_bsite.xvg is
# skipped, so it can be run again after an interruption.
# =====================================================================
set -u

R=/home/donmigcage/Research
MD=$R/05_MD_v2
APO=$R/06_MD_apo_v2
WHICH="${1:-all}"

# GROMACS 2021.6, the version every production run used.
# GMXRC expands LD_LIBRARY_PATH and MANPATH unguarded, so set -u is off here.
set +u; source /usr/local/gromacs/bin/GMXRC >/dev/null 2>&1; set -u
GMX="${GMX:-/usr/local/gromacs/bin/gmx}"
"$GMX" --version | grep -m1 "GROMACS version"

# resolve a named index group to its number
gi() { awk -v g="$2" 'BEGIN{c=0} /^\[/{n=$0; gsub(/[][]/,"",n);
       gsub(/^ +| +$/,"",n); if(n==g){print c; exit} c++}' "$1"; }

analyse_one () {
  local D=$1 LABEL=$2
  local OUT=$D/analysis

  if [ ! -f "$D/md_prod.tpr" ]; then
    echo "  SKIP $LABEL: no md_prod.tpr"; return
  fi

  # pick the trajectory. A run restarted with -noappend leaves parts, and
  # md_prod_full.xtc is the concatenation made at the time.
  local TRJ=""
  for c in md_prod_full.xtc md_prod.xtc; do
    [ -f "$D/$c" ] && { TRJ=$D/$c; break; }
  done
  if [ -z "$TRJ" ]; then
    echo "  SKIP $LABEL: no production trajectory"; return
  fi

  if [ -f "$OUT/rmsd_bsite.xvg" ]; then
    echo "  done already: $LABEL"; return
  fi

  echo "  == $LABEL"
  echo "     trajectory: $(basename "$TRJ") ($(du -h "$TRJ" | cut -f1))"
  mkdir -p "$OUT"
  cd "$D" || return

  printf 'q\n' | "$GMX" make_ndx -f md_prod.tpr -o "$OUT/d.ndx" >/dev/null 2>&1
  local NG; NG=$(grep -c '^\[' "$OUT/d.ndx")
  # Two custom groups. BindingSite is the window the argument rests on.
  # OmpA_BB is the structured domain, residues 78 to 186, every one of which
  # has pLDDT at or above 90 with a minimum of 90.69. See the note in the
  # header about why the whole-protein Backbone group is not used as the
  # headline figure.
  printf 'r 115-150 & a N CA C O\nname %d BindingSite\nr 78-186 & a N CA C O\nname %d OmpA_BB\nq\n' \
    "$NG" "$((NG + 1))" \
    | "$GMX" make_ndx -f md_prod.tpr -o "$OUT/an.ndx" > "$OUT/mkndx.log" 2>&1

  local P B S O Y
  P=$(gi "$OUT/an.ndx" Protein)
  B=$(gi "$OUT/an.ndx" Backbone)
  S=$(gi "$OUT/an.ndx" BindingSite)
  O=$(gi "$OUT/an.ndx" OmpA_BB)
  Y=$(gi "$OUT/an.ndx" System)
  echo "     groups Protein=$P Backbone=$B BindingSite=$S OmpA_BB=$O System=$Y"
  if [ -z "$S" ] || [ -z "$O" ] || [ -z "$P" ] || [ -z "$B" ] || [ -z "$Y" ]; then
    echo "     !! a required group did not resolve, see $OUT/mkndx.log"; return
  fi

  # three PBC passes, exactly as analyse_apo.sh
  printf '%s\n' "$Y" | "$GMX" trjconv -s md_prod.tpr -f "$TRJ" -n "$OUT/an.ndx" \
      -o "$OUT/whole.xtc" -pbc whole >/dev/null 2>&1
  printf '%s\n' "$Y" | "$GMX" trjconv -s md_prod.tpr -f "$OUT/whole.xtc" -n "$OUT/an.ndx" \
      -o "$OUT/nojump.xtc" -pbc nojump >/dev/null 2>&1
  rm -f "$OUT/whole.xtc"
  printf '%s\n%s\n' "$P" "$Y" | "$GMX" trjconv -s md_prod.tpr -f "$OUT/nojump.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/center.xtc" -center -ur compact >/dev/null 2>&1
  rm -f "$OUT/nojump.xtc"

  if [ ! -f "$OUT/center.xtc" ]; then
    echo "     !! PBC correction produced no center.xtc"; return
  fi

  printf '%s\n%s\n' "$B" "$B" | "$GMX" rms -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/rmsd_backbone.xvg" -tu ns >/dev/null 2>&1
  printf '%s\n%s\n' "$O" "$O" | "$GMX" rms -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/rmsd_ompa_bb.xvg" -tu ns >/dev/null 2>&1
  printf '%s\n%s\n' "$S" "$S" | "$GMX" rms -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/rmsd_bsite.xvg" -tu ns >/dev/null 2>&1
  printf '%s\n' "$P" | "$GMX" rmsf -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/rmsf.xvg" -res >/dev/null 2>&1
  printf '%s\n' "$P" | "$GMX" gyrate -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/gyrate.xvg" >/dev/null 2>&1

  # frame count and the two statistics, in Angstrom
  awk '!/^[@#]/ && NF>=2 { v=$2*10; n++; s+=v; ss+=v*v }
       END { if (n>0) { m=s/n;
         printf "     binding-site RMSD  n=%d  %.2f +/- %.2f A\n", n, m, sqrt(ss/n-m*m) } }' \
      "$OUT/rmsd_bsite.xvg"
  awk '!/^[@#]/ && NF>=2 { v=$2*10; n++; s+=v; ss+=v*v }
       END { if (n>0) { m=s/n;
         printf "     OmpA backbone      n=%d  %.2f +/- %.2f A   (78-186)\n", n, m, sqrt(ss/n-m*m) } }' \
      "$OUT/rmsd_ompa_bb.xvg"
  awk '!/^[@#]/ && NF>=2 { v=$2*10; n++; s+=v; ss+=v*v }
       END { if (n>0) { m=s/n;
         printf "     whole backbone     n=%d  %.2f +/- %.2f A   (1-195, termini dominate)\n", n, m, sqrt(ss/n-m*m) } }' \
      "$OUT/rmsd_backbone.xvg"
}

echo "########## v2 CONFORMATIONAL STABILITY  $(date) ##########"

if [ "$WHICH" = "all" ] || [ "$WHICH" = "apo" ]; then
  echo
  echo "=== apo, ligand-free, v2 ==="
  for REP in rep1 rep2 rep3; do
    analyse_one "$APO/$REP" "apo/$REP"
  done
fi

if [ "$WHICH" = "all" ] || [ "$WHICH" = "docked" ]; then
  echo
  echo "=== docked pose, v2 ==="
  for C in LOA22-B1 LOA22-B2 LOA22-B3; do
    for REP in rep1 rep2 rep3; do
      analyse_one "$MD/$C/$REP" "$C/$REP"
    done
  done
fi

echo
echo "########## COMPLETE $(date) ##########"
echo "Next: python3 summarize_conformational_stability.py"
