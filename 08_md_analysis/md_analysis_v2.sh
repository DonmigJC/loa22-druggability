#!/bin/bash
# =====================================================================
# Loa22 MD analysis v2 — PBC-safe, numeric group indices
# Run:  bash ~/Research/md_analysis_v2.sh
# =====================================================================
BASE=~/Research/04_Molecular_Dynamics

# return the 0-based index of a named group in an .ndx file
gidx() {
  awk -v g="$2" 'BEGIN{c=0}
       /^\[/{n=$0; gsub(/[][]/,"",n); gsub(/^ +| +$/,"",n);
             if(n==g){print c; exit} c++}' "$1"
}

for LEAD in Lead_1 Lead_4 Lead_5; do
  MD="$BASE/$LEAD/md_run"
  OUT="$MD/analysis_v2"
  echo "############### $LEAD  $(date) ###############"
  if [ ! -f "$MD/md_prod.xtc" ]; then echo "SKIP: no md_prod.xtc"; continue; fi
  rm -rf "$OUT"; mkdir -p "$OUT"
  cd "$MD" || continue

  # ---- [1/9] build index with a BindingSite group ------------------
  echo "[1/9] index groups"
  printf 'q\n' | gmx make_ndx -f md_prod.tpr -o "$OUT/d.ndx" >/dev/null 2>&1
  NG=$(grep -c '^\[' "$OUT/d.ndx")
  printf 'r 115-150 & a N CA C O\nname %d BindingSite\nq\n' "$NG" \
    | gmx make_ndx -f md_prod.tpr -o "$OUT/an.ndx" > "$OUT/mkndx.log" 2>&1

  echo "   --- groups available ---"
  grep '^\[' "$OUT/an.ndx" | nl -v0 -w3 -s'  '

  P=$(gidx "$OUT/an.ndx" Protein)
  B=$(gidx "$OUT/an.ndx" Backbone)
  L=$(gidx "$OUT/an.ndx" UNL)
  S=$(gidx "$OUT/an.ndx" BindingSite)
  Y=$(gidx "$OUT/an.ndx" System)
  echo "   Protein=$P  Backbone=$B  UNL=$L  BindingSite=$S  System=$Y"

  if [ -z "$L" ] || [ -z "$S" ]; then
    echo "   !! group lookup failed — see $OUT/mkndx.log"; continue
  fi

  # ---- [2/9] mindist on RAW trajectory (PBC-safe, definitive) ------
  echo "[2/9] mindist protein-ligand  <<<< DEFINITIVE"
  printf '%s\n%s\n' "$P" "$L" | gmx mindist -s md_prod.tpr -f md_prod.xtc \
      -n "$OUT/an.ndx" -od "$OUT/mindist.xvg" -on "$OUT/numcont.xvg" \
      -d 0.6 -tu ns 2>&1 | tail -3

  # ---- [3/9] mindist binding site vs ligand ------------------------
  echo "[3/9] mindist bindingsite-ligand"
  printf '%s\n%s\n' "$S" "$L" | gmx mindist -s md_prod.tpr -f md_prod.xtc \
      -n "$OUT/an.ndx" -od "$OUT/mindist_bsite.xvg" -tu ns 2>&1 | tail -3

  # ---- [4-6/9] three-pass PBC correction ---------------------------
  echo "[4/9] pbc whole"
  printf '%s\n' "$Y" | gmx trjconv -s md_prod.tpr -f md_prod.xtc \
      -n "$OUT/an.ndx" -o "$OUT/whole.xtc" -pbc whole 2>&1 | tail -2

  echo "[5/9] pbc nojump"
  printf '%s\n' "$Y" | gmx trjconv -s md_prod.tpr -f "$OUT/whole.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/nojump.xtc" -pbc nojump 2>&1 | tail -2
  rm -f "$OUT/whole.xtc"

  echo "[6/9] centre on protein"
  printf '%s\n%s\n' "$P" "$Y" | gmx trjconv -s md_prod.tpr -f "$OUT/nojump.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/center.xtc" -center -ur compact 2>&1 | tail -2
  rm -f "$OUT/nojump.xtc"

  # ---- [7/9] ligand RMSD relative to protein frame -----------------
  echo "[7/9] ligand RMSD (fit backbone, measure UNL)"
  printf '%s\n%s\n' "$B" "$L" | gmx rms -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/rmsd_ligand.xvg" -tu ns 2>&1 | tail -2

  # ---- [8/9] protein RMSD + RMSF -----------------------------------
  echo "[8/9] protein RMSD and RMSF"
  printf '%s\n%s\n' "$B" "$B" | gmx rms -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/rmsd_backbone.xvg" -tu ns 2>&1 | tail -2
  printf '%s\n%s\n' "$S" "$S" | gmx rms -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/rmsd_bsite.xvg" -tu ns 2>&1 | tail -2
  printf '%s\n' "$P" | gmx rmsf -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/rmsf.xvg" -res 2>&1 | tail -2

  # ---- [9/9] H-bonds and Rg ----------------------------------------
  echo "[9/9] H-bonds and Rg"
  printf '%s\n%s\n' "$P" "$L" | gmx hbond -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -num "$OUT/hbnum.xvg" -tu ns 2>&1 | tail -2
  printf '%s\n' "$P" | gmx gyrate -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/gyrate.xvg" 2>&1 | tail -2

  echo ">>> $LEAD complete $(date)"
  if [ -f "$OUT/mindist.xvg" ]; then
    echo ">>> mindist t=0 : $(awk '!/^[@#]/{print $2; exit}' "$OUT/mindist.xvg") nm"
    echo ">>> mindist t=end: $(awk '!/^[@#]/{v=$2} END{print v}' "$OUT/mindist.xvg") nm"
  fi
done
echo "########## ALL COMPLETE $(date) ##########"
