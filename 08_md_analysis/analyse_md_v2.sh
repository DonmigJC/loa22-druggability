#!/bin/bash
# =====================================================================
# MD STAGE 1 ANALYSIS
#
# The question: do the docked bridging poses survive 100 ns?
#
# Method notes - every one of these corrects a defect in the original
# analysis:
#
#   PBC handling.  Three passes: -pbc whole reassembles molecules split
#   across boundaries, -pbc nojump enforces frame-to-frame continuity,
#   then -center -ur compact repositions the protein.  The original
#   analysis applied none of these, which made ligands that crossed a
#   box face appear to fly 11.8 nm - the box width.
#
#   Ligand RMSD.  Fit on protein BACKBONE, measure the LIGAND.  The
#   original fitted the ligand onto itself, which mathematically removes
#   the translation and rotation being measured.
#
#   Distances.  gmx mindist reports distance to the NEAREST PERIODIC
#   COPY, which is ambiguous.  The definitive numbers here come from
#   plain geometry on PBC-corrected static frames (section 5).
#
# Run:  bash ~/Research/analyse_md_v2.sh
# =====================================================================
set -u
R=~/Research
MD=$R/05_MD_v2
LEADS=("LOA22-B1" "LOA22-B2" "LOA22-B3")
REP=${1:-1}

echo "########## MD STAGE $REP ANALYSIS  $(date) ##########"

gi() { awk -v g="$2" 'BEGIN{c=0} /^\[/{n=$0; gsub(/[][]/,"",n);
       gsub(/^ +| +$/,"",n); if(n==g){print c; exit} c++}' "$1"; }

for LEAD in "${LEADS[@]}"; do
  RUN=$MD/$LEAD/rep$REP
  OUT=$RUN/analysis
  echo
  echo "################ $LEAD  replicate $REP ################"
  if [ ! -f "$RUN/md_prod.xtc" ]; then
    echo "  trajectory not found - skipping"; continue
  fi
  mkdir -p "$OUT"; cd "$RUN" || continue

  # ---- 1. index groups -------------------------------------------
  echo "  [1/6] index groups"
  printf 'q\n' | gmx make_ndx -f md_prod.tpr -o "$OUT/d.ndx" >/dev/null 2>&1
  NG=$(grep -c '^\[' "$OUT/d.ndx")
  printf 'r 115-150 & a N CA C O\nname %d BindingSite\nr 121 142\nname %d KeyRes\nq\n' \
      "$NG" "$((NG+1))" \
    | gmx make_ndx -f md_prod.tpr -n "$OUT/d.ndx" -o "$OUT/an.ndx" \
      > "$OUT/mkndx.log" 2>&1

  P=$(gi "$OUT/an.ndx" Protein)
  B=$(gi "$OUT/an.ndx" Backbone)
  S=$(gi "$OUT/an.ndx" BindingSite)
  K=$(gi "$OUT/an.ndx" KeyRes)
  Y=$(gi "$OUT/an.ndx" System)
  L=$(gi "$OUT/an.ndx" "$LEAD")
  [ -z "$L" ] && L=$(gi "$OUT/an.ndx" Other)
  echo "     Protein=$P Backbone=$B BindingSite=$S KeyRes=$K Ligand=$L"
  if [ -z "$L" ] || [ -z "$K" ]; then
    echo "     !! group lookup failed - see $OUT/mkndx.log"; continue
  fi

  # ---- 2. PBC correction, three passes ----------------------------
  echo "  [2/6] PBC correction"
  if [ ! -f "$OUT/center.xtc" ]; then
    printf '%s\n' "$Y" | gmx trjconv -s md_prod.tpr -f md_prod.xtc \
        -n "$OUT/an.ndx" -o "$OUT/whole.xtc" -pbc whole >/dev/null 2>&1
    printf '%s\n' "$Y" | gmx trjconv -s md_prod.tpr -f "$OUT/whole.xtc" \
        -n "$OUT/an.ndx" -o "$OUT/nojump.xtc" -pbc nojump >/dev/null 2>&1
    rm -f "$OUT/whole.xtc"
    printf '%s\n%s\n' "$P" "$Y" | gmx trjconv -s md_prod.tpr \
        -f "$OUT/nojump.xtc" -n "$OUT/an.ndx" -o "$OUT/center.xtc" \
        -center -ur compact >/dev/null 2>&1
    rm -f "$OUT/nojump.xtc"
  fi
  echo "     done"

  # ---- 3. distances (PBC-safe, on the raw trajectory) -------------
  echo "  [3/6] contact distances"
  printf '%s\n%s\n' "$K" "$L" | gmx mindist -s md_prod.tpr -f md_prod.xtc \
      -n "$OUT/an.ndx" -od "$OUT/mindist_key.xvg" -tu ns >/dev/null 2>&1
  printf '%s\n%s\n' "$S" "$L" | gmx mindist -s md_prod.tpr -f md_prod.xtc \
      -n "$OUT/an.ndx" -od "$OUT/mindist_site.xvg" -tu ns >/dev/null 2>&1
  printf '%s\n%s\n' "$P" "$L" | gmx mindist -s md_prod.tpr -f md_prod.xtc \
      -n "$OUT/an.ndx" -od "$OUT/mindist_prot.xvg" -on "$OUT/numcont.xvg" \
      -d 0.6 -tu ns >/dev/null 2>&1

  # ---- 4. RMSD, RMSF, hydrogen bonds ------------------------------
  echo "  [4/6] RMSD, RMSF, hydrogen bonds"
  printf '%s\n%s\n' "$B" "$L" | gmx rms -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/rmsd_ligand.xvg" -tu ns >/dev/null 2>&1
  printf '%s\n%s\n' "$B" "$B" | gmx rms -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/rmsd_backbone.xvg" -tu ns >/dev/null 2>&1
  printf '%s\n%s\n' "$S" "$S" | gmx rms -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/rmsd_bsite.xvg" -tu ns >/dev/null 2>&1
  printf '%s\n' "$P" | gmx rmsf -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -o "$OUT/rmsf.xvg" -res >/dev/null 2>&1
  printf '%s\n%s\n' "$P" "$L" | gmx hbond -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -num "$OUT/hbnum.xvg" -tu ns >/dev/null 2>&1
  printf '%s\n%s\n' "$K" "$L" | gmx hbond -s md_prod.tpr -f "$OUT/center.xtc" \
      -n "$OUT/an.ndx" -num "$OUT/hbnum_key.xvg" -tu ns >/dev/null 2>&1

  # ---- 5. static frames, the definitive measurement ---------------
  echo "  [5/6] static frames"
  mkdir -p "$OUT/frames"
  for T in 0 25000 50000 75000 100000; do
    printf '%s\n' "$Y" | gmx trjconv -s md_prod.tpr -f "$OUT/center.xtc" \
        -n "$OUT/an.ndx" -o "$OUT/frames/${LEAD}_${T}ps.pdb" \
        -dump $T >/dev/null 2>&1
  done
  echo "     5 snapshots written"

  # ---- 6. per-system summary --------------------------------------
  echo "  [6/6] summary"
  ~/miniconda3/envs/reinvent4/bin/python - "$OUT" "$LEAD" <<'PY'
import sys, os, statistics as st
out, lead = sys.argv[1], sys.argv[2]
def series(f):
    v = []
    for l in open(os.path.join(out, f)):
        if l.startswith(("@", "#")): continue
        p = l.split()
        if len(p) >= 2:
            try: v.append(float(p[1]))
            except ValueError: pass
    return v
def stat(f, label, unit="nm"):
    if not os.path.exists(os.path.join(out, f)):
        print(f"     {label:<28} (missing)"); return
    v = series(f)
    if not v: return
    n = len(v)
    first = st.mean(v[:max(1, n//10)])
    last  = st.mean(v[-max(1, n//10):])
    print(f"     {label:<28} mean {st.mean(v):6.3f}  SD {st.pstdev(v):5.3f}  "
          f"max {max(v):6.3f}  first10% {first:5.3f} -> last10% {last:5.3f} {unit}")
stat("mindist_prot.xvg",  "ligand-protein distance")
stat("mindist_site.xvg",  "ligand-binding site")
stat("mindist_key.xvg",   "ligand-Asp121/Arg142")
stat("rmsd_ligand.xvg",   "ligand RMSD (fit backbone)")
stat("rmsd_bsite.xvg",    "binding-site RMSD")
stat("rmsd_backbone.xvg", "protein backbone RMSD")
stat("hbnum.xvg",         "H-bonds, whole protein", "")
stat("hbnum_key.xvg",     "H-bonds, key residues", "")
PY
done

# =====================================================================
echo
echo "############################################################"
echo "  DEFINITIVE MEASUREMENT — static frames, plain geometry"
echo "############################################################"
~/miniconda3/envs/reinvent4/bin/python - <<'PY'
import os, math, glob
R = os.path.expanduser("~/Research")
MD = os.path.join(R, "05_MD_v2")
KEY = {121:{"CB","CG","OD1","OD2"}, 142:{"CB","CG","CD","NE","CZ","NH1","NH2"}}
P0 = {61,64,65,66,67,68,69,70,71,72,73,74,75,88,89,92,93,96,97,115,116,117,
      119,163,165,169,171,172,183,185,186,187}
DOCKED = {"LOA22-B1": (3.43, 3.04), "LOA22-B2": (3.33, 3.13),
          "LOA22-B3": (3.12, 2.88)}

print(f"\n{'system':<22}{'to Asp121/Arg142':>18}{'to P_0':>10}   verdict")
print("-" * 72)
for lead in ("LOA22-B1", "LOA22-B2", "LOA22-B3"):
    d = os.path.join(MD, lead, "rep1", "analysis", "frames")
    if not os.path.isdir(d): continue
    dk0, dp0 = DOCKED[lead]
    print(f"{lead}  (docked pose)      {dk0:>8.2f} A{dp0:>10.2f} A   reference")
    for t in (0, 25000, 50000, 75000, 100000):
        f = os.path.join(d, f"{lead}_{t}ps.pdb")
        if not os.path.exists(f): continue
        key, p0, lig = [], [], []
        for l in open(f):
            if not l.startswith(("ATOM","HETATM")): continue
            nm = l[12:16].strip()
            if nm.startswith("H") or l[76:78].strip() == "H": continue
            try: rn = int(l[22:26])
            except ValueError: continue
            rs = l[17:20].strip()
            c = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            if rs in ("UNL","LIG") or rs == lead[:3].upper(): lig.append(c); continue
            if rn in KEY and nm in KEY[rn]: key.append(c)
            if rn in P0: p0.append(c)
        if not lig or not key:
            print(f"  {t//1000:>4} ns   (ligand or key residues not found)"); continue
        dk = min(math.dist(a,b) for a in lig for b in key)
        dp = min(math.dist(a,b) for a in lig for b in p0) if p0 else float('nan')
        v = "BRIDGING" if dk < 4.5 and dp < 4.5 else (
            "key only" if dk < 4.5 else ("P_0 only" if dp < 4.5 else "neither"))
        print(f"  {t//1000:>4} ns              {dk:>8.2f} A{dp:>10.2f} A   {v}")
    print()
print("van der Waals contact is 3.5-4.0 A; beyond ~5 A there is no")
print("direct interaction.  'BRIDGING' means both contacts are maintained.")
PY

echo
echo "########## ANALYSIS COMPLETE $(date) ##########"
