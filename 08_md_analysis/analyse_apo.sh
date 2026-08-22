#!/bin/bash
# =====================================================================
# APO CONTROL ANALYSIS
#
# Analyses the ligand-free Loa22 simulation and compares its pocket
# stability against the three ligand-bound (holo) runs.
#
# Answers: is the binding pocket stable BECAUSE a ligand is bound, or
# is the AlphaFold model simply rigid on this timescale?
#
# Run:  bash ~/Research/analyse_apo.sh
# =====================================================================
set -u
BASE=~/Research/04_Molecular_Dynamics
APO=$BASE/Apo/md_run
OUT=$APO/analysis

echo "########## APO ANALYSIS $(date) ##########"
cd "$APO" || exit 1
mkdir -p "$OUT"

# helper: resolve a named index group to its number
gi() { awk -v g="$2" 'BEGIN{c=0} /^\[/{n=$0; gsub(/[][]/,"",n);
       gsub(/^ +| +$/,"",n); if(n==g){print c; exit} c++}' "$1"; }

# ---------------------------------------------------------------- 1
echo
echo "=== 1. Building index groups ==="
printf 'q\n' | gmx make_ndx -f md_prod.tpr -o "$OUT/d.ndx" >/dev/null 2>&1
NG=$(grep -c '^\[' "$OUT/d.ndx")
printf 'r 115-150 & a N CA C O\nname %d BindingSite\nq\n' "$NG" \
  | gmx make_ndx -f md_prod.tpr -o "$OUT/an.ndx" > "$OUT/mkndx.log" 2>&1
grep '^\[' "$OUT/an.ndx" | nl -v0 -w3 -s'  ' | tail -6

P=$(gi "$OUT/an.ndx" Protein)
B=$(gi "$OUT/an.ndx" Backbone)
S=$(gi "$OUT/an.ndx" BindingSite)
Y=$(gi "$OUT/an.ndx" System)
echo "  Protein=$P  Backbone=$B  BindingSite=$S  System=$Y"
if [ -z "$S" ]; then echo "!! BindingSite not created - see $OUT/mkndx.log"; exit 1; fi

# ---------------------------------------------------------------- 2
echo
echo "=== 2. PBC correction (three passes) ==="
printf '%s\n' "$Y" | gmx trjconv -s md_prod.tpr -f md_prod.xtc -n "$OUT/an.ndx" \
    -o "$OUT/whole.xtc" -pbc whole 2>&1 | tail -2
printf '%s\n' "$Y" | gmx trjconv -s md_prod.tpr -f "$OUT/whole.xtc" -n "$OUT/an.ndx" \
    -o "$OUT/nojump.xtc" -pbc nojump 2>&1 | tail -2
rm -f "$OUT/whole.xtc"
printf '%s\n%s\n' "$P" "$Y" | gmx trjconv -s md_prod.tpr -f "$OUT/nojump.xtc" \
    -n "$OUT/an.ndx" -o "$OUT/center.xtc" -center -ur compact 2>&1 | tail -2
rm -f "$OUT/nojump.xtc"

# ---------------------------------------------------------------- 3
echo
echo "=== 3. RMSD, RMSF, radius of gyration ==="
printf '%s\n%s\n' "$B" "$B" | gmx rms -s md_prod.tpr -f "$OUT/center.xtc" \
    -n "$OUT/an.ndx" -o "$OUT/rmsd_backbone.xvg" -tu ns 2>&1 | tail -2
printf '%s\n%s\n' "$S" "$S" | gmx rms -s md_prod.tpr -f "$OUT/center.xtc" \
    -n "$OUT/an.ndx" -o "$OUT/rmsd_bsite.xvg" -tu ns 2>&1 | tail -2
printf '%s\n' "$P" | gmx rmsf -s md_prod.tpr -f "$OUT/center.xtc" \
    -n "$OUT/an.ndx" -o "$OUT/rmsf.xvg" -res 2>&1 | tail -2
printf '%s\n' "$P" | gmx gyrate -s md_prod.tpr -f "$OUT/center.xtc" \
    -n "$OUT/an.ndx" -o "$OUT/gyrate.xvg" 2>&1 | tail -2

# ---------------------------------------------------------------- 4
echo
echo "=== 4. Pocket geometry over time (does P_0 stay open?) ==="
printf '%s\n' "$Y" | gmx trjconv -s md_prod.tpr -f "$OUT/center.xtc" \
    -n "$OUT/an.ndx" -o "$OUT/apo_0ns.pdb"  -dump 0     >/dev/null 2>&1
printf '%s\n' "$Y" | gmx trjconv -s md_prod.tpr -f "$OUT/center.xtc" \
    -n "$OUT/an.ndx" -o "$OUT/apo_25ns.pdb" -dump 25000 >/dev/null 2>&1
printf '%s\n' "$Y" | gmx trjconv -s md_prod.tpr -f "$OUT/center.xtc" \
    -n "$OUT/an.ndx" -o "$OUT/apo_50ns.pdb" -dump 50000 >/dev/null 2>&1
echo "  snapshots written to $OUT/"

# ---------------------------------------------------------------- 5
echo
echo "=== 5. APO vs HOLO COMPARISON ==="
~/miniconda3/envs/reinvent4/bin/python - <<'PY'
import os, math, statistics as st

BASE = os.path.expanduser("~/Research/04_Molecular_Dynamics")

def stats(path):
    v = []
    for l in open(path):
        if l.startswith(("@", "#")):
            continue
        p = l.split()
        if len(p) >= 2:
            try:
                v.append(float(p[1]))
            except ValueError:
                pass
    if not v:
        return None
    return st.mean(v), st.pstdev(v), max(v)

rows = [("APO (no ligand)", f"{BASE}/Apo/md_run/analysis"),
        ("LOA22-1 (Lead_1)", f"{BASE}/Lead_1/md_run/analysis_v2"),
        ("LOA22-2 (Lead_4)", f"{BASE}/Lead_4/md_run/analysis_v2"),
        ("LOA22-3 (Lead_5)", f"{BASE}/Lead_5/md_run/analysis_v2")]

print(f"\n{'system':<20}{'backbone RMSD':>22}{'binding-site RMSD':>24}")
print(f"{'':<20}{'mean ± SD (nm)':>22}{'mean ± SD (nm)':>24}")
print("-" * 68)
for name, d in rows:
    bb = stats(os.path.join(d, "rmsd_backbone.xvg")) if \
         os.path.exists(os.path.join(d, "rmsd_backbone.xvg")) else None
    bs = None
    for fn in ("rmsd_bsite.xvg", "rmsd_bindingsite.xvg"):
        p = os.path.join(d, fn)
        if os.path.exists(p):
            bs = stats(p); break
    f1 = f"{bb[0]:.3f} ± {bb[1]:.3f}" if bb else "n/a"
    f2 = f"{bs[0]:.3f} ± {bs[1]:.3f}" if bs else "n/a"
    print(f"{name:<20}{f1:>22}{f2:>24}")

print("\nInterpretation:")
print("  If APO binding-site RMSD is similar to the holo values, the pocket")
print("  is intrinsically stable and ligand binding is not what holds it")
print("  together - the holo numbers cannot be attributed to the ligands.")
print("  If APO is markedly higher, the ligands are stabilising the pocket.")
print("\n  Note: apo ran under the CORRECTED CHARMM36 protocol (position")
print("  restraints on, force-switched vdW, DispCorr off, fixed seed).")
print("  Comparable behaviour indicates the protocol deviations in the")
print("  holo runs did not change the conclusions.")
PY

echo
echo "########## COMPLETE $(date) ##########"
