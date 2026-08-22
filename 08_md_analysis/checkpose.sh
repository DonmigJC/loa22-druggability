#!/bin/bash
# Ligand-position check using gmx mindist (applies minimum-image convention).
# Usage: bash ~/Research/checkpose.sh LOA22-B1 1 npt.gro
L=${1:-LOA22-B1}; R=${2:-1}; S=${3:-npt.gro}
D=~/Research/05_MD_v2/$L/rep$R
[ -d "$D" ] || { echo "  $L rep$R: directory does not exist yet"; exit 0; }
cd "$D" || exit 1
B="${S%.gro}"
[ -f "$S" ] || { echo "  $L rep$R: $S not written yet"; exit 0; }
[ -f "$B.tpr" ] || { echo "  $L rep$R: $B.tpr missing"; exit 0; }
gi() { awk -v g="$2" 'BEGIN{c=0} /^\[/{n=$0; gsub(/[][]/,"",n);
       gsub(/^ +| +$/,"",n); if(n==g){print c; exit} c++}' "$1"; }
NG=$(grep -c '^\[' index.ndx)
printf 'r 121 & a CB CG OD1 OD2\nname %d Asp\nr 142 & a CB CG CD NE CZ NH1 NH2\nname %d Arg\nq\n' \
  "$NG" "$((NG+1))" | gmx make_ndx -f "$S" -n index.ndx -o _chk.ndx >/dev/null 2>&1
A=$(gi _chk.ndx Asp); G=$(gi _chk.ndx Arg); U=$(gi _chk.ndx UNL)
for pair in "$A Asp121" "$G Arg142"; do
  set -- $pair
  printf '%s\n%s\n' "$1" "$U" | gmx mindist -f "$S" -s "$B.tpr" \
      -n _chk.ndx -od /tmp/_d.xvg >/dev/null 2>&1
  awk -v n="$L rep$R" -v r="$2" '!/^[@#]/{printf "  %-14s %s: %.2f A\n", n, r, $2*10; exit}
      END{if(NR==0) printf "  %-14s %s: no output\n", n, r}' /tmp/_d.xvg
done
