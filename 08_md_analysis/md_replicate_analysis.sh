#!/usr/bin/env bash
# =============================================================================
# Ligand-to-pharmacophore distance analysis for the docked-pose replicates,
# reproducing exactly what was done for rep1 and for the Arg142-pose replicates:
# gmx mindist between the ligand (UNL) and the pharmacophoric side chains,
# written in ns on the time axis.
#
# Index groups created here, matching the names already used in the project:
#   KeyRes      Asp121 and Arg142 together  (the pharmacophore)
#   Asp121only  Asp121 alone
#   Arg142only  Arg142 alone
#
# Usage:
#   ./md_replicate_analysis.sh /home/donmigcage/Research/05_MD_v2
#   ./md_replicate_analysis.sh /path/to/05_MD_v2 LOA22-B2 rep2
# =============================================================================
set -euo pipefail

MDROOT="${1:?usage: $0 <path to 05_MD_v2> [compound] [replicate]}"
ONLY_CMPD="${2:-}"
ONLY_REP="${3:-}"
GMX="${GMX:-gmx}"

# Residue numbers of the two pharmacophoric side chains in this model.
ASP=121
ARG=142

for CMPD in LOA22-B1 LOA22-B2 LOA22-B3; do
    [ -n "$ONLY_CMPD" ] && [ "$CMPD" != "$ONLY_CMPD" ] && continue
    for REP in rep1 rep2 rep3; do
        [ -n "$ONLY_REP" ] && [ "$REP" != "$ONLY_REP" ] && continue
        D="$MDROOT/$CMPD/$REP"
        [ -f "$D/md_prod.xtc" ] || { echo "SKIP $CMPD/$REP: no md_prod.xtc" >&2; continue; }

        echo "== $CMPD/$REP"
        mkdir -p "$D/analysis"
        cd "$D"

        # --- build the three residue groups once per replicate ---------------
        if [ ! -f analysis/split.ndx ]; then
            "$GMX" make_ndx -f md_prod.tpr -o analysis/split.ndx >/dev/null 2>&1 <<EOF
r ${ASP} | r ${ARG}
name 24 KeyRes
r ${ASP}
name 25 Asp121only
r ${ARG}
name 26 Arg142only
q
EOF
            echo "   index groups written, verify the names:"
            grep "^\[" analysis/split.ndx | tail -3
        fi

        # --- minimum heavy-atom distance, ligand to each target --------------
        # The two group names are supplied on stdin, reference group first.
        for PAIR in "UNL KeyRes:md_KeyRes" "UNL Asp121only:md_Asp121only" "UNL Arg142only:md_Arg142only"; do
            GROUPS="${PAIR%%:*}"; OUT="${PAIR##*:}"
            echo "$GROUPS" | "$GMX" mindist \
                -s md_prod.tpr -f md_prod.xtc -n analysis/split.ndx \
                -od "analysis/${OUT}.xvg" -tu ns >/dev/null 2>&1
            echo "   wrote analysis/${OUT}.xvg"
        done
    done
done

echo
echo "Done. Next: python3 summarize_docked_replicates.py $MDROOT"
