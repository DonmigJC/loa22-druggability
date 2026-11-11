#!/usr/bin/env python3
"""
Figure S1 - Needleman-Wunsch global alignment of the modelled Loa22 sequence
(L. interrogans serovar Grippotyphosa, AF-M4QAZ7-F1) against RefSeq
WP_004458759.1 (L. santarosai serovar Shermani).

Reads the EMBOSS Needle output directly, so the figure cannot drift from the
alignment it depicts. Self-verifying: aborts if the residues it is about to
highlight are not the residues the manuscript names.

Writes FigureS1.png at 300 dpi.
"""
import re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

NEEDLE = sys.argv[1] if len(sys.argv) > 1 else "emboss_needle.out"

# residues named in the manuscript, in query (L. interrogans) numbering
CDD = {80: "Y", 81: "K", 120: "T", 121: "D", 124: "G",
       128: "A", 135: "N", 138: "Y", 178: "D", 182: "R"}
PHARM = {119: "H", 121: "D", 142: "R"}

# ---------------------------------------------------------------- parse
top, mid, bot, stats = [], [], [], {}
lines = open(NEEDLE).read().split("\n")
for i, ln in enumerate(lines):
    m = re.match(r"^#\s+(Identity|Similarity|Gaps):\s+(.*)$", ln)
    if m:
        stats[m.group(1)] = m.group(2).strip()
    m = re.match(r"^EMBOSS_001\s+(\d+)\s+(\S+)\s+(\d+)\s*$", ln)
    if m:
        seq = m.group(2)
        if len(top) == len(bot):
            top.append(seq)
            mid.append(lines[i + 1][21:21 + len(seq)])
        else:
            bot.append(seq)

Q = "".join(top)
M = "".join(mid)
S = "".join(bot)
assert len(Q) == len(M) == len(S), "alignment rows differ in length"

# map alignment column -> query residue number
col2q, n = {}, 0
for c, ch in enumerate(Q):
    if ch != "-":
        n += 1
        col2q[c] = n
q2col = {v: k for k, v in col2q.items()}
qseq = Q.replace("-", "")
print("query length %d, alignment length %d" % (len(qseq), len(Q)))

# ---------------------------------------------------------------- verify
for num, aa in {**CDD, **PHARM}.items():
    got = qseq[num - 1]
    if got != aa:
        sys.exit("ABORTED: query residue %d is %s, manuscript says %s" % (num, got, aa))
print("verified: all 12 named residues match the manuscript")

gapcols = [c for c, ch in enumerate(Q) if ch == "-"]
if len(gapcols) != 1:
    sys.exit("ABORTED: expected exactly one gap in the query, found %d" % len(gapcols))
sgap = S[gapcols[0]]
if sgap != "S":
    sys.exit("ABORTED: the inserted residue is %s, manuscript says a serine" % sgap)
print("verified: single gap, inserted residue is Ser")

ident = M.count("|")
print("identities %d/%d" % (ident, len(Q)), stats)

# ---------------------------------------------------------------- draw
PER = 60
rows = [(i, min(i + PER, len(Q))) for i in range(0, len(Q), PER)]
BLOCK_H = 1.0
fig_h = 2.55 + BLOCK_H * len(rows)
fig, ax = plt.subplots(figsize=(9.5, fig_h))
ax.set_xlim(-4.5, PER + 4)
ax.set_ylim(0, fig_h)
ax.axis("off")

MONO = {"family": "DejaVu Sans Mono", "size": 8.2}
INK, MUT = "#1a1a1a", "#666666"
C_CDD, C_PHARM, C_GAP = "#F2A900", "#B22222", "#2E86C1"

ax.text(-0.5, fig_h - 0.35,
        "Figure S1. Global pairwise alignment, Loa22 (L. interrogans serovar "
        "Grippotyphosa) against WP_004458759.1 (L. santarosai serovar Shermani)",
        fontsize=10.5, fontweight="bold", color=INK, va="top")
ax.text(-0.5, fig_h - 0.68,
        "Needleman-Wunsch, EBLOSUM62, gap-open 10.0, gap-extend 0.5.   "
        "Identity %s   Similarity %s   Gaps %s"
        % (stats.get("Identity", "?"), stats.get("Similarity", "?"),
           stats.get("Gaps", "?")),
        fontsize=8.4, color=MUT, va="top")

y = fig_h - 1.15
for a, b in rows:
    qn = [col2q[c] for c in range(a, b) if c in col2q]
    qstart, qend = (qn[0], qn[-1]) if qn else ("", "")
    sstart = len([ch for ch in S[:a] if ch != "-"]) + 1
    send = len([ch for ch in S[:b] if ch != "-"])

    for c in range(a, b):
        x = c - a
        q = col2q.get(c)
        box = None
        if q in PHARM:
            box = C_PHARM
        elif q in CDD:
            box = C_CDD
        elif Q[c] == "-":
            box = C_GAP
        if box:
            ax.add_patch(Rectangle((x - 0.42, y - 0.60), 0.84, 0.90,
                                   facecolor=box, alpha=0.30,
                                   edgecolor=box, lw=0.6, zorder=0))
        ax.text(x, y, Q[c], fontdict=MONO, ha="center", va="center", color=INK, zorder=2)
        ax.text(x, y - 0.24, M[c], fontdict=MONO, ha="center", va="center",
                color=MUT if M[c] != "|" else "#999999", zorder=2)
        ax.text(x, y - 0.48, S[c], fontdict=MONO, ha="center", va="center", color=INK, zorder=2)

    ax.text(PER + 0.8, y, str(qend), fontsize=7.4, color=MUT, va="center")
    ax.text(PER + 0.8, y - 0.48, str(send), fontsize=7.4, color=MUT, va="center")
    ax.text(-1.6, y, str(qstart), fontsize=7.4, color=MUT, va="center", ha="right")
    ax.text(-1.6, y - 0.48, str(sstart), fontsize=7.4, color=MUT, va="center", ha="right")
    y -= BLOCK_H

leg = [(C_PHARM, "His119, Asp121, Arg142 (pharmacophoric and catalytic residues named in the text)"),
       (C_CDD, "the ten CDD-annotated peptidoglycan-binding residues of family cd07185"),
       (C_GAP, "single-residue gap: Ser43 of L. santarosai, the origin of the one-residue offset")]
ly = y + 0.10
for col, lab in leg:
    ax.add_patch(Rectangle((0, ly - 0.09), 1.1, 0.20, facecolor=col, alpha=0.30,
                           edgecolor=col, lw=0.6))
    ax.text(1.7, ly, lab, fontsize=8.0, color=INK, va="center")
    ly -= 0.30

plt.tight_layout()
plt.savefig("FigureS1.png", dpi=300, facecolor="white", bbox_inches="tight")
print("wrote FigureS1.png")
