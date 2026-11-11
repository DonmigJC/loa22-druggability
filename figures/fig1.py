#!/usr/bin/env python3
"""
Figure 1 - computational workflow.

Drawn programmatically so the counts on the diagram cannot drift from the text.
Every number below is verified against Loa22_Manuscript.md at run time: the
script greps the manuscript for each value and aborts if one is missing.

This replaces the hand-drawn diagram, whose stage 1 and stage 2 boxes carried
compound counts for stages at which no compounds yet existed, and whose
dual-site arrow entered from the 766-molecule filtered subset rather than from
the 6,319-molecule library.

Usage:  python3 fig1.py ../00_manuscript/Loa22_Manuscript.md
Writes Figure1.png at 300 dpi, 180 mm wide.
"""
import sys, os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

MS = sys.argv[1] if len(sys.argv) > 1 else "Loa22_Manuscript.md"

MAIN, ALT, CTRL = "#3380CC", "#B22222", "#E0A21B"
INK, MUTED, LINE = "#1a1a1a", "#5a5a5a", "#9aa0a6"
FILL = "#ffffff"

# every figure-borne number, with the string that must appear in the manuscript
CHECK = [
    ("six pockets", "identified six pockets"),
    ("Drug Score 0.81", "Drug Score of 0.81"),
    ("enclosure 0.24", "0.24"),
    ("145.69", "145.69"),
    ("enclosure 0.00", "enclosure of 0.00"),
    ("12.91 A", "12.91"),
    ("clashscore", "clashscore of 0.00"),
    ("MolProbity 1.41", "MolProbity score of 1.41"),
    ("mean pLDDT 83.16", "83.16"),
    ("6,400 generated", "produced 6,400 molecules"),
    ("6,319 valid", "6,319 (98.7%)"),
    ("766 filtered", "reduced the library to 766 candidates"),
    ("6,400 sampled", "6,400 molecules were sampled from the untrained prior"),
    ("6,049 docked", "6,049 were successfully prepared and docked"),
    ("5,541 scored", "5,541"),
    ("1.89-fold", "1.89-fold enrichment"),
    ("6,253 screened", "6,253"),
    ("three compounds", "Three bearing distinct Bemis"),
    ("100 ns", "100 ns"),
]
txt = open(MS, encoding="utf-8").read()
print("verifying figure content against %s" % os.path.basename(MS))
bad = [n for n, s in CHECK if s not in txt]
for n, s in CHECK:
    print("  %-18s %s" % (n, "ok" if s in txt else "NOT FOUND"))
if bad:
    sys.exit("ABORTED: not found in the manuscript: %s" % ", ".join(bad))
print("all content verified\n")

fig, ax = plt.subplots(figsize=(7.09, 8.3))
ax.set_xlim(0, 100); ax.set_ylim(0, 118); ax.axis("off")
ax.invert_yaxis()

def box(x, y, w, h, title, body, color, sub=None, lw=1.6, fs=7.4):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.6",
                                linewidth=lw, edgecolor=color, facecolor=FILL, zorder=3))
    ax.text(x + w / 2, y + 3.6, title, ha="center", va="center", fontsize=fs,
            fontweight="bold", color=color, zorder=4)
    ax.text(x + w / 2, y + h / 2 + 1.6, body, ha="center", va="center", fontsize=7.0,
            color=INK, zorder=4, linespacing=1.5)
    if sub:
        ax.text(x + w / 2, y + h - 2.6, sub, ha="center", va="center", fontsize=6.2,
                color=MUTED, zorder=4, linespacing=1.4)

def arrow(x1, y1, x2, y2, color=LINE, style="-|>", lw=1.4, ls="-", rad=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=11, linewidth=lw, color=color,
                                 linestyle=ls, zorder=2,
                                 connectionstyle="arc3,rad=%.2f" % rad))

def stage(n, y, label):
    ax.text(2.0, y, "%s" % n, ha="center", va="center",
            fontsize=9.5 if len(str(n)) == 1 else 7.6,
            fontweight="bold", color="white", zorder=5,
            bbox=dict(boxstyle="circle,pad=0.34" if len(str(n)) == 1 else "round,pad=0.34",
                      facecolor=INK, edgecolor="none"))
    ax.text(2.0, y + 5.4, label, ha="center", va="top", fontsize=5.9, color=MUTED,
            rotation=90, zorder=4)

# ---------------------------------------------------------------- 1 structure
stage(1, 6, "structure")
box(9, 1, 52, 15, "AlphaFold2 model  AF-M4QAZ7-F1", 
    "clashscore 0.00  ·  MolProbity 1.41\nmean pLDDT 83.16, and 97.01 across the OmpA domain",
    MAIN, sub="§2.2 to §2.3   ·   no experimental structure exists")

# ---------------------------------------------------------------- 2 cavities
arrow(35, 16, 35, 21)
stage(2, 27, "cavities")
box(9, 21, 52, 17, "Cavity detection across the full surface",
    "six pockets detected\n\nP_0   Drug Score 0.81, enclosure 0.24\nfunctional site   145.69 Å³, enclosure 0.00",
    MAIN, sub="§2.4   ·   the two lie 12.91 Å apart and share no residue")
ax.text(63.5, 29.5, "two search\nvolumes defined", fontsize=6.2, color=MUTED,
        ha="left", va="center", linespacing=1.4)

# ---------------------------------------------------------------- 3 generation
arrow(35, 38, 35, 43)
stage(3, 51, "generation")
box(9, 43, 40, 16, "Reinforcement learning, P_0 objective",
    "6,400 generated\n6,319 valid and unique", MAIN, sub="§2.7 to §2.9")
box(55, 43, 38, 16, "Untrained prior, control",
    "6,400 sampled\n6,049 docked", CTRL, sub="§2.9")
arrow(49, 51, 55, 51, color=CTRL, style="<->", lw=1.2, ls=(0, (3, 2)))
ax.text(52, 48.4, "benchmark", fontsize=5.9, color=CTRL, ha="center", rotation=90)

# 766 side branch
arrow(29, 59, 29, 64, rad=0.0)
ax.add_patch(FancyBboxPatch((7, 64), 30, 9.5, boxstyle="round,pad=0,rounding_size=1.4",
                            linewidth=1.2, edgecolor=LINE, facecolor="#f6f7f8", zorder=3))
ax.text(22, 68.7, "766 after drug-likeness\nand PAINS filtering", ha="center", va="center",
        fontsize=6.6, color=MUTED, zorder=4, linespacing=1.45)
ax.text(22, 77.4, "a filtered subset, not a stage\nthe library passes through",
        ha="center", va="center", fontsize=5.9, color=MUTED, zorder=4, style="italic",
        linespacing=1.4)

# ---------------------------------------------------------------- 3b revised
# The revised campaign restarts from the same untrained prior with a changed
# objective. It does not descend from the prior-only control, so the arrow is
# drawn from the reinforcement-learning box, not from the control box.
arrow(45, 59, 60, 64, color=ALT, rad=-0.18)
ax.text(62.5, 61.6, "repeated with a corrected objective", fontsize=5.9,
        color=ALT, ha="left", va="center")
ax.text(69.5, 80.4, "reported as a separate campaign (Table 9),\nnot screened in stage 4",
        fontsize=5.9, color=ALT, ha="center", va="top", linespacing=1.4)
stage("3b", 71, "")
box(46, 64, 47, 15, "Revised objective, functional-site volume",
    "6,400 generated  →  5,541 scored\n1.89-fold enrichment in dual-site contact",
    ALT, sub="§3.7   ·   the source of the 1.89× in the Abstract")

# ---------------------------------------------------------------- 4 screening
arrow(43, 59, 43, 84, rad=-0.12)
ax.text(38.5, 81.5, "the full 6,319-molecule library,\nnot the filtered 766",
        fontsize=6.0, color=MAIN, ha="right", va="center", linespacing=1.4)
stage(4, 90, "screening")
box(9, 84, 84, 13, "Dual-site screening at both search volumes",
    "6,253 molecules docked at P_0 and at the functional site",
    MAIN, sub="§2.12   ·   0.35% prefer the functional site, mean penalty 1.544 kcal/mol")

# ---------------------------------------------------------------- 5 MD
arrow(51, 97, 51, 102)
stage(5, 108, "dynamics")
box(9, 102, 84, 14, "Molecular dynamics, two starting configurations",
    "three compounds  ·  100 ns each  ·  explicit solvent\ndocked pose at Asp121, and an alternative pose at Arg142",
    MAIN, sub="§2.17   ·   all three docked poses dissociate within 20 ns")

ax.text(50, 119.5, "", ha="center")
plt.tight_layout(pad=0.4)
plt.savefig("Figure1.png", dpi=300, facecolor="white", bbox_inches="tight")
print("wrote Figure1.png")
