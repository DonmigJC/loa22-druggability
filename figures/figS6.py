#!/usr/bin/env python3
"""
Figure S6 - computational workflow.

Replaces the earlier FigureS6.tif, which carried section references from a
draft in which the Materials and Methods was section 2 and the Results was
section 3. In the submitted manuscript the Results is section 2 and the
Materials and Methods is section 4, so every reference on the old figure
pointed at the wrong place. It also cited a "Table 9" that does not exist, said
the 1.89-fold enrichment is the source of a figure in the Abstract when the
Abstract no longer carries it, and described the docked poses as dissociating
within 20 ns when the ligand leaves the 4.5 A contact threshold within 20 ps.

Every number below is checked against the manuscript before anything is drawn,
and the script writes no file if a check fails.

Run:  ~/miniconda3/envs/reinvent4/bin/python figS6.py [outdir]
Writes FigureS6.tif at 600 dpi, 180 mm wide, LZW, and FigureS6.png beside it.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
DPI = int(os.environ.get("EXPORT_DPI", "600"))
WIDTH_MM = 180.0

BLUE = "#1f5fa8"
AMBER = "#b8860b"
RED = "#b22222"
GREY = "#555555"
EDGE_GREY = "#999999"

# ----------------------------------------------------------------- the content
# Each box: (stage label, left margin, width, title, colour, lines)
# Lines are plain strings. Section references are to the submitted manuscript:
# Results is section 2, Materials and Methods is section 4.
STAGES = [
    dict(num="1", band="structure", color=BLUE, x=0.10, w=0.62,
         title="AlphaFold2 model  AF-M4QAZ7-F1",
         lines=["clashscore 0.00  .  MolProbity score 1.41",
                "mean pLDDT 83.16, and 97.01 across the OmpA domain",
                "Methods 4.2  .  no experimental structure exists"]),
    dict(num="2", band="cavities", color=BLUE, x=0.10, w=0.62,
         title="Cavity detection across the full surface",
         lines=["six pockets detected",
                "P₀  Drug Score 0.81, enclosure 0.24",
                "functional site in P₄  145.69 Å³, enclosure 0.00",
                "Methods 4.3  .  the P₀ centroid lies 12.91 Å from the",
                "centroid of Asp121 and Arg142, and they share no residue"]),
    dict(num="3", band="generation", color=BLUE, x=0.10, w=0.44,
         title="Reinforcement learning, P₀ objective",
         lines=["6,400 generated", "6,319 valid and unique", "Methods 4.4"]),
    dict(num="4", band="screening", color=BLUE, x=0.10, w=0.62,
         title="Dual-site screening at both search volumes",
         lines=["6,253 molecules docked at P₀ and at the functional site",
                "Methods 4.5  .  0.35% prefer the functional site,",
                "mean penalty 1.544 kcal/mol"]),
    dict(num="5", band="dynamics", color=BLUE, x=0.10, w=0.62,
         title="Molecular dynamics, three starting conditions",
         lines=["three compounds and the ligand-free control",
                "21 trajectories of 100 ns, 2,100 ns in total",
                "Methods 4.6  .  every docked pose leaves the 4.5 Å contact",
                "threshold within the first 20 ps"]),
]

SIDE_CONTROL = dict(color=AMBER, x=0.58, w=0.34,
                    title="Untrained prior, control",
                    lines=["6,400 sampled", "6,049 docked", "Methods 4.4"])

SIDE_REVISED = dict(color=RED, x=0.42, w=0.50,
                    title="Revised objective, functional-site volume",
                    lines=["6,400 generated, 5,541 scored at both volumes",
                           "1.89-fold enrichment in dual-site contact",
                           "Results 2.4  .  a separate campaign, reported in",
                           "Supplementary Table S8, not screened in stage 4"])

FILTER_NOTE = dict(color=GREY, x=0.06, w=0.32,
                   title="",
                   lines=["766 after drug-likeness", "and PAINS filtering"])

# ------------------------------------------------------------------ the checks
EXPECT = {
    "6,319 valid and unique": "6,319",
    "6,049 docked": "6,049",
    "6,253 molecules": "6,253",
    "5,541 scored": "5,541",
    "1.89-fold": "1.89",
    "0.35%": "0.35",
    "1.544 kcal/mol": "1.544",
    "145.69 Å³": "145.69",
    "12.91 Å": "12.91",
    "2,100 ns": "2,100",
    "20 ps": "20 ps",
}
allbox = " ".join(
    s["title"] + " " + " ".join(s["lines"])
    for s in STAGES + [SIDE_CONTROL, SIDE_REVISED, FILTER_NOTE])
missing = [k for k, v in EXPECT.items() if v not in allbox]
if missing:
    sys.exit("ABORTED: the figure does not carry %s" % missing)

# nothing from the draft numbering may survive
for bad in ("§2.", "§3.", "Table 9", "within 20 ns", "Abstract"):
    if bad in allbox:
        sys.exit("ABORTED: stale reference %r is still on the figure" % bad)
print("  content checks passed")

# ------------------------------------------------------------------- the figure
W_IN = WIDTH_MM / 25.4
H_IN = W_IN * 1.36
fig = plt.figure(figsize=(W_IN, H_IN))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

TITLE_PT = 9.0
BODY_PT = 8.0
SMALL_PT = 7.5


def box(x, w, ytop, d, title_pt=TITLE_PT):
    """Draw one rounded box, return its bottom y."""
    n = len(d["lines"]) + (1 if d["title"] else 0)
    h = 0.018 + n * 0.0205
    ax.add_patch(FancyBboxPatch(
        (x, ytop - h), w, h,
        boxstyle="round,pad=0.006,rounding_size=0.012",
        linewidth=1.1, edgecolor=d["color"], facecolor="white", zorder=2))
    y = ytop - 0.020
    if d["title"]:
        ax.text(x + w / 2, y, d["title"], ha="center", va="center",
                fontsize=title_pt, fontweight="bold", color=d["color"], zorder=3)
        y -= 0.0215
    for ln in d["lines"]:
        ax.text(x + w / 2, y, ln, ha="center", va="center",
                fontsize=BODY_PT, color="#111111", zorder=3)
        y -= 0.0205
    return ytop - h


def arrow(x, y0, y1, color="#777777", style="-|>", ls="-"):
    ax.add_patch(FancyArrowPatch(
        (x, y0), (x, y1), arrowstyle=style, mutation_scale=11,
        linewidth=1.0, color=color, linestyle=ls, zorder=1))


y = 0.975
bottoms = {}
for s in STAGES:
    ax.text(0.035, y - 0.021, s["num"], ha="center", va="center",
            fontsize=8.5, color="white", fontweight="bold", zorder=4,
            bbox=dict(boxstyle="circle,pad=0.32", fc="#222222", ec="none"))
    b = box(s["x"], s["w"], y, s)
    bottoms[s["num"]] = b
    ax.text(0.012, (y + b) / 2, s["band"], ha="center", va="center",
            fontsize=SMALL_PT, color=GREY, rotation=90)

    if s["num"] == "3":
        # the prior-only control sits beside the generation stage
        cb = box(SIDE_CONTROL["x"], SIDE_CONTROL["w"], y, SIDE_CONTROL)
        mid = (y + max(b, cb)) / 2
        ax.add_patch(FancyArrowPatch(
            (s["x"] + s["w"] + 0.004, mid), (SIDE_CONTROL["x"] - 0.004, mid),
            arrowstyle="<|-|>", mutation_scale=9, linewidth=1.0,
            color=AMBER, linestyle=(0, (4, 2)), zorder=1))
        ax.text((s["x"] + s["w"] + SIDE_CONTROL["x"]) / 2, mid - 0.016,
                "identical\npipeline", ha="center", va="top",
                fontsize=SMALL_PT, color=AMBER, linespacing=1.3)
        b = min(b, cb)

        # the filtered subset and the revised campaign
        y2 = b - 0.030
        fb = box(FILTER_NOTE["x"], FILTER_NOTE["w"], y2, FILTER_NOTE)
        ax.text(FILTER_NOTE["x"] + FILTER_NOTE["w"] / 2, fb - 0.012,
                "a filtered subset, not a stage the library passes through",
                ha="center", va="top", fontsize=SMALL_PT, color=GREY,
                style="italic")
        ax.text(FILTER_NOTE["x"] + FILTER_NOTE["w"] / 2, fb - 0.030,
                "stage 4 screens the full 6,319-molecule library",
                ha="center", va="top", fontsize=SMALL_PT, color=BLUE)
        arrow(FILTER_NOTE["x"] + FILTER_NOTE["w"] / 2, b, y2, GREY)

        rb = box(SIDE_REVISED["x"], SIDE_REVISED["w"], y2, SIDE_REVISED)
        ax.add_patch(FancyArrowPatch(
            (s["x"] + s["w"] * 0.72, b), (SIDE_REVISED["x"] + 0.10, y2),
            arrowstyle="-|>", mutation_scale=11, linewidth=1.0, color=RED,
            connectionstyle="arc3,rad=-0.25", zorder=1))
        ax.text(SIDE_REVISED["x"] + SIDE_REVISED["w"] / 2, y2 + 0.011,
                "repeated with a corrected objective", ha="center",
                va="bottom", fontsize=SMALL_PT, color=RED)
        b = min(fb - 0.050, rb)

    y = b - 0.034
    if s is not STAGES[-1]:
        arrow(s["x"] + s["w"] / 2 if s["num"] != "3" else 0.40, b, y + 0.002)

png = os.path.join(OUT, "FigureS6.png")
fig.savefig(png, dpi=DPI, facecolor="white")
plt.close(fig)
print("  wrote %s" % png)

from PIL import Image
Image.MAX_IMAGE_PIXELS = None
tif = os.path.join(OUT, "FigureS6.tif")
Image.open(png).convert("RGB").save(
    tif, format="TIFF", compression="tiff_lzw", dpi=(DPI, DPI))
print("  wrote %s" % tif)
print("  %.0f mm wide at %d dpi" % (WIDTH_MM, DPI))
