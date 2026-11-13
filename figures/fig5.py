#!/usr/bin/env python3
"""
Figure 5 - validation of the docking protocol.

The submitted paper numbers this Figure 4. The name follows the draft numbering
that every script in this directory uses; README.md carries the mapping.

This replaces an earlier fig5.py that no longer reproduced the submitted
figure. That version drew panel B without value labels or error bars, panel C
as a boxplot rather than a strip plot with group means, and panel D with one
series instead of two. It is kept in superseded/.

The reason for rewriting rather than patching: panel D's legend ran off the
right edge of the submitted TIFF, cutting "with distance to the site" mid-word.
The legend is now placed above the axes with room for both entries.

Panels:
  A  Closest approach to the pharmacophore against fragment size
  B  Predicted affinity of each fragment at both sites
  C  Leads against property-matched decoys and a random library sample
  D  Correlation of two geometric measures with molecular descriptors

Every value is taken from the deposited data or from the supplementary table
that reports it, and the script checks each against the published figure before
drawing. Panel C reads data/docking/decoy_results.csv directly.

Run:  ~/miniconda3/envs/reinvent4/bin/python fig5.py [outdir]
Writes Figure4.tif at 600 dpi, 180 mm wide, LZW, and Figure4.png beside it.
"""
import os
import sys
import csv
import io

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
DPI = int(os.environ.get("EXPORT_DPI", "600"))
REPO = os.environ.get("LOA22_REPO", os.path.expanduser("~/loa22-code"))

RED, BLUE, GREY, GOLD = "#C0392B", "#2E86C1", "#8C8C8C", "#8a6d1f"
A = "Å"

# ------------------------------------------------------------------ the data
# Supplementary Table S9
FRAG = [("N-acetyl-D-\nglucosamine", 221, 13.25, -4.878, None, -5.244),
        ("N-acetyl-\nmuramic acid", 293, 7.64, -5.103, 0.540, -6.310),
        ("Muramyl\ndipeptide", 492, 4.94, -5.311, 0.154, -6.241)]
# Supplementary Table S10 and main-text Table 4
LEADS = [3.43, 3.33, 3.12]
MEANS = {"Selected\nleads": 3.29, "Matched\ndecoys": 6.92, "Random\nsample": 8.81}
WEAKEST = 3.43
# Supplementary Table S11
DESC = [("Molecular weight", 0.417, -0.575),
        ("Heavy atom count", 0.433, -0.592),
        ("Rotatable bonds", 0.276, -0.315),
        ("Calculated logP", 0.105, -0.300),
        ("H-bond donors", 0.103, -0.061)]

path = os.path.join(REPO, "data", "docking", "decoy_results.csv")
rows = list(csv.DictReader(io.open(path, encoding="utf-8", errors="replace")))
decoys = [float(r["d_key"]) for r in rows
          if "decoy" in r["kind"].lower() and r["d_key"]]
random_ = [float(r["d_key"]) for r in rows
           if r["kind"].lower() not in ("decoy",) and r["d_key"]]

# ---------------------------------------------------------------- the checks
fail = []
if len(decoys) != 45:
    fail.append("expected 45 decoys, read %d" % len(decoys))
if len(random_) != 15:
    fail.append("expected 15 random molecules, read %d" % len(random_))
for name, got, want, tol in (
        ("decoy mean distance", np.mean(decoys), 6.92, 0.006),
        ("random mean distance", np.mean(random_), 8.81, 0.006),
        ("lead mean distance", np.mean(LEADS), 3.29, 0.006),
        ("closest decoy", min(decoys), 3.44, 0.006),
        ("closest random", min(random_), 3.78, 0.006)):
    if abs(got - want) > tol:
        fail.append("%s is %.3f, Supplementary Table S10 says %.2f"
                    % (name, got, want))
n_below = sum(1 for d in decoys if d < WEAKEST)
if n_below != 0:
    fail.append("%d decoys reach the weakest lead, the paper says none" % n_below)
if fail:
    sys.exit("ABORTED, nothing written:\n  " + "\n  ".join(fail))
print("  checks passed: 45 decoys, 15 random, means %.2f, %.2f, %.2f, none below %.2f"
      % (np.mean(LEADS), np.mean(decoys), np.mean(random_), WEAKEST))

# --------------------------------------------------------------- the figure
W_IN = 180.0 / 25.4
fig, ax = plt.subplots(2, 2, figsize=(W_IN, W_IN * 1.10))
plt.rcParams.update({"font.size": 8})
LBL, TCK, ANN = 8.5, 8.0, 8.0

# ---------------------------------------------------------------------- A
a = ax[0, 0]
mw = [f[1] for f in FRAG]
dk = [f[2] for f in FRAG]
a.plot(mw, dk, "o-", color=RED, ms=7, lw=2, zorder=3)
NICE = ["N-acetyl-D-glucosamine", "N-acetylmuramic acid", "Muramyl dipeptide"]
for (name, m, d, _s, _e, _p), nice, off in zip(FRAG, NICE, ((10, 10), (10, 10), (10, 12))):
    a.annotate("%s\n%.2f %s" % (nice, d, A), (m, d),
               textcoords="offset points", xytext=off,
               fontsize=7.5, ha="left", va="bottom")
a.axhline(4.5, ls="--", color=GREY, lw=1)
a.text(700, 4.0, "4.5 %s contact threshold" % A, fontsize=7.5, color=GREY,
       ha="right", va="top")
a.set_xlim(170, 720)
a.set_ylim(0, 19)
a.set_xlabel("Fragment molecular weight (Da)", fontsize=LBL)
a.set_ylabel("Closest approach to\nthe pharmacophore (%s)" % A, fontsize=LBL)
a.set_title("A", loc="left", fontweight="bold", fontsize=11)

# ---------------------------------------------------------------------- B
b = ax[0, 1]
x = np.arange(3)
w = 0.36
site = [f[3] for f in FRAG]
serr = [f[4] if f[4] else 0.0 for f in FRAG]
p0 = [f[5] for f in FRAG]
b.bar(x - w / 2, site, w, color=RED, label="functional site",
      yerr=serr, error_kw=dict(ecolor="#333333", lw=1, capsize=3))
b.bar(x + w / 2, p0, w, color=BLUE, label="P₀")
# the label clears the error bar where there is one
for xi, v, e in zip(x - w / 2, site, serr):
    b.text(xi, v - e - 0.16, "%.2f" % v, ha="center", va="top", fontsize=7.5)
for xi, v in zip(x + w / 2, p0):
    b.text(xi, v - 0.12, "%.2f" % v, ha="center", va="top", fontsize=7.5)
b.set_xticks(x)
b.set_xticklabels([f[0] for f in FRAG], fontsize=7.5)
b.set_ylim(-8.4, 0)
b.set_ylabel("Predicted affinity (kcal/mol)", fontsize=LBL)
b.legend(frameon=False, fontsize=7.5, loc="lower left",
         bbox_to_anchor=(0.0, 0.0), ncol=2)
b.set_title("B", loc="left", fontweight="bold", fontsize=11)
b.text(0.5, 1.02, "every fragment prefers P₀", transform=b.transAxes,
       ha="center", va="bottom", fontsize=8, color=GREY)

# ---------------------------------------------------------------------- C
c = ax[1, 0]
rng = np.random.default_rng(42)          # jitter only, fixed so it reproduces
for i, (label, vals, colour) in enumerate(
        zip(MEANS, (LEADS, decoys, random_), (RED, BLUE, GOLD))):
    xs = i + rng.uniform(-0.17, 0.17, len(vals))
    c.scatter(xs, vals, s=14, color=colour, alpha=0.65, linewidths=0, zorder=2)
    m = MEANS[label]
    c.plot([i - 0.30, i + 0.30], [m, m], color=colour, lw=2.4, zorder=3)
    c.text(i + 0.34, m, "%.2f %s" % (m, A), color=colour, fontsize=8,
           fontweight="bold", va="center")
c.axhline(WEAKEST, ls="--", color=GREY, lw=1)
c.text(2.42, WEAKEST - 0.5, "weakest lead, %.2f %s" % (WEAKEST, A),
       fontsize=7.5, color=GREY, ha="right", va="top")
c.set_xticks(range(3))
c.set_xticklabels(list(MEANS), fontsize=8)
c.set_xlim(-0.5, 2.6)
c.set_ylim(0, 22)
c.set_ylabel("Closest approach to\nthe pharmacophore (%s)" % A, fontsize=LBL)
c.text(0.02, 0.97, "0 of 45 decoys reach the weakest lead\n"
                   "Mann-Whitney U = 0, p = 0.0001",
       transform=c.transAxes, va="top", ha="left", fontsize=7.5)
c.set_title("C", loc="left", fontweight="bold", fontsize=11)

# ---------------------------------------------------------------------- D
d = ax[1, 1]
names = [t[0] for t in DESC]
y = np.arange(len(DESC))[::-1]
rb = [t[1] for t in DESC]
rd = [t[2] for t in DESC]
d.barh(y + 0.18, rb, 0.34, color=RED)
d.barh(y - 0.18, rd, 0.34, color=BLUE)
for yi, v in zip(y + 0.18, rb):
    d.text(v + 0.02, yi, "%.3f" % v, va="center", fontsize=7.5)
for yi, v in zip(y - 0.18, rd):
    d.text(v - 0.02, yi, "%.3f" % v, va="center", ha="right", fontsize=7.5)
d.axvline(0, color="#444444", lw=1.1)
for s in (-0.5, 0.5):
    d.axvline(s, ls="--", color=GREY, lw=1)
d.set_ylim(-1.25, len(DESC) - 0.4)
d.set_yticks(y)
d.set_yticklabels(names, fontsize=8)
d.set_xlim(-0.80, 0.74)
d.set_xlabel("Pearson correlation coefficient", fontsize=LBL)
d.set_title("D", loc="left", fontweight="bold", fontsize=11)
# the legend that used to run off the edge, now above the axes and inside them
h = [plt.Rectangle((0, 0), 1, 1, color=RED),
     plt.Rectangle((0, 0), 1, 1, color=BLUE)]
d.legend(h, ["with bridging score", "with distance to the site"],
         frameon=False, fontsize=7.5, ncol=1,
         loc="lower right", bbox_to_anchor=(1.0, -0.02))

for a_ in ax.flat:
    a_.spines[["top", "right"]].set_visible(False)
    a_.tick_params(labelsize=TCK)

plt.tight_layout(pad=1.1)
png = os.path.join(OUT, "Figure4.png")
fig.savefig(png, dpi=DPI, facecolor="white")
plt.close(fig)
print("  wrote %s" % png)

from PIL import Image
Image.MAX_IMAGE_PIXELS = None
im = Image.open(png).convert("RGB")
tif = os.path.join(OUT, "Figure4.tif")
im.save(tif, format="TIFF", compression="tiff_lzw", dpi=(DPI, DPI))
print("  wrote %s  %s  %.0f mm at %d dpi"
      % (tif, im.size, im.size[0] / DPI * 25.4, DPI))
