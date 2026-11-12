#!/usr/bin/env python3
"""
Figure 5 — validation of the docking protocol.

Panels:
  A  Closest approach to the pharmacophore for peptidoglycan fragments
  B  Predicted affinity of each fragment at both sites
  C  Leads vs property-matched decoys vs random library sample
  D  Correlation of bridging score with molecular descriptors

Fix relative to the first version: matplotlib renamed the boxplot keyword
`labels` to `tick_labels` in 3.9. This version sets tick labels explicitly,
which works across versions.

Run from ~/Research:
    ~/miniconda3/envs/reinvent4/bin/python fig5.py
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GREY, RED, BLUE = "#8C8C8C", "#C0392B", "#2E86C1"
fig, ax = plt.subplots(2, 2, figsize=(11, 8))

# ---------------------------------------------------------------- A
frag = ["GlcNAc\n221 Da", "MurNAc\n293 Da", "MDP\n492 Da"]
dk = [13.25, 7.64, 4.94]
ax[0, 0].plot(range(3), dk, "o-", color=RED, ms=10, lw=2)
ax[0, 0].set_xticks(range(3))
ax[0, 0].set_xticklabels(frag, fontsize=9)
ax[0, 0].axhline(4.5, ls="--", color=GREY, lw=1)
ax[0, 0].text(2.05, 4.9, "contact threshold", fontsize=7, color=GREY)
ax[0, 0].set_ylabel("Closest approach to pharmacophore (Å)")
ax[0, 0].set_title("A", loc="left", fontweight="bold")

# ---------------------------------------------------------------- B
site = [-4.878, -5.103, -5.311]
p0 = [-5.244, -6.310, -6.241]
x = np.arange(3)
w = 0.38
ax[0, 1].bar(x - w/2, site, w, color=RED, label="Functional site")
ax[0, 1].bar(x + w/2, p0, w, color=BLUE, label="P_0")
ax[0, 1].set_xticks(x)
ax[0, 1].set_xticklabels(frag, fontsize=9)
ax[0, 1].set_ylabel("Predicted affinity (kcal/mol)")
ax[0, 1].legend(frameon=False, fontsize=8)
ax[0, 1].invert_yaxis()
ax[0, 1].set_title("B", loc="left", fontweight="bold")

# ---------------------------------------------------------------- C
dec = pd.read_csv("decoy_validation/decoy_results.csv")

# The kind column may be named differently; fall back gracefully
kindcol = None
for c in ("kind", "group", "type", "category"):
    if c in dec.columns:
        kindcol = c
        break

distcol = None
for c in ("d_key", "d_pharm", "distance", "min_dist"):
    if c in dec.columns:
        distcol = c
        break

if distcol is None:
    raise SystemExit(f"No distance column found. Columns present: {list(dec.columns)}")

leads = [3.43, 3.33, 3.12]
if kindcol:
    decoys = dec[dec[kindcol].astype(str).str.contains("decoy", case=False,
                 na=False)][distcol].dropna().tolist()
    random = dec[dec[kindcol].astype(str).str.contains("random|baseline",
                 case=False, na=False)][distcol].dropna().tolist()
else:
    # no group column: assume first 45 are decoys, remainder random
    vals = dec[distcol].dropna().tolist()
    decoys, random = vals[:45], vals[45:]

groups = [leads, decoys, random]
labels = [f"Leads\n(n={len(leads)})",
          f"Matched decoys\n(n={len(decoys)})",
          f"Random\n(n={len(random)})"]

bp = ax[1, 0].boxplot(groups, patch_artist=True)
ax[1, 0].set_xticks(range(1, len(labels) + 1))
ax[1, 0].set_xticklabels(labels, fontsize=8)
for patch, colour in zip(bp["boxes"], [RED, GREY, "#D5D8DC"]):
    patch.set_facecolor(colour)
for median in bp["medians"]:
    median.set_color("black")
ax[1, 0].axhline(3.43, ls="--", color=RED, lw=1)
ax[1, 0].text(2.6, 3.6, "weakest lead", fontsize=7, color=RED)
ax[1, 0].set_ylabel("Closest approach to pharmacophore (Å)")
ax[1, 0].set_title("C", loc="left", fontweight="bold")

# ---------------------------------------------------------------- D
props = ["MW", "Heavy\natoms", "logP", "HBD", "Rotatable\nbonds"]
r = [0.417, 0.433, 0.105, 0.103, 0.276]
ax[1, 1].barh(props, r, color=BLUE)
ax[1, 1].axvline(0.5, ls="--", color=RED, lw=1.5)
ax[1, 1].text(0.51, 0.2, "|r| = 0.5", fontsize=7, color=RED, rotation=90)
ax[1, 1].set_xlabel("Pearson r with bridging score")
ax[1, 1].set_xlim(0, 0.6)
ax[1, 1].set_title("D", loc="left", fontweight="bold")

for a in ax.flat:
    a.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plt.savefig("Figure5.png", dpi=300)
print("wrote Figure5.png")
print(f"  decoys plotted: {len(decoys)}, random: {len(random)}")
