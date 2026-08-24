#!/usr/bin/env python3
"""
Figure 6 panel C - bridging occupancy per replicate, alternative pose set.

Source: main-text Table 17. Each replicate contributes 21 frames sampled at
5 ns intervals over 100 ns. A frame counts as bridging when the ligand is
simultaneously within contact of the pharmacophoric region and of P_0.

Colours match panels A and B of the same figure.

Writes Figure6C.png at 300 dpi.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Table 17, per replicate, out of 21 sampled frames
DATA = {"LOA22-B1": [11, 0, 20],
        "LOA22-B2": [21, 21, 21],
        "LOA22-B3": [18, 17, 20]}
COL = {"LOA22-B1": "#C0392B", "LOA22-B2": "#2E86C1", "LOA22-B3": "#27AE60"}
ORDER = ["LOA22-B1", "LOA22-B2", "LOA22-B3"]
NFRAMES = 21
INK, MUTED = "#1a1a1a", "#5a5a5a"

# totals must agree with Table 17
EXPECTED_TOTAL = {"LOA22-B1": 31, "LOA22-B2": 63, "LOA22-B3": 55}
for k, v in DATA.items():
    assert sum(v) == EXPECTED_TOTAL[k], "%s total %d does not match Table 17" % (k, sum(v))
print("totals verified against Table 17:",
      ", ".join("%s %d/63" % (k, sum(DATA[k])) for k in ORDER))

fig, ax = plt.subplots(figsize=(7.2, 4.4))
x = np.arange(3)
w = 0.26

for i, name in enumerate(ORDER):
    off = (i - 1) * w
    vals = DATA[name]
    ax.bar(x + off, vals, w, color=COL[name], zorder=3,
           label="%s  (%d/63)" % (name, sum(vals)))
    for xi, v in zip(x + off, vals):
        if v == 0:
            # a zero bar is invisible, so mark it explicitly
            ax.text(xi, 0.45, "0", ha="center", va="bottom", fontsize=10,
                    color=COL[name], fontweight="bold", zorder=4)
            ax.plot([xi - w/2.4, xi + w/2.4], [0.12, 0.12], color=COL[name],
                    lw=2.4, solid_capstyle="butt", zorder=4)
        else:
            ax.text(xi, v + 0.4, str(v), ha="center", va="bottom",
                    fontsize=9.5, color=INK, zorder=4)

ax.axhline(NFRAMES, ls="--", color="#888888", lw=1, zorder=2)
ax.set_xlim(-0.55, 2.95)
ax.text(2.92, NFRAMES + 0.3, "all sampled frames", ha="right", va="bottom",
        fontsize=8, color=MUTED)

ax.set_xticks(x)
ax.set_xticklabels(["Replicate 1", "Replicate 2", "Replicate 3"], fontsize=10,
                   color=INK)
ax.set_ylabel("Bridging frames (of %d)" % NFRAMES, fontsize=10, color=INK)
ax.set_ylim(0, 23.5)
ax.set_yticks([0, 5, 10, 15, 20])
ax.tick_params(axis="y", labelsize=9, colors=MUTED, length=3)
ax.tick_params(axis="x", length=0)
ax.spines[["top", "right"]].set_visible(False)
ax.spines["left"].set_color("#cccccc")
ax.spines["bottom"].set_color("#cccccc")
ax.grid(axis="y", color="#ececec", lw=0.7, zorder=0)
ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=9, ncol=3, loc="lower center",
          bbox_to_anchor=(0.5, -0.30), columnspacing=2.4)

plt.tight_layout()
plt.savefig("Figure6C.png", dpi=300, facecolor="white")
print("wrote Figure6C.png")
