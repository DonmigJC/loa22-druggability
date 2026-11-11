#!/usr/bin/env python3
"""
Figure 5 panel C (source Figure6C) - bridging occupancy per replicate,
alternative configuration, evaluated at EVERY frame.

SUPERSEDES fig6C.py, which is kept for provenance.

Why this exists. fig6C.py hardcoded 11, 0 and 20 bridging frames out of a
21-frame sample taken at 5 ns intervals, and asserted their totals against a
"Table 17" that no longer exists in the manuscript. The paper now evaluates
simultaneous contact at every frame, and Table 5 reports the resulting
percentages. A panel drawn on the retired 21-frame sample beside a table of
every-frame percentages reads as one measurement disagreeing with itself.

What a bridging frame is. The ligand is simultaneously within 4.5 A of the
pharmacophoric side chains (Asp121 and Arg142) and of P_0. The two distance
series per replicate are gmx mindist output:

    analysis/sc_KeyRes.xvg   ligand heavy atoms to the pharmacophoric side chains
    analysis/md_P0.xvg       ligand heavy atoms to the P_0 lining

Both carry 10,001 frames at 10 ps over 100 ns. Only the files written by the
10 November 2026 analysis pass are used. Several August files sit in the same
directories, were computed in the docked-pose rep1 working directory rather
than here, and two of them are truncated, so this script refuses any series
whose length does not match its partner.

Nothing is hardcoded. The occupancies are computed from the trajectories and
then checked against the values printed in Table 5. The script writes no file
if any of them disagrees.

Usage
    python3 fig6C_occupancy.py [MD_ROOT]
    EXPORT_DPI=600 python3 fig6C_occupancy.py      # for the 600 dpi figure set

Writes Figure6C.png, which assemble_fig6.py takes as panel C.
"""
import os
import sys
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MDROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    "~/Research/05_MD_v2")
DPI = float(os.environ.get("EXPORT_DPI", "300"))

ORDER = ["LOA22-B1", "LOA22-B2", "LOA22-B3"]
REPS = ["rep1_arg142pose", "rep2_arg142pose", "rep3_arg142pose"]
COL = {"LOA22-B1": "#C0392B", "LOA22-B2": "#2E86C1", "LOA22-B3": "#27AE60"}
INK, MUTED = "#1a1a1a", "#5a5a5a"
CUTOFF_NM = 0.45

# Table 5, "Simultaneous contact (% of frames)". The script aborts rather than
# write a figure that disagrees with the manuscript.
EXPECTED = {"LOA22-B1": (59.2, 52.3),
            "LOA22-B2": (99.7, 0.3),
            "LOA22-B3": (90.8, 7.0)}


def series(path):
    """Second column of a gmx .xvg, as a float array."""
    vals = []
    with open(path) as fh:
        for line in fh:
            if line[:1] in "@#":
                continue
            parts = line.split()
            if len(parts) >= 2:
                vals.append(float(parts[1]))
    return np.asarray(vals)


def occupancy(compound, rep):
    d = os.path.join(MDROOT, compound, rep, "analysis")
    key = os.path.join(d, "sc_KeyRes.xvg")
    p0 = os.path.join(d, "md_P0.xvg")
    for f in (key, p0):
        if not os.path.exists(f):
            sys.exit("ABORTED: missing %s" % f)
    a, b = series(key), series(p0)
    if a.size != b.size:
        sys.exit("ABORTED: %s %s row mismatch, sc_KeyRes %d against md_P0 %d. "
                 "One of these is a truncated August file."
                 % (compound, rep, a.size, b.size))
    if a.size < 9000:
        sys.exit("ABORTED: %s %s has only %d frames, expected about 10,001"
                 % (compound, rep, a.size))
    bridging = int(np.count_nonzero((a <= CUTOFF_NM) & (b <= CUTOFF_NM)))
    return 100.0 * bridging / a.size, bridging, a.size


print("computing bridging occupancy at every frame, cutoff %.2f nm" % CUTOFF_NM)
DATA, NFRAMES = {}, {}
for name in ORDER:
    pcts = []
    for rep in REPS:
        pct, b, n = occupancy(name, rep)
        pcts.append(pct)
        print("  %-10s %-18s %6d/%-6d = %6.2f%%" % (name, rep, b, n, pct))
    DATA[name] = pcts
    NFRAMES[name] = n

print()
print("verifying against Table 5")
for name in ORDER:
    v = DATA[name]
    mean = sum(v) / len(v)
    sd = math.sqrt(sum((x - mean) ** 2 for x in v) / (len(v) - 1))
    want_m, want_s = EXPECTED[name]
    ok = abs(mean - want_m) <= 0.05 and abs(sd - want_s) <= 0.05
    print("  %-10s %5.1f +/- %4.1f  vs Table 5 %5.1f +/- %4.1f  %s"
          % (name, mean, sd, want_m, want_s, "ok" if ok else "MISMATCH"))
    if not ok:
        sys.exit("ABORTED: %s does not reproduce Table 5" % name)
print("all gates passed")

# ---------------------------------------------------------------- plot
fig, ax = plt.subplots(figsize=(7.2, 4.4))
x = np.arange(3)
w = 0.26

for i, name in enumerate(ORDER):
    off = (i - 1) * w
    vals = DATA[name]
    mean = sum(vals) / len(vals)
    sd = math.sqrt(sum((q - mean) ** 2 for q in vals) / (len(vals) - 1))
    ax.bar(x + off, vals, w, color=COL[name], zorder=3,
           label="%s  (%.1f %s %.1f%%)" % (name, mean, u"±", sd))
    for xi, v in zip(x + off, vals):
        if v < 1.0:
            # a bar this short is invisible, so mark it and give the value
            ax.text(xi, 2.2, "%.1f" % v, ha="center", va="bottom", fontsize=9.5,
                    color=COL[name], fontweight="bold", zorder=4)
            ax.plot([xi - w / 2.4, xi + w / 2.4], [0.6, 0.6], color=COL[name],
                    lw=2.4, solid_capstyle="butt", zorder=4)
        else:
            ax.text(xi, v + 1.6, "%.1f" % v, ha="center", va="bottom",
                    fontsize=9.5, color=INK, zorder=4)

ax.axhline(100, ls="--", color="#888888", lw=1, zorder=2)
ax.set_xlim(-0.55, 2.95)
ax.text(2.92, 101.2, "every frame", ha="right", va="bottom",
        fontsize=8, color=MUTED)

ax.set_xticks(x)
ax.set_xticklabels(["Replicate 1", "Replicate 2", "Replicate 3"], fontsize=10,
                   color=INK)
ax.set_ylabel("Bridging occupancy (% of frames)", fontsize=10, color=INK)
ax.set_ylim(0, 112)
ax.set_yticks([0, 25, 50, 75, 100])
ax.tick_params(axis="y", labelsize=9, colors=MUTED, length=3)
ax.tick_params(axis="x", length=0)
ax.spines[["top", "right"]].set_visible(False)
ax.spines["left"].set_color("#cccccc")
ax.spines["bottom"].set_color("#cccccc")
ax.grid(axis="y", color="#ececec", lw=0.7, zorder=0)
ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=9, ncol=3, loc="lower center",
          bbox_to_anchor=(0.5, -0.30), columnspacing=1.8)

plt.tight_layout()
plt.savefig("Figure6C.png", dpi=DPI, facecolor="white")
print()
print("wrote Figure6C.png at %d dpi (%d frames per replicate)"
      % (DPI, NFRAMES[ORDER[0]]))
