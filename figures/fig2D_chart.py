#!/usr/bin/env python3
"""
Figure 2 panel D - quantitative comparison of all six detected cavities.

Replaces the attempted cross-section render. Enclosure is a number between 0
and 1, and the claim the panel has to carry is that the functional site sits
at exactly 0.00 while the druggable cavity sits at 0.24. A shaded surface can
only ever hint at that. A chart states it.

All values read directly from the DoGSiteScorer output, not retyped:
    DoGSiteScorer Results/*_desc.txt
"""
import csv, glob, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, RED, GREY = "#3380CC", "#B22222", "#8C8C8C"
INK, MUTED = "#1a1a1a", "#5a5a5a"

# ---- read the descriptor table (single source of truth)
f = glob.glob(os.path.join("DoGSiteScorer Results", "*_desc.txt"))
if not f:
    raise SystemExit("desc.txt not found - run from the Research root")
rows = list(csv.reader(open(f[0]), delimiter="\t"))
idx = {h: i for i, h in enumerate(rows[0])}
P = []
for r in rows[1:]:
    if not r:
        continue
    P.append(dict(name=r[idx["name"]],
                  vol=float(r[idx["volume"]]),
                  enc=float(r[idx["enclosure"]]),
                  ds=float(r[idx["drugScore"]])))
P.sort(key=lambda p: -p["vol"])

def style(name):
    if name == "P_0":
        return BLUE, "P_0  druggable cavity"
    if name == "P_4":
        return RED, "P_4  functional site"
    return GREY, name

labels = [style(p["name"])[1] for p in P]
cols = [style(p["name"])[0] for p in P]
y = range(len(P))[::-1]

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0))

# ---- left: volume
a = ax[0]
a.axvspan(300, 450, color="#000000", alpha=0.055, lw=0, zorder=0)
a.text(375, 1.55, "drug-sized ligand", ha="center", va="center",
       rotation=90, fontsize=7.5, color=MUTED)

a.barh(list(y), [p["vol"] for p in P], color=cols, height=0.62, zorder=2)
for yy, p in zip(y, P):
    a.text(p["vol"] + 18, yy, f'{p["vol"]:.2f}', va="center",
           fontsize=8.5, color=INK)
a.set_xlim(0, 1080)
a.set_xlabel("Cavity volume (Å$^3$)   ·   shaded band = volume a 300-500 Da ligand requires",
             fontsize=9.5, color=INK)
a.set_title("Volume", loc="left", fontsize=10, color=INK, fontweight="bold")

# ---- right: enclosure
b = ax[1]
b.barh(list(y), [p["enc"] for p in P], color=cols, height=0.62, zorder=2)
for yy, p in zip(y, P):
    if p["enc"] == 0.0:
        b.text(0.008, yy, "0.00  fully solvent-exposed", va="center",
               fontsize=8.5, color=RED, fontweight="bold")
    else:
        b.text(p["enc"] + 0.008, yy, f'{p["enc"]:.2f}', va="center",
               fontsize=8.5, color=INK)
b.set_xlim(0, 0.42)
b.set_xlabel("Enclosure  (0 = planar surface, 1 = fully buried)",
             fontsize=9.5, color=INK)
b.set_title("Enclosure", loc="left", fontsize=10, color=INK, fontweight="bold")

for a_ in ax:
    a_.set_yticks(list(y))
    a_.set_yticklabels(labels, fontsize=9, color=INK)
    a_.spines[["top", "right", "left"]].set_visible(False)
    a_.spines["bottom"].set_color("#cccccc")
    a_.tick_params(axis="x", labelsize=8.5, colors=MUTED, length=3)
    a_.tick_params(axis="y", length=0)
    a_.grid(axis="x", color="#e8e8e8", lw=0.7, zorder=0)
    a_.set_axisbelow(True)

plt.tight_layout()
plt.savefig("fig2D_chart.png", dpi=300, facecolor="white")
print("wrote fig2D_chart.png")
for p in P:
    print(f'  {p["name"]}  vol {p["vol"]:>7.2f}  enc {p["enc"]:.2f}  drugScore {p["ds"]:.3f}')
