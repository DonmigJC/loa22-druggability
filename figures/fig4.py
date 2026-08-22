import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = pd.read_csv("dual_site/dual_site_results.csv").dropna(subset=["aff_site","aff_p0","bridge"])
s = pd.read_csv("pharmacophore_split/armB_split.csv").dropna(subset=["d_asp","d_arg"])
GREY, RED, BLUE = "#8C8C8C", "#C0392B", "#2E86C1"
fig, ax = plt.subplots(2, 2, figsize=(11, 9))

# A: affinity vs affinity
ax[0,0].scatter(d.aff_p0, d.aff_site, s=3, alpha=.25, color=GREY, edgecolors="none")
lim = [-11, -2]
ax[0,0].plot(lim, lim, "k--", lw=1)
ax[0,0].set_xlim(lim); ax[0,0].set_ylim(lim)
ax[0,0].set_xlabel("Affinity at P_0 (kcal/mol)")
ax[0,0].set_ylabel("Affinity at functional site (kcal/mol)")
n_pref = (d.aff_site < d.aff_p0).sum()
ax[0,0].text(.05,.93,f"{n_pref} of {len(d)} ({100*n_pref/len(d):.2f}%)\nprefer functional site",
             transform=ax[0,0].transAxes, fontsize=9, va="top")
ax[0,0].set_title("A", loc="left", fontweight="bold")

# B: penalty distribution
pen = d.aff_site - d.aff_p0
ax[0,1].hist(pen, bins=60, color=BLUE, alpha=.75)
ax[0,1].axvline(pen.mean(), color=RED, lw=2, label=f"mean {pen.mean():.3f}")
ax[0,1].axvline(0, color="k", ls="--", lw=1)
ax[0,1].set_xlabel("Affinity penalty at functional site (kcal/mol)")
ax[0,1].set_ylabel("Molecules"); ax[0,1].legend(frameon=False)
ax[0,1].set_title("B", loc="left", fontweight="bold")

# C: bridging distribution
ax[1,0].hist(d.bridge, bins=60, color=GREY, alpha=.8)
for t,lab in [(0.5,"0.5"),(0.7,"0.7"),(0.8,"0.8")]:
    ax[1,0].axvline(t, ls="--", color=RED, lw=1)
    ax[1,0].text(t, ax[1,0].get_ylim()[1]*.9, lab, fontsize=8, color=RED)
ax[1,0].set_xlabel("Bridging score"); ax[1,0].set_ylabel("Molecules")
ax[1,0].set_yscale("log")
ax[1,0].set_title("C", loc="left", fontweight="bold")

# D: Asp121 vs Arg142 distance
ax[1,1].scatter(s.d_asp, s.d_arg, s=3, alpha=.25, color=GREY, edgecolors="none")
ax[1,1].axvline(4.5, ls="--", color=RED, lw=1)
ax[1,1].axhline(4.5, ls="--", color=RED, lw=1)
ax[1,1].set_xlabel("Distance to Asp121 (Å)")
ax[1,1].set_ylabel("Distance to Arg142 (Å)")
both = ((s.d_asp<4.5)&(s.d_arg<4.5)).sum()
ax[1,1].text(.55,.93,f"both contacted:\n{both} of {len(s)}",
             transform=ax[1,1].transAxes, fontsize=9, va="top", color=RED)
ax[1,1].set_title("D", loc="left", fontweight="bold")

for a in ax.flat: a.spines[["top","right"]].set_visible(False)
plt.tight_layout(); plt.savefig("Figure4.png", dpi=300)
print("wrote Figure4.png")
