import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

GREY, RED = "#8C8C8C", "#C0392B"
fig, ax = plt.subplots(2, 2, figsize=(11, 8))

# --- A: learning curves
b = pd.read_csv("loa22_generation_1.csv"); b = b[b.SMILES_state==1]
c = pd.read_csv("arm_c/armC_generation_1.csv"); c = c[c.SMILES_state==1]
ax[0,0].plot(b.groupby("step").Score.mean(), color=GREY, lw=2, label="Initial campaign")
ax[0,0].plot(c.groupby("step").Score.mean(), color=RED, lw=2, label="Revised objective")
ax[0,0].set_xlabel("Training step"); ax[0,0].set_ylabel("Mean composite score")
ax[0,0].legend(frameon=False); ax[0,0].set_title("A", loc="left", fontweight="bold")

# --- B: affinity distributions
p = pd.read_csv("prior_control/prior_docked.csv").dropna(subset=["affinity"])
b["kcal"] = -12.0*b["Vina_Loa22 (raw)"]
bins = np.linspace(-11, -2, 60)
ax[0,1].hist(p.affinity, bins=bins, color=GREY, alpha=.65, density=True, label=f"Untrained prior (n={len(p)})")
ax[0,1].hist(b.kcal, bins=bins, color=RED, alpha=.65, density=True, label=f"Trained agent (n={len(b)})")
ax[0,1].axvline(-6.862, ls="--", color="k", lw=1.5, label="Doxycycline")
ax[0,1].set_xlabel("Predicted affinity (kcal/mol)"); ax[0,1].set_ylabel("Density")
ax[0,1].legend(frameon=False, fontsize=8); ax[0,1].set_title("B", loc="left", fontweight="bold")

# --- C: threshold fractions
th = [-6.862, -7.5, -8.0, -8.5, -9.0]
pf = [100*(p.affinity<=t).mean() for t in th]
bf = [100*(b.kcal<=t).mean() for t in th]
x = np.arange(len(th)); w=0.38
ax[1,0].bar(x-w/2, pf, w, color=GREY, label="Prior")
ax[1,0].bar(x+w/2, bf, w, color=RED, label="Trained")
ax[1,0].set_xticks(x); ax[1,0].set_xticklabels([f"≤{t}" for t in th], fontsize=8)
ax[1,0].set_ylabel("Molecules (%)"); ax[1,0].set_xlabel("Affinity threshold (kcal/mol)")
ax[1,0].legend(frameon=False); ax[1,0].set_title("C", loc="left", fontweight="bold")

# --- D: QED
ax[1,1].hist(p.QED.dropna(), bins=40, color=GREY, alpha=.65, density=True, label="Prior")
ax[1,1].hist(b["QED (raw)"], bins=40, color=RED, alpha=.65, density=True, label="Trained")
ax[1,1].set_xlabel("QED"); ax[1,1].set_ylabel("Density")
ax[1,1].legend(frameon=False); ax[1,1].set_title("D", loc="left", fontweight="bold")

for a in ax.flat: a.spines[["top","right"]].set_visible(False)
plt.tight_layout(); plt.savefig("Figure3.png", dpi=300)
print("wrote Figure3.png")
