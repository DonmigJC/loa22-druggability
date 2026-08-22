#!/usr/bin/env python3
"""
Prior-only control analysis.

Compares the docking-score distribution of molecules sampled from the
UNTRAINED prior against the RL-optimised library, to quantify what
reinforcement learning actually contributed.

Run:  ~/miniconda3/envs/reinvent4/bin/python ~/Research/compare_prior.py
"""
import os
import pandas as pd
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import FilterCatalog

RDLogger.DisableLog("rdApp.*")
R = os.path.expanduser("~/Research")
DOXY = -6.862          # corrected reference, 3 seeds, exhaustiveness 32
FILTER = -7.5          # library filtering threshold

# ---------------------------------------------------------------- load
prior = pd.read_csv(os.path.join(R, "prior_control", "prior_docked.csv"))
prior = prior.dropna(subset=["affinity"])

rl = pd.read_csv(os.path.join(R, "loa22_generation_1.csv"))
rl = rl.dropna(subset=["SMILES"])
rl = rl[rl["SMILES_state"] == 1].copy()
rl["affinity"] = -12.0 * rl["Vina_Loa22 (raw)"]      # normalised -> kcal/mol

print("=" * 66)
print("PRIOR-ONLY CONTROL vs RL-OPTIMISED LIBRARY")
print("=" * 66)
print(f"\nprior sampled & docked : {len(prior)}")
print(f"RL library (valid)     : {len(rl)}")

# ------------------------------------------------------- distributions
print("\n" + "-" * 66)
print(f"{'metric':<34}{'PRIOR':>15}{'RL':>15}")
print("-" * 66)


def row(label, fn, fmt="{:>15.3f}"):
    print(f"{label:<34}" + fmt.format(fn(prior.affinity))
          + fmt.format(fn(rl.affinity)))


row("mean affinity (kcal/mol)", np.mean)
row("median", np.median)
row("std dev", np.std)
row("best (most negative)", np.min)
row("25th percentile", lambda x: np.percentile(x, 25))
row("10th percentile", lambda x: np.percentile(x, 10))

print("-" * 66)
for label, thr in [("beat doxycycline (-6.862)", DOXY),
                   ("pass filter (<= -7.5)", FILTER),
                   ("strong binders (<= -8.5)", -8.5),
                   ("very strong (<= -9.0)", -9.0)]:
    p = (prior.affinity <= thr).sum()
    r = (rl.affinity <= thr).sum()
    print(f"{label:<34}{p:>7} ({100*p/len(prior):>5.2f}%)"
          f"{r:>7} ({100*r/len(rl):>5.2f}%)")

# ------------------------------------------------------- enrichment
print("-" * 66)
pd_ = (prior.affinity <= DOXY).mean()
rd_ = (rl.affinity <= DOXY).mean()
pf_ = (prior.affinity <= FILTER).mean()
rf_ = (rl.affinity <= FILTER).mean()
print(f"{'enrichment vs doxycycline':<34}{'':>15}{rd_/pd_ if pd_ else float('inf'):>14.2f}x")
print(f"{'enrichment vs -7.5 filter':<34}{'':>15}{rf_/pf_ if pf_ else float('inf'):>14.2f}x")
print(f"{'mean affinity shift':<34}{'':>15}"
      f"{rl.affinity.mean()-prior.affinity.mean():>14.3f} kcal/mol")

# ------------------------------------------------------- significance
try:
    from scipy.stats import mannwhitneyu
    u, p = mannwhitneyu(rl.affinity, prior.affinity, alternative="less")
    n1, n2 = len(rl), len(prior)
    cliffs = 2 * u / (n1 * n2) - 1          # effect size, -1 to +1
    print("-" * 66)
    print(f"Mann-Whitney U      : {u:.0f}")
    print(f"p-value             : {p:.3e}")
    print(f"Cliff's delta       : {abs(cliffs):.3f}"
          "   (0.15 small, 0.33 medium, 0.47 large)")
except ImportError:
    print("\n(scipy unavailable - skipping significance test)")

# ------------------------------------------------------- quality check
print("-" * 66)
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
cat = FilterCatalog.FilterCatalog(params)

for name, df, qcol in [("PRIOR", prior, "QED"), ("RL", rl, "QED (raw)")]:
    smi = df["SMILES"]
    mols = [Chem.MolFromSmiles(s) for s in smi]
    ok = [m for m in mols if m]
    uniq = len({Chem.MolToSmiles(m) for m in ok})
    pains = sum(1 for m in ok if cat.HasMatch(m))
    print(f"{name:<6} unique {uniq:>5}/{len(ok):<5}  "
          f"PAINS {pains:>4} ({100*pains/len(ok):>5.2f}%)  "
          f"mean QED {df[qcol].mean():.3f}")

# ------------------------------------------------------- plot
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))
    bins = np.linspace(-11, 0, 60)
    ax.hist(prior.affinity, bins=bins, alpha=0.6, label=f"Prior (n={len(prior)})",
            color="#9E9E9E", edgecolor="white", density=True)
    ax.hist(rl.affinity, bins=bins, alpha=0.6, label=f"RL-optimised (n={len(rl)})",
            color="#E91E63", edgecolor="white", density=True)
    ax.axvline(DOXY, color="orange", ls="--", lw=2,
               label=f"Doxycycline ({DOXY})")
    ax.axvline(FILTER, color="red", ls="--", lw=2, label="Filter (-7.5)")
    ax.set_xlabel("AutoDock Vina score (kcal/mol)")
    ax.set_ylabel("Density")
    ax.set_title("Reinforcement learning shifts the affinity distribution",
                 fontweight="bold")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    out = os.path.join(R, "prior_control", "prior_vs_rl.png")
    plt.savefig(out, dpi=300)
    print(f"\nplot written: {out}")
except Exception as e:
    print(f"\n(plot skipped: {e})")

print("=" * 66)
