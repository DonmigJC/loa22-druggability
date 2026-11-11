#!/usr/bin/env python3
"""
Correlation of predicted binding affinity with molecular weight and drug-likeness
in the reinforcement-learning library (Loa22 druggability study).

Supports the statement in Discussion section 3.2 that geometric-mean aggregation
disadvantaged high-affinity molecules because those molecules are larger and
less drug-like. Produces Supplementary Table S15.

Input : 09_data/loa22_generation_1.csv   (REINVENT 4 staged-RL output, 6,400 rows)
Output: SupplementaryTableS15.csv
        affinity_property_correlation_results.txt

Note on the affinity column. REINVENT stores the transformed reward in
"Vina_Loa22 (raw)". The manuscript's reward transform is
    reward = clip(affinity, -12, 0) / -12
so predicted affinity in kcal/mol is recovered as  affinity = -12 * reward.
This reconstruction reproduces every figure quoted in the manuscript:
best -10.660 kcal/mol, mean -6.608 over all valid molecules, mean -6.654 over
docked molecules only, n = 6,319 valid and n = 6,275 docked, mean QED 0.691.

Usage:  python affinity_property_correlation.py ../09_data/loa22_generation_1.csv
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

CSV = sys.argv[1] if len(sys.argv) > 1 else "../09_data/loa22_generation_1.csv"
REWARD_CLIP = -12.0


def fisher_ci(r, n, alpha=0.05):
    """Two-sided confidence interval for a Pearson r via Fisher z."""
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    crit = stats.norm.ppf(1 - alpha / 2)
    return np.tanh(z - crit * se), np.tanh(z + crit * se)


def correlate(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)
    pr, pp = stats.pearsonr(x, y)
    sr, sp = stats.spearmanr(x, y)
    lo, hi = fisher_ci(pr, n)
    return dict(n=n, pearson_r=pr, pearson_ci_low=lo, pearson_ci_high=hi,
                pearson_p=pp, spearman_rho=sr, spearman_p=sp)


def main():
    df = pd.read_csv(CSV)

    valid = df[df["SMILES_state"] == 1].copy()
    valid["affinity"] = REWARD_CLIP * valid["Vina_Loa22 (raw)"]
    docked = valid[valid["Vina_Loa22 (raw)"] > 0].copy()

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("Reconstruction check against the manuscript")
    out("-" * 78)
    out(f"  valid molecules (SMILES_state == 1)      {len(valid):>10d}   manuscript 6,319")
    out(f"  docked molecules (reward > 0)            {len(docked):>10d}   manuscript 6,275")
    out(f"  best predicted affinity (kcal/mol)       {valid['affinity'].min():>10.3f}   manuscript -10.660")
    out(f"  mean affinity, all valid (kcal/mol)      {valid['affinity'].mean():>10.3f}   manuscript  -6.608")
    out(f"  mean affinity, docked only (kcal/mol)    {docked['affinity'].mean():>10.3f}   manuscript  -6.654")
    out(f"  mean QED, all valid                      {valid['QED (raw)'].mean():>10.3f}   manuscript   0.691")
    out()

    rows = []
    out("Library-wide correlations (docked molecules only)")
    out("-" * 78)
    pairs = [("Predicted affinity", "Molecular weight", "affinity", "MW (raw)"),
             ("Predicted affinity", "QED", "affinity", "QED (raw)"),
             ("Molecular weight", "QED", "MW (raw)", "QED (raw)")]
    for label_x, label_y, cx, cy in pairs:
        r = correlate(docked[cx].to_numpy(float), docked[cy].to_numpy(float))
        rows.append(dict(comparison=f"{label_x} vs {label_y}", scope="full docked library", **r))
        out(f"  {label_x} vs {label_y:<18s} n = {r['n']:<6d} "
            f"Pearson r = {r['pearson_r']:+.3f} "
            f"[{r['pearson_ci_low']:+.3f}, {r['pearson_ci_high']:+.3f}], "
            f"Spearman rho = {r['spearman_rho']:+.3f}")
    out()

    out("High-affinity tail against the remainder of the library")
    out("-" * 78)
    tail_rows = []
    for threshold in (-8.5, -9.0):
        hi = docked[docked["affinity"] <= threshold]
        rest = docked[docked["affinity"] > threshold]
        mw_p = stats.mannwhitneyu(hi["MW (raw)"], rest["MW (raw)"]).pvalue
        qed_p = stats.mannwhitneyu(hi["QED (raw)"], rest["QED (raw)"]).pvalue
        out(f"  affinity <= {threshold} kcal/mol   n = {len(hi)} "
            f"({100 * len(hi) / len(docked):.2f}% of the docked library)")
        out(f"      mean molecular weight   {hi['MW (raw)'].mean():7.2f} Da  against "
            f"{rest['MW (raw)'].mean():7.2f} Da   Mann-Whitney p = {mw_p:.3g}")
        out(f"      mean QED                {hi['QED (raw)'].mean():7.3f}     against "
            f"{rest['QED (raw)'].mean():7.3f}       Mann-Whitney p = {qed_p:.3g}")
        tail_rows.append(dict(threshold_kcal_mol=threshold, n_tail=len(hi),
                              pct_of_docked_library=100 * len(hi) / len(docked),
                              mean_MW_tail=hi["MW (raw)"].mean(),
                              mean_MW_remainder=rest["MW (raw)"].mean(),
                              MW_mannwhitney_p=mw_p,
                              mean_QED_tail=hi["QED (raw)"].mean(),
                              mean_QED_remainder=rest["QED (raw)"].mean(),
                              QED_mannwhitney_p=qed_p))

    pd.concat([pd.DataFrame(rows), pd.DataFrame(tail_rows)], axis=0,
              ignore_index=True).to_csv("SupplementaryTableS15.csv", index=False)
    with open("affinity_property_correlation_results.txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\nWrote SupplementaryTableS15.csv and affinity_property_correlation_results.txt")


if __name__ == "__main__":
    main()
