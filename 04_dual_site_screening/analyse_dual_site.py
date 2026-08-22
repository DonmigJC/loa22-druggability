#!/usr/bin/env python3
"""
DUAL-SITE RESULTS ANALYSIS AND LEAD SELECTION
=============================================

Ranks the library on a composite of bridging geometry and binding
affinity, benchmarks the original three leads, and selects a new
candidate set for molecular dynamics.

Run:  ~/miniconda3/envs/reinvent4/bin/python analyse_dual_site.py
"""
import os
import pandas as pd
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, QED, rdMolDescriptors, FilterCatalog

RDLogger.DisableLog("rdApp.*")

R = os.path.expanduser("~/Research")
CSV = os.path.join(R, "dual_site", "dual_site_results.csv")
OUT = os.path.join(R, "dual_site")

ORIGINAL = {
    "LOA22-1": "CC(C)CNC(=O)Nc1ccc2nnc(-c3cc(C(F)(F)F)cc(C(F)(F)F)c3)n2n1",
    "LOA22-2": "CC(CC(=O)N1CCc2ccc(C(F)(F)F)cc2C1)c1ccc(O)cc1",
    "LOA22-3": "CC(=O)Nc1ccc(F)c(C(=O)NS(=O)(=O)c2ccccc2Cl)c1",
}

# ------------------------------------------------------------------ load
d = pd.read_csv(CSV)
print("=" * 74)
print("DUAL-SITE SCREENING RESULTS")
print("=" * 74)
print(f"rows in file            : {len(d)}")
d = d.dropna(subset=["aff_site", "aff_p0", "d_key", "d_p0"])
print(f"successfully docked     : {len(d)}")

# ------------------------------------------------------- site comparison
print()
print("-" * 74)
print("SITE COMPARISON")
print("-" * 74)
print(f"mean affinity, P_0             : {d.aff_p0.mean():7.3f} kcal/mol")
print(f"mean affinity, functional site : {d.aff_site.mean():7.3f} kcal/mol")
print(f"mean penalty at functional site: {(d.aff_site - d.aff_p0).mean():7.3f} kcal/mol")
n_pref = (d.aff_site < d.aff_p0).sum()
print(f"molecules preferring the site  : {n_pref} ({100*n_pref/len(d):.2f}%)")
print(f"best P_0 affinity              : {d.aff_p0.min():7.3f}")
print(f"best functional-site affinity  : {d.aff_site.min():7.3f}")

# ------------------------------------------------------- bridging
print()
print("-" * 74)
print("BRIDGING")
print("-" * 74)
for t in (0.3, 0.5, 0.7, 0.8, 0.9):
    n = (d.bridge > t).sum()
    print(f"  bridge > {t:.1f} : {n:5d}  ({100*n/len(d):5.2f}%)")
print(f"\n  closest approach to Asp121/Arg142 : {d.d_key.min():.2f} A")
print(f"  median approach                   : {d.d_key.median():.2f} A")

# ------------------------------------------------------- composite score
# normalise both terms to 0-1 then take the geometric mean, so a molecule
# must satisfy BOTH geometry and affinity
aff_norm = ((-d.aff_site) - (-d.aff_site).min()) / \
           ((-d.aff_site).max() - (-d.aff_site).min())
d["composite"] = np.sqrt(d.bridge.clip(0, 1) * aff_norm.clip(0, 1))

# ------------------------------------------------------- original leads
print()
print("-" * 74)
print("ORIGINAL LEADS — where they rank")
print("-" * 74)
print(f"{'lead':<10}{'bridge':>8}{'rank':>7}{'aff_site':>10}{'aff_p0':>9}"
      f"{'d_key':>7}{'d_p0':>7}{'composite':>11}{'rank':>7}")
for name, smi in ORIGINAL.items():
    m = d[d.SMILES == smi]
    if not len(m):
        print(f"{name:<10} not found in results")
        continue
    r = m.iloc[0]
    br = int((d.bridge > r.bridge).sum()) + 1
    cr = int((d.composite > r.composite).sum()) + 1
    print(f"{name:<10}{r.bridge:>8.3f}{br:>7}{r.aff_site:>10.3f}{r.aff_p0:>9.3f}"
          f"{r.d_key:>7.2f}{r.d_p0:>7.2f}{r.composite:>11.3f}{cr:>7}")

# ------------------------------------------------------- drug-likeness
print()
print("-" * 74)
print("TOP 25 BY COMPOSITE — with drug-likeness filters")
print("-" * 74)
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
cat = FilterCatalog.FilterCatalog(params)

top = d.sort_values("composite", ascending=False).head(300).copy()
recs = []
for _, r in top.iterrows():
    m = Chem.MolFromSmiles(r.SMILES)
    if m is None:
        continue
    if cat.HasMatch(m):
        continue
    q = QED.qed(m)
    mw = Descriptors.MolWt(m)
    if q < 0.4 or mw > 600:
        continue
    recs.append(dict(SMILES=r.SMILES, bridge=r.bridge, aff_site=r.aff_site,
                     aff_p0=r.aff_p0, d_key=r.d_key, d_p0=r.d_p0,
                     composite=r.composite, QED=q, MW=mw,
                     logP=Descriptors.MolLogP(m),
                     HA=m.GetNumHeavyAtoms(),
                     TPSA=Descriptors.TPSA(m),
                     LE=abs(r.aff_site) / m.GetNumHeavyAtoms()))

sel = pd.DataFrame(recs)
print(f"passing PAINS, QED>=0.4, MW<=600 : {len(sel)} of top 300\n")
print(f"{'#':>3}{'brdg':>7}{'site':>8}{'P_0':>8}{'dkey':>6}{'dp0':>6}"
      f"{'QED':>6}{'MW':>7}{'LE':>6}  SMILES")
print("-" * 100)
for i, (_, r) in enumerate(sel.head(25).iterrows(), 1):
    print(f"{i:>3}{r.bridge:>7.3f}{r.aff_site:>8.3f}{r.aff_p0:>8.3f}"
          f"{r.d_key:>6.2f}{r.d_p0:>6.2f}{r.QED:>6.3f}{r.MW:>7.1f}{r.LE:>6.3f}"
          f"  {r.SMILES[:46]}")

sel.to_csv(os.path.join(OUT, "bridging_candidates.csv"), index=False)
d.sort_values("composite", ascending=False).to_csv(
    os.path.join(OUT, "dual_site_ranked.csv"), index=False)
print(f"\nwritten: bridging_candidates.csv ({len(sel)} molecules)")
print(f"written: dual_site_ranked.csv (all {len(d)}, ranked)")

# ------------------------------------------------------- MD proposal
print()
print("=" * 74)
print("PROPOSED CANDIDATES FOR MOLECULAR DYNAMICS")
print("=" * 74)
if len(sel) >= 3:
    picks = []
    # 1. best composite
    picks.append(("best overall", sel.iloc[0]))
    # 2. closest approach to the pharmacophore
    picks.append(("closest to Asp121/Arg142",
                  sel.sort_values("d_key").iloc[0]))
    # 3. strongest functional-site binder
    picks.append(("strongest site affinity",
                  sel.sort_values("aff_site").iloc[0]))
    seen = set()
    for label, r in picks:
        if r.SMILES in seen:
            continue
        seen.add(r.SMILES)
        print(f"\n{label}")
        print(f"  SMILES   {r.SMILES}")
        print(f"  bridge {r.bridge:.3f}   site {r.aff_site:.3f}   "
              f"P_0 {r.aff_p0:.3f}   d_key {r.d_key:.2f} A   d_p0 {r.d_p0:.2f} A")
        print(f"  QED {r.QED:.3f}   MW {r.MW:.1f}   logP {r.logP:.2f}   "
              f"TPSA {r.TPSA:.1f}   LE {r.LE:.3f}")

print()
print("Benchmark — LOA22-3 (current best, discovered by chance):")
m = d[d.SMILES == ORIGINAL["LOA22-3"]]
if len(m):
    r = m.iloc[0]
    print(f"  bridge {r.bridge:.3f}   site {r.aff_site:.3f}   "
          f"d_key {r.d_key:.2f} A   d_p0 {r.d_p0:.2f} A")
print("=" * 74)
