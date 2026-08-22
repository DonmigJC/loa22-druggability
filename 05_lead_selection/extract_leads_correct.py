#!/usr/bin/env python3
# File: ~/Research/extract_leads_correct.py
# Uses exact column names from loa22_generation_1.csv

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit.Chem import rdMolDescriptors, FilterCatalog
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("=" * 60)
print("LOA22 LEAD EXTRACTION — REINVENT4 RESULTS")
print("=" * 60)

# ── Load CSV ────────────────────────────────────────────────
df = pd.read_csv('loa22_generation_1.csv')
print(f"Total rows loaded: {len(df)}")

# Remove failed SMILES
df = df.dropna(subset=['SMILES'])
df = df[df['SMILES_state'] == 1].copy()
print(f"Valid SMILES: {len(df)}")

# ── Key columns from YOUR CSV ───────────────────────────────
# Vina_Loa22 (raw) = actual Vina affinity in kcal/mol (no conversion needed)
# QED (raw)        = actual QED value (0-1)
# MW (raw)         = actual molecular weight in Da
# Score            = composite REINVENT4 score (0-1)

# Rename for convenience
df['Affinity_kcal_mol'] = df['Vina_Loa22 (raw)']
df['QED_value']         = df['QED (raw)']
df['MW_value']          = df['MW (raw)']

# ── Summary statistics ──────────────────────────────────────
print("\n=== RAW VINA AFFINITY DISTRIBUTION ===")
print(df['Affinity_kcal_mol'].describe().round(3))

print("\n=== QED DISTRIBUTION ===")
print(df['QED_value'].describe().round(3))

print("\n=== COMPOSITE SCORE DISTRIBUTION ===")
print(df['Score'].describe().round(4))

# How many molecules beat doxycycline?
doxycy_affinity = -6.918
better_than_doxy = (df['Affinity_kcal_mol'] <= doxycy_affinity).sum()
print(f"\nMolecules beating doxycycline ({doxycy_affinity} kcal/mol): "
      f"{better_than_doxy}")

# ── PAINS filter ────────────────────────────────────────────
print("\nApplying PAINS structural alert filter...")
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(
    FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
catalog = FilterCatalog.FilterCatalog(params)

pains_pass = []
for smi in df['SMILES']:
    mol = Chem.MolFromSmiles(smi)
    if mol and not catalog.HasMatch(mol):
        pains_pass.append(True)
    else:
        pains_pass.append(False)

df['Passes_PAINS'] = pains_pass
print(f"Molecules passing PAINS: {sum(pains_pass)} / {len(df)}")

# ── Apply quality filters ───────────────────────────────────
print("\nApplying quality filters...")

filtered = df[
    (df['Affinity_kcal_mol'] <= -7.5) &   # Beats doxycycline threshold
    (df['QED_value'] >= 0.5)             & # Drug-like
    (df['MW_value'] <= 600)              & # Reasonable size
    (df['Passes_PAINS'] == True)           # No structural alerts
].copy()

print(f"After all filters: {len(filtered)} candidates remain")

# ── Top 20 by binding affinity ──────────────────────────────
top_20 = filtered.sort_values(
    'Affinity_kcal_mol', ascending=True
).drop_duplicates(subset='SMILES').head(20)

print("\n=== TOP 20 LEAD CANDIDATES ===")
print(f"{'SMILES':<55} {'Affinity':>10} {'QED':>6} {'MW':>7} {'Score':>7}")
print("-" * 85)
for _, row in top_20.iterrows():
    smi = row['SMILES'][:52] + "..." if len(row['SMILES']) > 52 else row['SMILES']
    print(f"{smi:<55} {row['Affinity_kcal_mol']:>10.3f} "
          f"{row['QED_value']:>6.3f} {row['MW_value']:>7.1f} "
          f"{row['Score']:>7.4f}")

# ── Save results ────────────────────────────────────────────
top_20.to_csv('top_20_leads_Loa22.csv', index=False)
filtered.to_csv('all_filtered_leads_Loa22.csv', index=False)
print(f"\nSaved: top_20_leads_Loa22.csv")
print(f"Saved: all_filtered_leads_Loa22.csv ({len(filtered)} molecules)")

# ── Draw 2D structures ──────────────────────────────────────
print("Drawing 2D structures...")
mols, legends = [], []

for _, row in top_20.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        mols.append(mol)
        legends.append(
            f"{row['Affinity_kcal_mol']:.2f} kcal/mol\n"
            f"QED: {row['QED_value']:.2f}  "
            f"MW: {row['MW_value']:.0f}"
        )

if mols:
    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=4,
        subImgSize=(450, 350),
        legends=legends
    )
    img.save('Top_20_Loa22_Leads.png')
    print("Saved: Top_20_Loa22_Leads.png")

# ── Score progression plot ──────────────────────────────────
print("Plotting score progression across training steps...")
step_avg = df.groupby('step')['Score'].mean()

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].plot(step_avg.index, step_avg.values,
             color='#2196F3', linewidth=1.5)
axes[0].set_xlabel('Training Step', fontsize=12)
axes[0].set_ylabel('Mean Composite Score', fontsize=12)
axes[0].set_title('REINVENT4 Learning Curve — Loa22', fontweight='bold')
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

axes[1].hist(df['Affinity_kcal_mol'].dropna(), bins=50,
             color='#E91E63', edgecolor='white', alpha=0.85)
axes[1].axvline(x=doxycy_affinity, color='orange',
                linestyle='--', linewidth=2,
                label=f'Doxycycline ({doxycy_affinity})')
axes[1].axvline(x=-7.5, color='red',
                linestyle='--', linewidth=2,
                label='Filter threshold (-7.5)')
axes[1].set_xlabel('Vina Affinity (kcal/mol)', fontsize=12)
axes[1].set_ylabel('Count', fontsize=12)
axes[1].set_title('Binding Affinity Distribution', fontweight='bold')
axes[1].legend()
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('loa22_results_summary.png', dpi=300, bbox_inches='tight')
print("Saved: loa22_results_summary.png")

print("\n" + "=" * 60)
print("DONE. Paste the Top 20 table here for next steps.")
print("=" * 60)