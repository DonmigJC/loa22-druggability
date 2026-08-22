import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, FilterCatalog
RDLogger.DisableLog("rdApp.*")

df = pd.read_csv("loa22_generation_1.csv").dropna(subset=["SMILES"])
v = df[df["SMILES_state"] == 1].copy()
v["A"] = v["Vina_Loa22 (raw)"]; v["Q"] = v["QED (raw)"]; v["M"] = v["MW (raw)"]
v["c"] = [Chem.MolToSmiles(m) if (m := Chem.MolFromSmiles(s)) else None
          for s in v["SMILES"]]
v = v.dropna(subset=["c"])

print(f"total generated   : {len(df)}")
print(f"valid             : {len(v)}")
print(f"UNIQUE molecules  : {v['c'].nunique()}")
print(f"duplicates        : {len(v)-v['c'].nunique()}")

for lbl, t in [("old -6.918", -6.918), ("new -6.862", -6.862)]:
    s = v[v["A"] <= t]
    print(f"beat {lbl}: {len(s)} rows ({100*len(s)/len(v):.2f}%), "
          f"{s['c'].nunique()} unique")

p = FilterCatalog.FilterCatalogParams()
p.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
cat = FilterCatalog.FilterCatalog(p)
s = v.copy(); print(f"\ncascade start     : {len(s)}")
s = s[s["A"] <= -7.5];  print(f"  affinity <=-7.5 : {len(s)}")
s = s[s["Q"] >= 0.5];   print(f"  QED >= 0.5      : {len(s)}")
s = s[s["M"] <= 600];   print(f"  MW <= 600       : {len(s)}")
s = s[[not cat.HasMatch(Chem.MolFromSmiles(x)) for x in s["SMILES"]]]
print(f"  PAINS-free      : {len(s)}   <- the '766'")
print(f"  UNIQUE of those : {s['c'].nunique()}")
s.to_csv("all_filtered_leads_Loa22_REGEN.csv", index=False)

print("\nLead      Vina    HA      LE    logP     pKd     LLE")
for n, (smi, a) in {
  "LOA22-1": ("CC(C)CNC(=O)Nc1ccc2nnc(-c3cc(C(F)(F)F)cc(C(F)(F)F)c3)n2n1", -9.099),
  "LOA22-2": ("CC(CC(=O)N1CCc2ccc(C(F)(F)F)cc2C1)c1ccc(O)cc1", -8.955),
  "LOA22-3": ("CC(=O)Nc1ccc(F)c(C(=O)NS(=O)(=O)c2ccccc2Cl)c1", -8.945),
}.items():
    m = Chem.MolFromSmiles(smi); ha = m.GetNumHeavyAtoms()
    lp = Descriptors.MolLogP(m); pk = abs(a)/1.364
    print(f"{n:<10}{a:7.3f}{ha:5d}{abs(a)/ha:8.3f}{lp:8.3f}{pk:8.3f}{pk-lp:8.3f}")
