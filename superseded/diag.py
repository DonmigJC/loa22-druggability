import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import FilterCatalog
RDLogger.DisableLog("rdApp.*")

df = pd.read_csv("loa22_generation_1.csv").dropna(subset=["SMILES"])
v = df[df["SMILES_state"] == 1].copy()
print("columns:", list(df.columns))
print("\nVina_Loa22 (raw) describe:")
print(v["Vina_Loa22 (raw)"].describe().round(4))

v["kcal"] = -12.0 * v["Vina_Loa22 (raw)"]
print("\nconverted (-12 x raw):")
print(v["kcal"].describe().round(4))

v["c"] = [Chem.MolToSmiles(m) if (m := Chem.MolFromSmiles(s)) else None for s in v["SMILES"]]
v = v.dropna(subset=["c"])
print(f"\nvalid {len(v)}   unique {v['c'].nunique()}   dupes {len(v)-v['c'].nunique()}")

for lbl, t in [("old -6.918", -6.918), ("new -6.862", -6.862)]:
    s = v[v["kcal"] <= t]
    print(f"beat {lbl}: {len(s)} ({100*len(s)/len(v):.2f}%), {s['c'].nunique()} unique")

p = FilterCatalog.FilterCatalogParams()
p.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
cat = FilterCatalog.FilterCatalog(p)
s = v.copy(); print(f"\ncascade start   : {len(s)}")
s = s[s["kcal"] <= -7.5]; print(f"  affinity<=-7.5: {len(s)}")
s = s[s["QED (raw)"] >= 0.5]; print(f"  QED>=0.5      : {len(s)}")
s = s[s["MW (raw)"] <= 600]; print(f"  MW<=600       : {len(s)}")
s = s[[not cat.HasMatch(Chem.MolFromSmiles(x)) for x in s["SMILES"]]]
print(f"  PAINS-free    : {len(s)}   <- compare to 766")
print(f"  UNIQUE        : {s['c'].nunique()}")
s.to_csv("all_filtered_leads_Loa22_REGEN.csv", index=False)
print("\ntop 5 by affinity (compare to -9.080, -9.040, -8.984...):")
print(s.sort_values("kcal").drop_duplicates("c").head(5)[["kcal","QED (raw)","MW (raw)"]].to_string())
