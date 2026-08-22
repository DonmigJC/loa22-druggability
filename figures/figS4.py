#!/usr/bin/env python3
"""
Figure S4 — two-dimensional structures of the three selected compounds,
with Bemis-Murcko scaffolds shown separately.

Fix relative to the first version: RDKit's MolsToGridImage now returns a
PIL PngImageFile rather than an object carrying a .data attribute. This
version handles both return types.

Also prints the InChIKey, molecular formula and scaffold for each compound,
which are required for Supplementary Table S5.

Run from ~/Research:
    ~/miniconda3/envs/reinvent4/bin/python figS4.py
"""
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

SMILES = {
    "LOA22-B1": "O=C(Nc1ccc(O)c(C(=O)O)c1)c1cccc(C(F)(F)F)c1",
    "LOA22-B2": "CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1",
    "LOA22-B3": "COC(=O)c1sc(NC(=O)COc2ccc(Cl)cc2F)nc1C",
}

CLASS = {
    "LOA22-B1": "benzamide-salicylate",
    "LOA22-B2": "pyridine-oxadiazole",
    "LOA22-B3": "aminothiophene-carboxylate",
}


def save(img, path):
    """Handle both RDKit return types across versions."""
    if hasattr(img, "data"):                 # older: IPython image object
        with open(path, "wb") as fh:
            fh.write(img.data)
    elif hasattr(img, "save"):               # newer: PIL image
        img.save(path)
    else:
        with open(path, "wb") as fh:
            fh.write(img)
    print(f"  wrote {path}")


# ---------------------------------------------------------------- molecules
mols, names, scaffolds = [], [], []
for name, smi in SMILES.items():
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise SystemExit(f"Failed to parse SMILES for {name}")
    AllChem.Compute2DCoords(m)
    mols.append(m)
    names.append(f"{name}\n{CLASS[name]}")
    s = MurckoScaffold.GetScaffoldForMol(m)
    AllChem.Compute2DCoords(s)
    scaffolds.append(s)

# ---------------------------------------------------------------- panel A
img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(420, 380),
                           legends=names)
save(img, "FigureS4A_compounds.png")

# ---------------------------------------------------------------- panel B
img = Draw.MolsToGridImage(scaffolds, molsPerRow=3, subImgSize=(420, 300),
                           legends=[n.split("\n")[0] + " scaffold"
                                    for n in names])
save(img, "FigureS4B_scaffolds.png")

# ---------------------------------------------------------------- table data
print()
print("Supplementary Table S5 values:")
print(f"{'ID':<12}{'Formula':<16}{'MW':>8}  {'InChIKey':<30}")
print("-" * 70)
for name, smi in SMILES.items():
    m = Chem.MolFromSmiles(smi)
    print(f"{name:<12}{rdMolDescriptors.CalcMolFormula(m):<16}"
          f"{Descriptors.MolWt(m):>8.2f}  {Chem.MolToInchiKey(m):<30}")

print()
print("Murcko scaffolds (SMILES):")
for name, smi in SMILES.items():
    m = Chem.MolFromSmiles(smi)
    sc = MurckoScaffold.GetScaffoldForMol(m)
    print(f"  {name}: {Chem.MolToSmiles(sc)}")
