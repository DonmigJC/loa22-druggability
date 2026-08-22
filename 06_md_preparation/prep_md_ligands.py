import os
from rdkit import Chem
from rdkit.Chem import AllChem

md_dir = '/home/donmigcage/Research/04_Molecular_Dynamics'
os.makedirs(md_dir, exist_ok=True)

candidates = {
    'Lead_1': 'CC(C)CNC(=O)Nc1ccc2nnc(-c3cc(C(F)(F)F)cc(C(F)(F)F)c3)n2n1',
    'Lead_4': 'CC(CC(=O)N1CCc2ccc(C(F)(F)F)cc2C1)c1ccc(O)cc1',
    'Lead_5': 'CC(=O)Nc1ccc(F)c(C(=O)NS(=O)(=O)c2ccccc2Cl)c1'
}

for name, smiles in candidates.items():
    print(f"Generating 3D coordinates for {name}...")
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        raise ValueError(f"3D embedding failed for {name}")
    AllChem.MMFFOptimizeMolecule(mol)

    ligand_dir = os.path.join(md_dir, name)
    os.makedirs(ligand_dir, exist_ok=True)
    sdf_path = os.path.join(ligand_dir, f"{name}.sdf")

    with Chem.SDWriter(sdf_path) as writer:
        writer.write(mol)

    print(f"  Saved: {sdf_path}")

print(f"\nSuccess! 3D SDF files saved to {md_dir}")
