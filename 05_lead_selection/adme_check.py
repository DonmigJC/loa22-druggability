import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, QED

candidates = {
    'Lead_1 (Best Affinity)': 'CC(C)CNC(=O)Nc1ccc2nnc(-c3cc(C(F)(F)F)cc(C(F)(F)F)c3)n2n1',
    'Lead_4 (Highest QED)': 'CC(CC(=O)N1CCc2ccc(C(F)(F)F)cc2C1)c1ccc(O)cc1',
    'Lead_5 (High QED/Clean)': 'CC(=O)Nc1ccc(F)c(C(=O)NS(=O)(=O)c2ccccc2Cl)c1'
}

print(f"{'Name':<25} {'MW':>7} {'logP':>7} {'HBD':>5} {'HBA':>5} {'TPSA':>7} {'RotB':>5} {'QED':>6}")
print('-' * 80)

for name, smi in candidates.items():
    mol = Chem.MolFromSmiles(smi)
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = Descriptors.TPSA(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    qed  = QED.qed(mol)
    print(f"{name:<25} {mw:>7.1f} {logp:>7.2f} {hbd:>5} {hba:>5} {tpsa:>7.1f} {rotb:>5} {qed:>6.3f}")
