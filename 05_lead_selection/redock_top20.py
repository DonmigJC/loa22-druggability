import pandas as pd
import subprocess
import os
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

RECEPTOR = "/home/donmigcage/Research/loa22_receptor.pdbqt"
CENTER = (2.30, 7.09, 2.02)
BOX_SIZE = (25, 25, 25)
EXHAUSTIVENESS = 32   # Proper validation exhaustiveness

def smiles_to_pdbqt(smiles, out_path):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:
        return False
    AllChem.MMFFOptimizeMolecule(mol)
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    if not mol_setups:
        return False
    pdbqt_string, is_ok, _ = PDBQTWriterLegacy.write_string(mol_setups[0])
    if not is_ok:
        return False
    with open(out_path, 'w') as f:
        f.write(pdbqt_string)
    return True

def run_vina(ligand_pdbqt, out_pdbqt):
    cmd = [
        "vina",
        "--receptor", RECEPTOR,
        "--ligand", ligand_pdbqt,
        "--center_x", str(CENTER[0]),
        "--center_y", str(CENTER[1]),
        "--center_z", str(CENTER[2]),
        "--size_x", str(BOX_SIZE[0]),
        "--size_y", str(BOX_SIZE[1]),
        "--size_z", str(BOX_SIZE[2]),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", "9",
        "--seed", "12345",
        "--out", out_pdbqt,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    for line in result.stdout.splitlines():
        parts = line.split()
        if parts and parts[0] == "1":
            try:
                return float(parts[1])
            except:
                return None
    return None

# Load top 20
df = pd.read_csv('top_20_leads_Loa22.csv')
print(f"Re-docking {len(df)} candidates with exhaustiveness={EXHAUSTIVENESS}")
print("This will take approximately 10-20 minutes...\n")

results = []
for i, (_, row) in enumerate(df.iterrows(), 1):
    smi = row['SMILES']
    training_affinity = -12.0 * row['Vina_Loa22 (raw)']
    
    print(f"[{i:2d}/20] Docking... ", end='', flush=True)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        lig_path = os.path.join(tmpdir, 'lig.pdbqt')
        out_path = os.path.join(tmpdir, 'out.pdbqt')
        
        if not smiles_to_pdbqt(smi, lig_path):
            print("FAILED (SMILES preparation)")
            results.append({'SMILES': smi,
                           'Training_Affinity': training_affinity,
                           'Validated_Affinity': None})
            continue
        
        affinity = run_vina(lig_path, out_path)
        
        if affinity:
            delta = affinity - training_affinity
            print(f"Training: {training_affinity:.3f}  "
                  f"Validated: {affinity:.3f}  "
                  f"Delta: {delta:+.3f} kcal/mol")
        else:
            print("FAILED (Vina error)")
        
        results.append({
            'SMILES': smi,
            'QED': row['QED (raw)'],
            'MW': row['MW (raw)'],
            'Training_Affinity': training_affinity,
            'Validated_Affinity': affinity
        })

df_results = pd.DataFrame(results)
df_results = df_results.dropna(subset=['Validated_Affinity'])
df_results = df_results.sort_values('Validated_Affinity')

print("\n=== VALIDATED RANKING (exhaustiveness=32) ===")
print(f"{'#':<3} {'Validated':>10} {'Training':>10} "
      f"{'Delta':>8} {'QED':>6} {'MW':>7}")
print("-" * 55)
for i, (_, r) in enumerate(df_results.iterrows(), 1):
    delta = r['Validated_Affinity'] - r['Training_Affinity']
    print(f"{i:<3} {r['Validated_Affinity']:>10.3f} "
          f"{r['Training_Affinity']:>10.3f} "
          f"{delta:>+8.3f} {r['QED']:>6.3f} {r['MW']:>7.1f}")

df_results.to_csv('top20_validated_docking.csv', index=False)
print("\nSaved: top20_validated_docking.csv")
print("\nPaste the VALIDATED RANKING table — "
      "top 3 from this list go into molecular dynamics.")
