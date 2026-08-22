#!/usr/bin/env python3
import sys
import os
import json
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

from rdkit import Chem
from rdkit.Chem import AllChem

RECEPTOR = "/home/donmigcage/Research/loa22_receptor.pdbqt"
CENTER = (2.30, 7.09, 2.02)
BOX_SIZE = (25, 25, 25)

NUM_WORKERS = 2
EXHAUSTIVENESS = 4


def smiles_to_pdbqt(smiles, out_path):
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        mol = Chem.AddHs(mol)

        if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:
            if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) == -1:
                return False

        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:
            pass

        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)
        if not mol_setups:
            return False

        pdbqt_string, is_ok, _err = PDBQTWriterLegacy.write_string(mol_setups[0])
        if not is_ok:
            return False

        with open(out_path, "w") as f:
            f.write(pdbqt_string)
        return True
    except Exception:
        return False


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
        "--num_modes", "1",
        "--out", out_pdbqt,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return None

    for line in result.stdout.splitlines():
        parts = line.split()
        if parts and parts[0] == "1":
            try:
                return float(parts[1])
            except (IndexError, ValueError):
                return None
    return None


def affinity_to_score(affinity):
    if affinity is None:
        return 0.0
    clipped = max(min(affinity, 0.0), -12.0)
    return clipped / -12.0


def score_smiles(smiles):
    with tempfile.TemporaryDirectory() as tmpdir:
        ligand_pdbqt = os.path.join(tmpdir, "ligand.pdbqt")
        out_pdbqt = os.path.join(tmpdir, "out.pdbqt")

        if not smiles_to_pdbqt(smiles, ligand_pdbqt):
            return 0.0

        affinity = run_vina(ligand_pdbqt, out_pdbqt)
        return affinity_to_score(affinity)


def main():
    smilies = [line.strip() for line in sys.stdin if line.strip()]

    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        scores = list(executor.map(score_smiles, smilies))

    output = {"version": 1, "payload": {"vina_score": scores}}
    print(json.dumps(output))


if __name__ == "__main__":
    main()
