#!/usr/bin/env python3
"""Dock a set of SMILES against Loa22 P_0 using the same protocol as the
RL scoring bridge (vina_scorer.py): ETKDG default embedding seed 42,
MMFF optimisation, Meeko PDBQT, Vina exhaustiveness 4, num_modes 1.
Usage: dock_prior.py <input.csv> <output.csv>
"""
import sys, os, csv, subprocess, tempfile
from concurrent.futures import ProcessPoolExecutor
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, QED

RDLogger.DisableLog("rdApp.*")

RECEPTOR = "/home/donmigcage/Research/loa22_receptor.pdbqt"
CENTER = (2.30, 7.09, 2.02)
BOX = (25, 25, 25)
EXHAUSTIVENESS = 4
WORKERS = 4
SEED = 42


def prep(smiles, out_path):
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:
            if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) == -1:
                return None
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:
            pass
        setups = MoleculePreparation().prepare(mol)
        if not setups:
            return None
        s, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return None
        open(out_path, "w").write(s)
        return mol
    except Exception:
        return None


def dock(smiles):
    with tempfile.TemporaryDirectory() as td:
        lig = os.path.join(td, "l.pdbqt")
        mol = prep(smiles, lig)
        if mol is None:
            return (smiles, None, None, None)
        cmd = ["vina", "--receptor", RECEPTOR, "--ligand", lig,
               "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]),
               "--center_z", str(CENTER[2]),
               "--size_x", str(BOX[0]), "--size_y", str(BOX[1]),
               "--size_z", str(BOX[2]),
               "--exhaustiveness", str(EXHAUSTIVENESS),
               "--num_modes", "1", "--seed", str(SEED),
               "--out", os.path.join(td, "o.pdbqt")]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        except subprocess.TimeoutExpired:
            return (smiles, None, None, None)
        aff = None
        for line in r.stdout.splitlines():
            p = line.split()
            if p and p[0] == "1":
                try:
                    aff = float(p[1])
                except (IndexError, ValueError):
                    pass
                break
        try:
            plain = Chem.MolFromSmiles(smiles)
            q = QED.qed(plain); mw = Descriptors.MolWt(plain)
        except Exception:
            q = mw = None
        return (smiles, aff, q, mw)


def main():
    inp, outp = sys.argv[1], sys.argv[2]
    smiles = []
    with open(inp) as f:
        rd = csv.DictReader(f)
        col = "SMILES" if "SMILES" in rd.fieldnames else rd.fieldnames[0]
        for row in rd:
            s = (row.get(col) or "").strip()
            if s:
                smiles.append(s)
    print(f"loaded {len(smiles)} SMILES; docking with {WORKERS} workers")

    done = 0
    with open(outp, "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["SMILES", "affinity", "QED", "MW"])
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            for res in ex.map(dock, smiles, chunksize=8):
                w.writerow(res)
                done += 1
                if done % 200 == 0:
                    fo.flush()
                    print(f"  {done}/{len(smiles)}", flush=True)
    print("wrote", outp)


if __name__ == "__main__":
    main()
