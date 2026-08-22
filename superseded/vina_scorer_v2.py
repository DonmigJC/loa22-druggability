#!/usr/bin/env python3
"""
Loa22 scoring bridge v2 — ARM C
================================
Differences from vina_scorer.py (Arm B):

  1. REVERSE-SIGMOID affinity transform instead of linear clip.
     Concentrates the reward gradient between -7.5 and -10 kcal/mol
     so the agent cannot coast on mediocre binding.

  2. PHARMACOPHORE CONSTRAINT.  Parses the docked pose and measures the
     minimum distance from ligand heavy atoms to the Asp121 and Arg142
     side chains.  Contact scores high; binding elsewhere in the box
     scores near zero.

  3. Exhaustiveness 8 (was 4) with a FIXED SEED, reducing reward noise.

The two terms are combined into one value so that REINVENT's weighted
geometric mean reproduces the intended four-component weighting:

    intended : vina 0.45 | pharmacophore 0.25 | QED 0.20 | MW 0.10
    returned : combined = vina_t^(0.45/0.70) * pharm^(0.25/0.70)
    TOML     : combined weight 0.70, QED 0.20, MW 0.10

which expands to exactly vina_t^0.45 * pharm^0.25 * QED^0.20 * MW^0.10.
"""
import sys
import os
import json
import math
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

# ----------------------------------------------------------- configuration
RECEPTOR = "/home/donmigcage/Research/loa22_receptor.pdbqt"
CENTER = (2.30, 7.09, 2.02)        # unchanged from Arm B for comparability
BOX_SIZE = (25, 25, 25)            # unchanged from Arm B for comparability

NUM_WORKERS = 5                    # leave threads for the concurrent MD run
EXHAUSTIVENESS = 8                 # was 4 in Arm B; halves search noise
VINA_SEED = 42                     # was unseeded in Arm B

# affinity transform: logistic, midpoint -8.5, steepness 0.85
AFF_MID, AFF_STEEP = -8.5, 0.85
# pharmacophore transform: logistic on min distance in angstroms
PHARM_MID, PHARM_STEEP = 5.0, 0.8

W_VINA, W_PHARM = 0.45, 0.25
W_TOTAL = W_VINA + W_PHARM

# key residues and the side-chain atoms that define the pharmacophore
KEY_RESIDUES = {
    121: {"CB", "CG", "OD1", "OD2"},                       # Asp121
    142: {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},     # Arg142
}


# ----------------------------------------------------------- receptor setup
def load_key_atoms(path):
    """Extract coordinates of the pharmacophoric side-chain atoms."""
    coords = []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        if resnum not in KEY_RESIDUES:
            continue
        name = line[12:16].strip()
        if name not in KEY_RESIDUES[resnum]:
            continue
        coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return coords


KEY_ATOMS = load_key_atoms(RECEPTOR)
if len(KEY_ATOMS) < 8:
    sys.stderr.write(f"WARNING: only {len(KEY_ATOMS)} key atoms found "
                     f"(expected 11). Check residue numbering.\n")


# ----------------------------------------------------------- transforms
def affinity_transform(aff):
    """Reverse sigmoid. -6 -> ~0.05, -8.5 -> 0.50, -9 -> ~0.64, -10.5 -> ~0.91"""
    if aff is None or aff >= 0:
        return 0.0
    try:
        return 1.0 / (1.0 + math.exp((aff - AFF_MID) / AFF_STEEP))
    except OverflowError:
        return 0.0


def pharmacophore_transform(dist):
    """Logistic on distance. 3.5 A -> ~0.87, 5 A -> 0.50, 8 A -> ~0.02"""
    if dist is None:
        return 0.0
    try:
        return 1.0 / (1.0 + math.exp((dist - PHARM_MID) / PHARM_STEEP))
    except OverflowError:
        return 0.0


def min_distance_to_key(pose_path):
    """Minimum distance (A) from any ligand heavy atom to any key atom."""
    if not KEY_ATOMS or not os.path.exists(pose_path):
        return None
    best = None
    for line in open(pose_path):
        if line.startswith("ENDMDL"):
            break                                   # first pose only
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[77:79].strip() in ("HD", "H"):      # skip hydrogens
            continue
        try:
            lx, ly, lz = (float(line[30:38]), float(line[38:46]),
                          float(line[46:54]))
        except ValueError:
            continue
        for kx, ky, kz in KEY_ATOMS:
            d2 = (lx - kx) ** 2 + (ly - ky) ** 2 + (lz - kz) ** 2
            if best is None or d2 < best:
                best = d2
    return math.sqrt(best) if best is not None else None


# ----------------------------------------------------------- pipeline
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
        setups = MoleculePreparation().prepare(mol)
        if not setups:
            return False
        s, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return False
        open(out_path, "w").write(s)
        return True
    except Exception:
        return False


def run_vina(ligand_pdbqt, out_pdbqt):
    cmd = [
        "vina",
        "--receptor", RECEPTOR, "--ligand", ligand_pdbqt,
        "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]),
        "--center_z", str(CENTER[2]),
        "--size_x", str(BOX_SIZE[0]), "--size_y", str(BOX_SIZE[1]),
        "--size_z", str(BOX_SIZE[2]),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", "1", "--seed", str(VINA_SEED),
        "--out", out_pdbqt,
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    except subprocess.TimeoutExpired:
        return None
    for line in r.stdout.splitlines():
        p = line.split()
        if p and p[0] == "1":
            try:
                return float(p[1])
            except (IndexError, ValueError):
                return None
    return None


def score_smiles(smiles):
    with tempfile.TemporaryDirectory() as td:
        lig = os.path.join(td, "l.pdbqt")
        out = os.path.join(td, "o.pdbqt")
        if not smiles_to_pdbqt(smiles, lig):
            return 0.0
        aff = run_vina(lig, out)
        if aff is None:
            return 0.0
        v = affinity_transform(aff)
        p = pharmacophore_transform(min_distance_to_key(out))
        if v <= 0.0 or p <= 0.0:
            return 0.0
        return (v ** (W_VINA / W_TOTAL)) * (p ** (W_PHARM / W_TOTAL))


def main():
    smiles = [ln.strip() for ln in sys.stdin if ln.strip()]
    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as ex:
        scores = list(ex.map(score_smiles, smiles))
    print(json.dumps({"version": 1, "payload": {"combined_score": scores}}))


if __name__ == "__main__":
    main()
