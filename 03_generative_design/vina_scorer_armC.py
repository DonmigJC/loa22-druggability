#!/usr/bin/env python3
"""
ARM C SCORING BRIDGE — BRIDGING OBJECTIVE AT THE FUNCTIONAL SITE
===============================================================

The third arm of the generative ablation.

  Arm A   untrained prior, no optimisation                    [done]
  Arm B   RL with naive scoring at P_0                        [done]
  Arm C   RL with corrected scoring at the functional site    [this]

What Arm B did, and why it failed
---------------------------------
Arm B docked into P_0 - which we now know does not contain the
pharmacophore - and scored with a LINEAR affinity transform inside a
geometric mean:

    reward = clip(affinity, -12, 0) / -12

Improving from -6 to -9 kcal/mol raised that term by only 0.25, while QED
moving from 0.4 to 0.9 raised its term by 0.5.  QED is computable
directly from the SMILES string; docking affinity is not.  The agent
therefore optimised the tractable component and abandoned the
high-affinity tail: relative to the untrained prior it produced 2.3x
fewer molecules below -8.5 kcal/mol.

Nothing in Arm B's scoring asked whether a molecule touched Asp121 or
Arg142.  It could not have selected for that.

What Arm C changes
------------------
  1. TARGET SITE.  Docking is performed in the functional-site box
     centred on the CDD-annotated peptidoglycan residues, not P_0.

  2. AFFINITY TRANSFORM.  A reverse sigmoid centred at -6.0 kcal/mol
     replaces the linear ramp, concentrating the reward gradient across
     the range where functional-site affinities actually fall.  A
     molecule at -5.0 now earns almost nothing.

  3. BRIDGING TERM.  A new component measures simultaneous contact with
     the pharmacophore AND with P_0, using the same metric as the
     dual-site library screen so the results are directly comparable.

Weighting:  affinity 0.40 | bridging 0.30 | QED 0.20 | MW 0.10

This script returns the affinity and bridging terms combined, so that
REINVENT's weighted geometric mean reproduces those four weights:

    returned = aff^(0.40/0.70) * bridge^(0.30/0.70)     weight 0.70
    QED                                                  weight 0.20
    MW                                                   weight 0.10
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

R = "/home/donmigcage/Research"
RECEPTOR = os.path.join(R, "loa22_receptor.pdbqt")

# functional-site box, identical to the dual-site library screen
CENTER = (7.381, -2.952, -10.291)
BOX = (27.6, 27.6, 27.6)

EXHAUSTIVENESS = 8
VINA_SEED = 42
# The molecular dynamics runs hold the GPU and most CPU threads.  Keep
# this low while they are running; raise it to 5-6 afterwards.
NUM_WORKERS = int(os.environ.get("ARMC_WORKERS", "3"))
TIMEOUT = 240

# transforms
AFF_MID, AFF_STEEP = -6.0, 0.5          # reverse sigmoid on affinity
CONTACT_MID, CONTACT_STEEP = 4.5, 0.8   # sigmoid on contact distance

W_AFF, W_BRIDGE = 0.40, 0.30
W_TOTAL = W_AFF + W_BRIDGE

KEY = {121: {"CB", "CG", "OD1", "OD2"},
       142: {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}}
P0_RES = {61, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 88, 89, 92,
          93, 96, 97, 115, 116, 117, 119, 163, 165, 169, 171, 172, 183,
          185, 186, 187}


def load_reference():
    key, p0 = [], []
    for l in open(RECEPTOR):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        try:
            rn = int(l[22:26])
        except ValueError:
            continue
        nm = l[12:16].strip()
        if nm.startswith("H"):
            continue
        c = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if rn in KEY and nm in KEY[rn]:
            key.append(c)
        if rn in P0_RES:
            p0.append(c)
    return key, p0


KEY_ATOMS, P0_ATOMS = load_reference()


def affinity_term(a):
    """Reverse sigmoid: -5 -> 0.12, -6 -> 0.50, -7 -> 0.88, -8 -> 0.98"""
    if a is None or a >= 0:
        return 0.0
    try:
        return 1.0 / (1.0 + math.exp((a - AFF_MID) / AFF_STEEP))
    except OverflowError:
        return 0.0


def contact_term(d):
    """Sigmoid on distance: 3.5 A -> 0.78, 4.5 -> 0.50, 6.5 -> 0.08"""
    if d is None:
        return 0.0
    try:
        return 1.0 / (1.0 + math.exp((d - CONTACT_MID) / CONTACT_STEEP))
    except OverflowError:
        return 0.0


def pose_min_distances(path):
    coords = []
    if not os.path.exists(path):
        return None, None
    for l in open(path):
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("HD", "H"):
            try:
                coords.append((float(l[30:38]), float(l[38:46]),
                               float(l[46:54])))
            except ValueError:
                pass
    if not coords:
        return None, None
    dk = min(math.dist(a, b) for a in coords for b in KEY_ATOMS)
    dp = min(math.dist(a, b) for a in coords for b in P0_ATOMS)
    return dk, dp


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


def run_vina(lig, out):
    cmd = ["vina", "--receptor", RECEPTOR, "--ligand", lig,
           "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]),
           "--center_z", str(CENTER[2]),
           "--size_x", str(BOX[0]), "--size_y", str(BOX[1]),
           "--size_z", str(BOX[2]),
           "--exhaustiveness", str(EXHAUSTIVENESS),
           "--num_modes", "1", "--seed", str(VINA_SEED),
           "--cpu", "1", "--out", out]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
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
        a = affinity_term(aff)
        dk, dp = pose_min_distances(out)
        bridge = math.sqrt(contact_term(dk) * contact_term(dp))
        if a <= 0.0 or bridge <= 0.0:
            return 0.0
        return (a ** (W_AFF / W_TOTAL)) * (bridge ** (W_BRIDGE / W_TOTAL))


def main():
    smiles = [l.strip() for l in sys.stdin if l.strip()]
    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as ex:
        scores = list(ex.map(score_smiles, smiles))
    print(json.dumps({"version": 1, "payload": {"bridge_score": scores}}))


if __name__ == "__main__":
    main()
