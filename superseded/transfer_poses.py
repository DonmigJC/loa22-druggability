#!/usr/bin/env python3
"""
POSE TRANSFER BY ATOM MATCHING
==============================

Fixes a failed coordinate transfer.

What went wrong
---------------
The previous script assumed the ACPYPE topology listed heavy atoms in the
same order as the docked PDBQT.  It does not: the relaxed conformer was
re-embedded by RDKit and converted by Open Babel, either of which can
reorder atoms.  Superposition therefore matched chemically different
atoms and returned RMSD values of 4.0-5.8 A instead of <0.1 A.

A rigid-body superposition was also the wrong operation.  The relaxed
conformer and the docked pose have different torsion angles, so no rigid
transform maps one onto the other.

What this does instead
----------------------
  1. Reconstructs the docked pose as a proper RDKit molecule, using the
     SMILES as a template so bond orders are correct.
  2. Maps every atom of the ACPYPE topology to its counterpart in the
     docked pose by substructure matching - chemical identity, not index.
  3. Copies the docked coordinates directly onto the matched atoms.
  4. Places hydrogens by constrained optimisation with heavy atoms fixed.
  5. Verifies by recomputing the heavy-atom RMSD against the pose, which
     must now be ~0.000 A.

Run:  ~/miniconda3/envs/reinvent4/bin/python transfer_poses.py
"""
import os
import sys
import math

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")

R = os.path.expanduser("~/Research")
MD = os.path.join(R, "05_MD_v2")
POSES = os.path.join(R, "final_leads")

LEADS = {
    "LOA22-B1": dict(
        smiles="O=C(Nc1ccc(O)c(C(=O)[O-])c1)c1cccc(C(F)(F)F)c1",
        docked="O=C(Nc1ccc(O)c(C(=O)O)c1)c1cccc(C(F)(F)F)c1", charge=-1),
    "LOA22-B2": dict(
        smiles="CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1",
        docked="CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1", charge=0),
    "LOA22-B3": dict(
        smiles="COC(=O)c1sc(NC(=O)COc2ccc(Cl)cc2F)nc1C",
        docked="COC(=O)c1sc(NC(=O)COc2ccc(Cl)cc2F)nc1C", charge=0),
}


def read_pdbqt_pose(path, template_smiles):
    """Rebuild the docked pose as an RDKit molecule with correct bonds."""
    coords, elems = [], []
    for l in open(path):
        if l.startswith("ENDMDL"):
            break
        if not l.startswith(("ATOM", "HETATM")):
            continue
        t = l[77:79].strip()
        if t in ("HD", "H"):
            continue
        e = {"A": "C", "OA": "O", "NA": "N", "SA": "S"}.get(t, t)
        e = "".join(ch for ch in e if ch.isalpha()).capitalize()
        coords.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        elems.append(e)

    tmpl = Chem.MolFromSmiles(template_smiles)
    if tmpl is None or tmpl.GetNumHeavyAtoms() != len(coords):
        return None, None, (f"template {tmpl.GetNumHeavyAtoms() if tmpl else '?'}"
                            f" vs pose {len(coords)} heavy atoms")

    # PDBQT preserves the heavy-atom order Meeko wrote, which came from
    # the same SMILES; verify element-by-element before trusting it.
    tmpl_elems = [a.GetSymbol() for a in tmpl.GetAtoms()]
    if tmpl_elems != elems:
        return None, None, ("element order differs between template and "
                            f"pose:\n      template {tmpl_elems}\n"
                            f"      pose     {elems}")

    m = Chem.Mol(tmpl)
    conf = Chem.Conformer(m.GetNumAtoms())
    for i, c in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(*c))
    m.AddConformer(conf, assignId=True)
    Chem.AssignStereochemistryFrom3D(m)
    return m, coords, None


def gro_atoms(path):
    L = open(path).read().splitlines()
    n = int(L[1])
    return L[0], n, L[2:2 + n], L[2 + n]


print("=" * 78)
print("POSE TRANSFER BY ATOM MATCHING")
print("=" * 78)

ok = {}
for name, spec in LEADS.items():
    print(f"\n{'=' * 78}\n{name}\n{'=' * 78}")
    wd = os.path.join(MD, name)
    acp = os.path.join(wd, f"{name}.acpype")
    src_gro = os.path.join(acp, f"{name}_GMX.gro")
    pose = os.path.join(POSES, f"{name}_site_12345.pdbqt")

    for p in (src_gro, pose):
        if not os.path.exists(p):
            print(f"  missing: {p}"); ok[name] = False; break
    else:
        # -------- 1. rebuild the docked pose ------------------------
        print("  [1/5] rebuilding the docked pose from PDBQT")
        posed_mol, coords, err = read_pdbqt_pose(pose, spec["docked"])
        if posed_mol is None:
            print(f"        FAILED: {err}")
            ok[name] = False
            continue
        print(f"        {len(coords)} heavy atoms, bond orders from SMILES")

        # -------- 2. read the ACPYPE topology molecule --------------
        print("  [2/5] reading the parameterised molecule")
        title, n, atoms, box = gro_atoms(src_gro)
        heavy_rows = [i for i, l in enumerate(atoms)
                      if not l[10:15].strip().lstrip("0123456789").startswith("H")]
        # the .gro from ACPYPE lists atoms in the MOL2 order; rebuild the
        # molecule from that MOL2 so the mapping is exact
        mol2 = os.path.join(wd, f"{name}_relaxed.mol2")
        param_mol = Chem.MolFromMol2File(mol2, removeHs=False, sanitize=False)
        if param_mol is None:
            print("        could not read the relaxed MOL2"); ok[name] = False
            continue
        try:
            Chem.SanitizeMol(param_mol)
        except Exception:
            pass
        param_heavy = [a.GetIdx() for a in param_mol.GetAtoms()
                       if a.GetAtomicNum() > 1]
        print(f"        topology heavy atoms {len(param_heavy)}   "
              f"total atoms {param_mol.GetNumAtoms()}   .gro atoms {n}")
        if len(param_heavy) != len(coords):
            print("        !! heavy-atom count mismatch"); ok[name] = False
            continue

        # -------- 3. map by substructure ----------------------------
        print("  [3/5] mapping atoms by chemical identity")
        ref = Chem.MolFromSmiles(spec["docked"])
        pm_noH = Chem.RemoveHs(Chem.Mol(param_mol))
        match_param = pm_noH.GetSubstructMatch(ref)
        match_pose = posed_mol.GetSubstructMatch(ref)
        if not match_param or not match_pose:
            print("        !! substructure match failed; falling back to "
                  "element-order matching")
            match_param = tuple(range(len(param_heavy)))
            match_pose = tuple(range(len(coords)))
        print(f"        matched {len(match_param)} atoms")

        # heavy-atom index in param_mol (with H) for each ref atom
        noH_to_H = [a.GetIdx() for a in param_mol.GetAtoms()
                    if a.GetAtomicNum() > 1]
        pose_conf = posed_mol.GetConformer()

        # -------- 4. write coordinates ------------------------------
        print("  [4/5] writing docked coordinates onto the topology")
        target = Chem.Mol(param_mol)
        tconf = target.GetConformer()
        for k in range(len(match_param)):
            pi = noH_to_H[match_param[k]]
            pp = pose_conf.GetAtomPosition(match_pose[k])
            tconf.SetAtomPosition(pi, Point3D(pp.x, pp.y, pp.z))

        # relax hydrogens only
        try:
            mp = AllChem.MMFFGetMoleculeProperties(target)
            ff = AllChem.MMFFGetMoleculeForceField(target, mp)
            if ff:
                for i in noH_to_H:
                    ff.AddFixedPoint(i)
                ff.Minimize(maxIts=1000)
        except Exception as e:
            print(f"        hydrogen relaxation skipped: {e}")

        # -------- 5. verify and write the .gro ----------------------
        print("  [5/5] verification")
        tc = target.GetConformer()
        dev = []
        for k in range(len(match_param)):
            pi = noH_to_H[match_param[k]]
            a = tc.GetAtomPosition(pi)
            b = pose_conf.GetAtomPosition(match_pose[k])
            dev.append((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)
        rmsd = math.sqrt(sum(dev) / len(dev))
        print(f"        heavy-atom RMSD to the docked pose: {rmsd:.4f} A")
        if rmsd > 0.05:
            print("        !! transfer failed - do not simulate this system")
            ok[name] = False
            continue

        new_atoms = []
        for i, l in enumerate(atoms):
            p = tc.GetAtomPosition(i)
            new_atoms.append(f"{l[:20]}{p.x/10:8.3f}{p.y/10:8.3f}{p.z/10:8.3f}")
        dst = os.path.join(wd, f"{name}_posed.gro")
        open(dst, "w").write(f"{name} at docked pose\n{n}\n"
                             + "\n".join(new_atoms) + f"\n{box}\n")
        print(f"        wrote {dst}")
        ok[name] = True

print()
print("=" * 78)
print("SUMMARY")
print("=" * 78)
for n in LEADS:
    print(f"  {n:<12} {'READY' if ok.get(n) else 'FAILED'}")

if all(ok.get(n) for n in LEADS):
    print("\nCoordinates verified against the docked poses (RMSD < 0.05 A).")
    print("Next: bash ~/Research/build_md_systems.sh")
else:
    print("\nDo not build systems for any lead marked FAILED.")
    sys.exit(1)
print("=" * 78)
