#!/usr/bin/env python3
"""
POSE TRANSFER BY CONNECTIVITY MATCHING
======================================

Third attempt, and the correct approach.

History of the two failures
---------------------------
  Attempt 1  Rigid-body superposition assuming matching atom order.
             RMSD 4.0-5.8 A - the orders did not match, so chemically
             different atoms were superimposed.

  Attempt 2  Element-order comparison as a guard.  Correctly REFUSED to
             proceed, revealing that Meeko permutes atoms when it writes
             PDBQT.  Same composition, different sequence.

Why this works
--------------
Atom order is an artefact of whichever program wrote the file.  Molecular
identity is not.  This script therefore ignores order completely:

  1. Read the docked pose coordinates and elements from the PDBQT.
  2. Infer bonds from interatomic distances and covalent radii, giving a
     connectivity graph.
  3. Use RDKit's AssignBondOrdersFromTemplate to impose correct bond
     orders from the SMILES.
  4. Graph-match the parameterised molecule onto the pose molecule with
     GetSubstructMatch, which solves the atom correspondence exactly.
  5. Copy coordinates through that mapping.
  6. Verify: heavy-atom RMSD must be ~0.000 A.

Run:  ~/miniconda3/envs/reinvent4/bin/python transfer_poses_v3.py
"""
import os
import sys
import math

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem import rdDetermineBonds
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")

R = os.path.expanduser("~/Research")
MD = os.path.join(R, "05_MD_v2")
POSES = os.path.join(R, "final_leads")

LEADS = {
    "LOA22-B1": dict(
        sim="O=C(Nc1ccc(O)c(C(=O)[O-])c1)c1cccc(C(F)(F)F)c1",
        docked="O=C(Nc1ccc(O)c(C(=O)O)c1)c1cccc(C(F)(F)F)c1", charge=-1),
    "LOA22-B2": dict(
        sim="CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1",
        docked="CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1", charge=0),
    "LOA22-B3": dict(
        sim="COC(=O)c1sc(NC(=O)COc2ccc(Cl)cc2F)nc1C",
        docked="COC(=O)c1sc(NC(=O)COc2ccc(Cl)cc2F)nc1C", charge=0),
}

PDBQT_ELEM = {"A": "C", "OA": "O", "NA": "N", "SA": "S", "HD": "H"}


def read_pose(path):
    """Heavy-atom elements and coordinates from the top-ranked pose."""
    el, xyz = [], []
    for l in open(path):
        if l.startswith("ENDMDL"):
            break
        if not l.startswith(("ATOM", "HETATM")):
            continue
        t = l[77:79].strip()
        e = PDBQT_ELEM.get(t, t)
        if e == "H":
            continue
        e = e.capitalize() if len(e) > 1 else e.upper()
        el.append(e)
        xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return el, xyz


def pose_to_mol(el, xyz, template_smiles):
    """Build an RDKit molecule from coordinates, then impose bond orders."""
    rw = Chem.RWMol()
    conf = Chem.Conformer(len(el))
    for i, (e, c) in enumerate(zip(el, xyz)):
        rw.AddAtom(Chem.Atom(e))
        conf.SetAtomPosition(i, Point3D(*c))
    m = rw.GetMol()
    m.AddConformer(conf, assignId=True)

    # infer connectivity from geometry
    try:
        rdDetermineBonds.DetermineConnectivity(m)
    except Exception as e:
        return None, f"connectivity determination failed: {e}"

    tmpl = Chem.MolFromSmiles(template_smiles)
    if tmpl is None:
        return None, "template SMILES failed to parse"
    try:
        m = AllChem.AssignBondOrdersFromTemplate(tmpl, m)
    except Exception as e:
        return None, f"bond-order assignment failed: {e}"
    return m, None


def gro_read(path):
    L = open(path).read().splitlines()
    n = int(L[1])
    return L[0], n, L[2:2 + n], L[2 + n]


print("=" * 78)
print("POSE TRANSFER BY CONNECTIVITY MATCHING")
print("=" * 78)

ok = {}
for name, spec in LEADS.items():
    print(f"\n{'=' * 78}\n{name}\n{'=' * 78}")
    wd = os.path.join(MD, name)
    src_gro = os.path.join(wd, f"{name}.acpype", f"{name}_GMX.gro")
    mol2 = os.path.join(wd, f"{name}_relaxed.mol2")
    pose = os.path.join(POSES, f"{name}_site_12345.pdbqt")

    missing = [p for p in (src_gro, mol2, pose) if not os.path.exists(p)]
    if missing:
        for p in missing:
            print(f"  missing: {p}")
        ok[name] = False
        continue

    # ---- 1. docked pose as a molecule ------------------------------
    print("  [1/5] reading the docked pose and inferring connectivity")
    el, xyz = read_pose(pose)
    print(f"        {len(el)} heavy atoms")
    pmol, err = pose_to_mol(el, xyz, spec["docked"])
    if pmol is None:
        print(f"        FAILED: {err}")
        ok[name] = False
        continue
    print(f"        bond orders assigned from the SMILES template")

    # ---- 2. parameterised molecule ---------------------------------
    print("  [2/5] reading the parameterised molecule")
    pm = Chem.MolFromMol2File(mol2, removeHs=False, sanitize=True)
    if pm is None:
        pm = Chem.MolFromMol2File(mol2, removeHs=False, sanitize=False)
        if pm is not None:
            try:
                Chem.SanitizeMol(pm)
            except Exception:
                pass
    if pm is None:
        print("        could not read the relaxed MOL2")
        ok[name] = False
        continue
    heavy_all = [a.GetIdx() for a in pm.GetAtoms() if a.GetAtomicNum() > 1]
    print(f"        {pm.GetNumAtoms()} atoms, {len(heavy_all)} heavy")
    if len(heavy_all) != len(el):
        print("        !! heavy-atom count mismatch")
        ok[name] = False
        continue

    # ---- 3. graph match --------------------------------------------
    print("  [3/5] matching atoms by connectivity")
    pm_noH = Chem.RemoveHs(Chem.Mol(pm))
    # map: index in pm_noH -> index in pmol
    match = pm_noH.GetSubstructMatch(pmol, useChirality=False)
    if not match:
        match = pmol.GetSubstructMatch(pm_noH, useChirality=False)
        if match:
            inv = [0] * len(match)
            for i, j in enumerate(match):
                inv[j] = i
            match = tuple(inv)
    if not match or len(match) != len(el):
        print("        !! graph match failed")
        print("        this means the two molecules differ in connectivity")
        ok[name] = False
        continue
    print(f"        matched all {len(match)} heavy atoms")

    # ---- 4. copy coordinates ---------------------------------------
    print("  [4/5] copying docked coordinates")
    # pm_noH index k corresponds to heavy_all[k] in pm
    target = Chem.Mol(pm)
    tconf = target.GetConformer()
    pconf = pmol.GetConformer()
    for k, j in enumerate(match):
        p = pconf.GetAtomPosition(j)
        tconf.SetAtomPosition(heavy_all[k], Point3D(p.x, p.y, p.z))

    try:
        mp = AllChem.MMFFGetMoleculeProperties(target)
        ff = AllChem.MMFFGetMoleculeForceField(target, mp)
        if ff:
            for i in heavy_all:
                ff.AddFixedPoint(i)
            ff.Minimize(maxIts=2000)
            print("        hydrogens relaxed with heavy atoms fixed")
    except Exception as e:
        print(f"        hydrogen relaxation skipped: {e}")

    # ---- 5. verify and write ---------------------------------------
    print("  [5/5] verification")
    dev = 0.0
    for k, j in enumerate(match):
        a = tconf.GetAtomPosition(heavy_all[k])
        b = pconf.GetAtomPosition(j)
        dev += (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
    rmsd = math.sqrt(dev / len(match))
    print(f"        heavy-atom RMSD to the docked pose: {rmsd:.4f} A")
    if rmsd > 0.05:
        print("        !! transfer failed - not writing")
        ok[name] = False
        continue

    title, n, atoms, box = gro_read(src_gro)
    if n != target.GetNumAtoms():
        print(f"        !! .gro has {n} atoms, molecule has "
              f"{target.GetNumAtoms()}")
        ok[name] = False
        continue
    new = []
    for i, l in enumerate(atoms):
        p = tconf.GetAtomPosition(i)
        new.append(f"{l[:20]}{p.x/10:8.3f}{p.y/10:8.3f}{p.z/10:8.3f}")
    dst = os.path.join(wd, f"{name}_posed.gro")
    open(dst, "w").write(f"{name} at docked pose\n{n}\n"
                         + "\n".join(new) + f"\n{box}\n")
    print(f"        wrote {dst}")

    pdb = os.path.join(wd, f"{name}_posed.pdb")
    Chem.MolToPDBFile(target, pdb)
    print(f"        wrote {pdb}  (for visual inspection)")
    ok[name] = True

print()
print("=" * 78)
print("SUMMARY")
print("=" * 78)
for n in LEADS:
    print(f"  {n:<12} {'READY' if ok.get(n) else 'FAILED'}")

if all(ok.get(n) for n in LEADS):
    print("\nAll poses transferred and verified (RMSD < 0.05 A).")
    print("\nBefore building, confirm visually in PyMOL:")
    print("  pymol ~/Research/loa22_alphafold2_raw.pdb \\")
    print("        ~/Research/05_MD_v2/LOA22-B1/LOA22-B1_posed.pdb")
    print("  then:  show sticks, resi 121+142 ; show sticks, not polymer")
    print("\nNext: bash ~/Research/build_md_systems.sh")
else:
    print("\nDo not build systems for any lead marked FAILED.")
    sys.exit(1)
print("=" * 78)
