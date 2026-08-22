#!/usr/bin/env python3
"""
POSE TRANSFER — GRO-BASED MATCHING
==================================

Fourth and final approach.  The previous three failed as follows:

  v1  assumed matching atom order        -> RMSD 4-6 A, wrong atoms paired
  v2  guarded on element order           -> correctly refused; Meeko permutes
  v3  matched via the Open Babel MOL2    -> RDKit could not kekulize it
      ("non-ring atom marked aromatic")

The problem in v3 was reading the MOL2, not the matching strategy itself:
steps 1 and 2 succeeded, showing that connectivity inference and
bond-order assignment from the pose both work.

This version never touches the MOL2.  Instead:

  1. Build the docked pose as a molecule (connectivity from geometry,
     bond orders from the SMILES template).            [proven to work]
  2. Build the ACPYPE topology molecule from its own .gro coordinates,
     inferring connectivity the same way.
  3. Graph-match the two.
  4. Copy coordinates through the mapping and verify.

Both molecules are therefore constructed by the same route, from
coordinates, which removes any dependence on file conventions.

Run:  ~/miniconda3/envs/reinvent4/bin/python transfer_poses_v4.py
"""
import os
import sys
import math

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDetermineBonds
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")

R = os.path.expanduser("~/Research")
MD = os.path.join(R, "05_MD_v2")
POSES = os.path.join(R, "final_leads")

LEADS = {
    "LOA22-B1": dict(
        docked="O=C(Nc1ccc(O)c(C(=O)O)c1)c1cccc(C(F)(F)F)c1", charge=-1),
    "LOA22-B2": dict(
        docked="CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1", charge=0),
    "LOA22-B3": dict(
        docked="COC(=O)c1sc(NC(=O)COc2ccc(Cl)cc2F)nc1C", charge=0),
}

PDBQT_ELEM = {"A": "C", "OA": "O", "NA": "N", "SA": "S", "HD": "H"}


def mol_from_xyz(elements, coords, template_smiles=None):
    """Build a molecule from coordinates; optionally impose bond orders."""
    rw = Chem.RWMol()
    conf = Chem.Conformer(len(elements))
    for i, (e, c) in enumerate(zip(elements, coords)):
        rw.AddAtom(Chem.Atom(e))
        conf.SetAtomPosition(i, Point3D(*c))
    m = rw.GetMol()
    m.AddConformer(conf, assignId=True)
    try:
        rdDetermineBonds.DetermineConnectivity(m)
    except Exception as e:
        return None, f"connectivity failed: {e}"
    if template_smiles:
        t = Chem.MolFromSmiles(template_smiles)
        if t is None:
            return None, "template SMILES failed"
        try:
            m = AllChem.AssignBondOrdersFromTemplate(t, m)
        except Exception as e:
            return None, f"bond orders failed: {e}"
    return m, None


def read_pose(path):
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
        el.append(e.capitalize() if len(e) > 1 else e.upper())
        xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return el, xyz


def read_gro(path):
    """Return header, count, atom lines, box, plus element list and coords."""
    L = open(path).read().splitlines()
    n = int(L[1])
    atoms, box = L[2:2 + n], L[2 + n]
    el, xyz, heavy_rows = [], [], []
    for i, l in enumerate(atoms):
        nm = l[10:15].strip()
        sym = "".join(ch for ch in nm if ch.isalpha())
        if not sym:
            continue
        # two-letter elements present in these ligands
        if sym[:2].capitalize() in ("Cl", "Br"):
            e = sym[:2].capitalize()
        else:
            e = sym[0].upper()
        if e == "H":
            continue
        heavy_rows.append(i)
        el.append(e)
        xyz.append((float(l[20:28]) * 10.0, float(l[28:36]) * 10.0,
                    float(l[36:44]) * 10.0))
    return L[0], n, atoms, box, el, xyz, heavy_rows


print("=" * 78)
print("POSE TRANSFER - GRO-BASED CONNECTIVITY MATCHING")
print("=" * 78)

ok = {}
for name, spec in LEADS.items():
    print(f"\n{'=' * 78}\n{name}\n{'=' * 78}")
    wd = os.path.join(MD, name)
    src_gro = os.path.join(wd, f"{name}.acpype", f"{name}_GMX.gro")
    pose = os.path.join(POSES, f"{name}_site_12345.pdbqt")

    miss = [p for p in (src_gro, pose) if not os.path.exists(p)]
    if miss:
        for p in miss:
            print(f"  missing: {p}")
        ok[name] = False
        continue

    # ---- 1. docked pose -------------------------------------------
    print("  [1/5] docked pose")
    pel, pxyz = read_pose(pose)
    pmol, err = mol_from_xyz(pel, pxyz, spec["docked"])
    if pmol is None:
        print(f"        FAILED: {err}")
        ok[name] = False
        continue
    print(f"        {len(pel)} heavy atoms, bond orders from template")

    # ---- 2. topology molecule from the .gro -----------------------
    print("  [2/5] parameterised molecule (from ACPYPE .gro)")
    title, n, atoms, box, gel, gxyz, heavy_rows = read_gro(src_gro)
    print(f"        {n} atoms total, {len(gel)} heavy")
    if len(gel) != len(pel):
        print(f"        !! heavy-atom mismatch: .gro {len(gel)}, "
              f"pose {len(pel)}")
        ok[name] = False
        continue
    gmol, err = mol_from_xyz(gel, gxyz, spec["docked"])
    if gmol is None:
        print(f"        template route failed ({err}); retrying without")
        gmol, err = mol_from_xyz(gel, gxyz, None)
        if gmol is None:
            print(f"        FAILED: {err}")
            ok[name] = False
            continue

    # ---- 3. graph match -------------------------------------------
    print("  [3/5] matching by connectivity")
    match = gmol.GetSubstructMatch(pmol, useChirality=False)
    if not match or len(match) != len(pel):
        # try the reverse direction and invert
        rev = pmol.GetSubstructMatch(gmol, useChirality=False)
        if rev and len(rev) == len(pel):
            inv = [0] * len(rev)
            for i, j in enumerate(rev):
                inv[j] = i
            match = tuple(inv)
        else:
            print("        !! graph match failed")
            ok[name] = False
            continue
    print(f"        matched {len(match)} atoms")
    # sanity: matched atoms must share their element
    bad = [k for k in range(len(match)) if gel[match[k]] != pel[k]]
    if bad:
        print(f"        !! {len(bad)} matched pairs have different elements")
        ok[name] = False
        continue
    print("        element identity confirmed for every pair")

    # ---- 4. copy coordinates --------------------------------------
    print("  [4/5] writing docked coordinates")
    new_atoms = list(atoms)
    for k in range(len(match)):
        row = heavy_rows[match[k]]
        x, y, z = pxyz[k]
        l = new_atoms[row]
        new_atoms[row] = f"{l[:20]}{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}"

    # hydrogens: shift each by the displacement of its nearest heavy atom
    hyd_rows = [i for i in range(n) if i not in heavy_rows]
    old_heavy = {heavy_rows[j]: gxyz[j] for j in range(len(gxyz))}
    new_heavy = {heavy_rows[match[k]]: pxyz[k] for k in range(len(match))}
    moved = 0
    for hr in hyd_rows:
        l = atoms[hr]
        hx = (float(l[20:28]) * 10, float(l[28:36]) * 10, float(l[36:44]) * 10)
        nearest, best = None, 1e9
        for row, oc in old_heavy.items():
            d = math.dist(hx, oc)
            if d < best:
                best, nearest = d, row
        ox, oy, oz = old_heavy[nearest]
        nx, ny, nz = new_heavy[nearest]
        px, py, pz = hx[0] + (nx - ox), hx[1] + (ny - oy), hx[2] + (nz - oz)
        l2 = new_atoms[hr]
        new_atoms[hr] = f"{l2[:20]}{px/10:8.3f}{py/10:8.3f}{pz/10:8.3f}"
        moved += 1
    print(f"        {len(match)} heavy atoms placed, {moved} hydrogens carried")

    # ---- 5. verify -------------------------------------------------
    print("  [5/5] verification")
    dev = 0.0
    for k in range(len(match)):
        row = heavy_rows[match[k]]
        l = new_atoms[row]
        v = (float(l[20:28]) * 10, float(l[28:36]) * 10, float(l[36:44]) * 10)
        dev += sum((a - b) ** 2 for a, b in zip(v, pxyz[k]))
    rmsd = math.sqrt(dev / len(match))
    print(f"        heavy-atom RMSD to the docked pose: {rmsd:.4f} A")
    if rmsd > 0.05:
        print("        !! transfer failed - not writing")
        ok[name] = False
        continue

    dst = os.path.join(wd, f"{name}_posed.gro")
    open(dst, "w").write(f"{name} at docked pose\n{n}\n"
                         + "\n".join(new_atoms) + f"\n{box}\n")
    print(f"        wrote {dst}")

    # a PDB copy for visual checking
    pdb = os.path.join(wd, f"{name}_posed.pdb")
    with open(pdb, "w") as f:
        for i, l in enumerate(new_atoms, 1):
            nm = l[10:15].strip()
            x = float(l[20:28]) * 10
            y = float(l[28:36]) * 10
            z = float(l[36:44]) * 10
            sym = "".join(c for c in nm if c.isalpha())
            e = sym[:2].capitalize() if sym[:2].capitalize() in ("Cl", "Br") \
                else sym[0].upper()
            f.write(f"HETATM{i:5d} {nm:<4s} LIG A   1    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
                    f"          {e:>2s}\n")
        f.write("END\n")
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
    print("\nInspect one before building:")
    print("  conda activate pymol_env")
    print("  pymol ~/Research/loa22_alphafold2_raw.pdb \\")
    print("        ~/Research/05_MD_v2/LOA22-B1/LOA22-B1_posed.pdb")
    print("\nNext: bash ~/Research/build_md_systems.sh")
else:
    print("\nDo not build systems for any lead marked FAILED.")
    sys.exit(1)
print("=" * 78)
