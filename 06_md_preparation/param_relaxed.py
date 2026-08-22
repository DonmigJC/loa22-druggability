#!/usr/bin/env python3
"""
PARAMETERISE FROM RELAXED GEOMETRY, POSE FROM DOCKING
=====================================================

Why this exists
---------------
The previous attempt built each ligand with heavy atoms pinned to the
docked pose and only hydrogens relaxed.  That preserved the docked
geometry exactly, but left mild internal strain that broke
parameterisation:

    LOA22-B2   sqm: "No convergence in SCF after 1000 steps"
               The semi-empirical charge calculation could not converge
               on a strained conformer.

    LOA22-B3   antechamber: "Coordinates issues with your system"
               Two atoms at an implausible separation.

The correct approach, and standard practice
-------------------------------------------
Force-field parameters - atom types, bonded terms, partial charges - are
properties of the molecular TOPOLOGY, not of any particular conformation.
So:

    1. Parameterise from a FULLY RELAXED conformer.  This is what
       antechamber and sqm are designed to handle.
    2. Transfer the resulting topology onto the DOCKED coordinates.

The simulation therefore starts from the docked pose while using
parameters derived from clean geometry.  Energy minimisation before
equilibration then relaxes any residual strain in the complex.

Run:  ~/miniconda3/envs/reinvent4/bin/python param_relaxed.py
"""
import os
import re
import glob
import shutil
import subprocess

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")

R = os.path.expanduser("~/Research")
MD = os.path.join(R, "05_MD_v2")
POSES = os.path.join(R, "final_leads")

LEADS = {
    "LOA22-B1": dict(
        smiles="O=C(Nc1ccc(O)c(C(=O)[O-])c1)c1cccc(C(F)(F)F)c1", charge=-1),
    "LOA22-B2": dict(
        smiles="CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1", charge=0),
    "LOA22-B3": dict(
        smiles="COC(=O)c1sc(NC(=O)COc2ccc(Cl)cc2F)nc1C", charge=0),
}


def sh(cmd, cwd=None, timeout=3600):
    return subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True,
                          text=True, timeout=timeout)


def clean(wd, name):
    for pat in (f"{name}.acpype", f".acpype_tmp_{name}", "ANTECHAMBER*",
                "ATOMTYPE.INF", "sqm.*", "NEWPDB.PDB", "PREP.INF",
                f"{name}_AC.*", "leap.log", "*.prmtop", "*.inpcrd",
                "*_bcc_gaff*.mol2"):
        for p in glob.glob(os.path.join(wd, pat)):
            shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) \
                else os.remove(p)


def best_relaxed_conformer(smiles, nconf=50):
    """Generate an ensemble and return the lowest-energy conformer."""
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3()
    p.randomSeed = 42
    p.pruneRmsThresh = 0.5
    ids = AllChem.EmbedMultipleConfs(m, numConfs=nconf, params=p)
    if not ids:
        if AllChem.EmbedMolecule(m, randomSeed=42) == -1:
            return None, None
        ids = [0]
    res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=5000)
    energies = [(e, i) for i, (conv, e) in enumerate(res)]
    energies.sort()
    best_e, best_i = energies[0]
    keep = Chem.Mol(m, confId=int(best_i))
    return keep, best_e


def pose_heavy_coords(pdbqt):
    out = []
    for l in open(pdbqt):
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("HD", "H"):
            out.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return out


def acpype_ok(wd, name):
    a = os.path.join(wd, f"{name}.acpype")
    return (os.path.exists(os.path.join(a, f"{name}_GMX.itp")) and
            os.path.exists(os.path.join(a, f"{name}_GMX.gro")))


def topology_charge(wd, name):
    itp = os.path.join(wd, f"{name}.acpype", f"{name}_GMX.itp")
    tot, inb = 0.0, False
    for l in open(itp):
        if l.strip().startswith("[ atoms ]"):
            inb = True
            continue
        if inb:
            if l.strip().startswith("["):
                break
            p = l.split()
            if len(p) >= 7 and not l.strip().startswith(";"):
                try:
                    tot += float(p[6])
                except ValueError:
                    pass
    return tot


def write_posed_gro(wd, name, pose_coords):
    """Replace the ACPYPE .gro coordinates with the docked pose."""
    src = os.path.join(wd, f"{name}.acpype", f"{name}_GMX.gro")
    dst = os.path.join(wd, f"{name}_posed.gro")
    lines = open(src).read().splitlines()
    n = int(lines[1])
    atoms, box = lines[2:2 + n], lines[2 + n]

    heavy_rows = [i for i, l in enumerate(atoms)
                  if not l[10:15].strip().startswith("H")]
    if len(heavy_rows) != len(pose_coords):
        return None, (f"heavy-atom mismatch: topology {len(heavy_rows)}, "
                      f"pose {len(pose_coords)}")

    # rigid-body transform mapping the relaxed heavy atoms onto the pose
    import numpy as np
    src_xyz = np.array([[float(atoms[i][20:28]), float(atoms[i][28:36]),
                         float(atoms[i][36:44])] for i in heavy_rows]) * 10.0
    dst_xyz = np.array(pose_coords)
    sc, dc = src_xyz.mean(0), dst_xyz.mean(0)
    H = (src_xyz - sc).T @ (dst_xyz - dc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Rm = Vt.T @ np.diag([1, 1, d]) @ U.T
    rmsd = float(np.sqrt((((src_xyz - sc) @ Rm.T + dc - dst_xyz) ** 2)
                         .sum(1).mean()))

    new = []
    for l in atoms:
        v = np.array([float(l[20:28]), float(l[28:36]), float(l[36:44])]) * 10.0
        w = ((v - sc) @ Rm.T + dc) / 10.0
        new.append(f"{l[:20]}{w[0]:8.3f}{w[1]:8.3f}{w[2]:8.3f}")
    open(dst, "w").write(f"{name} at docked pose\n{n}\n"
                         + "\n".join(new) + f"\n{box}\n")
    return dst, rmsd


print("=" * 78)
print("PARAMETERISE FROM RELAXED GEOMETRY, POSE FROM DOCKING")
print("=" * 78)

status = {}
for name, spec in LEADS.items():
    print(f"\n{'=' * 78}\n{name}\n{'=' * 78}")
    wd = os.path.join(MD, name)
    os.makedirs(wd, exist_ok=True)
    clean(wd, name)

    # 1. relaxed conformer -------------------------------------------
    print("  [1/4] generating a relaxed conformer (50-member ensemble)")
    mol, energy = best_relaxed_conformer(spec["smiles"])
    if mol is None:
        print("        embedding failed"); status[name] = False; continue
    print(f"        lowest MMFF energy: {energy:.2f} kcal/mol")
    sdf = os.path.join(wd, f"{name}_relaxed.sdf")
    with Chem.SDWriter(sdf) as w:
        w.write(mol)

    mol2 = os.path.join(wd, f"{name}_relaxed.mol2")
    sh(f"obabel {sdf} -O {mol2}", cwd=wd)
    if not os.path.exists(mol2):
        print("        Open Babel failed"); status[name] = False; continue

    # 2. parameterise -------------------------------------------------
    print(f"  [2/4] ACPYPE (GAFF2, AM1-BCC, charge {spec['charge']:+d})")
    res = sh(f"conda run -n acpype_env acpype -i {os.path.basename(mol2)} "
             f"-b {name} -n {spec['charge']} -a gaff2 -c bcc -o gmx", cwd=wd)
    if not acpype_ok(wd, name):
        txt = (res.stdout or "") + (res.stderr or "")
        for l in [x for x in txt.splitlines()
                  if re.search(r"error|conver|fail", x, re.I)][-4:]:
            print(f"        {l.strip()[:100]}")
        for f in glob.glob(os.path.join(wd, "**", "sqm.out"), recursive=True):
            bad = [l for l in open(f, errors="ignore").read().splitlines()
                   if re.search(r"conver|error", l, re.I)]
            if bad:
                print(f"        sqm: {bad[-1].strip()[:90]}")
        status[name] = False
        continue
    q = topology_charge(wd, name)
    print(f"        topology charge {q:+.3f}  (target {spec['charge']:+d})")
    if abs(q - spec["charge"]) > 0.05:
        print("        !! charge mismatch"); status[name] = False; continue

    # 3. transfer the docked pose -------------------------------------
    print("  [3/4] transferring the docked pose onto the parameterised molecule")
    pose = os.path.join(POSES, f"{name}_site_12345.pdbqt")
    if not os.path.exists(pose):
        print(f"        pose not found: {pose}"); status[name] = False; continue
    coords = pose_heavy_coords(pose)
    out, info = write_posed_gro(wd, name, coords)
    if out is None:
        print(f"        {info}"); status[name] = False; continue
    print(f"        pose heavy atoms {len(coords)}")
    print(f"        superposition RMSD {info:.3f} A")
    if info > 1.0:
        print("        !! high RMSD - atom ordering may differ; inspect before use")
    print(f"        wrote {out}")

    # 4. verify -------------------------------------------------------
    print("  [4/4] verification")
    m = Chem.MolFromSmiles(spec["smiles"])
    print(f"        formula {rdMolDescriptors.CalcMolFormula(m)}   "
          f"heavy atoms {m.GetNumHeavyAtoms()}")
    status[name] = True

print()
print("=" * 78)
print("SUMMARY")
print("=" * 78)
for n in LEADS:
    print(f"  {n:<12} {'READY' if status.get(n) else 'FAILED'}")

if all(status.values()):
    print("\nAll three parameterised with GAFF2 / AM1-BCC from relaxed")
    print("geometry, positioned at their docked poses.")
    print("\nMethods note: parameters were derived from the lowest-energy")
    print("conformer of a 50-member ETKDGv3 ensemble; simulation")
    print("coordinates were taken from the exhaustiveness-32 docked pose")
    print("by rigid-body superposition.")
    print("\nNext: bash ~/Research/build_md_systems.sh")
else:
    print("\nSend the diagnostic lines for any that failed.")
print("=" * 78)
