#!/usr/bin/env python3
"""
MD LIGAND PREPARATION v2 — FROM DOCKED POSES
============================================

Prepares the three bridging leads for molecular dynamics, correcting
every defect identified in the original preparation.

WHAT WENT WRONG BEFORE, AND WHAT THIS FIXES

  1. Complexes were built by CONCATENATING a protein .gro with an
     independently embedded ligand .gro (build_complexes.py).  The ligand
     position came from RDKit's coordinate origin, not from docking.
     -> This script extracts coordinates from the exhaustiveness-32
        DOCKED POSE and preserves them exactly.

  2. All ligands were modelled NEUTRAL, including one whose acyl
     sulfonamide should have carried a negative charge at pH 7.4.
     -> Protonation is assigned explicitly per compound, with the
        reasoning recorded, and the net charge is passed to ACPYPE.

  3. Stereochemistry was undefined (isomeric_smiles = false during
     generation), so RDKit chose enantiomers arbitrarily.
     -> Stereocentres are detected and reported.  None of the three
        current leads has one; if that changes, it is flagged loudly.

  4. One ligand atom was manually displaced 0.3 nm to escape a steric
     clash, and this was never disclosed.
     -> Clashes are detected and reported BEFORE simulation.  A clashing
        pose is rejected rather than edited.

Run:  ~/miniconda3/envs/reinvent4/bin/python prep_md_leads.py
"""
import os
import sys
import math
import shutil
import subprocess

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

R = os.path.expanduser("~/Research")
POSES = os.path.join(R, "final_leads")
OUT = os.path.join(R, "05_MD_v2")
os.makedirs(OUT, exist_ok=True)


def find_receptor():
    """Locate the AlphaFold PDB; it may live on the Windows side."""
    candidates = [
        os.path.join(R, "loa22_alphafold2_raw.pdb"),
        "/mnt/c/Users/donmi/Desktop/Research/loa22_alphafold2_raw.pdb",
        os.path.join(R, "04_Molecular_Dynamics", "loa22_alphafold2_raw.pdb"),
        os.path.join(R, "04_Molecular_Dynamics", "loa22_clean.pdb"),
        os.path.join(R, "loa22_clean.pdb"),
    ]
    for c in candidates:
        if os.path.exists(c):
            # keep a local copy so later steps have a stable path
            local = os.path.join(R, "loa22_alphafold2_raw.pdb")
            if c != local:
                try:
                    shutil.copy(c, local)
                    print(f"copied receptor from {c}")
                except Exception:
                    return c
            return local
    print("!! receptor PDB not found. Searched:")
    for c in candidates:
        print(f"   {c}")
    sys.exit(1)


RECEPTOR_PDB = find_receptor()

# ---------------------------------------------------------------------
# Ligand definitions.
#
# `neutral_smiles` is the form as generated and docked.
# `sim_smiles`     is the form to SIMULATE, at pH 7.4.
# `charge`         is the resulting net formal charge.
#
# Protonation reasoning is recorded for the manuscript.
# ---------------------------------------------------------------------
LEADS = {
    "LOA22-B1": dict(
        neutral_smiles="O=C(Nc1ccc(O)c(C(=O)O)c1)c1cccc(C(F)(F)F)c1",
        sim_smiles="O=C(Nc1ccc(O)c(C(=O)[O-])c1)c1cccc(C(F)(F)F)c1",
        charge=-1,
        reason="Benzoic acid moiety, pKa ~4.2. Fully deprotonated at "
               "pH 7.4. The phenol (pKa ~10) remains protonated.",
    ),
    "LOA22-B2": dict(
        neutral_smiles="CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1",
        sim_smiles="CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1",
        charge=0,
        reason="Two amides (non-ionisable), a pyridine (conjugate acid "
               "pKa ~5, therefore neutral at pH 7.4) and a 1,3,4-"
               "oxadiazole (non-basic). Neutral.",
    ),
    "LOA22-B3": dict(
        neutral_smiles="COC(=O)c1sc(NC(=O)COc2ccc(Cl)cc2F)nc1C",
        sim_smiles="COC(=O)c1sc(NC(=O)COc2ccc(Cl)cc2F)nc1C",
        charge=0,
        reason="Methyl ester, amide and aryl ether; the 2-aminothiazole "
               "nitrogen is weakly basic (pKa ~5). Neutral.",
    ),
}


def sh(cmd, cwd=None, check=True):
    r = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    if check and r.returncode != 0:
        print(f"    command failed: {cmd}")
        print(f"    {r.stderr[-600:]}")
    return r


def pose_coords(pdbqt):
    """Heavy-atom coordinates of the top-ranked pose, in input order."""
    coords, names = [], []
    for l in open(pdbqt):
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")):
            t = l[77:79].strip()
            if t in ("HD", "H"):
                continue
            coords.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            names.append(l[12:16].strip())
    return coords, names


def receptor_heavy():
    out = []
    for l in open(RECEPTOR_PDB):
        if l.startswith("ATOM") and not l[12:16].strip().startswith("H"):
            out.append((l[17:20].strip(), int(l[22:26]), l[12:16].strip(),
                        (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return out


REC = receptor_heavy()
print("=" * 78)
print("MD LIGAND PREPARATION v2 — FROM DOCKED POSES")
print("=" * 78)
print(f"receptor heavy atoms: {len(REC)}")
print(f"pose directory      : {POSES}")
print(f"output directory    : {OUT}\n")

ok_all = True
summary = []

for name, spec in LEADS.items():
    print("=" * 78)
    print(name)
    print("=" * 78)
    wd = os.path.join(OUT, name)
    os.makedirs(wd, exist_ok=True)

    # ---------------- 1. chemistry checks --------------------------
    m_neu = Chem.MolFromSmiles(spec["neutral_smiles"])
    m_sim = Chem.MolFromSmiles(spec["sim_smiles"])
    if m_neu is None or m_sim is None:
        print("  !! SMILES failed to parse"); ok_all = False; continue

    q = Chem.GetFormalCharge(m_sim)
    print(f"  docked form   : {spec['neutral_smiles']}")
    print(f"  simulated form: {spec['sim_smiles']}")
    print(f"  protonation   : {spec['reason']}")
    print(f"  net charge    : {q}  (expected {spec['charge']})")
    if q != spec["charge"]:
        print("  !! charge mismatch — check the SMILES"); ok_all = False; continue

    stereo = Chem.FindMolChiralCenters(m_sim, includeUnassigned=True,
                                       useLegacyImplementation=False)
    if stereo:
        print(f"  !! STEREOCENTRES PRESENT: {stereo}")
        print("     Undefined stereochemistry must be resolved before "
              "simulation — see the LOA22-2 precedent.")
        ok_all = False
    else:
        print(f"  stereocentres : none")

    print(f"  formula {rdMolDescriptors.CalcMolFormula(m_sim)}   "
          f"MW {Descriptors.MolWt(m_sim):.2f}   "
          f"heavy atoms {m_sim.GetNumHeavyAtoms()}")

    # ---------------- 2. recover the docked pose -------------------
    pose = os.path.join(POSES, f"{name}_site_12345.pdbqt")
    if not os.path.exists(pose):
        print(f"  !! docked pose not found: {pose}"); ok_all = False; continue
    coords, _ = pose_coords(pose)
    n_heavy = m_neu.GetNumHeavyAtoms()
    print(f"\n  docked pose   : {os.path.basename(pose)}")
    print(f"  pose heavy atoms {len(coords)}   molecule heavy atoms {n_heavy}")
    if len(coords) != n_heavy:
        print("  !! atom-count mismatch — cannot map the pose reliably")
        ok_all = False; continue

    # ---------------- 3. clash detection ---------------------------
    worst, worst_desc = 99.0, ""
    for c in coords:
        for rs, rn, an, rc in REC:
            d = math.dist(c, rc)
            if d < worst:
                worst, worst_desc = d, f"{rs}{rn}:{an}"
    print(f"  closest receptor contact: {worst:.2f} A  ({worst_desc})")
    if worst < 2.0:
        print("  !! STERIC CLASH — this pose must not be simulated.")
        print("     Do NOT displace atoms manually (the LOA22-2 precedent).")
        print("     Use a different pose or a different seed instead.")
        ok_all = False; continue
    elif worst < 2.5:
        print("     tight but tolerable; energy minimisation should resolve it")
    else:
        print("     clash-free")

    # ---------------- 4. build the 3D molecule at the pose ---------
    mh = Chem.AddHs(m_sim)
    if AllChem.EmbedMolecule(mh, randomSeed=42) == -1:
        print("  !! embedding failed"); ok_all = False; continue
    AllChem.MMFFOptimizeMolecule(mh, maxIters=2000)

    # transfer docked coordinates onto the heavy atoms, in order
    conf = mh.GetConformer()
    heavy_idx = [a.GetIdx() for a in mh.GetAtoms() if a.GetAtomicNum() > 1]
    if len(heavy_idx) != len(coords):
        print(f"  !! heavy-atom count changed on protonation "
              f"({len(heavy_idx)} vs {len(coords)})")
        ok_all = False; continue
    from rdkit.Geometry import Point3D
    for idx, c in zip(heavy_idx, coords):
        conf.SetAtomPosition(idx, Point3D(*c))
    # re-optimise hydrogens only, keeping heavy atoms fixed
    ff = AllChem.MMFFGetMoleculeForceField(
        mh, AllChem.MMFFGetMoleculeProperties(mh))
    if ff:
        for idx in heavy_idx:
            ff.AddFixedPoint(idx)
        ff.Minimize(maxIts=500)

    sdf = os.path.join(wd, f"{name}.sdf")
    with Chem.SDWriter(sdf) as w:
        w.write(mh)
    print(f"\n  wrote {sdf}")

    # ---------------- 5. parameterise with ACPYPE ------------------
    mol2 = os.path.join(wd, f"{name}.mol2")
    r = sh(f"obabel {sdf} -O {mol2} 2>&1", check=False)
    if not os.path.exists(mol2):
        print("  !! Open Babel conversion failed"); ok_all = False; continue

    print(f"  running ACPYPE (GAFF2, AM1-BCC, net charge {q}) ...")
    r = sh(f"conda run -n acpype_env acpype -i {os.path.basename(mol2)} "
           f"-b {name} -n {q} -a gaff2 -c bcc -o gmx", cwd=wd, check=False)
    acp = os.path.join(wd, f"{name}.acpype")
    itp = os.path.join(acp, f"{name}_GMX.itp")
    gro = os.path.join(acp, f"{name}_GMX.gro")
    if not (os.path.exists(itp) and os.path.exists(gro)):
        print("  !! ACPYPE did not produce GROMACS files")
        print(f"     {r.stdout[-500:] if r.stdout else ''}")
        ok_all = False; continue

    # verify the topology charge sums to the intended value
    tot, inblock = 0.0, False
    for l in open(itp):
        if l.strip().startswith("[ atoms ]"):
            inblock = True; continue
        if inblock:
            if l.strip().startswith("["):
                break
            p = l.split()
            if len(p) >= 7 and not l.strip().startswith(";"):
                try:
                    tot += float(p[6])
                except ValueError:
                    pass
    print(f"  topology total charge: {tot:+.3f}  (target {q:+d})")
    if abs(tot - q) > 0.05:
        print("  !! charge mismatch in the generated topology")
        ok_all = False; continue

    print(f"  parameterisation complete: {itp}")
    summary.append(dict(name=name, charge=q, itp=itp, gro=gro,
                        sdf=sdf, clash=worst))
    print()

# =====================================================================
print("=" * 78)
print("SUMMARY")
print("=" * 78)
print(f"{'lead':<12}{'charge':>8}{'closest contact':>18}  status")
print("-" * 56)
for s in summary:
    print(f"{s['name']:<12}{s['charge']:>+8d}{s['clash']:>15.2f} A  ready")
for name in LEADS:
    if not any(s["name"] == name for s in summary):
        print(f"{name:<12}{'—':>8}{'—':>18}  FAILED")

if ok_all and len(summary) == len(LEADS):
    print("\nAll ligands prepared from docked poses.")
    print("Next: bash ~/Research/build_md_systems.sh")
else:
    print("\n!! Some ligands failed. Resolve before proceeding —")
    print("   do not edit coordinates by hand to force a pose through.")
    sys.exit(1)
