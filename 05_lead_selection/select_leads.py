#!/usr/bin/env python3
"""
FINAL LEAD SELECTION AND POSE GENERATION
========================================

Selects bridging leads on balanced multi-criteria grounds, then re-docks
them at exhaustiveness 32 with three seeds, SAVING the poses so they can
seed the molecular dynamics.

This fixes the flaw where MD complexes were built from RDKit-embedded
coordinates rather than docked poses.

Run:  ~/miniconda3/envs/reinvent4/bin/python select_leads.py
"""
import os
import subprocess
import tempfile
import math
import statistics as st
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, QED, rdMolDescriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")

R = os.path.expanduser("~/Research")
RECEPTOR = os.path.join(R, "loa22_receptor.pdbqt")
OUT = os.path.join(R, "final_leads")
os.makedirs(OUT, exist_ok=True)

BOX_SITE = dict(cx=7.381, cy=-2.952, cz=-10.291, s=27.6)
BOX_P0 = dict(cx=2.30, cy=7.09, cz=2.02, s=25.0)

KEY = {121: {"CB", "CG", "OD1", "OD2"},
       142: {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}}
P0_RES = {61,64,65,66,67,68,69,70,71,72,73,74,75,88,89,92,93,96,97,115,116,
          117,119,163,165,169,171,172,183,185,186,187}


def load_ref():
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


KEY_ATOMS, P0_ATOMS = load_ref()

# ------------------------------------------------------------- selection
print("=" * 76)
print("MULTI-CRITERIA LEAD SELECTION")
print("=" * 76)

d = pd.read_csv(os.path.join(R, "dual_site", "bridging_candidates.csv"))
print(f"candidates passing drug-likeness filters: {len(d)}")

# hard requirements
q = d[(d.bridge >= 0.75) &        # genuine dual-site contact
      (d.d_key <= 4.0) &          # within vdW range of the pharmacophore
      (d.aff_site <= -5.5) &      # meaningful binding at the functional site
      (d.QED >= 0.55) &           # drug-like
      (d.MW <= 500) &             # Lipinski
      (d.logP <= 5.0) &           # Lipinski
      (d.LE >= 0.20)].copy()      # efficient for its size
print(f"passing all hard criteria               : {len(q)}")
print("\n  bridge >= 0.75 | d_key <= 4.0 A | aff_site <= -5.5 | "
      "QED >= 0.55 | MW <= 500 | logP <= 5 | LE >= 0.20")

if len(q) < 3:
    print("\n!! relaxing QED to 0.45 and bridge to 0.70")
    q = d[(d.bridge >= 0.70) & (d.d_key <= 4.2) & (d.aff_site <= -5.5) &
          (d.QED >= 0.45) & (d.MW <= 550) & (d.logP <= 5.0)].copy()
    print(f"   now passing: {len(q)}")

# balanced score: all four criteria normalised, geometric mean
for col, invert in (("bridge", False), ("aff_site", True),
                    ("QED", False), ("LE", False)):
    v = -q[col] if invert else q[col]
    q[col + "_n"] = (v - v.min()) / (v.max() - v.min() + 1e-9)
q["balanced"] = (q.bridge_n * q.aff_site_n * q.QED_n * q.LE_n) ** 0.25
q = q.sort_values("balanced", ascending=False)

print(f"\n{'#':>3}{'bal':>7}{'brdg':>7}{'site':>8}{'dkey':>6}{'QED':>7}"
      f"{'MW':>7}{'logP':>6}{'LE':>6}  SMILES")
print("-" * 104)
for i, (_, r) in enumerate(q.head(12).iterrows(), 1):
    print(f"{i:>3}{r.balanced:>7.3f}{r.bridge:>7.3f}{r.aff_site:>8.3f}"
          f"{r.d_key:>6.2f}{r.QED:>7.3f}{r.MW:>7.1f}{r.logP:>6.2f}{r.LE:>6.3f}"
          f"  {r.SMILES[:44]}")

# pick three with distinct scaffolds
picks, scaffolds = [], set()
for _, r in q.iterrows():
    m = Chem.MolFromSmiles(r.SMILES)
    if m is None:
        continue
    from rdkit.Chem.Scaffolds import MurckoScaffold
    sc = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    if sc in scaffolds:
        continue
    scaffolds.add(sc)
    picks.append(r)
    if len(picks) == 3:
        break

print(f"\nselected {len(picks)} leads with distinct Murcko scaffolds")

# ------------------------------------------------------------- docking
def prep(smi, path):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    if AllChem.EmbedMolecule(m, randomSeed=42) == -1:
        return False
    AllChem.MMFFOptimizeMolecule(m)
    s, ok, _ = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(m)[0])
    if not ok:
        return False
    open(path, "w").write(s)
    return True


def dock(lig, out, box, seed):
    cmd = ["vina", "--receptor", RECEPTOR, "--ligand", lig,
           "--center_x", str(box["cx"]), "--center_y", str(box["cy"]),
           "--center_z", str(box["cz"]),
           "--size_x", str(box["s"]), "--size_y", str(box["s"]),
           "--size_z", str(box["s"]),
           "--exhaustiveness", "32", "--num_modes", "9",
           "--seed", str(seed), "--out", out]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    for line in r.stdout.splitlines():
        p = line.split()
        if p and p[0] == "1":
            try:
                return float(p[1])
            except (IndexError, ValueError):
                return None
    return None


def contacts(pose):
    coords = []
    for l in open(pose):
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("HD", "H"):
            coords.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    dk = min(math.dist(a, b) for a in coords for b in KEY_ATOMS)
    dp = min(math.dist(a, b) for a in coords for b in P0_ATOMS)
    return dk, dp


print()
print("=" * 76)
print("RE-DOCKING AT EXHAUSTIVENESS 32 — POSES SAVED FOR MD")
print("=" * 76)

records = []
for i, r in enumerate(picks, 1):
    name = f"LOA22-B{i}"
    print(f"\n{name}")
    print(f"  SMILES {r.SMILES}")
    lig = os.path.join(OUT, f"{name}_input.pdbqt")
    if not prep(r.SMILES, lig):
        print("  preparation failed")
        continue
    vals, poses = [], []
    for sd in (12345, 67890, 24680):
        po = os.path.join(OUT, f"{name}_site_{sd}.pdbqt")
        v = dock(lig, po, BOX_SITE, sd)
        if v is not None:
            vals.append(v)
            poses.append(po)
    if not vals:
        print("  docking failed")
        continue
    dk, dp = contacts(poses[0])
    m = Chem.MolFromSmiles(r.SMILES)
    print(f"  functional site : {st.mean(vals):.3f} +/- "
          f"{st.pstdev(vals) if len(vals)>1 else 0:.3f} kcal/mol")
    print(f"  d(Asp121/Arg142): {dk:.2f} A     d(P_0): {dp:.2f} A")
    print(f"  QED {QED.qed(m):.3f}  MW {Descriptors.MolWt(m):.1f}  "
          f"logP {Descriptors.MolLogP(m):.2f}  "
          f"TPSA {Descriptors.TPSA(m):.1f}  "
          f"HBD/HBA {rdMolDescriptors.CalcNumHBD(m)}/{rdMolDescriptors.CalcNumHBA(m)}")
    print(f"  best pose -> {poses[0]}")

    # also dock at P_0 for the comparison table
    pp = os.path.join(OUT, f"{name}_p0_12345.pdbqt")
    vp = dock(lig, pp, BOX_P0, 12345)
    print(f"  P_0 (reference) : {vp:.3f} kcal/mol" if vp else "  P_0 failed")

    records.append(dict(name=name, SMILES=r.SMILES,
                        aff_site=st.mean(vals),
                        aff_site_sd=st.pstdev(vals) if len(vals) > 1 else 0,
                        aff_p0=vp, d_key=dk, d_p0=dp,
                        bridge=r.bridge, QED=QED.qed(m),
                        MW=Descriptors.MolWt(m),
                        logP=Descriptors.MolLogP(m),
                        TPSA=Descriptors.TPSA(m),
                        HA=m.GetNumHeavyAtoms(),
                        LE=st.mean(vals) / -m.GetNumHeavyAtoms(),
                        pose=poses[0]))

if records:
    out_csv = os.path.join(OUT, "final_leads.csv")
    pd.DataFrame(records).to_csv(out_csv, index=False)
    print(f"\nwritten: {out_csv}")

print()
print("=" * 76)
print("COMPARISON WITH THE ORIGINAL LEADS")
print("=" * 76)
print(f"{'compound':<12}{'d_key':>8}{'aff_site':>10}{'bridge':>9}")
print("-" * 40)
for r in records:
    print(f"{r['name']:<12}{r['d_key']:>8.2f}{r['aff_site']:>10.3f}{r['bridge']:>9.3f}")
print("-" * 40)
print(f"{'LOA22-1':<12}{6.71:>8.2f}{-6.436:>10.3f}{0.028:>9.3f}")
print(f"{'LOA22-2':<12}{7.35:>8.2f}{-6.416:>10.3f}{0.012:>9.3f}")
print(f"{'LOA22-3':<12}{7.41:>8.2f}{-5.705:>10.3f}{0.016:>9.3f}")
print("\nvan der Waals contact is 3.5-4.0 A; beyond ~5 A there is no")
print("direct interaction with the pharmacophoric residues.")
print("=" * 76)
