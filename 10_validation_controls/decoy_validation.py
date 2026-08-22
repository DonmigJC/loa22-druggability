#!/usr/bin/env python3
"""
DECOY VALIDATION AND DISCRIMINATION TESTING
===========================================

Tests whether the docking protocol and the bridging metric discriminate,
or whether they simply score everything favourably.

Three analyses:

  A. PROPERTY CORRELATION (no docking required)
     Is the bridging score measuring molecular geometry, or is it just
     rewarding large, floppy molecules that can reach further?  If bridge
     score correlates strongly with molecular weight or rotatable-bond
     count, the metric measures size rather than complementarity.

  B. PROPERTY-MATCHED DECOYS
     For each lead, selects molecules matched on molecular weight, logP,
     hydrogen-bond donors and acceptors, and rotatable bonds, but
     topologically dissimilar (Tanimoto < 0.3).  Docks them at both sites
     at exhaustiveness 32, identically to the leads.  If the leads are
     indistinguishable from their property-matched decoys, the selection
     carries no information.

  C. RANDOM BASELINE
     A random sample from the library, for context.

A note on what this can and cannot establish: with no experimentally
confirmed binders, a classical enrichment or ROC analysis is impossible.
Selecting "actives" by the same docking score being validated would be
circular.  These analyses avoid that by testing whether the metric is
confounded by molecular properties, and whether property-matched
molecules behave differently.

Run:  ~/miniconda3/envs/reinvent4/bin/python decoy_validation.py
"""
import os
import math
import subprocess
import tempfile
import statistics as st
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")

R = os.path.expanduser("~/Research")
RECEPTOR = os.path.join(R, "loa22_receptor.pdbqt")
OUT = os.path.join(R, "decoy_validation")
os.makedirs(OUT, exist_ok=True)

BOX_SITE = dict(cx=7.381, cy=-2.952, cz=-10.291, s=27.6)
BOX_P0 = dict(cx=2.30, cy=7.09, cz=2.02, s=25.0)

KEY = {121: {"CB", "CG", "OD1", "OD2"},
       142: {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}}
P0_RES = {61,64,65,66,67,68,69,70,71,72,73,74,75,88,89,92,93,96,97,115,116,
          117,119,163,165,169,171,172,183,185,186,187}

LEADS = {
    "LOA22-B1": "O=C(Nc1ccc(O)c(C(=O)O)c1)c1cccc(C(F)(F)F)c1",
    "LOA22-B2": "CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1",
    "LOA22-B3": "COC(=O)c1sc(NC(=O)COc2ccc(Cl)cc2F)nc1C",
}

N_DECOYS_PER_LEAD = 15
N_RANDOM = 15
SEED = 12345


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
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def props(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return dict(mol=m,
                MW=Descriptors.MolWt(m),
                logP=Descriptors.MolLogP(m),
                HBD=rdMolDescriptors.CalcNumHBD(m),
                HBA=rdMolDescriptors.CalcNumHBA(m),
                rotb=rdMolDescriptors.CalcNumRotatableBonds(m),
                TPSA=Descriptors.TPSA(m),
                HA=m.GetNumHeavyAtoms(),
                fp=FPGEN.GetFingerprint(m))


# ==================================================================
print("=" * 78)
print("A. IS THE BRIDGING METRIC MEASURING GEOMETRY OR SIZE?")
print("=" * 78)
print("""If bridge score is driven by molecular weight or flexibility, it is
rewarding molecules that can physically reach further rather than
molecules that are geometrically complementary to both sites.
""")

d = pd.read_csv(os.path.join(R, "dual_site", "dual_site_results.csv"))
d = d.dropna(subset=["aff_site", "aff_p0", "d_key", "d_p0", "bridge"])
print(f"library with complete data: {len(d)}")

print("\ncomputing descriptors ...")
rows = []
for smi in d.SMILES:
    p = props(smi)
    rows.append((p["MW"], p["logP"], p["HBD"], p["HBA"], p["rotb"],
                 p["TPSA"], p["HA"]) if p else (np.nan,) * 7)
desc = pd.DataFrame(rows, columns=["MW", "logP", "HBD", "HBA",
                                   "rotb", "TPSA", "HA"], index=d.index)
d = pd.concat([d, desc], axis=1).dropna(subset=["MW"])

print(f"\n{'property':<10}{'r with bridge':>16}{'r with aff_site':>18}"
      f"{'r with d_key':>15}")
print("-" * 60)
for c in ["MW", "HA", "logP", "HBD", "HBA", "rotb", "TPSA"]:
    r1 = d[c].corr(d.bridge)
    r2 = d[c].corr(d.aff_site)
    r3 = d[c].corr(d.d_key)
    print(f"{c:<10}{r1:>16.3f}{r2:>18.3f}{r3:>15.3f}")

print("""
Reading this: r is the Pearson correlation coefficient, from -1 to +1.
Values near 0 mean the property does not explain the score.  |r| > 0.5
would indicate the metric is substantially confounded by that property.
Note that affinity is expected to correlate with size - that is a known
property of empirical scoring functions - so the informative column is
the bridging score.""")

# ==================================================================
print()
print("=" * 78)
print("B. SELECTING PROPERTY-MATCHED DECOYS")
print("=" * 78)
print("""Decoys are matched on molecular weight, logP, hydrogen-bond counts and
flexibility, but constrained to be topologically dissimilar (Morgan
fingerprint Tanimoto < 0.3).  If a lead cannot be distinguished from
molecules with the same physicochemical profile, its selection carries
no information beyond those properties.
""")

lead_props = {n: props(s) for n, s in LEADS.items()}
selected = []
used = set(LEADS.values())

for name, lp in lead_props.items():
    print(f"\n{name}  MW {lp['MW']:.1f}  logP {lp['logP']:.2f}  "
          f"HBD/HBA {lp['HBD']}/{lp['HBA']}  rotb {lp['rotb']}")
    cand = d[(d.MW.between(lp["MW"] - 30, lp["MW"] + 30)) &
             (d.logP.between(lp["logP"] - 0.7, lp["logP"] + 0.7)) &
             (d.HBD.between(lp["HBD"] - 1, lp["HBD"] + 1)) &
             (d.HBA.between(lp["HBA"] - 2, lp["HBA"] + 2)) &
             (d.rotb.between(lp["rotb"] - 2, lp["rotb"] + 2))]
    cand = cand[~cand.SMILES.isin(used)]
    print(f"  property-matched in library: {len(cand)}")

    picked = 0
    for _, r in cand.sample(frac=1.0, random_state=SEED).iterrows():
        p = props(r.SMILES)
        if p is None:
            continue
        tan = DataStructs.TanimotoSimilarity(lp["fp"], p["fp"])
        if tan >= 0.30:
            continue
        selected.append(dict(smiles=r.SMILES, matched_to=name,
                             tanimoto=tan, kind="decoy"))
        used.add(r.SMILES)
        picked += 1
        if picked >= N_DECOYS_PER_LEAD:
            break
    print(f"  dissimilar decoys selected : {picked}")

rnd = d[~d.SMILES.isin(used)].sample(n=min(N_RANDOM, len(d)),
                                     random_state=SEED)
for _, r in rnd.iterrows():
    selected.append(dict(smiles=r.SMILES, matched_to="—",
                         tanimoto=np.nan, kind="random"))
print(f"\nrandom baseline molecules: {len(rnd)}")
print(f"total to dock: {len(selected)}")

# ==================================================================
def prepare(smiles, path):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False
    mh = Chem.AddHs(m)
    if AllChem.EmbedMolecule(mh, randomSeed=42) == -1:
        return False
    try:
        AllChem.MMFFOptimizeMolecule(mh)
    except Exception:
        pass
    s = MoleculePreparation().prepare(mh)
    if not s:
        return False
    txt, ok, _ = PDBQTWriterLegacy.write_string(s[0])
    if not ok:
        return False
    open(path, "w").write(txt)
    return True


def dock(lig, out, box):
    cmd = ["vina", "--receptor", RECEPTOR, "--ligand", lig,
           "--center_x", str(box["cx"]), "--center_y", str(box["cy"]),
           "--center_z", str(box["cz"]),
           "--size_x", str(box["s"]), "--size_y", str(box["s"]),
           "--size_z", str(box["s"]),
           "--exhaustiveness", "32", "--num_modes", "9",
           "--seed", str(SEED), "--out", out]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
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


def mind(pose, ref):
    c = []
    for l in open(pose):
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("HD", "H"):
            c.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    if not c:
        return None
    return min(math.dist(a, b) for a in c for b in ref)


def run_one(rec):
    with tempfile.TemporaryDirectory() as td:
        lig = os.path.join(td, "l.pdbqt")
        if not prepare(rec["smiles"], lig):
            return None
        ps = os.path.join(td, "s.pdbqt")
        a_s = dock(lig, ps, BOX_SITE)
        pp = os.path.join(td, "p.pdbqt")
        a_p = dock(lig, pp, BOX_P0)
        if a_s is None or a_p is None:
            return None
        dk = mind(ps, KEY_ATOMS)
        dp = mind(ps, P0_ATOMS)
        out = dict(rec)
        out.update(aff_site=a_s, aff_p0=a_p, d_key=dk, d_p0=dp)
        return out


print()
print("=" * 78)
print("C. DOCKING AT EXHAUSTIVENESS 32 (identical protocol to the leads)")
print("=" * 78)
print(f"docking {len(selected)} molecules at two sites ...")

done = []
with ProcessPoolExecutor(max_workers=5) as ex:
    for i, res in enumerate(ex.map(run_one, selected, chunksize=1), 1):
        if res:
            done.append(res)
        if i % 10 == 0:
            print(f"  {i}/{len(selected)}", flush=True)

res = pd.DataFrame(done)
res.to_csv(os.path.join(OUT, "decoy_results.csv"), index=False)
print(f"\nsuccessfully docked: {len(res)}")

# ==================================================================
print()
print("=" * 78)
print("D. LEADS versus DECOYS")
print("=" * 78)

LEAD_DATA = [
    dict(name="LOA22-B1", aff_site=-6.335, aff_p0=-7.465, d_key=3.43, d_p0=3.04),
    dict(name="LOA22-B2", aff_site=-5.970, aff_p0=-7.706, d_key=3.33, d_p0=3.13),
    dict(name="LOA22-B3", aff_site=-5.625, aff_p0=-7.044, d_key=3.12, d_p0=2.88),
]
ld = pd.DataFrame(LEAD_DATA)

dec = res[res.kind == "decoy"]
ran = res[res.kind == "random"]

def summarise(label, frame):
    if not len(frame):
        print(f"{label:<24} (none)")
        return
    print(f"{label:<24}{frame.aff_site.mean():>11.3f}"
          f"{frame.aff_p0.mean():>11.3f}{frame.d_key.mean():>10.2f}"
          f"{frame.d_key.min():>10.2f}")

print(f"{'group':<24}{'aff_site':>11}{'aff_P_0':>11}{'d_key':>10}{'best d_key':>10}")
print("-" * 66)
summarise("LEADS (n=3)", ld)
summarise(f"property-matched decoys (n={len(dec)})", dec)
summarise(f"random library (n={len(ran)})", ran)

if len(dec):
    print()
    print("Discrimination:")
    n_closer = (dec.d_key < ld.d_key.max()).sum()
    print(f"  decoys reaching within {ld.d_key.max():.2f} A of the "
          f"pharmacophore : {n_closer}/{len(dec)} "
          f"({100*n_closer/len(dec):.1f}%)")
    n_better = (dec.aff_site < ld.aff_site.max()).sum()
    print(f"  decoys binding better than the weakest lead        : "
          f"{n_better}/{len(dec)} ({100*n_better/len(dec):.1f}%)")
    both = ((dec.d_key < 4.5) & (dec.aff_site < -5.6)).sum()
    print(f"  decoys satisfying BOTH criteria                    : "
          f"{both}/{len(dec)} ({100*both/len(dec):.1f}%)")

    try:
        from scipy.stats import mannwhitneyu
        u, p = mannwhitneyu(ld.d_key, dec.d_key, alternative="less")
        print(f"\n  Mann-Whitney (lead d_key < decoy d_key): "
              f"U={u:.0f}, p={p:.4f}")
        print("  (n=3 leads, so this is indicative rather than conclusive)")
    except ImportError:
        pass

print()
print("=" * 78)
print("INTERPRETATION")
print("=" * 78)
print("""
The bridging criterion is validated if property-matched decoys rarely
reach the pharmacophore.  A low percentage means the leads were selected
for geometry that similar molecules do not share.

A high percentage would mean the criterion is satisfied by any molecule
of that size and polarity, and the leads carry no specific information.

The random baseline shows what an unselected molecule achieves, giving
the scale against which both other groups should be read.
""")
print(f"written: {os.path.join(OUT, 'decoy_results.csv')}")
print("=" * 78)
