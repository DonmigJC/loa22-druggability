#!/usr/bin/env python3
"""
DUAL-SITE RE-DOCKING
====================
Re-docks the full generated library at BOTH sites and computes a
bridging metric identifying molecules that engage the druggable pocket
(P_0) and the functional peptidoglycan site simultaneously.

For every molecule it records:
    aff_p0     affinity in the P_0 box
    aff_site   affinity in the functional-site box
    d_key      closest approach to Asp121/Arg142   (functional-site pose)
    d_p0       closest approach to P_0 residues    (functional-site pose)
    bridge     0-1 score, high when BOTH contacts are made

Usage
    python dual_site_dock.py            # full library, resumable
    python dual_site_dock.py --test 50  # quick 50-molecule check first
"""
import os
import sys
import csv
import math
import argparse
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

R = os.path.expanduser("~/Research")
RECEPTOR = os.path.join(R, "loa22_receptor.pdbqt")
OUTDIR = os.path.join(R, "dual_site")
OUTCSV = os.path.join(OUTDIR, "dual_site_results.csv")

# ---- boxes -----------------------------------------------------------
BOX_P0 = dict(cx=2.30, cy=7.09, cz=2.02, s=25.0)          # as originally used
BOX_SITE = dict(cx=7.381, cy=-2.952, cz=-10.291, s=27.6)  # functional site

EXHAUSTIVENESS = 8
SEED = 42
WORKERS = 5
TIMEOUT = 180

# ---- pharmacophore and pocket definitions ---------------------------
KEY = {121: {"CB", "CG", "OD1", "OD2"},
       142: {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}}
P0_RES = {61, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 88, 89, 92,
          93, 96, 97, 115, 116, 117, 119, 163, 165, 169, 171, 172, 183,
          185, 186, 187}

# bridging transform parameters (angstroms)
KEY_MID, KEY_STEEP = 4.5, 0.8
P0_MID, P0_STEEP = 4.5, 0.8


def load_reference_atoms():
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


KEY_ATOMS, P0_ATOMS = load_reference_atoms()


def sigmoid_contact(d, mid, steep):
    """1.0 at close contact, 0.0 when far."""
    if d is None:
        return 0.0
    try:
        return 1.0 / (1.0 + math.exp((d - mid) / steep))
    except OverflowError:
        return 0.0


def pose_coords(path):
    out = []
    if not os.path.exists(path):
        return out
    for l in open(path):
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("HD", "H"):
            try:
                out.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError:
                pass
    return out


def min_dist(pose, ref):
    if not pose or not ref:
        return None
    return min(math.dist(a, b) for a in pose for b in ref)


def prepare(smiles, path):
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return False
        m = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m, randomSeed=42) == -1:
            if AllChem.EmbedMolecule(m, AllChem.ETKDGv3()) == -1:
                return False
        try:
            AllChem.MMFFOptimizeMolecule(m)
        except Exception:
            pass
        st = MoleculePreparation().prepare(m)
        if not st:
            return False
        s, ok, _ = PDBQTWriterLegacy.write_string(st[0])
        if not ok:
            return False
        open(path, "w").write(s)
        return True
    except Exception:
        return False


def dock(lig, out, box):
    cmd = ["vina", "--receptor", RECEPTOR, "--ligand", lig,
           "--center_x", str(box["cx"]), "--center_y", str(box["cy"]),
           "--center_z", str(box["cz"]),
           "--size_x", str(box["s"]), "--size_y", str(box["s"]),
           "--size_z", str(box["s"]),
           "--exhaustiveness", str(EXHAUSTIVENESS),
           "--num_modes", "1", "--seed", str(SEED),
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


def process(smiles):
    """Dock at both sites; return the full record."""
    with tempfile.TemporaryDirectory() as td:
        lig = os.path.join(td, "l.pdbqt")
        if not prepare(smiles, lig):
            return dict(SMILES=smiles, aff_p0=None, aff_site=None,
                        d_key=None, d_p0=None, bridge=0.0)

        o1 = os.path.join(td, "p0.pdbqt")
        a_p0 = dock(lig, o1, BOX_P0)

        o2 = os.path.join(td, "site.pdbqt")
        a_site = dock(lig, o2, BOX_SITE)

        pose = pose_coords(o2)
        dk = min_dist(pose, KEY_ATOMS)
        dp = min_dist(pose, P0_ATOMS)

        # bridging: geometric mean of both contact terms -> a molecule must
        # touch BOTH to score; touching only one gives near zero
        b = math.sqrt(sigmoid_contact(dk, KEY_MID, KEY_STEEP) *
                      sigmoid_contact(dp, P0_MID, P0_STEEP))

        return dict(SMILES=smiles, aff_p0=a_p0, aff_site=a_site,
                    d_key=dk, d_p0=dp, bridge=b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", type=int, default=0,
                    help="dock only N molecules as a check")
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    print(f"reference atoms: {len(KEY_ATOMS)} pharmacophore, {len(P0_ATOMS)} P_0")
    if len(KEY_ATOMS) < 8:
        sys.exit("!! pharmacophore atoms not found - check residue numbering")

    # --- load library --------------------------------------------------
    import pandas as pd
    df = pd.read_csv(os.path.join(R, "loa22_generation_1.csv"))
    df = df[df.SMILES_state == 1].dropna(subset=["SMILES"])
    smiles = df["SMILES"].tolist()
    if args.test:
        smiles = smiles[:args.test]
    print(f"molecules to dock: {len(smiles)}")

    # --- resume support -------------------------------------------------
    done = set()
    if os.path.exists(OUTCSV) and not args.test:
        with open(OUTCSV) as f:
            for row in csv.DictReader(f):
                done.add(row["SMILES"])
        print(f"already done: {len(done)} (resuming)")
    todo = [s for s in smiles if s not in done]
    print(f"remaining: {len(todo)}")
    if not todo:
        print("nothing to do")
        return

    cols = ["SMILES", "aff_p0", "aff_site", "d_key", "d_p0", "bridge"]
    new = not os.path.exists(OUTCSV) or args.test
    mode = "w" if new else "a"
    path = OUTCSV if not args.test else OUTCSV.replace(".csv", "_test.csv")

    n = 0
    with open(path, mode, newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols)
        if new:
            w.writeheader()
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            for rec in ex.map(process, todo, chunksize=4):
                w.writerow(rec)
                n += 1
                if n % 100 == 0:
                    fo.flush()
                    print(f"  {n}/{len(todo)}", flush=True)

    print(f"\nwrote {path}")

    # --- summary --------------------------------------------------------
    d = pd.read_csv(path).dropna(subset=["aff_site"])
    print(f"\ndocked successfully: {len(d)}")
    print(f"mean affinity  P_0 {d.aff_p0.mean():.3f}   "
          f"functional site {d.aff_site.mean():.3f}")
    print(f"molecules preferring the functional site: "
          f"{(d.aff_site < d.aff_p0).sum()} "
          f"({100*(d.aff_site < d.aff_p0).mean():.1f}%)")
    print(f"\nbridging candidates (bridge > 0.5): {(d.bridge > 0.5).sum()}")
    print(f"                    (bridge > 0.3): {(d.bridge > 0.3).sum()}")

    top = d.sort_values("bridge", ascending=False).head(15)
    print(f"\n{'bridge':>7}{'aff_site':>10}{'aff_p0':>9}"
          f"{'d_key':>8}{'d_p0':>8}  SMILES")
    print("-" * 90)
    for _, r in top.iterrows():
        print(f"{r.bridge:>7.3f}{r.aff_site:>10.3f}{r.aff_p0:>9.3f}"
              f"{r.d_key:>8.2f}{r.d_p0:>8.2f}  {r.SMILES[:44]}")

    ref = "CC(=O)Nc1ccc(F)c(C(=O)NS(=O)(=O)c2ccccc2Cl)c1"   # LOA22-3
    m = d[d.SMILES == ref]
    if len(m):
        r = m.iloc[0]
        rank = (d.bridge > r.bridge).sum() + 1
        print(f"\nLOA22-3 bridge {r.bridge:.3f}, rank {rank} of {len(d)}")


if __name__ == "__main__":
    main()
