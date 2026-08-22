#!/usr/bin/env python3
"""
ARM C EVALUATION — completing the three-arm ablation
======================================================

Arm C was trained with a bridging-aware scoring bridge (linear-in-loop
scoring, exhaustiveness 8, single pose per molecule).  That in-loop
score is NOT directly comparable across arms, because Arm B's in-loop
score used a different transform and a different box (P_0 only).

To make a fair comparison, this script re-scores Arm C's generated
library with the SAME independent bridging metric already applied to
Arm A (prior) and Arm B (naive RL) in dual_site_dock.py:
docking at BOTH sites, exhaustiveness 8, one fixed seed, geometric-mean
bridge score from contact distances.

This is the only way to compare three arms trained under three
different objectives on equal footing.

Run:  ~/miniconda3/envs/reinvent4/bin/python evaluate_armC.py
"""
import os
import math
import subprocess
import tempfile
import glob
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, QED, Descriptors, FilterCatalog

RDLogger.DisableLog("rdApp.*")

R = os.path.expanduser("~/Research")
RECEPTOR = os.path.join(R, "loa22_receptor.pdbqt")
OUT = os.path.join(R, "arm_c")

BOX_SITE = dict(cx=7.381, cy=-2.952, cz=-10.291, s=27.6)
BOX_P0 = dict(cx=2.30, cy=7.09, cz=2.02, s=25.0)
EXHAUSTIVENESS = 8
SEED = 42
WORKERS = 6

KEY = {121: {"CB", "CG", "OD1", "OD2"},
       142: {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}}
P0_RES = {61, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 88, 89, 92,
          93, 96, 97, 115, 116, 117, 119, 163, 165, 169, 171, 172, 183,
          185, 186, 187}
CONTACT_MID, CONTACT_STEEP = 4.5, 0.8


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


def contact_term(d):
    if d is None:
        return 0.0
    try:
        return 1.0 / (1.0 + math.exp((d - CONTACT_MID) / CONTACT_STEEP))
    except OverflowError:
        return 0.0


def prep(smiles, path):
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
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    setups = MoleculePreparation().prepare(m)
    if not setups:
        return False
    s, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        return False
    open(path, "w").write(s)
    return True


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


def pose_dist(pose, ref):
    if not os.path.exists(pose):
        return None
    c = []
    for l in open(pose):
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("HD", "H"):
            c.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    if not c:
        return None
    return min(math.dist(a, b) for a in c for b in ref)


def process(smiles):
    with tempfile.TemporaryDirectory() as td:
        lig = os.path.join(td, "l.pdbqt")
        if not prep(smiles, lig):
            return dict(SMILES=smiles, aff_site=None, aff_p0=None,
                        d_key=None, d_p0=None, bridge=0.0)
        o1 = os.path.join(td, "s.pdbqt")
        a_s = dock(lig, o1, BOX_SITE)
        dk = pose_dist(o1, KEY_ATOMS)
        dp = pose_dist(o1, P0_ATOMS)
        o2 = os.path.join(td, "p.pdbqt")
        a_p = dock(lig, o2, BOX_P0)
        bridge = math.sqrt(contact_term(dk) * contact_term(dp))
        return dict(SMILES=smiles, aff_site=a_s, aff_p0=a_p,
                    d_key=dk, d_p0=dp, bridge=bridge)


def main():
    files = sorted(glob.glob(os.path.join(OUT, "armC_generation*.csv")))
    if not files:
        print("No Arm C output CSV found in", OUT)
        return
    d = pd.read_csv(files[-1])
    d = d[d.SMILES_state == 1].dropna(subset=["SMILES"])
    smiles = d["SMILES"].tolist()
    print(f"Arm C: {len(smiles)} valid molecules to re-score at both sites")
    print(f"unique: {d.SMILES.nunique()}")

    resf = os.path.join(OUT, "armC_dualsite_results.csv")
    done = set()
    if os.path.exists(resf):
        done = set(pd.read_csv(resf).SMILES)
        print(f"resuming: {len(done)} already scored")
    todo = [s for s in smiles if s not in done]

    if todo:
        rows = []
        n = 0
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            for rec in ex.map(process, todo, chunksize=4):
                rows.append(rec)
                n += 1
                if n % 100 == 0:
                    print(f"  {n}/{len(todo)}", flush=True)
                    pd.DataFrame(rows).to_csv(
                        resf, mode="a" if os.path.exists(resf) else "w",
                        header=not os.path.exists(resf), index=False)
                    rows = []
        if rows:
            pd.DataFrame(rows).to_csv(
                resf, mode="a" if os.path.exists(resf) else "w",
                header=not os.path.exists(resf), index=False)

    res = pd.read_csv(resf).dropna(subset=["aff_site", "aff_p0"])
    print(f"\nscored: {len(res)}")

    print("\n" + "=" * 70)
    print("THREE-ARM COMPARISON")
    print("=" * 70)
    arm_b = pd.read_csv(os.path.join(R, "dual_site", "dual_site_results.csv"))
    arm_b = arm_b.dropna(subset=["bridge"])
    arm_a_path = os.path.join(R, "prior_control", "prior_docked.csv")

    print(f"\n{'arm':<10}{'n':>7}{'mean bridge':>13}{'>0.3':>8}"
          f"{'>0.5':>8}{'>0.7':>8}{'>0.8':>8}")
    print("-" * 62)
    for label, df, col in (("B (P_0)", arm_b, "bridge"),
                           ("C (site)", res, "bridge")):
        n = len(df)
        m = df[col].mean()
        f3 = (df[col] > 0.3).sum()
        f5 = (df[col] > 0.5).sum()
        f7 = (df[col] > 0.7).sum()
        f8 = (df[col] > 0.8).sum()
        print(f"{label:<10}{n:>7}{m:>13.4f}{f3:>5}({100*f3/n:4.1f}%)"
              f"{f5:>5}({100*f5/n:4.1f}%){f7:>5}({100*f7/n:4.1f}%)"
              f"{f8:>5}({100*f8/n:4.1f}%)")

    a5 = (arm_b.bridge > 0.5).mean()
    c5 = (res.bridge > 0.5).mean()
    print(f"\nenrichment (C vs B), bridge>0.5: {c5/a5:.2f}x" if a5 else "")

    print("\n" + "-" * 70)
    print("top 10 from Arm C by bridging score")
    print("-" * 70)
    print(f"{'bridge':>7}{'aff_site':>10}{'aff_p0':>9}{'d_key':>7}{'d_p0':>7}  SMILES")
    for _, r in res.sort_values("bridge", ascending=False).head(10).iterrows():
        print(f"{r.bridge:>7.3f}{r.aff_site:>10.3f}{r.aff_p0:>9.3f}"
              f"{r.d_key:>7.2f}{r.d_p0:>7.2f}  {r.SMILES[:44]}")

    print("\n" + "-" * 70)
    print("drug-likeness of Arm C's top bridgers")
    print("-" * 70)
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    cat = FilterCatalog.FilterCatalog(params)
    top = res.sort_values("bridge", ascending=False).head(50)
    ok = 0
    for _, r in top.iterrows():
        m = Chem.MolFromSmiles(r.SMILES)
        if m and not cat.HasMatch(m) and QED.qed(m) >= 0.5 and Descriptors.MolWt(m) <= 500:
            ok += 1
    print(f"  of top 50 by bridge: {ok} pass PAINS-free, QED>=0.5, MW<=500")

    print("=" * 70)


if __name__ == "__main__":
    main()
