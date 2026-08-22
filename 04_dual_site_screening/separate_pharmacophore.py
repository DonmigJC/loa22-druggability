#!/usr/bin/env python3
"""
SEPARATED PHARMACOPHORE ANALYSIS
================================

Why this exists
---------------
Every bridging figure in this project so far - the 22.5% for Arm B, the
42.6% for Arm C, the 1.91x enrichment, the decoy discrimination, the
lead selection criteria - was computed against a COMBINED Asp121+Arg142
atom group.  Both gmx mindist and the docking scorer return the shortest
distance to EITHER residue.

Because Arg142 sits consistently at ~0.25 nm and Asp121 at ~0.72 nm,
every one of those numbers is in practice measuring Arg142 alone while
being described as measuring both.

The comparisons remain valid (all arms were measured identically), but
the description does not match the measurement.  This script recomputes
everything with the two residues separated.

Three questions
---------------
  1. How far is each library molecule from Asp121, and from Arg142,
     measured independently?
  2. Does ANY molecule in the 6,253-compound library reach both?
     If none does, that is itself a finding: it would mean the two
     residues are too far apart for a drug-sized molecule to span,
     reinforcing the protein-polymer interface argument.
  3. What are the corrected Arm B vs Arm C figures?

Method
------
Re-docks each molecule at the functional site and measures distances to
Asp121 side-chain atoms and Arg142 side-chain atoms separately, storing
both.  Uses the same box, exhaustiveness and seed as the original
dual-site screen so the numbers are directly comparable.

Run:  python separate_pharmacophore.py           # both arms, full
      python separate_pharmacophore.py --test 50 # quick check first
"""
import os
import sys
import csv
import math
import argparse
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

R = os.path.expanduser("~/Research")
RECEPTOR = os.path.join(R, "loa22_receptor.pdbqt")
OUTDIR = os.path.join(R, "pharmacophore_split")
os.makedirs(OUTDIR, exist_ok=True)

BOX = dict(cx=7.381, cy=-2.952, cz=-10.291, s=27.6)
EXHAUSTIVENESS = 8
SEED = 42
WORKERS = 6
TIMEOUT = 180

# side-chain atoms only - backbone is common to all residues and would
# not distinguish a specific interaction
ASP121 = {"CB", "CG", "OD1", "OD2"}
ARG142 = {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}
P0_RES = {61, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 88, 89, 92,
          93, 96, 97, 115, 116, 117, 119, 163, 165, 169, 171, 172, 183,
          185, 186, 187}

CONTACT = 4.5     # angstroms


def load_reference():
    asp, arg, p0 = [], [], []
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
        if rn == 121 and nm in ASP121:
            asp.append(c)
        if rn == 142 and nm in ARG142:
            arg.append(c)
        if rn in P0_RES:
            p0.append(c)
    return asp, arg, p0


ASP_ATOMS, ARG_ATOMS, P0_ATOMS = load_reference()


def prep(smiles, path):
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
        s = MoleculePreparation().prepare(m)
        if not s:
            return False
        txt, ok, _ = PDBQTWriterLegacy.write_string(s[0])
        if not ok:
            return False
        open(path, "w").write(txt)
        return True
    except Exception:
        return False


def dock(lig, out):
    cmd = ["vina", "--receptor", RECEPTOR, "--ligand", lig,
           "--center_x", str(BOX["cx"]), "--center_y", str(BOX["cy"]),
           "--center_z", str(BOX["cz"]),
           "--size_x", str(BOX["s"]), "--size_y", str(BOX["s"]),
           "--size_z", str(BOX["s"]),
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


def pose_coords(path):
    c = []
    if not os.path.exists(path):
        return c
    for l in open(path):
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("HD", "H"):
            try:
                c.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError:
                pass
    return c


def mind(pose, ref):
    if not pose or not ref:
        return None
    return min(math.dist(a, b) for a in pose for b in ref)


def process(smiles):
    with tempfile.TemporaryDirectory() as td:
        lig = os.path.join(td, "l.pdbqt")
        out = os.path.join(td, "o.pdbqt")
        if not prep(smiles, lig):
            return dict(SMILES=smiles, aff=None, d_asp=None,
                        d_arg=None, d_p0=None)
        aff = dock(lig, out)
        if aff is None:
            return dict(SMILES=smiles, aff=None, d_asp=None,
                        d_arg=None, d_p0=None)
        pose = pose_coords(out)
        return dict(SMILES=smiles, aff=aff,
                    d_asp=mind(pose, ASP_ATOMS),
                    d_arg=mind(pose, ARG_ATOMS),
                    d_p0=mind(pose, P0_ATOMS))


def run_arm(name, smiles, outfile):
    print(f"\n{'=' * 74}\n{name}: {len(smiles)} molecules\n{'=' * 74}")
    done = set()
    if os.path.exists(outfile):
        done = set(pd.read_csv(outfile).SMILES)
        print(f"  resuming, {len(done)} already done")
    todo = [s for s in smiles if s not in done]
    if not todo:
        print("  nothing to do")
        return
    cols = ["SMILES", "aff", "d_asp", "d_arg", "d_p0"]
    new = not os.path.exists(outfile)
    n = 0
    with open(outfile, "a" if not new else "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        if new:
            w.writeheader()
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            for rec in ex.map(process, todo, chunksize=4):
                w.writerow(rec)
                n += 1
                if n % 200 == 0:
                    f.flush()
                    print(f"  {n}/{len(todo)}", flush=True)


def report(path, label):
    if not os.path.exists(path):
        print(f"{label}: no data")
        return None
    d = pd.read_csv(path).dropna(subset=["d_asp", "d_arg", "d_p0"])
    n = len(d)
    print(f"\n{'-' * 74}")
    print(f"{label}   n = {n}")
    print(f"{'-' * 74}")
    print(f"  distance to Asp121 : mean {d.d_asp.mean():6.2f}  "
          f"median {d.d_asp.median():6.2f}  min {d.d_asp.min():6.2f} A")
    print(f"  distance to Arg142 : mean {d.d_arg.mean():6.2f}  "
          f"median {d.d_arg.median():6.2f}  min {d.d_arg.min():6.2f} A")
    print(f"  distance to P_0    : mean {d.d_p0.mean():6.2f}  "
          f"median {d.d_p0.median():6.2f}  min {d.d_p0.min():6.2f} A")
    print()
    a = (d.d_asp < CONTACT).sum()
    g = (d.d_arg < CONTACT).sum()
    both = ((d.d_asp < CONTACT) & (d.d_arg < CONTACT)).sum()
    p0 = (d.d_p0 < CONTACT).sum()
    print(f"  contacting Asp121        : {a:5d}  ({100*a/n:5.2f}%)")
    print(f"  contacting Arg142        : {g:5d}  ({100*g/n:5.2f}%)")
    print(f"  contacting BOTH residues : {both:5d}  ({100*both/n:5.2f}%)  <<<")
    print(f"  contacting P_0           : {p0:5d}  ({100*p0/n:5.2f}%)")
    print()
    arg_p0 = ((d.d_arg < CONTACT) & (d.d_p0 < CONTACT)).sum()
    all3 = ((d.d_asp < CONTACT) & (d.d_arg < CONTACT) & (d.d_p0 < CONTACT)).sum()
    print(f"  Arg142 + P_0 (the bridge as measured) : {arg_p0:5d}  "
          f"({100*arg_p0/n:5.2f}%)")
    print(f"  Asp121 + Arg142 + P_0 (full span)     : {all3:5d}  "
          f"({100*all3/n:5.2f}%)")
    if both:
        print(f"\n  closest dual-residue binders:")
        t = d[(d.d_asp < CONTACT) & (d.d_arg < CONTACT)].nsmallest(5, "d_asp")
        for _, r in t.iterrows():
            print(f"    Asp {r.d_asp:5.2f}  Arg {r.d_arg:5.2f}  P_0 {r.d_p0:5.2f}  "
                  f"aff {r.aff:7.3f}  {r.SMILES[:40]}")
    else:
        print(f"\n  NO molecule contacts both residues.")
        print(f"  closest approach to Asp121 among Arg142 binders:")
        t = d[d.d_arg < CONTACT].nsmallest(5, "d_asp")
        for _, r in t.iterrows():
            print(f"    Asp {r.d_asp:5.2f}  Arg {r.d_arg:5.2f}  P_0 {r.d_p0:5.2f}  "
                  f"aff {r.aff:7.3f}  {r.SMILES[:40]}")
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", type=int, default=0)
    args = ap.parse_args()

    print("=" * 74)
    print("SEPARATED PHARMACOPHORE ANALYSIS")
    print("=" * 74)
    print(f"Asp121 side-chain atoms : {len(ASP_ATOMS)}")
    print(f"Arg142 side-chain atoms : {len(ARG_ATOMS)}")
    print(f"P_0 atoms               : {len(P0_ATOMS)}")
    sep = min(math.dist(a, b) for a in ASP_ATOMS for b in ARG_ATOMS)
    print(f"\nAsp121 to Arg142, closest side-chain atoms: {sep:.2f} A")
    print("  (a ligand must span at least this to touch both)")
    if len(ASP_ATOMS) < 4 or len(ARG_ATOMS) < 7:
        sys.exit("!! reference atoms not found - check residue numbering")

    # Arm B
    b = pd.read_csv(os.path.join(R, "loa22_generation_1.csv"))
    b = b[b.SMILES_state == 1].dropna(subset=["SMILES"])
    bs = b.SMILES.tolist()

    # Arm C
    cpath = os.path.join(R, "arm_c", "armC_generation_1.csv")
    cs = []
    if os.path.exists(cpath):
        c = pd.read_csv(cpath)
        c = c[c.SMILES_state == 1].dropna(subset=["SMILES"])
        cs = c.SMILES.tolist()

    if args.test:
        bs, cs = bs[:args.test], cs[:args.test]
        suffix = "_test"
    else:
        suffix = ""

    fb = os.path.join(OUTDIR, f"armB_split{suffix}.csv")
    fc = os.path.join(OUTDIR, f"armC_split{suffix}.csv")

    run_arm("ARM B (naive RL at P_0)", bs, fb)
    if cs:
        run_arm("ARM C (bridging objective)", cs, fc)

    print("\n" + "=" * 74)
    print("RESULTS")
    print("=" * 74)
    db = report(fb, "ARM B")
    dc = report(fc, "ARM C") if cs else None

    if db is not None and dc is not None:
        print("\n" + "=" * 74)
        print("CORRECTED ARM COMPARISON")
        print("=" * 74)
        for label, fn in (("Arg142 contact", lambda d: (d.d_arg < CONTACT).mean()),
                          ("Arg142 + P_0",
                           lambda d: ((d.d_arg < CONTACT) & (d.d_p0 < CONTACT)).mean()),
                          ("Asp121 contact", lambda d: (d.d_asp < CONTACT).mean()),
                          ("both residues",
                           lambda d: ((d.d_asp < CONTACT) & (d.d_arg < CONTACT)).mean())):
            vb, vc = fn(db), fn(dc)
            enr = f"{vc/vb:.2f}x" if vb > 0 else "n/a"
            print(f"  {label:<20} B {100*vb:6.2f}%   C {100*vc:6.2f}%   "
                  f"enrichment {enr}")

    print("\n" + "=" * 74)


if __name__ == "__main__":
    main()
