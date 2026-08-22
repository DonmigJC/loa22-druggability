#!/usr/bin/env python3
"""
Where is the ligand, really?

Reads static PDB snapshots dumped from the PBC-corrected trajectory and
measures, with plain coordinate geometry, the ligand's closest approach to:

  - Asp121 / Arg142 side chains  (the pharmacophore)
  - the 10 CDD-annotated peptidoglycan-binding residues
  - P_0, the druggable but non-functional pocket 12.9 A away

No trajectory tools, no periodicity, no structural fitting - so the
answer is unambiguous.
"""
import os
import math

KEY = {121: {"CB", "CG", "OD1", "OD2"},
       142: {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}}
CDD = [80, 81, 120, 121, 124, 128, 135, 138, 178, 182]
P0 = [61, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 88, 89, 92, 93,
      96, 97, 115, 116, 117, 119, 163, 165, 169, 171, 172, 183, 185, 186, 187]

BASE = os.path.expanduser("~/Research/04_Molecular_Dynamics")
LEADS = {"Lead_1": "LOA22-1", "Lead_4": "LOA22-2", "Lead_5": "LOA22-3"}
TIMES = [0, 25000, 50000]


def closest(lig, group):
    if not group:
        return None
    return min(math.dist(a, b) for a in lig for b in group) / 10.0   # A -> nm


print("=" * 72)
print("LIGAND POSITION IN STATIC FRAMES (nm, closest heavy-atom approach)")
print("=" * 72)
print(f"{'system':<22}{'Asp121/Arg142':>14}{'CDD site':>10}{'P_0':>8}   verdict")
print("-" * 72)

summary = {}
for folder, label in LEADS.items():
    for t in TIMES:
        f = os.path.join(BASE, folder, "md_run", "analysis_v2",
                         "frames", f"{folder}_{t}ps.pdb")
        if not os.path.exists(f):
            print(f"{label} {t//1000:>2}ns  frame missing")
            continue
        key, cdd, p0, lig = [], [], [], []
        for l in open(f):
            if not l.startswith(("ATOM", "HETATM")):
                continue
            nm = l[12:16].strip()
            if nm.startswith("H") or l[76:78].strip() == "H":
                continue
            try:
                rn = int(l[22:26])
            except ValueError:
                continue
            rs = l[17:20].strip()
            c = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            if rs in ("UNL", "LIG"):
                lig.append(c)
                continue
            if rn in KEY and nm in KEY[rn]:
                key.append(c)
            if rn in CDD and nm not in ("N", "CA", "C", "O"):
                cdd.append(c)
            if rn in P0:
                p0.append(c)

        if not lig:
            print(f"{label} {t//1000:>2}ns  no ligand atoms found")
            continue

        dk, dc, dp = closest(lig, key), closest(lig, cdd), closest(lig, p0)
        if dk is not None and dk < 0.45:
            v = "AT PHARMACOPHORE"
        elif dp is not None and dk is not None and dp < dk:
            v = "nearer P_0"
        else:
            v = "elsewhere"
        summary.setdefault(label, []).append((dk, dc, dp, v))
        print(f"{label + ' ' + str(t // 1000) + 'ns':<22}"
              f"{dk:>14.3f}{dc:>10.3f}{dp:>8.3f}   {v}")

print("-" * 72)
print("\nSUMMARY")
for label, rows in summary.items():
    n_ph = sum(1 for r in rows if r[3] == "AT PHARMACOPHORE")
    mk = sum(r[0] for r in rows) / len(rows)
    mp = sum(r[2] for r in rows) / len(rows)
    pref = "pharmacophore" if mk < mp else "P_0"
    print(f"  {label}: {n_ph}/{len(rows)} frames at pharmacophore   "
          f"mean {mk:.3f} nm to key vs {mp:.3f} nm to P_0   -> prefers {pref}")

print("\nReference: van der Waals contact is 0.15-0.35 nm.")
print("Frames are 0, 25 and 50 ns from the PBC-corrected trajectory.")
print("=" * 72)
