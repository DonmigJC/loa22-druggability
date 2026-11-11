#!/usr/bin/env python3
"""Recompute the paper's central structural numbers from the AlphaFold model.

Claims checked (manuscript v21):
  2.1 / 3.1 / abstract / significance : P0 to pharmacophore centroid 12.91 A
  2.1 / 3.1                           : nearest P0 atom to nearest
                                        pharmacophoric side-chain atom 6.95 A
  2.1                                 : Asp121 to Arg142 closest side-chain
                                        atoms 5.38 A
  2.1                                 : Asp121 side-chain centroid 16.27 A from
                                        the P0 docking box centre, box
                                        half-width 12.5 A, so not enclosed
  2.1                                 : P4 Drug Score 0.38 "ranks third of six"
  2.1 / 3.1                           : two of the ten CDD residues fall
                                        outside every detected cavity
"""
import numpy as np
import os

R = "/home/donmigcage/Research"
PDB = os.path.join(R, "loa22_alphafold2_raw.pdb")
P0PDB = os.path.join(R, "P_0_residues.pdb")

KEY = {121: {"CB", "CG", "OD1", "OD2"},
       142: {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}}
BACKBONE = {"N", "CA", "C", "O", "OXT"}

CDD = [80, 81, 120, 121, 124, 128, 135, 138, 178, 182]
POCKET_LINING = {
    0: [61, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 88, 89, 92, 93, 96,
        97, 115, 116, 117, 119, 163, 165, 169, 171, 172, 183, 185, 186, 187],
    1: [48, 49, 50, 52, 54, 55, 56, 58, 59, 60, 63, 103, 107],
    2: [85, 86, 87, 90, 91, 145, 148, 149, 152, 153],
    3: [110, 111, 112, 113, 158, 159, 161, 191, 192, 193, 194],
    4: [80, 81, 82, 120, 121, 123, 124, 125, 132, 133, 134, 135, 138, 141, 142, 182],
    5: [53, 54, 57, 58, 61, 112, 187, 188, 189, 190],
}
DRUGSCORE = {0: 0.81, 1: 0.58, 2: 0.36, 3: 0.30, 4: 0.38, 5: 0.20}
BOX_P0_CENTRE = np.array([2.30, 7.09, 2.02])
BOX_P0_HALF = 12.5


def read_pdb(path):
    out = []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        resn = line[17:20].strip()
        try:
            resi = int(line[22:26])
        except ValueError:
            continue
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        elem = line[76:78].strip() or name[:1]
        out.append(dict(name=name, resn=resn, resi=resi, xyz=xyz, elem=elem))
    return out


def check(label, measured, claimed, tol=0.02, unit="A"):
    ok = abs(measured - claimed) <= tol
    print("  %-58s measured %8.3f %s   paper %7.2f   %s"
          % (label, measured, unit, claimed, "OK" if ok else "MISMATCH"))
    return ok


prot = read_pdb(PDB)
p0 = read_pdb(P0PDB)
print("AlphaFold model: %d atoms, %d residues" % (len(prot), len({a["resi"] for a in prot})))
print("P_0_residues.pdb: %d atoms" % len(p0))
print()

# --- pharmacophore atom sets ----------------------------------------------
key_atoms = [a for a in prot if a["resi"] in KEY and a["name"] in KEY[a["resi"]]]
sc_atoms = [a for a in prot if a["resi"] in KEY and a["name"] not in BACKBONE
            and a["elem"] != "H"]
print("pharmacophoric atoms, KEY definition  : %d  (expected 11)" % len(key_atoms))
print("all side-chain heavy atoms of 121+142 : %d" % len(sc_atoms))
print()

key_xyz = np.array([a["xyz"] for a in key_atoms])
sc_xyz = np.array([a["xyz"] for a in sc_atoms])
p0_xyz = np.array([a["xyz"] for a in p0])
p0_cen = p0_xyz.mean(axis=0)

print("P0 centroid from P_0_residues.pdb : (%.2f, %.2f, %.2f)   TableS3 says (3.10, 4.78, 2.72)"
      % tuple(p0_cen))
print()

print("CENTROID SEPARATION, the abstract and significance headline")
for label, xyz in (("KEY atoms (CB,CG,OD1,OD2 + CB..NH2)", key_xyz),
                   ("all side-chain heavy atoms of 121+142", sc_xyz)):
    cen = xyz.mean(axis=0)
    check("P0 centroid to pharmacophore centroid, %s" % label,
          float(np.linalg.norm(p0_cen - cen)), 12.91, tol=0.05)
print()

print("NEAREST-ATOM SEPARATION")
d = np.sqrt(((p0_xyz[:, None, :] - key_xyz[None, :, :]) ** 2).sum(-1))
check("nearest P0 atom to nearest KEY atom", float(d.min()), 6.95, tol=0.05)
d2 = np.sqrt(((p0_xyz[:, None, :] - sc_xyz[None, :, :]) ** 2).sum(-1))
check("nearest P0 atom to nearest side-chain heavy atom", float(d2.min()), 6.95, tol=0.05)
print()

print("Asp121 TO Arg142")
a121 = np.array([a["xyz"] for a in key_atoms if a["resi"] == 121])
a142 = np.array([a["xyz"] for a in key_atoms if a["resi"] == 142])
d3 = np.sqrt(((a121[:, None, :] - a142[None, :, :]) ** 2).sum(-1))
check("closest Asp121 to Arg142 side-chain atoms", float(d3.min()), 5.38, tol=0.05)
print()

print("DOCKING BOX CONTAINMENT")
c121 = a121.mean(axis=0)
c142 = a142.mean(axis=0)
check("Asp121 side-chain centroid to P0 box centre",
      float(np.linalg.norm(c121 - BOX_P0_CENTRE)), 16.27, tol=0.05)
print("    box half-width %.1f A, so Asp121 centroid %s enclosed"
      % (BOX_P0_HALF,
         "IS" if np.all(np.abs(c121 - BOX_P0_CENTRE) <= BOX_P0_HALF) else "is NOT"))
print("    Arg142 centroid distance %.2f A, per-axis offset (%.2f, %.2f, %.2f)"
      % (float(np.linalg.norm(c142 - BOX_P0_CENTRE)),
         *(np.abs(c142 - BOX_P0_CENTRE))))
print("    paper says the box 'placed the Arg142 centroid at its boundary'")
print()

print("DRUG SCORE RANKING OF P4")
order = sorted(DRUGSCORE.items(), key=lambda kv: -kv[1])
print("    descending: " + ", ".join("P%d=%.2f" % (k, v) for k, v in order))
rank = [k for k, _ in order].index(4) + 1
print("    P4 Drug Score 0.38 ranks %d of 6.   paper says 'third of the six'   %s"
      % (rank, "OK" if rank == 3 else "MISMATCH"))
print()

print("CDD RESIDUES OUTSIDE EVERY CAVITY")
inside = set()
for pk, res in POCKET_LINING.items():
    inside |= set(res) & set(CDD)
outside = [r for r in CDD if r not in inside]
print("    the ten annotated: %s" % CDD)
print("    P4 contains %d of them: %s"
      % (len(set(POCKET_LINING[4]) & set(CDD)),
         sorted(set(POCKET_LINING[4]) & set(CDD))))
print("    outside every cavity: %s  (%d)   paper says two   %s"
      % (outside, len(outside), "OK" if len(outside) == 2 else "MISMATCH"))
print("    P0 contains %d of the ten   paper says none   %s"
      % (len(set(POCKET_LINING[0]) & set(CDD)),
         "OK" if not (set(POCKET_LINING[0]) & set(CDD)) else "MISMATCH"))
