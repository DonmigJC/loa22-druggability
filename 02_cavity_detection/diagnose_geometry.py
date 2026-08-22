#!/usr/bin/env python3
"""
Diagnose two questions:
  A. Is the docking grid box correctly centred on Asp121/Arg142?
  B. Did the MD complexes start from the DOCKED pose, or from an
     arbitrary RDKit-embedded position?
"""
import os, math, glob

R = os.path.expanduser("~/Research")
MD = os.path.join(R, "04_Molecular_Dynamics")
BOX_CENTER = (2.30, 7.09, 2.02)          # angstroms
BOX_SIZE = 25.0

KEY = {121: {"CB", "CG", "OD1", "OD2"},
       142: {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}}


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


# ---------------------------------------------------------------- A
print("=" * 68)
print("A.  DOCKING BOX vs KEY RESIDUES")
print("=" * 68)

key_atoms, per_res = [], {121: [], 142: []}
rec = os.path.join(R, "loa22_receptor.pdbqt")
for line in open(rec):
    if not line.startswith(("ATOM", "HETATM")):
        continue
    try:
        rn = int(line[22:26])
    except ValueError:
        continue
    if rn in KEY and line[12:16].strip() in KEY[rn]:
        c = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        key_atoms.append(c); per_res[rn].append(c)

print(f"key side-chain atoms found : {len(key_atoms)}  "
      f"(Asp121 {len(per_res[121])}, Arg142 {len(per_res[142])})")
if not key_atoms:
    print("!! none found - residue numbering problem"); raise SystemExit

for rn in (121, 142):
    if per_res[rn]:
        cen = [sum(a[i] for a in per_res[rn]) / len(per_res[rn]) for i in range(3)]
        d = dist(cen, BOX_CENTER)
        inside = all(abs(cen[i] - BOX_CENTER[i]) <= BOX_SIZE / 2 for i in range(3))
        print(f"  residue {rn} centroid  ({cen[0]:7.2f},{cen[1]:7.2f},{cen[2]:7.2f})"
              f"   {d:6.2f} A from box centre   inside box: {inside}")

allc = [sum(a[i] for a in key_atoms) / len(key_atoms) for i in range(3)]
print(f"  combined key centroid ({allc[0]:7.2f},{allc[1]:7.2f},{allc[2]:7.2f})"
      f"   {dist(allc, BOX_CENTER):6.2f} A from box centre")
print(f"\n  IDEAL box centre would be ({allc[0]:.2f}, {allc[1]:.2f}, {allc[2]:.2f})")

# pocket file: centroid vs bounding-box midpoint
p0 = os.path.join(R, "P_0_residues.pdb")
if os.path.exists(p0):
    pts = [(float(l[30:38]), float(l[38:46]), float(l[46:54]))
           for l in open(p0) if l.startswith(("ATOM", "HETATM"))]
    if pts:
        cen = [sum(p[i] for p in pts) / len(pts) for i in range(3)]
        mid = [(min(p[i] for p in pts) + max(p[i] for p in pts)) / 2 for i in range(3)]
        ext = [max(p[i] for p in pts) - min(p[i] for p in pts) for i in range(3)]
        print(f"\n  P_0 atoms              : {len(pts)}")
        print(f"  P_0 true centroid      : ({cen[0]:7.2f},{cen[1]:7.2f},{cen[2]:7.2f})")
        print(f"  P_0 bbox midpoint      : ({mid[0]:7.2f},{mid[1]:7.2f},{mid[2]:7.2f})")
        print(f"  centroid vs midpoint   : {dist(cen, mid):.2f} A apart")
        print(f"  box centre used        : {BOX_CENTER}")
        print(f"    -> matches centroid  : {dist(cen, BOX_CENTER):.2f} A")
        print(f"    -> matches midpoint  : {dist(mid, BOX_CENTER):.2f} A")
        print(f"  P_0 bbox dimensions    : "
              f"{ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} A  (box is {BOX_SIZE} A)")
        if max(ext) > BOX_SIZE:
            print("  !! pocket is LARGER than the search box on at least one axis")

# ---------------------------------------------------------------- B
print()
print("=" * 68)
print("B.  MD STARTING LIGAND POSITION")
print("=" * 68)
print("A .gro file stores nanometres; converted to angstroms below.\n")

for lead in ("Lead_1", "Lead_4", "Lead_5"):
    gro = os.path.join(MD, lead, "md_run", f"complex_{lead}.gro")
    if not os.path.exists(gro):
        print(f"{lead}: complex .gro not found"); continue
    lines = open(gro).read().splitlines()
    n = int(lines[1].strip())
    lig = []
    for l in lines[2:2 + n]:
        if l[5:10].strip() in ("UNL", lead):
            lig.append((float(l[20:28]) * 10, float(l[28:36]) * 10,
                        float(l[36:44]) * 10))
    if not lig:
        print(f"{lead}: no UNL atoms found in complex .gro"); continue
    cen = [sum(a[i] for a in lig) / len(lig) for i in range(3)]
    dmin = min(dist(a, k) for a in lig for k in key_atoms)
    print(f"{lead}: {len(lig)} ligand atoms")
    print(f"   centroid            ({cen[0]:7.2f},{cen[1]:7.2f},{cen[2]:7.2f}) A")
    print(f"   dist to box centre  {dist(cen, BOX_CENTER):6.2f} A")
    print(f"   MIN dist to key res {dmin:6.2f} A")
    print(f"   -> {'AT the pharmacophore' if dmin < 5 else 'NOT at the pharmacophore'}")

# ---------------------------------------------------------------- C
print()
print("=" * 68)
print("C.  A DOCKED POSE FOR COMPARISON (doxycycline)")
print("=" * 68)
for f in sorted(glob.glob(os.path.join(R, "prior_control", "doxy_docked_*.pdbqt")))[:1]:
    pose = []
    for l in open(f):
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("HD", "H"):
            pose.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    if pose:
        cen = [sum(a[i] for a in pose) / len(pose) for i in range(3)]
        dmin = min(dist(a, k) for a in pose for k in key_atoms)
        print(f"{os.path.basename(f)}: {len(pose)} heavy atoms")
        print(f"   centroid            ({cen[0]:7.2f},{cen[1]:7.2f},{cen[2]:7.2f}) A")
        print(f"   dist to box centre  {dist(cen, BOX_CENTER):6.2f} A")
        print(f"   MIN dist to key res {dmin:6.2f} A")
print("=" * 68)
