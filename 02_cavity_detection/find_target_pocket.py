#!/usr/bin/env python3
"""
Find the pocket that actually contains Asp121 and Arg142.

1. Check every DoGSiteScorer pocket already computed.
2. Characterise the site geometrically from the receptor itself.
3. Propose a corrected docking box.
"""
import os, glob, math, zipfile, io

R = os.path.expanduser("~/Research")
DESK = "/mnt/c/Users/donmi/Desktop/Research"
KEY = {121: {"CB", "CG", "OD1", "OD2"},
       142: {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}}
CDD = [80, 81, 120, 121, 124, 128, 135, 138, 178, 182]


def d(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def cen(p):
    return [sum(q[i] for q in p) / len(p) for i in range(3)]


# ------------------------------------------------ receptor
atoms, key_atoms, cdd_atoms = [], [], []
for l in open(os.path.join(R, "loa22_receptor.pdbqt")):
    if not l.startswith(("ATOM", "HETATM")):
        continue
    try:
        rn = int(l[22:26])
    except ValueError:
        continue
    nm = l[12:16].strip()
    c = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    heavy = l[77:79].strip() not in ("HD", "H")
    atoms.append((rn, nm, c, heavy))
    if rn in KEY and nm in KEY[rn]:
        key_atoms.append(c)
    if rn in CDD and nm not in ("N", "CA", "C", "O") and heavy:
        cdd_atoms.append(c)

kc = cen(key_atoms)
print("=" * 70)
print("TARGET SITE")
print("=" * 70)
print(f"Asp121+Arg142 side-chain centroid : "
      f"({kc[0]:.2f}, {kc[1]:.2f}, {kc[2]:.2f})")
if cdd_atoms:
    cc = cen(cdd_atoms)
    print(f"CDD 10-residue site centroid      : "
          f"({cc[0]:.2f}, {cc[1]:.2f}, {cc[2]:.2f})")
    print(f"  separation from key centroid    : {d(cc, kc):.2f} A")

# ------------------------------------------------ existing pockets
print()
print("=" * 70)
print("EXISTING DoGSiteScorer POCKETS")
print("=" * 70)

cands = []
for pat in ("P_*.pdb", "pocket*.pdb", "*residues*.pdb"):
    cands += glob.glob(os.path.join(R, pat)) + glob.glob(os.path.join(DESK, pat))
for z in glob.glob(os.path.join(DESK, "*.zip")):
    try:
        with zipfile.ZipFile(z) as zf:
            for n in zf.namelist():
                if n.endswith((".pdb", ".txt")) and ("P_" in n or "pocket" in n.lower()):
                    cands.append(f"{z}::{n}")
    except Exception:
        pass

if not cands:
    print("No pocket files found beyond P_0.")
    print("-> Re-run DoGSiteScorer, or download the full results archive.")
else:
    seen = set()
    for c in sorted(set(cands)):
        try:
            if "::" in c:
                zp, inner = c.split("::")
                with zipfile.ZipFile(zp) as zf:
                    txt = zf.read(inner).decode("utf-8", "ignore")
            else:
                txt = open(c, errors="ignore").read()
        except Exception:
            continue
        pts, res = [], set()
        for l in txt.splitlines():
            if l.startswith(("ATOM", "HETATM")):
                try:
                    pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
                    res.add(int(l[22:26]))
                except ValueError:
                    pass
        if not pts or frozenset(res) in seen:
            continue
        seen.add(frozenset(res))
        pcen = cen(pts)
        mind = min(d(p, a) for p in pts for a in key_atoms)
        has = [r for r in (121, 142) if r in res]
        ncdd = len([r for r in CDD if r in res])
        flag = "  <<< CONTAINS TARGET" if has else ""
        print(f"\n{os.path.basename(c)}")
        print(f"  atoms {len(pts):4d}  residues {len(res):3d}")
        print(f"  centroid ({pcen[0]:7.2f},{pcen[1]:7.2f},{pcen[2]:7.2f})"
              f"   {d(pcen, kc):6.2f} A from target")
        print(f"  contains 121/142: {has if has else 'neither'}   "
              f"CDD residues: {ncdd}/10   min dist {mind:.2f} A{flag}")

# ------------------------------------------------ characterise the site
print()
print("=" * 70)
print("GEOMETRY AT THE TARGET SITE")
print("=" * 70)

for rad in (8, 10, 12):
    near = {rn for rn, nm, c, h in atoms
            if h and min(d(c, a) for a in key_atoms) <= rad}
    print(f"  within {rad:2d} A : {len(near):3d} residues")

near12 = sorted({rn for rn, nm, c, h in atoms
                 if h and min(d(c, a) for a in key_atoms) <= 12})
print(f"\n  residues lining the site (12 A): {near12}")
print(f"  CDD residues present: {[r for r in CDD if r in near12]}")

# crude enclosure estimate: protein density in a shell around the centroid
shell = [c for rn, nm, c, h in atoms if h and 3 <= d(c, kc) <= 10]
print(f"\n  heavy atoms in 3-10 A shell around site centroid: {len(shell)}")
print("  (a buried cavity typically shows 150+; an open groove far fewer)")

# ------------------------------------------------ proposed box
print()
print("=" * 70)
print("PROPOSED CORRECTED DOCKING BOX")
print("=" * 70)
site_atoms = [c for rn, nm, c, h in atoms
              if h and min(d(c, a) for a in key_atoms) <= 10]
sc = cen(site_atoms)
print(f"  centre : ({sc[0]:.2f}, {sc[1]:.2f}, {sc[2]:.2f})")
print(f"  size   : 22 x 22 x 22 A")
print(f"  check  : Asp121 {min(d(sc, a) for a in key_atoms[:4]):.2f} A, "
      f"Arg142 {min(d(sc, a) for a in key_atoms[4:]):.2f} A from centre")
print(f"           both well inside an 11 A half-width: "
      f"{max(d(sc, a) for a in key_atoms) < 11}")
print()
print("  vina flags:")
print(f"    --center_x {sc[0]:.2f} --center_y {sc[1]:.2f} --center_z {sc[2]:.2f} \\")
print( "    --size_x 22 --size_y 22 --size_z 22")
print("=" * 70)
