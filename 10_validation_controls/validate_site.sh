#!/bin/bash
# =====================================================================
# FUNCTIONAL SITE VALIDATION (v2)
#   - extracts P_4 from the DoGSiteScorer archive
#   - defines the box from the CDD-annotated functional site, not from
#     P_4's 145 A^3 sub-cavity (which is smaller than the ligands)
#   - re-docks the three leads + doxycycline, reporting pharmacophore contact
#
# Run:  bash ~/Research/validate_site.sh
# =====================================================================
set -u
R=~/Research
W=$R/site_validation
PY=~/miniconda3/envs/reinvent4/bin/python
mkdir -p "$W"; cd "$W" || exit 1

echo "########## FUNCTIONAL SITE VALIDATION $(date) ##########"

# ---------------------------------------------------------------- 1
echo
echo "=== 1. Extracting pocket files from the DoGSiteScorer archive ==="
$PY - <<'PY'
import glob, zipfile, os
desk = "/mnt/c/Users/donmi/Desktop/Research"
out = os.path.expanduser("~/Research/site_validation")
n = 0
for z in glob.glob(os.path.join(desk, "*.zip")):
    try:
        zf = zipfile.ZipFile(z)
    except Exception:
        continue
    for name in zf.namelist():
        if name.endswith(".pdb") and "_res" in name:
            tgt = os.path.join(out, os.path.basename(name))
            open(tgt, "wb").write(zf.read(name))
            n += 1
print(f"  extracted {n} pocket residue files")
for f in sorted(glob.glob(os.path.join(out, "*_res.pdb"))):
    res = {int(l[22:26]) for l in open(f) if l.startswith(("ATOM", "HETATM"))}
    tag = "  <-- functional site" if (121 in res or 142 in res) else ""
    print(f"    {os.path.basename(f)[-14:]:>14}  {len(res):3d} residues{tag}")
PY

# ---------------------------------------------------------------- 2
echo
echo "=== 2. Defining the docking box from the FUNCTIONAL SITE ==="
echo "    (CDD-annotated peptidoglycan residues + Asp121/Arg142)"
$PY - <<'PY' | tee box_site.txt
import os, math
R = os.path.expanduser("~/Research")
CDD = [80, 81, 120, 121, 124, 128, 135, 138, 178, 182]
KEY = {121, 142}
pts, seen = [], set()
for l in open(os.path.join(R, "loa22_receptor.pdbqt")):
    if not l.startswith(("ATOM", "HETATM")):
        continue
    try:
        rn = int(l[22:26])
    except ValueError:
        continue
    if rn not in CDD and rn not in KEY:
        continue
    if l[77:79].strip() in ("HD", "H"):
        continue
    pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    seen.add(rn)
cen = [sum(p[i] for p in pts) / len(pts) for i in range(3)]
ext = [max(p[i] for p in pts) - min(p[i] for p in pts) for i in range(3)]
size = max(24.0, min(30.0, max(ext) + 6.0))
print(f"# functional-site residues found: {sorted(seen)}")
print(f"# heavy atoms: {len(pts)}")
print(f"# site span: {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} A")
print(f"center_x = {cen[0]:.3f}")
print(f"center_y = {cen[1]:.3f}")
print(f"center_z = {cen[2]:.3f}")
print(f"size_x = {size:.1f}")
print(f"size_y = {size:.1f}")
print(f"size_z = {size:.1f}")
PY

CX=$(awk -F'= *' '/^center_x/{print $2}' box_site.txt)
CY=$(awk -F'= *' '/^center_y/{print $2}' box_site.txt)
CZ=$(awk -F'= *' '/^center_z/{print $2}' box_site.txt)
SZ=$(awk -F'= *' '/^size_x/{print $2}' box_site.txt)
if [ -z "$CX" ]; then echo "!! box definition failed"; exit 1; fi
echo "    box centre ($CX, $CY, $CZ)  size ${SZ} A"

# ---------------------------------------------------------------- 3
echo
echo "=== 3. Re-docking into the functional site (exh 32, 3 seeds) ==="
$PY - "$CX" "$CY" "$CZ" "$SZ" <<'PY'
import sys, os, subprocess, tempfile, math, statistics as st
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
RDLogger.DisableLog("rdApp.*")

CX, CY, CZ, SZ = map(float, sys.argv[1:5])
R = os.path.expanduser("~/Research")
REC = os.path.join(R, "loa22_receptor.pdbqt")
OUT = os.path.join(R, "site_validation")

KEY = {121: {"CB", "CG", "OD1", "OD2"},
       142: {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}}
key = []
for l in open(REC):
    if l.startswith(("ATOM", "HETATM")):
        try:
            rn = int(l[22:26])
        except ValueError:
            continue
        if rn in KEY and l[12:16].strip() in KEY[rn]:
            key.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))

MOLS = {
 "LOA22-1": "CC(C)CNC(=O)Nc1ccc2nnc(-c3cc(C(F)(F)F)cc(C(F)(F)F)c3)n2n1",
 "LOA22-2": "CC(CC(=O)N1CCc2ccc(C(F)(F)F)cc2C1)c1ccc(O)cc1",
 "LOA22-3": "CC(=O)Nc1ccc(F)c(C(=O)NS(=O)(=O)c2ccccc2Cl)c1",
}
doxy = os.path.join(R, "prior_control", "doxycycline_correct.pdbqt")


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


def dock(lig, out, seed):
    cmd = ["vina", "--receptor", REC, "--ligand", lig,
           "--center_x", str(CX), "--center_y", str(CY), "--center_z", str(CZ),
           "--size_x", str(SZ), "--size_y", str(SZ), "--size_z", str(SZ),
           "--exhaustiveness", "32", "--num_modes", "9", "--seed", str(seed),
           "--out", out]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    for line in r.stdout.splitlines():
        p = line.split()
        if p and p[0] == "1":
            try:
                return float(p[1])
            except (IndexError, ValueError):
                return None
    return None


def mindist(pose):
    best = None
    for l in open(pose):
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("HD", "H"):
            c = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            for k in key:
                d = math.dist(c, k)
                if best is None or d < best:
                    best = d
    return best


OLD = {"LOA22-1": -9.099, "LOA22-2": -8.955, "LOA22-3": -8.945,
       "doxycycline": -6.862}
print(f"\n{'compound':<12}{'site':>9}{'SD':>7}{'P_0':>9}{'delta':>8}"
      f"{'minDist':>9}  engages")
print("-" * 64)
for name, smi in list(MOLS.items()) + [("doxycycline", None)]:
    with tempfile.TemporaryDirectory() as td:
        lig = os.path.join(td, "l.pdbqt")
        if smi is None:
            if not os.path.exists(doxy):
                print(f"{name:<12} (pdbqt missing)"); continue
            lig = doxy
        elif not prep(smi, lig):
            print(f"{name:<12} prep failed"); continue
        vals, poses = [], []
        for sd in (12345, 67890, 24680):
            po = os.path.join(OUT, f"{name}_site_{sd}.pdbqt")
            v = dock(lig, po, sd)
            if v is not None:
                vals.append(v); poses.append(po)
        if not vals:
            print(f"{name:<12} docking failed"); continue
        m = st.mean(vals)
        s = st.pstdev(vals) if len(vals) > 1 else 0.0
        md = mindist(poses[0])
        old = OLD.get(name)
        print(f"{name:<12}{m:>9.3f}{s:>7.3f}{old:>9.3f}{m - old:>8.3f}"
              f"{md:>9.2f}  {'YES' if md and md < 4.5 else 'no'}")
print("\nminDist = closest heavy atom to Asp121/Arg142 side chains")
print("engages = within 4.5 A (direct pharmacophore contact)")
print("\nPoses saved to ~/Research/site_validation/ for inspection.")
PY

echo
echo "########## COMPLETE $(date) ##########"
