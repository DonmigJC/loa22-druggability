#!/usr/bin/env python3
"""
Supplementary Table S3 - complete DoGSiteScorer output for all six cavities.

Nothing is transcribed. Descriptors come from DoGSiteScorer Results/*_desc.txt;
lining residues, centroids and distances are computed from the *_P_n_res.pdb
files against the Asp121/Arg142 side-chain centroid.

Writes:  TableS3.md   paste into Word
         TableS3.csv  deposit with the supplement
Run from the Research root.
"""
import csv, glob, os, re
import numpy as np

CDD = [80, 81, 120, 121, 124, 128, 135, 138, 178, 182]   # cd07185 annotated
PHARM = (121, 142)

d = glob.glob(os.path.join("DoGSiteScorer Results", "*_desc.txt"))
if not d:
    raise SystemExit("desc.txt not found - run from the Research root")
rows = list(csv.reader(open(d[0]), delimiter="\t"))
ix = {h: i for i, h in enumerate(rows[0])}
D = {}
for r in rows[1:]:
    if r:
        D[int(r[ix["name"]].split("_")[1])] = r

atoms = {}
for l in open("loa22_alphafold2_raw.pdb"):
    if l.startswith("ATOM"):
        atoms.setdefault(int(l[22:26]), []).append(
            (float(l[30:38]), float(l[38:46]), float(l[46:54]), l[12:16].strip()))
BB = {"N", "CA", "C", "O", "OXT"}
pc = np.array([a[:3] for r in PHARM for a in atoms[r] if a[3] not in BB]).mean(0)

POC = {}
for f in glob.glob("site_validation/*_P_*_res.pdb"):
    k = int(re.search(r"_P_(\d)_res", f).group(1))
    pts, res = [], set()
    for l in open(f):
        if l.startswith(("ATOM", "HETATM")):
            pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            res.add(int(l[22:26]))
    c = np.array(pts).mean(0)
    POC[k] = dict(cen=c, res=sorted(res), natoms=len(pts),
                  dist=float(np.linalg.norm(c - pc)))

def g(k, col):
    return D[k][ix[col]].strip()

HDR = ["Pocket", "Drug Score", "Simple Score", "Volume (Å³)", "Surface (Å²)",
       "Depth (Å)", "Enclosure", "Hydrophobicity ratio", "Lining atoms",
       "Lining residues", "H-bond acceptors", "H-bond donors",
       "Apolar / polar / +ve / −ve AA", "Centroid (x, y, z)",
       "Centroid to pharmacophore (Å)", "CDD residues contained",
       "Contains Asp121", "Contains Arg142"]

table = []
for k in sorted(POC):
    r = POC[k]
    n_cdd = sum(1 for c in CDD if c in r["res"])
    table.append([
        "P_%d" % k, "%.2f" % float(g(k, "drugScore")), g(k, "simpleScore"),
        g(k, "volume"), g(k, "surface"), g(k, "depth"), g(k, "enclosure"),
        g(k, "hydrophobicity"), str(r["natoms"]), str(len(r["res"])),
        g(k, "accept"), g(k, "donor"),
        "%s / %s / %s / %s" % (g(k, "apolarAA"), g(k, "polarAA"),
                               g(k, "posAA"), g(k, "negAA")),
        "(%.2f, %.2f, %.2f)" % tuple(r["cen"]), "%.2f" % r["dist"],
        "%d / 10" % n_cdd,
        "Yes" if PHARM[0] in r["res"] else "No",
        "Yes" if PHARM[1] in r["res"] else "No"])

with open("TableS3.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(HDR); w.writerows(table)

# transposed markdown: descriptors down the side, pockets across
with open("TableS3.md", "w") as fh:
    fh.write("| Descriptor | " + " | ".join(r[0] for r in table) + " |\n")
    fh.write("|" + "---|" * (len(table) + 1) + "\n")
    for i, h in enumerate(HDR[1:], start=1):
        fh.write("| %s | %s |\n" % (h, " | ".join(r[i] for r in table)))
    fh.write("\n**Lining residues**\n\n")
    for k in sorted(POC):
        fh.write("- **P_%d** (%d): %s\n" % (k, len(POC[k]["res"]),
                 ", ".join(str(x) for x in POC[k]["res"])))

print("wrote TableS3.md and TableS3.csv\n")
print(open("TableS3.md").read())
