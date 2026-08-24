#!/usr/bin/env python3
"""
Figure S4 - structures of the three selected compounds and their scaffolds.

Drawn from the SMILES rather than reusing earlier renders, so the depicted
structures cannot drift from the strings deposited in Supplementary Table S5.
Panel B shows Bemis-Murcko scaffolds computed with RDKit, not hand-drawn.

Verified against main-text Table 10 at run time: molecular weight, QED and
Crippen logP must match the published values or the script aborts.

Run:  conda activate reinvent4
      cd ~/Research
      python3 figS4.py
Writes FigureS4.png at 300 dpi.
"""
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, QED, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold
from PIL import Image, ImageDraw, ImageFont
import io, os
import matplotlib

SMILES = {
    "LOA22-B1": "O=C(Nc1ccc(O)c(C(=O)O)c1)c1cccc(C(F)(F)F)c1",
    "LOA22-B2": "CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1",
    "LOA22-B3": "COC(=O)c1sc(NC(=O)COc2ccc(Cl)cc2F)nc1C",
}
# scaffold class as it should read. B3 was published as
# "aminothiophene-carboxylate"; the ring is a 1,3-thiazole, not a thiophene.
CLASS = {
    "LOA22-B1": "benzamide–salicylate",
    "LOA22-B2": "pyridine–oxadiazole",
    "LOA22-B3": "aminothiazole–carboxylate",
}
# main-text Table 10, used as an assertion
TABLE10 = {"LOA22-B1": (325.2, 0.755, 3.36),
           "LOA22-B2": (323.3, 0.762, 2.34),
           "LOA22-B3": (358.8, 0.831, 3.05)}
TOL = 0.06

CELL = (900, 470)
DPI = 300


def font(sz, bold=False):
    n = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    p = os.path.join(os.path.dirname(matplotlib.__file__),
                     "mpl-data", "fonts", "ttf", n)
    return ImageFont.truetype(p, sz) if os.path.exists(p) else ImageFont.load_default()


def draw(mol, size):
    AllChem.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    o = d.drawOptions()
    o.bondLineWidth = 3
    o.addStereoAnnotation = True
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText())).convert("RGB")


names = list(SMILES)
mols, scafs = {}, {}
print("VERIFICATION against Table 10")
ok = True
for n in names:
    m = Chem.MolFromSmiles(SMILES[n])
    if m is None:
        sys.exit("ERROR: %s has an invalid SMILES" % n)
    mw, qed, lp = Descriptors.MolWt(m), QED.qed(m), Crippen.MolLogP(m)
    e = TABLE10[n]
    good = abs(mw - e[0]) < 0.1 and abs(qed - e[1]) < TOL and abs(lp - e[2]) < TOL
    ok &= good
    print("  %-10s MW %7.2f (%.1f)  QED %.3f (%.3f)  logP %5.2f (%.2f)  %s"
          % (n, mw, e[0], qed, e[1], lp, e[2], "OK" if good else "MISMATCH"))
    mols[n] = m
    scafs[n] = MurckoScaffold.GetScaffoldForMol(m)
if not ok:
    sys.exit("ABORTED: structures do not match Table 10. Nothing written.")

# sanity check the corrected scaffold name
if not mols["LOA22-B3"].HasSubstructMatch(Chem.MolFromSmarts("c1cscn1")):
    sys.exit("ABORTED: LOA22-B3 does not contain a thiazole ring")
print("  LOA22-B3 thiazole ring confirmed (not thiophene)")

f_lbl, f_name, f_sub = font(46, True), font(34, True), font(28)
lbl_h, cap_h, pad = 70, 84, 26
W = CELL[0] * 3 + pad * 4
H = lbl_h + CELL[1] + cap_h + lbl_h + CELL[1] + cap_h + pad
canvas = Image.new("RGB", (W, H), "white")
dr = ImageDraw.Draw(canvas)

dr.text((pad, pad // 2), "A", fill="black", font=f_lbl)
y = lbl_h
for i, n in enumerate(names):
    x = pad + i * (CELL[0] + pad)
    canvas.paste(draw(mols[n], CELL), (x, y))
    dr.text((x + CELL[0] // 2, y + CELL[1] + 6), n, fill="black",
            font=f_name, anchor="ma")
    dr.text((x + CELL[0] // 2, y + CELL[1] + 46), CLASS[n], fill="#4a4a4a",
            font=f_sub, anchor="ma")

y2 = lbl_h + CELL[1] + cap_h
dr.text((pad, y2 - 4), "B", fill="black", font=f_lbl)
y3 = y2 + lbl_h
for i, n in enumerate(names):
    x = pad + i * (CELL[0] + pad)
    canvas.paste(draw(scafs[n], CELL), (x, y3))
    dr.text((x + CELL[0] // 2, y3 + CELL[1] + 6),
            "%s scaffold" % n, fill="black", font=f_name, anchor="ma")

canvas.save("FigureS4.png", dpi=(DPI, DPI))
print("\nwrote FigureS4.png  %d x %d px" % (canvas.width, canvas.height))
for n in names:
    print("  %-10s Murcko %s" % (n, Chem.MolToSmiles(scafs[n])))
