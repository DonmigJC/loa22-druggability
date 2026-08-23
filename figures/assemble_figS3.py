#!/usr/bin/env python3
"""
Assemble Figure S3: three orientations side by side, with a colour legend
listing each pocket's volume, enclosure and distance to the pharmacophore.

Every number is read from the DoGSiteScorer descriptor table and the pocket
residue files, not retyped. Run from the Research root.
"""
import csv, glob, os, re
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import matplotlib

DPI, W_MM = 300, 180
W = int(W_MM / 25.4 * DPI)
MARGIN, GAP = int(0.012 * W), int(0.020 * W)
OUT = "FigureS3.png"
VIEWS = ["FigureS3_v1.png", "FigureS3_v2.png", "FigureS3_v3.png"]
ANGLES = ["0°", "60°", "120°"]
HEX = {0: "#3380CC", 1: "#E69F00", 2: "#009E73",
       3: "#CC79A7", 4: "#B22222", 5: "#6A3D9A"}

def font(sz, bold=False):
    n = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    p = os.path.join(os.path.dirname(matplotlib.__file__),
                     "mpl-data", "fonts", "ttf", n)
    return ImageFont.truetype(p, sz) if os.path.exists(p) else ImageFont.load_default()

def autocrop(im, tol=248, pad=10):
    m = im.convert("L").point(lambda v: 255 if v < tol else 0)
    bb = m.getbbox()
    if not bb: return im
    l, t, r, b = bb
    return im.crop((max(0, l-pad), max(0, t-pad),
                    min(im.width, r+pad), min(im.height, b+pad)))

# ---- descriptors straight from the DoGSiteScorer output
d = glob.glob(os.path.join("DoGSiteScorer Results", "*_desc.txt"))
if not d: raise SystemExit("desc.txt not found")
rows = list(csv.reader(open(d[0]), delimiter="\t"))
ix = {h: i for i, h in enumerate(rows[0])}
DESC = {}
for r in rows[1:]:
    if not r: continue
    k = int(r[ix["name"]].split("_")[1])
    DESC[k] = (float(r[ix["volume"]]), float(r[ix["enclosure"]]))

# ---- distances, recomputed from the pocket residue files
atoms = {}
for l in open("loa22_alphafold2_raw.pdb"):
    if l.startswith("ATOM"):
        atoms.setdefault(int(l[22:26]), []).append(
            (float(l[30:38]), float(l[38:46]), float(l[46:54]), l[12:16].strip()))
BB = {"N", "CA", "C", "O", "OXT"}
pc = np.array([a[:3] for r in (121, 142) for a in atoms[r] if a[3] not in BB]).mean(0)
DIST = {}
for f in glob.glob("site_validation/*_P_*_res.pdb"):
    k = int(re.search(r"_P_(\d)_res", f).group(1))
    c = np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54]))
                  for l in open(f) if l.startswith(("ATOM", "HETATM"))]).mean(0)
    DIST[k] = float(np.linalg.norm(c - pc))

panels = [autocrop(Image.open(v).convert("RGB")) for v in VIEWS]
cell = (W - 2*MARGIN - 2*GAP) // 3
fitted, letter_h = [], int(34 * 1.5)
for im in panels:
    s = cell / im.width
    fitted.append(im.resize((cell, int(im.height * s)), Image.LANCZOS))
row_h = max(im.height for im in fitted)

f_letter, f_head, f_body = font(34, True), font(21, True), font(21)
line_h = int(21 * 1.75)
legend_h = line_h * 8
H = MARGIN + letter_h + row_h + GAP + legend_h + MARGIN
canvas = Image.new("RGB", (W, H), "white")
dr = ImageDraw.Draw(canvas)

for i, im in enumerate(fitted):
    x = MARGIN + i * (cell + GAP)
    dr.text((x, MARGIN), "ABC"[i], fill="black", font=f_letter)
    dr.text((x + int(cell*0.13), MARGIN + 6), "rotated %s" % ANGLES[i],
            fill="#4a4a4a", font=f_body)
    canvas.paste(im, (x, MARGIN + letter_h))

y = MARGIN + letter_h + row_h + GAP
cols = [MARGIN + 60, MARGIN + 190, MARGIN + 430, MARGIN + 640, MARGIN + 900]
for c, h in zip(cols, ["Pocket", "Volume (Å³)", "Enclosure",
                       "Centroid to pharmacophore", ""]):
    dr.text((c, y), h, fill="black", font=f_head)
y += line_h
for k in range(6):
    sw = 30
    dr.rectangle([MARGIN, y + 4, MARGIN + sw, y + 4 + sw], fill=HEX[k],
                 outline="#666666")
    note = ""
    if k == 0: note = "highest Drug Score"
    if k == 4: note = "contains Asp121 and Arg142"
    bold = k in (0, 4)
    fb = f_head if bold else f_body
    dr.text((cols[0], y), "P_%d" % k, fill="black", font=fb)
    dr.text((cols[1], y), "%.2f" % DESC[k][0], fill="black", font=fb)
    dr.text((cols[2], y), "%.2f" % DESC[k][1], fill="black", font=fb)
    dr.text((cols[3], y), "%.2f Å" % DIST[k], fill="black", font=fb)
    dr.text((cols[4], y), note, fill="#4a4a4a", font=f_body)
    y += line_h

canvas.save(OUT, dpi=(DPI, DPI))
print("wrote %s  %d x %d px" % (OUT, canvas.width, canvas.height))
for k in range(6):
    print("  P_%d  vol %8.2f  encl %.2f  dist %6.2f A" % (k, DESC[k][0], DESC[k][1], DIST[k]))
