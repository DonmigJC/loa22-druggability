#!/usr/bin/env python3
"""
Assemble Figure 6: panels A and B (existing wide image) on top,
panels C and D side by side beneath. Panel letters A and B are already
baked into Figure6AB.png, so only C and D are drawn here.
Output: Figure6.png, 180 mm at 300 dpi.
"""
import os
from PIL import Image, ImageDraw, ImageFont
import matplotlib

DPI, W_MM = 300, 180
W = int(W_MM / 25.4 * DPI)
MARGIN, GAP = int(0.012 * W), int(0.026 * W)
AB, C, D, OUT = "Figure6AB.png", "Figure6C.png", "Figure6D.png", "Figure6.png"

def font(sz):
    p = os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data",
                     "fonts", "ttf", "DejaVuSans-Bold.ttf")
    return ImageFont.truetype(p, sz) if os.path.exists(p) else ImageFont.load_default()

def autocrop(im, tol=248, pad=10):
    m = im.convert("L").point(lambda v: 255 if v < tol else 0)
    bb = m.getbbox()
    if not bb: return im
    l, t, r, b = bb
    return im.crop((max(0, l-pad), max(0, t-pad),
                    min(im.width, r+pad), min(im.height, b+pad)))

for f in (AB, C, D):
    if not os.path.exists(f): raise SystemExit("missing %s" % f)

ab = Image.open(AB).convert("RGB")
ab = ab.resize((W - 2*MARGIN, int(ab.height * (W - 2*MARGIN) / ab.width)), Image.LANCZOS)

cell = (W - 2*MARGIN - GAP) // 2
c = Image.open(C).convert("RGB")
c = c.resize((cell, int(c.height * cell / c.width)), Image.LANCZOS)
d = autocrop(Image.open(D).convert("RGB"))
d = d.resize((cell, int(d.height * cell / d.width)), Image.LANCZOS)
row2_h = max(c.height, d.height)

lh = 70
H = MARGIN + ab.height + GAP + lh + row2_h + MARGIN
canvas = Image.new("RGB", (W, H), "white")
dr = ImageDraw.Draw(canvas)
f_lbl = font(46)

canvas.paste(ab, (MARGIN, MARGIN))
y = MARGIN + ab.height + GAP
dr.text((MARGIN, y), "C", fill="black", font=f_lbl)
dr.text((MARGIN + cell + GAP, y), "D", fill="black", font=f_lbl)
canvas.paste(c, (MARGIN, y + lh + (row2_h - c.height) // 2))
canvas.paste(d, (MARGIN + cell + GAP, y + lh + (row2_h - d.height) // 2))

canvas.save(OUT, dpi=(DPI, DPI))
print("wrote %s  %d x %d px" % (OUT, canvas.width, canvas.height))
