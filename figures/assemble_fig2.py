#!/usr/bin/env python3
"""
Assemble Figure 2 from its four panels.

Layout: A, B, C across the top at equal height, D full width beneath.
Output: Figure2.png, 180 mm wide at 300 dpi (2126 px), white background.

Each PyMOL panel is auto-cropped to its content first. The three renders share
one camera but frame the protein differently, so without cropping the molecule
appears at three different apparent sizes. Cropping to content and matching
height fixes that. Panel D is a chart and is used as-is.

Panel letters are drawn here, so they sit outside the artwork rather than on it.
"""
import os
from PIL import Image, ImageDraw, ImageFont
import matplotlib

DPI = 300
W_MM = 180
W = int(W_MM / 25.4 * DPI)          # 2126 px
GAP = int(0.030 * W)                # gap between panels
MARGIN = int(0.012 * W)
LETTER_PT = 46

PANELS = ["fig2A_t25.png", "fig2B.png", "fig2C.png"]
PANEL_D = "fig2D_chart.png"
OUT = "Figure2.png"

def load_font(size):
    for name in ("DejaVuSans-Bold.ttf", "DejaVuSans.ttf"):
        p = os.path.join(os.path.dirname(matplotlib.__file__),
                         "mpl-data", "fonts", "ttf", name)
        if os.path.exists(p):
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()

def autocrop(im, tol=248, pad=12):
    """Trim near-white margin, keeping a small pad."""
    g = im.convert("L")
    mask = g.point(lambda v: 255 if v < tol else 0)
    bb = mask.getbbox()
    if not bb:
        return im
    l, t, r, b = bb
    l = max(0, l - pad); t = max(0, t - pad)
    r = min(im.width, r + pad); b = min(im.height, b + pad)
    return im.crop((l, t, r, b))

for f in PANELS + [PANEL_D]:
    if not os.path.exists(f):
        raise SystemExit("missing %s - run from the Research root" % f)

# ---- top row: crop, then scale each to a common height
tops = [autocrop(Image.open(f).convert("RGB")) for f in PANELS]
cell_w = (W - 2 * MARGIN - 2 * GAP) // 3
# common height = the tallest panel once each is fitted to cell_w
heights = [int(im.height * cell_w / im.width) for im in tops]
row_h = max(heights)
fitted = []
for im in tops:
    s = min(cell_w / im.width, row_h / im.height)
    fitted.append(im.resize((max(1, int(im.width * s)),
                             max(1, int(im.height * s))), Image.LANCZOS))

# ---- panel D: full content width
d = Image.open(PANEL_D).convert("RGB")
d_w = W - 2 * MARGIN
d = d.resize((d_w, int(d.height * d_w / d.width)), Image.LANCZOS)

letter_h = int(LETTER_PT * 1.5)
H = MARGIN + letter_h + row_h + GAP + letter_h + d.height + MARGIN
canvas = Image.new("RGB", (W, H), "white")
draw = ImageDraw.Draw(canvas)
font = load_font(LETTER_PT)

y_letter = MARGIN
y_row = MARGIN + letter_h
for i, im in enumerate(fitted):
    x_cell = MARGIN + i * (cell_w + GAP)
    draw.text((x_cell, y_letter), "ABC"[i], fill="black", font=font)
    # centre the panel horizontally in its cell, top-align in the row
    canvas.paste(im, (x_cell + (cell_w - im.width) // 2, y_row))

y_dl = y_row + row_h + GAP
draw.text((MARGIN, y_dl), "D", fill="black", font=font)
canvas.paste(d, (MARGIN, y_dl + letter_h))

canvas.save(OUT, dpi=(DPI, DPI))
print("wrote %s  %d x %d px  (%d mm wide at %d dpi)"
      % (OUT, canvas.width, canvas.height, W_MM, DPI))
print("  top row cell %d px wide, row height %d px" % (cell_w, row_h))
for f, im in zip(PANELS, fitted):
    print("    %-18s -> %d x %d" % (f, im.width, im.height))
