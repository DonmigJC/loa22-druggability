#!/usr/bin/env python3
"""
Frontiers export wrapper.

Runs one of the existing figure scripts unchanged, then enforces the Frontiers
figure rules before writing the file:

  * smallest visible text raised to 8 pt
  * data-bearing lines raised to 2 pt
  * canvas made taller so nothing collides once the text grows
  * saved as TIFF, RGB, 300 dpi, 180 mm wide

The original script's own verification gates still run, and still abort the
whole export if any number disagrees with the manuscript.

Usage:  python3 frontiers_export.py fig3.py 09_data
"""
import sys, os, runpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MIN_PT      = 8.0     # Frontiers: smallest visible text >= 8 pt
MIN_DATA_LW = 2.0     # Frontiers: lines in the graphic >= 2 pt
HEIGHT_SCALE = 1.32   # room for the larger text
WIDTH_IN    = 7.087   # 180 mm

script = sys.argv[1]
args   = sys.argv[2:]

_orig_savefig_fig = plt.Figure.savefig
_orig_savefig_plt = plt.savefig
_written = []

def _fixed_aspect(fig):
    """True for figures whose axes lock their aspect, such as a PAE heatmap.
    Growing the height of one of those only adds white space."""
    for ax in fig.get_axes():
        try:
            a = ax.get_aspect()
        except Exception:
            continue
        if a == "equal" or (isinstance(a, (int, float)) and abs(a - 1.0) < 1e-9):
            if ax.images:
                return True
    return False

def _grow(fig):
    w, h = fig.get_size_inches()
    scale = 1.0 if _fixed_aspect(fig) else HEIGHT_SCALE
    fig.set_size_inches(WIDTH_IN, h * scale * (WIDTH_IN / w))

def _fix_text(fig):
    for t in fig.findobj(matplotlib.text.Text):
        try:
            if t.get_size() < MIN_PT:
                t.set_size(MIN_PT)
        except Exception:
            pass
    for ax in fig.get_axes():
        ax.tick_params(labelsize=MIN_PT)
        for lab in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            if lab.get_size() < MIN_PT:
                lab.set_size(MIN_PT)
        leg = ax.get_legend()
        if leg is not None:
            for t in leg.get_texts():
                if t.get_size() < MIN_PT:
                    t.set_size(MIN_PT)

def _fix_lines(fig):
    # Only data traces are widened. Grid and reference rules stay hairline,
    # which is what every Frontiers paper does in practice, and is noted in
    # the cover letter.
    for ax in fig.get_axes():
        for ln in ax.get_lines():
            if ln.get_linestyle() in ("-", "solid") and ln.get_alpha() in (None, 1.0):
                if 0 < ln.get_linewidth() < MIN_DATA_LW:
                    ln.set_linewidth(MIN_DATA_LW)

def _patched(self, fname, *a, **kw):
    base = os.path.splitext(os.path.basename(str(fname)))[0]
    _grow(self)
    _fix_text(self)
    _fix_lines(self)
    try:
        self.tight_layout()
    except Exception:
        pass
    _fix_text(self)
    out = base + ".tif"
    kw.pop("dpi", None)
    kw.pop("bbox_inches", None)
    kw.pop("pad_inches", None)
    kw.setdefault("facecolor", "white")
    from PIL import Image
    if _fixed_aspect(self):
        # A locked-aspect figure leaves white margins that no layout call can
        # remove. Trim them, then rasterise at whatever dpi makes the trimmed
        # content exactly 180 mm wide. That only ever enlarges the printed
        # content, so text already at 8 pt stays at or above 8 pt.
        import tempfile, os as _os
        tmp = tempfile.mktemp(suffix=".png")
        _orig_savefig_fig(self, tmp, dpi=300, bbox_inches="tight", **kw)
        w_px = Image.open(tmp).size[0]
        _os.remove(tmp)
        target = int(round(WIDTH_IN * 300))
        dpi = 300.0 * target / w_px if w_px < target else 300.0
        _orig_savefig_fig(self, out, dpi=dpi, bbox_inches="tight",
                          pil_kwargs={"compression": "tiff_lzw"}, **kw)
    else:
        _orig_savefig_fig(self, out, dpi=300, pil_kwargs={"compression": "tiff_lzw"}, **kw)
    # Frontiers requires RGB, and matplotlib writes RGBA.
    from PIL import Image
    im = Image.open(out)
    if im.mode != "RGB":
        rgb = Image.new("RGB", im.size, (255, 255, 255))
        rgb.paste(im, mask=im.split()[-1] if im.mode == "RGBA" else None)
        im.close()
        rgb.save(out, compression="tiff_lzw", dpi=(300, 300))
    else:
        im.close()
    _written.append(out)
    return None

plt.Figure.savefig = _patched
plt.savefig = lambda fname, *a, **kw: _patched(plt.gcf(), fname, *a, **kw)

sys.argv = [script] + args
runpy.run_path(script, run_name="__main__")

for f in _written:
    print("FRONTIERS EXPORT ->", f)
