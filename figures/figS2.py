import json, numpy as np, glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
f = glob.glob("*predicted_aligned_error*.json") + glob.glob("*pae*.json")
if not f: raise SystemExit("PAE JSON not found — download from AlphaFold DB")
d = json.load(open(f[0]))
pae = np.array(d[0]["predicted_aligned_error"] if isinstance(d, list)
               else d["predicted_aligned_error"])
plt.figure(figsize=(6,5))
plt.imshow(pae, cmap="Greens_r", vmin=0, vmax=30)
plt.colorbar(label="Expected position error (Å)")
plt.xlabel("Scored residue"); plt.ylabel("Aligned residue")
plt.axhline(47, color="r", lw=.8); plt.axvline(47, color="r", lw=.8)
plt.tight_layout(); plt.savefig("FigureS2.png", dpi=300)
