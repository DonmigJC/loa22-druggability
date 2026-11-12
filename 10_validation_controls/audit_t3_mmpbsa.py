#!/usr/bin/env python3
"""Recompute Table 3 from pharmacophore_split/armB_split.csv, and the
section 2.6 free energies from the gmx_MMPBSA output files.

Two questions about the free energies:
  1. The paper reports every across-replicate deviation as a sample standard
     deviation (ddof=1), which section 4.6 states. Earlier drafts had the
     end-state totals on the population estimator while the per-residue term
     next to them used the sample one, so this check exists to keep the two
     from drifting apart again. Section 4.6 promises "mean +/- standard
     deviation with the unit
     of replication stated".
  2. LOA22-B1 and LOA22-B3 have only one replicate each with MMPBSA output, so
     their "+/-" cannot be a cross-replicate SD. What is it?
"""
import csv
import os
import re

import numpy as np

R = "/home/donmigcage/Research"


def load(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


print("=" * 98)
print("TABLE 3 recomputed from pharmacophore_split/armB_split.csv")
print("=" * 98)
rows = load(os.path.join(R, "pharmacophore_split", "armB_split.csv"))
print("  rows: %d" % len(rows))
good = [r for r in rows
        if num(r["d_asp"]) is not None and num(r["d_arg"]) is not None
        and num(r["d_p0"]) is not None]
print("  rows with all three distances: %d      paper Table 3 says n = 6,222  %s"
      % (len(good), "OK" if len(good) == 6222 else "MISMATCH"))

asp = np.array([num(r["d_asp"]) for r in good])
arg = np.array([num(r["d_arg"]) for r in good])
p0 = np.array([num(r["d_p0"]) for r in good])
aff = np.array([num(r["aff"]) if num(r["aff"]) is not None else np.nan for r in good])


def chk(label, got, claim, tol):
    print("  %-56s got %9.4f   paper %8.2f   %s"
          % (label, got, claim, "OK" if abs(got - claim) <= tol else "MISMATCH"))


def chk_i(label, got, claim):
    print("  %-56s got %9d   paper %8d   %s"
          % (label, got, claim, "OK" if got == claim else "MISMATCH"))


chk("mean distance to Asp121 (A)", float(asp.mean()), 8.96, 0.02)
chk("mean distance to Arg142 (A)", float(arg.mean()), 9.49, 0.02)
chk("minimum distance to Asp121 (A)", float(asp.min()), 3.04, 0.02)
chk("minimum distance to Arg142 (A)", float(arg.min()), 2.92, 0.02)
n_asp = int((asp <= 4.5).sum())
n_arg = int((arg <= 4.5).sum())
chk_i("molecules within 4.5 A of Asp121", n_asp, 1211)
print("  %-56s got %9.4f%%  paper   19.46%%" % ("as a percentage", 100.0 * n_asp / len(good)))
chk_i("molecules within 4.5 A of Arg142", n_arg, 1)
both = int(((asp <= 4.5) & (arg <= 4.5)).sum())
chk_i("molecules contacting both residues", both, 1)
tri = int(((asp <= 4.5) & (arg <= 4.5) & (p0 <= 4.5)).sum())
chk_i("contacting both residues and P0", tri, 0)

idx = np.where((asp <= 4.5) & (arg <= 4.5))[0]
if idx.size:
    i = int(idx[0])
    print("  the single molecule reaching both: aff %.3f, SMILES %s"
          % (aff[i], good[i]["SMILES"]))
    print("    paper calls it a three-heavy-atom fragment of affinity -2.563")

print()
print("=" * 98)
print("SECTION 2.6 free energies, LOA22-B2 alternative pose, three replicates")
print("=" * 98)


def parse_results(path):
    """Return {'GB': delta_total, 'PB': delta_total} from FINAL_RESULTS.dat."""
    out = {}
    model = None
    for line in open(path, encoding="utf-8", errors="replace"):
        s = line.strip()
        if re.search(r"GENERALIZED BORN", s, re.I):
            model = "GB"
        elif re.search(r"POISSON BOLTZMANN", s, re.I):
            model = "PB"
        m = re.match(r"^[Δ\?]?TOTAL\s+(-?\d+\.\d+)", s)
        if m and model and model not in out:
            pass
        if s.startswith("ΔTOTAL") or s.startswith("DELTA TOTAL"):
            v = re.findall(r"-?\d+\.\d+", s)
            if v and model:
                out.setdefault(model, float(v[0]))
    return out


gb, pb = [], []
for rep in ("rep1_arg142pose", "rep2_arg142pose", "rep3_arg142pose"):
    p = os.path.join(R, "05_MD_v2", "LOA22-B2", rep, "mmpbsa", "FINAL_RESULTS.dat")
    if not os.path.isfile(p):
        print("  %s MISSING" % rep)
        continue
    d = parse_results(p)
    print("  %-18s GB %8s   PB %8s"
          % (rep, d.get("GB", "?"), d.get("PB", "?")))
    if "GB" in d:
        gb.append(d["GB"])
    if "PB" in d:
        pb.append(d["PB"])

for name, vals, claim_m, claim_s in (("GB", gb, -27.63, 4.82), ("PB", pb, -17.73, 1.39)):
    if len(vals) < 2:
        print("  %s: only %d value(s), cannot compute an SD" % (name, len(vals)))
        continue
    a = np.array(vals)
    print("  %s mean %.3f   population SD (ddof=0) %.3f   sample SD (ddof=1) %.3f"
          % (name, a.mean(), a.std(ddof=0), a.std(ddof=1)))
    print("      paper states %.2f +/- %.2f  ->  matches %s"
          % (claim_m, claim_s,
             "ddof=0 (population)" if abs(a.std(ddof=0) - claim_s) < 0.05
             else ("ddof=1 (sample)" if abs(a.std(ddof=1) - claim_s) < 0.05
                   else "NEITHER")))

print()
print("=" * 98)
print("Arg142 per-residue van der Waals term, and what the +/- means")
print("=" * 98)
print("  MMPBSA output present for:")
for c in ("LOA22-B1", "LOA22-B2", "LOA22-B3"):
    reps = [r for r in ("rep1_arg142pose", "rep2_arg142pose", "rep3_arg142pose")
            if os.path.isfile(os.path.join(R, "05_MD_v2", c, r, "mmpbsa",
                                           "FINAL_DECOMP.dat"))]
    print("    %-10s %d replicate(s): %s" % (c, len(reps), ", ".join(reps) or "none"))
print()
print("  The paper reports Arg142 vdW as -2.32 +/- 0.69 (B2), -1.15 +/- 0.22 (B3)")
print("  and -0.71 +/- 0.20 (B1). B1 and B3 have one replicate each, so their")
print("  +/- cannot be a cross-replicate SD. Check what it is before the")
print("  revision restates them, because section 4.6 promises the unit of")
print("  replication is stated.")
