#!/usr/bin/env python3
"""
Rebuild Supplementary Table S12 from the v2 trajectories.

Table S12 as submitted reports the v1 simulations in every row. This script
reads the v2 analysis output written by analyse_conformational_stability_v2.sh
and reports, per system, the mean of each replicate and the mean and standard
deviation across replicates.

Two statistics are printed and they answer different questions:

  within-trajectory    mean and SD of the RMSD time series inside one run.
                       This is what Table S12 reports now. Its SD describes
                       fluctuation over time, not reproducibility.
  across-replicate     mean and SD of the three replicate means. Its SD
                       describes how far independent trajectories disagree,
                       which is what the three-replica requirement is for.

The across-replicate figure is the one to put in the table, with the number of
replicates stated.

The backbone column is reported for residues 78 to 186, the structured OmpA
domain, not for all 195 residues. The default whole-protein selection is
dominated by the two disordered termini, which have mean pLDDT 46.8 at the N
terminus and fall to 31.8 at the C terminus, and on the trajectories Table S12
reports it runs 14.60 to 20.67 A against 1.18 to 2.80 A for the domain. It also
inverts the ranking: LOA22-B1 has the largest whole-protein value and the
smallest domain value. Both are printed, domain first.

Self-check: if the v1 files are still present the script recomputes their
statistics and compares them against the four value pairs printed in the
submitted Table S12. That confirms the statistics code here matches whatever
produced the published numbers before any v2 figure is trusted.

Usage:  python3 summarize_conformational_stability.py [research_root]
"""
import os
import sys
import math

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/home/donmigcage/Research"

APO = os.path.join(ROOT, "06_MD_apo_v2")
MD = os.path.join(ROOT, "05_MD_v2")
V1 = os.path.join(ROOT, "04_Molecular_Dynamics")

COMPOUNDS = ["LOA22-B1", "LOA22-B2", "LOA22-B3"]
REPS = ["rep1", "rep2", "rep3"]

# Table S12 exactly as submitted, for the self-check.
# label -> (backbone mean, backbone sd, bsite mean, bsite sd, v1 analysis dir)
SUBMITTED = [
    ("Apo (ligand-free)", 17.62, 4.50, 1.21, 0.26,
     os.path.join(V1, "Apo", "md_run", "analysis")),
    ("LOA22-B1", 20.67, 6.52, 1.09, 0.20,
     os.path.join(V1, "Lead_1", "md_run", "analysis_v2")),
    ("LOA22-B2", 14.60, 2.43, 1.32, 0.27,
     os.path.join(V1, "Lead_4", "md_run", "analysis_v2")),
    ("LOA22-B3", 16.22, 5.16, 2.43, 0.59,
     os.path.join(V1, "Lead_5", "md_run", "analysis_v2")),
]

# The same v1 trajectories measured over the structured domain, residues 78 to
# 186, on 5 November 2026. Printed beside the submitted values to show how much
# of the submitted backbone column is the disordered termini.
V1_OMPA = {
    "Apo (ligand-free)": (1.18, 0.22),
    "LOA22-B1": (1.70, 0.22),
    "LOA22-B2": (2.80, 0.58),
    "LOA22-B3": (2.68, 0.57),
}


def read_xvg(path):
    """Return the second column of an .xvg, converted from nm to Angstrom."""
    vals = []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line or line[0] in "@#&":
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        vals.append(float(parts[1]) * 10.0)
                    except ValueError:
                        continue
    except OSError:
        return None
    return vals or None


def mean_sd(vals):
    n = len(vals)
    m = sum(vals) / n
    if n < 2:
        return m, 0.0
    var = sum((v - m) ** 2 for v in vals) / n
    return m, math.sqrt(var)


def sample_sd(vals):
    """Sample SD, n-1, for the across-replicate figure."""
    n = len(vals)
    if n < 2:
        return 0.0
    m = sum(vals) / n
    return math.sqrt(sum((v - m) ** 2 for v in vals) / (n - 1))


def series(d, name):
    return read_xvg(os.path.join(d, "analysis", name))


# ------------------------------------------------------------------ self-check
print("=" * 74)
print("SELF-CHECK: recomputing the submitted Table S12 from the v1 files")
print("=" * 74)
checked = failed = 0
for label, bb_m, bb_s, bs_m, bs_s, d in SUBMITTED:
    bb = read_xvg(os.path.join(d, "rmsd_backbone.xvg"))
    bs = read_xvg(os.path.join(d, "rmsd_bsite.xvg"))
    if bb is None or bs is None:
        print("  %-20s v1 files not present, cannot check" % label)
        continue
    gm, gs = mean_sd(bb)
    sm, ss = mean_sd(bs)
    ok = (abs(gm - bb_m) < 0.015 and abs(gs - bb_s) < 0.015
          and abs(sm - bs_m) < 0.015 and abs(ss - bs_s) < 0.015)
    checked += 1
    if not ok:
        failed += 1
    print("  %-20s n=%-6d backbone %.2f+/-%.2f (table %.2f+/-%.2f)  "
          "site %.2f+/-%.2f (table %.2f+/-%.2f)  %s"
          % (label, len(bs), gm, gs, bb_m, bb_s, sm, ss, bs_m, bs_s,
             "ok" if ok else "MISMATCH"))
if checked and not failed:
    print("\n  All %d rows reproduce. The statistics code here matches the one that"
          % checked)
    print("  produced the submitted table, and every submitted row is v1 data.")
elif failed:
    print("\n  %d row(s) did not reproduce. Resolve that before trusting anything"
          % failed)
    print("  below, because the two pipelines then differ.")

print()
print("  How much of the submitted backbone column is the disordered termini:")
print("    %-20s %-22s %-22s" % ("system", "whole 1-195 (as submitted)",
                                 "OmpA 78-186 (measured)"))
for label, bb_m, bb_s, _, _, _ in SUBMITTED:
    om = V1_OMPA.get(label)
    if om:
        print("    %-20s %6.2f +/- %-14.2f %6.2f +/- %-14.2f"
              % (label, bb_m, bb_s, om[0], om[1]))
print("    The whole-protein column also inverts the ranking: LOA22-B1 is the")
print("    largest of the three on it and the smallest on the domain.")

# ------------------------------------------------------------------ v2 table
systems = [("Apo (ligand-free)", [os.path.join(APO, r) for r in REPS])]
for c in COMPOUNDS:
    systems.append((c, [os.path.join(MD, c, r) for r in REPS]))

print()
print("=" * 74)
print("v2 TRAJECTORIES, PER REPLICATE")
print("=" * 74)
print("  %-20s %-6s %-7s %-21s %-21s %-21s"
      % ("system", "rep", "frames", "OmpA backbone 78-186",
         "binding-site 115-150", "whole backbone 1-195"))

table = {}
for label, dirs in systems:
    rows = []
    for d, rep in zip(dirs, REPS):
        od = series(d, "rmsd_ompa_bb.xvg")
        bs = series(d, "rmsd_bsite.xvg")
        bb = series(d, "rmsd_backbone.xvg")
        if od is None or bs is None:
            print("  %-20s %-6s %s" % (label, rep, "not analysed yet"))
            continue
        om, os_ = mean_sd(od)
        sm, ss = mean_sd(bs)
        if bb is None:
            gm, gs = float("nan"), float("nan")
        else:
            gm, gs = mean_sd(bb)
        rows.append((rep, len(bs), om, os_, sm, ss, gm, gs))
        print("  %-20s %-6s %-7d %7.2f +/- %-9.2f %7.2f +/- %-9.2f %7.2f +/- %-9.2f"
              % (label, rep, len(bs), om, os_, sm, ss, gm, gs))
    table[label] = rows

print()
print("=" * 74)
print("REBUILT TABLE S12, across-replicate mean and SD")
print("=" * 74)
print("  %-20s %-5s %-24s %-24s" % ("System", "n",
                                    "OmpA backbone RMSD (A)",
                                    "Binding-site RMSD (A)"))
ready = True
for label, _ in systems:
    rows = table.get(label, [])
    if len(rows) < 3:
        print("  %-20s %-5s %s" % (label, len(rows),
                                   "incomplete, %d of 3 replicates" % len(rows)))
        ready = False
        continue
    obm = [r[2] for r in rows]
    bsm = [r[4] for r in rows]
    print("  %-20s %-5d %9.2f +/- %-12.2f %9.2f +/- %-12.2f"
          % (label, len(rows), sum(obm) / len(obm), sample_sd(obm),
             sum(bsm) / len(bsm), sample_sd(bsm)))

if ready:
    print()
    print("  Whole-protein backbone, for the record, not for the table:")
    for label, _ in systems:
        rows = table[label]
        gm = [r[6] for r in rows]
        print("    %-20s %9.2f +/- %-12.2f"
              % (label, sum(gm) / len(gm), sample_sd(gm)))

if ready:
    bs_all = {lab: [r[4] for r in table[lab]] for lab, _ in systems}
    apo = sum(bs_all["Apo (ligand-free)"]) / 3.0
    holo = [sum(bs_all[c]) / 3.0 for c in COMPOUNDS]
    print()
    print("  For the sentence in section 2.7:")
    print("    apo binding-site RMSD        %.2f A" % apo)
    print("    ligand-bound range           %.2f to %.2f A" % (min(holo), max(holo)))
    print("    apo inside that range        %s"
          % ("yes" if min(holo) <= apo <= max(holo) else "NO, the claim changes"))
    print()
    print("  Every figure above comes from v2 trajectories only, three replicates")
    print("  per system, so it is reproducible from Methods section 4.6 as written.")
else:
    print()
    print("  Run analyse_conformational_stability_v2.sh once the simulations finish.")
