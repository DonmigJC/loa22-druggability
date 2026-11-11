#!/usr/bin/env python3
"""Recompute Tables 1, 2, 3 and the section 2.5 decoy claims from the
deposited CSV files. Every value the paper states that is computable from
these files is checked."""
import csv
import math
import os

import numpy as np

R = "/home/donmigcage/Research"


def load(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh))


def f(row, key):
    v = row.get(key, "")
    if v in ("", "NA", "nan", None):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def chk(label, got, claim, tol, unit=""):
    if got is None:
        print("  %-62s NOT COMPUTABLE" % label)
        return
    ok = abs(got - claim) <= tol
    print("  %-62s got %10.4f%s  paper %9.3f  %s"
          % (label, got, unit, claim, "OK" if ok else "MISMATCH"))


def chk_int(label, got, claim):
    print("  %-62s got %10d    paper %9d  %s"
          % (label, got, claim, "OK" if got == claim else "MISMATCH"))


print("=" * 100)
print("TABLE 2 and section 2.3: dual-site screen")
print("=" * 100)
ds = load(os.path.join(R, "dual_site", "dual_site_results.csv"))
print("  rows in dual_site_results.csv: %d" % len(ds))
both = [r for r in ds if f(r, "aff_p0") is not None and f(r, "aff_site") is not None]
chk_int("molecules docked at both volumes", len(both), 6253)

p0 = np.array([f(r, "aff_p0") for r in both])
site = np.array([f(r, "aff_site") for r in both])
chk("mean affinity P0 (kcal/mol)", float(p0.mean()), -6.699, 0.002)
chk("mean affinity functional site (kcal/mol)", float(site.mean()), -5.155, 0.002)
chk("mean penalty at functional site (kcal/mol)", float((site - p0).mean()), 1.544, 0.002)
chk("best affinity P0 (kcal/mol)", float(p0.min()), -10.670, 0.002)
chk("best affinity functional site (kcal/mol)", float(site.min()), -8.567, 0.002)
pref = int((site < p0).sum())
chk_int("molecules preferring the functional site", pref, 22)
print("  %-62s got %9.4f%%   paper      0.35%%" % ("as a percentage", 100.0 * pref / len(both)))

br = np.array([f(r, "bridge") for r in both if f(r, "bridge") is not None])
chk_int("bridging score > 0.5", int((br > 0.5).sum()), 1409)
chk_int("bridging score > 0.8", int((br > 0.8).sum()), 39)
print("  %-62s got %9.4f%%   paper     22.53%%" % ("fraction > 0.5", 100.0 * (br > 0.5).mean()))
print("  %-62s got %9.4f%%   paper      0.62%%" % ("fraction > 0.8", 100.0 * (br > 0.8).mean()))

dk = np.array([f(r, "d_key") for r in both if f(r, "d_key") is not None])
chk("closest approach anywhere in the library (A)", float(dk.min()), 2.92, 0.02)
chk("median pharmacophore distance (A)", float(np.median(dk)), 7.33, 0.02)

# thermodynamics of the penalty
dG = float((site - p0).mean())
for T, claim in ((298.0, 13.0), (310.0, 12.0)):
    fold = math.exp(dG / (0.0019872 * T))
    print("  %-62s got %10.2f    paper about %d-fold"
          % ("Kd ratio at %.0f K from %.3f kcal/mol" % (T, dG), fold, claim))

print()
print("=" * 100)
print("TABLE 3 and section 2.3: residue-resolved contact")
print("=" * 100)
rk = load(os.path.join(R, "dual_site", "dual_site_ranked.csv"))
print("  rows in dual_site_ranked.csv: %d" % len(rk))
print("  NOTE: Table 3 reports n = 6,222 with separate Asp121 and Arg142")
print("  distances. Neither CSV here carries per-residue distances, only the")
print("  combined d_key, so Table 3 cannot be recomputed from the deposited")
print("  files. Flag for the data-availability statement.")

print()
print("=" * 100)
print("SECTION 2.5: decoy and random controls")
print("=" * 100)
dc = load(os.path.join(R, "decoy_validation", "decoy_results.csv"))
kinds = {}
for r in dc:
    kinds.setdefault(r.get("kind", "?"), []).append(r)
print("  rows by kind: " + ", ".join("%s=%d" % (k, len(v)) for k, v in sorted(kinds.items())))

LEADS = [3.43, 3.33, 3.12]
WEAKEST = 3.43
for kind, claim_mean in (("decoy", 6.92), ("random", 8.81)):
    rows = kinds.get(kind, [])
    if not rows:
        continue
    d = np.array([f(r, "d_key") for r in rows if f(r, "d_key") is not None])
    chk("mean closest approach, %s (A)" % kind, float(d.mean()), claim_mean, 0.02)
    chk_int("n %s" % kind, len(d), 45 if kind == "decoy" else 15)

dec = np.array([f(r, "d_key") for r in kinds.get("decoy", []) if f(r, "d_key") is not None])
n_within = int((dec < WEAKEST).sum())
chk_int("decoys within 3.43 A of the pharmacophore", n_within, 0)
print("  %-62s got %10.3f    paper      3.29" % ("mean of the three leads (A)", float(np.mean(LEADS))))

# 42.2 percent of decoys bound the functional site better than the weakest lead
WEAKEST_AFF_SITE = -5.625   # LOA22-B3, the weakest of the three on site affinity
dsite = np.array([f(r, "aff_site") for r in kinds.get("decoy", []) if f(r, "aff_site") is not None])
n_better = int((dsite < WEAKEST_AFF_SITE).sum())
print("  %-62s got %2d of %d = %.1f%%   paper 42.2%%"
      % ("decoys beating the weakest lead on site affinity", n_better, len(dsite),
         100.0 * n_better / len(dsite)))

try:
    from scipy.stats import mannwhitneyu
    u, p = mannwhitneyu(LEADS, dec, alternative="less")
    print("  %-62s got U=%.1f p=%.5f  paper U=0 p=0.0001" % ("Mann-Whitney, leads vs decoys", u, p))
except ImportError:
    print("  scipy not available, Mann-Whitney not recomputed")

print()
print("=" * 100)
print("TABLE 1: prior-only control")
print("=" * 100)
pr = load(os.path.join(R, "prior_control", "prior_docked.csv"))
print("  rows in prior_docked.csv: %d" % len(pr))
aff = np.array([f(r, "affinity") for r in pr if f(r, "affinity") is not None])
qed = np.array([f(r, "QED") for r in pr if f(r, "QED") is not None])
chk_int("prior molecules with an affinity", len(aff), 6049)
chk("prior mean affinity (kcal/mol)", float(aff.mean()), -6.556, 0.002)
chk("prior median affinity (kcal/mol)", float(np.median(aff)), -6.632, 0.002)
chk("prior best single molecule (kcal/mol)", float(aff.min()), -12.160, 0.002)
chk("prior 25th percentile (kcal/mol)", float(np.percentile(aff, 25)), -7.318, 0.01)
chk("prior 10th percentile (kcal/mol)", float(np.percentile(aff, 10)), -7.909, 0.01)
chk("prior mean QED", float(qed.mean()), 0.618, 0.001)
DOXY = -6.862
print("  %-62s got %9.3f%%   paper    41.08%%"
      % ("prior beating doxycycline", 100.0 * (aff <= DOXY).mean()))
print("  %-62s got %9.3f%%   paper    19.53%%"
      % ("prior passing the -7.5 filter", 100.0 * (aff <= -7.5).mean()))
print("  %-62s got %9.3f%%   paper     3.42%%"
      % ("prior <= -8.5", 100.0 * (aff <= -8.5).mean()))
print("  %-62s got %9.3f%%   paper     1.16%%"
      % ("prior <= -9.0", 100.0 * (aff <= -9.0).mean()))
