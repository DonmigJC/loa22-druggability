#!/usr/bin/env python3
"""
BRIDGING OCCUPANCY BY REPLICATE
===============================

Measures, for one replicate, how often each lead maintains simultaneous
contact with both the pharmacophore (Asp121/Arg142) and the druggable
pocket P_0.

Method: plain coordinate geometry on PBC-corrected static frames sampled
every 5 ns.  No trajectory tools, no structural fitting, no minimum-image
convention - so there is no ambiguity in the numbers.

Handles the PDB column overflow that occurs past 99,999 atoms (these
systems have ~117,000), which broke the earlier embedded version.

Usage:
    python bridging_by_rep.py 1
    python bridging_by_rep.py 2
    python bridging_by_rep.py 3
"""
import os
import sys
import math
import subprocess

REP = sys.argv[1] if len(sys.argv) > 1 else "1"

R = os.path.expanduser("~/Research")
MD = os.path.join(R, "05_MD_v2")

KEY = {121: {"CB", "CG", "OD1", "OD2"},
       142: {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}}
P0 = {61, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 88, 89, 92, 93,
      96, 97, 115, 116, 117, 119, 163, 165, 169, 171, 172, 183, 185, 186, 187}
SKIP = {"SOL", "HOH", "WAT", "NA", "CL", "TIP3", "K"}

DOCKED = {"LOA22-B1": (3.43, 3.04),
          "LOA22-B2": (3.33, 3.13),
          "LOA22-B3": (3.12, 2.88)}

CONTACT = 4.5   # angstroms; beyond this there is no direct interaction


def parse_frame(path):
    """Read a PDB frame, tolerating column overflow past 99,999 atoms."""
    key, p0, lig = [], [], []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        rs = l[17:20].strip()
        if rs in SKIP:
            continue
        nm = l[12:16].strip()
        if nm.startswith("H") or l[76:78].strip() == "H":
            continue
        try:
            rn = int(l[22:26])
            c = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        except ValueError:
            continue          # column overflow - skip this line safely
        if rs in ("UNL", "LIG"):
            lig.append(c)
            continue
        if rn in KEY and nm in KEY[rn]:
            key.append(c)
        if rn in P0:
            p0.append(c)
    return key, p0, lig


print("=" * 74)
print(f"BRIDGING OCCUPANCY — REPLICATE {REP}")
print("=" * 74)
print("Contact threshold 4.5 A.  van der Waals contact is 3.5-4.0 A.")
print("BRIDGING = simultaneous contact with BOTH sites.\n")

summary = {}

for lead in ("LOA22-B1", "LOA22-B2", "LOA22-B3"):
    run = os.path.join(MD, lead, f"rep{REP}")
    out = os.path.join(run, "analysis")
    fr = os.path.join(out, "frames5")
    ndx = os.path.join(out, "an.ndx")
    ctr = os.path.join(out, "center.xtc")
    tpr = os.path.join(run, "md_prod.tpr")

    print(f"=== {lead}  replicate {REP} ===")
    if not os.path.exists(ctr):
        print(f"  PBC-corrected trajectory missing: {ctr}")
        print("  Run:  bash ~/Research/analyse_md_v2.sh " + REP + "\n")
        continue
    os.makedirs(fr, exist_ok=True)

    dk0, dp0 = DOCKED[lead]
    print(f"  {'time':>8}{'to key':>10}{'to P_0':>10}   verdict")
    print(f"  {'docked':>8}{dk0:>10.2f}{dp0:>10.2f}   reference")

    n_bridge = n_key = n_p0 = n_frames = 0
    dk_vals, dp_vals = [], []

    for t in range(0, 100001, 5000):
        f = os.path.join(fr, f"{t}.pdb")
        if not os.path.exists(f):
            subprocess.run(
                ["gmx", "trjconv", "-s", tpr, "-f", ctr, "-n", ndx,
                 "-o", f, "-dump", str(t)],
                input="0\n", text=True, capture_output=True, cwd=run)
        if not os.path.exists(f):
            continue
        key, p0, lig = parse_frame(f)
        if not lig or not key:
            print(f"  {t//1000:>6} ns   (parse failed)")
            continue
        dk = min(math.dist(a, b) for a in lig for b in key)
        dp = min(math.dist(a, b) for a in lig for b in p0) if p0 else 99.0
        dk_vals.append(dk)
        dp_vals.append(dp)
        n_frames += 1
        if dk < CONTACT:
            n_key += 1
        if dp < CONTACT:
            n_p0 += 1
        if dk < CONTACT and dp < CONTACT:
            v = "BRIDGING"
            n_bridge += 1
        elif dk < CONTACT:
            v = "key only"
        elif dp < CONTACT:
            v = "P_0 only"
        else:
            v = "neither"
        print(f"  {t//1000:>6} ns{dk:>10.2f}{dp:>10.2f}   {v}")

    if n_frames:
        mean_dk = sum(dk_vals) / len(dk_vals)
        mean_dp = sum(dp_vals) / len(dp_vals)
        sd_dk = math.sqrt(sum((x - mean_dk) ** 2 for x in dk_vals) / len(dk_vals))
        print(f"\n  bridging      : {n_bridge}/{n_frames} frames "
              f"({100*n_bridge/n_frames:.0f}%)")
        print(f"  key contact   : {n_key}/{n_frames}")
        print(f"  P_0 contact   : {n_p0}/{n_frames}")
        print(f"  mean d(key)   : {mean_dk:.2f} +/- {sd_dk:.2f} A")
        print(f"  mean d(P_0)   : {mean_dp:.2f} A")
        summary[lead] = dict(bridge=n_bridge, key=n_key, p0=n_p0,
                             n=n_frames, dk=mean_dk, sd=sd_dk, dp=mean_dp)
    print()

if summary:
    print("=" * 74)
    print(f"SUMMARY — REPLICATE {REP}")
    print("=" * 74)
    print(f"{'lead':<12}{'bridging':>12}{'key':>8}{'P_0':>8}"
          f"{'mean d(key)':>14}{'mean d(P_0)':>13}")
    print("-" * 68)
    for lead, s in summary.items():
        print(f"{lead:<12}{s['bridge']}/{s['n']:<10}{s['key']:>8}{s['p0']:>8}"
              f"{s['dk']:>10.2f} +/-{s['sd']:>4.2f}{s['dp']:>13.2f}")
    print()
    print("Replicate 1 for comparison:")
    print("  LOA22-B1    11/21        21      21      4.71          2.26")
    print("  LOA22-B2    21/21        21      21      3.36          3.05")
    print("  LOA22-B3    18/21        21      21      4.04          2.91")
print("=" * 74)
