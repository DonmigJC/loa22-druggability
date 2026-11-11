#!/usr/bin/env python3
"""
Summarize the docked-pose replicates: when each ligand left the functional site,
reported per replicate and as a mean across replicates.

Reads the gmx mindist output written by md_replicate_analysis.sh:
    <MDROOT>/<compound>/<rep>/analysis/md_KeyRes.xvg

Departure is defined exactly as the manuscript defines contact, by the 4.5 A
heavy-atom threshold (0.45 nm). Two departure times are reported, because they
answer different questions:

  first crossing  the first frame above 0.45 nm, which can be a brief excursion
  final departure the last crossing after which the ligand never returns below
                  0.45 nm for the remainder of the trajectory

The manuscript's statement that each compound "left the functional site" refers
to the final departure. Report that one, and keep the first crossing in the
supplement so a reviewer can see the difference.

Usage:
    python3 summarize_docked_replicates.py /home/donmigcage/Research/05_MD_v2
    python3 summarize_docked_replicates.py /path/to/05_MD_v2 --threshold 0.45
"""
import argparse
import glob
import os
import sys

import numpy as np

COMPOUNDS = ["LOA22-B1", "LOA22-B2", "LOA22-B3"]
REPLICATES = ["rep1", "rep2", "rep3"]


def read_xvg(path):
    """Return (time_ns, distance_nm) from a GROMACS .xvg file."""
    t, d = [], []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line[:1] in ("#", "@"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    t.append(float(parts[0]))
                    d.append(float(parts[1]))
                except ValueError:
                    continue
    return np.asarray(t), np.asarray(d)


def departures(t, d, threshold):
    """First crossing and final departure, in the time units of t."""
    above = d > threshold
    if not above.any():
        return None, None
    first = float(t[np.argmax(above)])
    # final departure: start of the last contiguous run of "above"
    below_idx = np.flatnonzero(~above)
    if below_idx.size == 0:
        final = float(t[0])
    elif below_idx[-1] == len(d) - 1:
        final = None          # ligand is in contact at the final frame
    else:
        final = float(t[below_idx[-1] + 1])
    return first, final


def fmt(x, unit="ns"):
    if x is None:
        return "no departure"
    if x < 1.0:
        return f"{x * 1000:.0f} ps"
    return f"{x:.1f} {unit}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mdroot")
    ap.add_argument("--threshold", type=float, default=0.45,
                    help="contact threshold in nm (default 0.45, i.e. 4.5 A)")
    ap.add_argument("--metric", default="md_KeyRes",
                    help="xvg basename to read (md_KeyRes, md_Asp121only, md_Arg142only)")
    args = ap.parse_args()

    print(f"Docked-pose replicates, departure from the functional site")
    print(f"contact threshold {args.threshold} nm, metric {args.metric}")
    print("=" * 78)

    missing = []
    table = {}
    for cmpd in COMPOUNDS:
        finals, occupancies = [], []
        for rep in REPLICATES:
            path = os.path.join(args.mdroot, cmpd, rep, "analysis", f"{args.metric}.xvg")
            if not os.path.isfile(path):
                missing.append(f"{cmpd}/{rep}")
                print(f"{cmpd:10s} {rep:5s}  MISSING  {path}")
                continue
            t, d = read_xvg(path)
            first, final = departures(t, d, args.threshold)
            occ = 100.0 * float((d <= args.threshold).mean())
            occupancies.append(occ)
            if final is not None:
                finals.append(final)
            print(f"{cmpd:10s} {rep:5s}  length {t[-1]:6.1f} ns   "
                  f"first crossing {fmt(first):>10s}   "
                  f"final departure {fmt(final):>12s}   "
                  f"contact {occ:5.1f}% of frames")
        if finals:
            m, s = float(np.mean(finals)), float(np.std(finals, ddof=1)) if len(finals) > 1 else 0.0
            table[cmpd] = (m, s, len(finals), float(np.mean(occupancies)))
            print(f"{cmpd:10s} ---->  final departure {m:.2f} +/- {s:.2f} ns "
                  f"across n = {len(finals)} replicates, "
                  f"mean contact {np.mean(occupancies):.1f}% of frames")
        print("-" * 78)

    if table:
        print("\nSentence-ready summary for the Results section")
        print("=" * 78)
        for cmpd, (m, s, n, occ) in table.items():
            if n > 1:
                print(f"  {cmpd} left the functional site at {m:.1f} +/- {s:.1f} ns "
                      f"across {n} replicates.")
            else:
                print(f"  {cmpd}: only {n} replicate analyzed, no standard deviation.")

    if missing:
        print("\nReplicates not yet available: " + ", ".join(missing))
        print("Run md_docked_replicates.sh, then md_replicate_analysis.sh, then this script again.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
