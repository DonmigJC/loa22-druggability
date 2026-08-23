#!/usr/bin/env python3
"""
Figure S5 - per-residue energy decomposition, LOA22-B2, alternative (Arg142) pose.

PLOTS: Generalized Born model, DELTAS block, Total Decomposition Contribution
       (TDC) section, VAN DER WAALS component.
       Error bars are standard errors of the mean over 101 frames, as reported
       by gmx_MMPBSA in FINAL_DECOMP.dat.

WHY VAN DER WAALS AND NOT TOTAL:
  The van der Waals term is identical between the Generalized Born and
  Poisson-Boltzmann solvation models, whereas the TOTAL changes sign between
  them (Arg142: GB -2.736, PB +0.198). van der Waals is therefore the
  solvation-model-invariant quantity. See Discussion 4.4.

WHY FINAL_DECOMP.dat AND NOT DECOMP_v2.csv:
  The v2 rerun changed only the &pb namelist (inp 2->1, radiopt 0->1).
  The &gb namelist was untouched, so the GB decomposition is identical in
  FINAL_DECOMP.dat and DECOMP_v2.dat. Either file gives the same GB numbers;
  FINAL_DECOMP.dat is used because it carries the SD/SEM columns directly.

Run from the Research root:   python3 figS5.py
Outputs: FigureS5.png (rep1, matches Table 20)
         FigureS5_3replicate.png (all three replicates)
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------- constants
BASE = "05_MD_v2/LOA22-B2"
REPS = ["rep1_arg142pose", "rep2_arg142pose", "rep3_arg142pose"]
PRIMARY = "rep1_arg142pose"          # the replicate Table 20 reports
HIGHLIGHT = "R:A:ARG:142"
TOP_N = 15                            # residues shown in the main panel

# FINAL_DECOMP.dat column layout: residue, then 6 terms x (Avg, SD, SEM)
TERMS = ["Internal", "van der Waals", "Electrostatic",
         "Polar Solvation", "Non-Polar Solv.", "TOTAL"]
COL = {t: (1 + 3 * i, 2 + 3 * i, 3 + 3 * i) for i, t in enumerate(TERMS)}

# Values independently confirmed against Table 20. The script refuses to write
# a figure if the parse drifts from these. Do not "fix" by loosening tolerance.
EXPECT = {"van der Waals": (-1.715, 0.255), "TOTAL": (-2.736, 2.361)}
TOL = 0.005

GREY, RED, BLUE, GREEN = "#8C8C8C", "#C0392B", "#2E86C1", "#27AE60"


# -------------------------------------------------------------------- parsing
def parse_decomp(path, model="GB", block="DELTAS", section="TDC"):
    """Return {residue: {term: (avg, sd, sem)}} for one model/block/section.

    FINAL_DECOMP.dat is organised as:
        Energy Decomposition Analysis ...: Generalized Born model
          Complex: / Receptor: / Ligand: / DELTAS:
            Total Energy Decomposition:      <- TDC
            Sidechain Energy Decomposition:  <- SDC (absent when idecomp=2)
            Backbone Energy Decomposition:   <- BDC
        ... then the whole thing repeats for the Poisson Boltzmann model.
    Reading the wrong block is the classic failure here: the Complex block
    returns roughly -163 kcal/mol for this system.
    """
    cur_model = cur_block = cur_sect = None
    out = {}
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if "Generalized Born model" in s:
                cur_model = "GB"; continue
            if "Poisson Boltzmann model" in s:
                cur_model = "PB"; continue
            if s.startswith("DELTAS:"):
                cur_block = "DELTAS"; continue
            if s.startswith(("Complex:", "Receptor:", "Ligand:")):
                cur_block = s[:-1]; continue
            if s.startswith("Total Energy Decomposition"):
                cur_sect = "TDC"; continue
            if s.startswith("Sidechain Energy Decomposition"):
                cur_sect = "SDC"; continue
            if s.startswith("Backbone Energy Decomposition"):
                cur_sect = "BDC"; continue
            if not s.startswith(("R:", "L:")):
                continue
            if (cur_model, cur_block, cur_sect) != (model, block, section):
                continue
            f = s.split(",")
            out[f[0]] = {t: tuple(float(f[c]) for c in COL[t]) for t in TERMS}
    return out


def label(res):
    """'R:A:ARG:142' -> 'Arg142'."""
    p = res.split(":")
    return p[2].capitalize() + p[3]


# ------------------------------------------------------------------ load data
def main():
    if not os.path.isdir(BASE):
        sys.exit(f"ERROR: '{BASE}' not found. Run this from the Research root.")

    data = {}
    for rep in REPS:
        p = os.path.join(BASE, rep, "mmpbsa", "FINAL_DECOMP.dat")
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing {p}")
        data[rep] = parse_decomp(p)

    # ------------------------------------------------------------ verification
    print("VERIFICATION - %s, GB/DELTAS/TDC, %s" % (PRIMARY, HIGHLIGHT))
    ok = True
    for term, (exp_avg, exp_sem) in EXPECT.items():
        avg, sd, sem = data[PRIMARY][HIGHLIGHT][term]
        good = abs(avg - exp_avg) < TOL and abs(sem - exp_sem) < TOL
        ok &= good
        print("  %-16s %8.3f +/- %-6.3f  expected %8.3f +/- %-6.3f  %s"
              % (term, avg, sem, exp_avg, exp_sem, "OK" if good else "MISMATCH"))
    if not ok:
        sys.exit("ABORTED: parsed values do not match Table 20. Nothing written.")

    # cross-replicate mean of the van der Waals term (Discussion 4.4)
    vals = [data[r][HIGHLIGHT]["van der Waals"][0] for r in REPS]
    m, s = float(np.mean(vals)), float(np.std(vals, ddof=1))
    print("  cross-replicate vdW mean %.3f +/- %.3f  (manuscript -2.32 +/- 0.70)"
          % (m, s))

    # ------------------------------------------------- panel data (rep1, vdW)
    rec = {k: v for k, v in data[PRIMARY].items() if k.startswith("R:")}
    lig = [k for k in data[PRIMARY] if k.startswith("L:")]
    print("\n  %d receptor residues within 6 A; %d ligand self-term entr%s"
          % (len(rec), len(lig), "y" if len(lig) == 1 else "ies"))
    print("  (gmx_MMPBSA lists the ligand as a residue: %d + %d = %d rows total)"
          % (len(rec), len(lig), len(rec) + len(lig)))

    order = sorted(rec, key=lambda r: rec[r]["van der Waals"][0])[:TOP_N]
    names = [label(r) for r in order]
    avgs = np.array([rec[r]["van der Waals"][0] for r in order])
    sems = np.array([rec[r]["van der Waals"][2] for r in order])
    cols = [RED if r == HIGHLIGHT else GREY for r in order]

    print("\n  Top %d by van der Waals (rep1):" % TOP_N)
    for r in order:
        a, _, e = rec[r]["van der Waals"]
        print("    %-10s %7.3f +/- %.3f%s"
              % (label(r), a, e, "   <- highlighted" if r == HIGHLIGHT else ""))

    # ------------------------------------------------------ Figure S5 (rep1)
    fig, ax = plt.subplots(figsize=(7.0, 5.6))
    y = np.arange(len(order))[::-1]
    ax.barh(y, avgs, xerr=sems, color=cols, height=0.72,
            error_kw=dict(ecolor="#404040", elinewidth=1.0, capsize=2.5))
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("van der Waals contribution (kcal/mol)", fontsize=10)
    ax.axvline(0, color="k", lw=0.8)
    ax.invert_xaxis()                      # more negative = further right
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=9)
    ax.set_title("Generalized Born, DELTAS/TDC - van der Waals component",
                 fontsize=10, loc="left")
    a, _, e = rec[HIGHLIGHT]["van der Waals"]
    ax.annotate("Arg142  %.2f $\\pm$ %.2f" % (a, e),
                xy=(a, y[order.index(HIGHLIGHT)]),
                xytext=(0.97, 0.06), textcoords="axes fraction",
                ha="right", fontsize=9, color=RED)
    fig.tight_layout()
    fig.savefig("FigureS5.png", dpi=300)
    plt.close(fig)
    print("\n  wrote FigureS5.png")

    # ------------------------------------- Figure S5 alternative: 3 replicates
    common = sorted(
        set.intersection(*[{k for k in data[r] if k.startswith("R:")} for r in REPS]),
        key=lambda r: np.mean([data[rp][r]["van der Waals"][0] for rp in REPS]))[:TOP_N]
    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    yy = np.arange(len(common))[::-1]
    h = 0.26
    for i, (rep, c) in enumerate(zip(REPS, [RED, BLUE, GREEN])):
        a = [data[rep][r]["van der Waals"][0] for r in common]
        e = [data[rep][r]["van der Waals"][2] for r in common]
        ax.barh(yy + (1 - i) * h, a, h, xerr=e, color=c,
                label="Replicate %d" % (i + 1),
                error_kw=dict(ecolor="#404040", elinewidth=0.8, capsize=2))
    ax.set_yticks(yy)
    ax.set_yticklabels([label(r) for r in common], fontsize=9)
    ax.set_xlabel("van der Waals contribution (kcal/mol)", fontsize=10)
    ax.axvline(0, color="k", lw=0.8)
    ax.invert_xaxis()
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.set_title("Generalized Born, DELTAS/TDC - van der Waals, three replicates",
                 fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig("FigureS5_3replicate.png", dpi=300)
    plt.close(fig)
    print("  wrote FigureS5_3replicate.png")


if __name__ == "__main__":
    main()
