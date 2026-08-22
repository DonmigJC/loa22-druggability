#!/usr/bin/env python3
"""
Extract per-residue AlphaFold confidence (pLDDT) by structural region.

AlphaFold stores the pLDDT score in the B-factor column of its PDB output.
This script reads the CA atom of each residue and reports mean and minimum
confidence for the structural regions used in the manuscript.

Corresponds to Methods section 2.2 and Results Table 1.

Usage:
    python extract_plddt.py [structure.pdb]

Default input: loa22_alphafold2_raw.pdb
"""
import sys
import statistics as st

path = sys.argv[1] if len(sys.argv) > 1 else "loa22_alphafold2_raw.pdb"

# Read one confidence value per residue, taken from its CA atom
plddt = {}
try:
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            resnum = int(line[22:26])
            resname = line[17:20].strip()
            plddt[resnum] = (resname, float(line[60:66]))
except FileNotFoundError:
    sys.exit(f"Structure file not found: {path}")

if not plddt:
    sys.exit(f"No CA atoms found in {path}")

# Regions reported in the manuscript
regions = [
    (1, 195, "Full model"),
    (1, 20, "Signal peptide"),
    (1, 77, "N-terminal region"),
    (78, 186, "OmpA domain"),
    (115, 150, "Analysis window"),
    (187, 195, "C-terminal tail"),
]

print(f"{'Region':<22}{'Residues':>12}{'n':>5}{'Mean':>9}{'Min':>9}")
print("-" * 57)

for lo, hi, label in regions:
    vals = [plddt[i][1] for i in plddt if lo <= i <= hi]
    if vals:
        print(f"{label:<22}{f'{lo}-{hi}':>12}{len(vals):>5}"
              f"{st.mean(vals):>9.2f}{min(vals):>9.2f}")

# Residues of mechanistic interest
print()
print("Residues of interest:")
for r in (119, 121, 142):
    if r in plddt:
        name, value = plddt[r]
        print(f"  {name}{r}: {value:.2f}")

# Partition at residue 47, the boundary between the low-confidence
# N-terminal region and the remainder (Results Table 2)
print()
print("Partition at residue 47:")
for lo, hi, label in [(1, 47, "Residues 1-47"), (48, 195, "Residues 48-195")]:
    vals = [plddt[i][1] for i in plddt if lo <= i <= hi]
    if vals:
        print(f"  {label:<18} n={len(vals):>3}  mean pLDDT {st.mean(vals):.2f}")
