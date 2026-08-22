#!/bin/bash
# =====================================================================
# JOB A (CPU) — Prior-only control + corrected doxycycline baseline
#   1. Fetch real doxycycline from PubChem, verify, dock 3x seeded
#   2. Sample 6400 SMILES from the UNTRAINED prior
#   3. Dock them identically to the RL run
# Run:  bash ~/Research/run_prior_control.sh
# =====================================================================
set -u
R=~/Research
WORK=$R/prior_control
mkdir -p "$WORK"; cd "$WORK" || exit 1
CONDA=~/miniconda3/envs/reinvent4/bin/python

echo "########## JOB A START $(date) ##########"

# ---------------------------------------------------------------
# PART 1 — doxycycline
# ---------------------------------------------------------------
echo; echo "===== PART 1: doxycycline baseline ====="
curl -s "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/doxycycline/SDF?record_type=3d" \
     -o doxy_pubchem.sdf
if [ ! -s doxy_pubchem.sdf ]; then
  echo "!! PubChem fetch failed — skipping Part 1"
else
  $CONDA - <<'PY'
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy
m = Chem.MolFromMolFile("doxy_pubchem.sdf", removeHs=False)
if m is None:
    print("!! could not parse SDF"); raise SystemExit
f = rdMolDescriptors.CalcMolFormula(m)
print("formula:", f, " MW: %.2f" % Descriptors.MolWt(m))
print("stereocentres:", Chem.FindMolChiralCenters(m, includeUnassigned=True,
                                                  useLegacyImplementation=False))
print("isomeric SMILES:", Chem.MolToSmiles(m))
if "C22H24N2O8" not in f:
    print("!! WARNING: formula is not C22H24N2O8")
m = Chem.AddHs(m, addCoords=True)
setups = MoleculePreparation().prepare(m)
s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
if ok:
    open("doxycycline_correct.pdbqt", "w").write(s)
    print("wrote doxycycline_correct.pdbqt")
else:
    print("!! PDBQT failed:", err)
PY

  if [ -f doxycycline_correct.pdbqt ]; then
    echo "--- docking 3 seeded replicates, exhaustiveness 32 ---"
    for SEED in 12345 67890 24680; do
      vina --receptor $R/loa22_receptor.pdbqt --ligand doxycycline_correct.pdbqt \
           --center_x 2.30 --center_y 7.09 --center_z 2.02 \
           --size_x 25 --size_y 25 --size_z 25 \
           --exhaustiveness 32 --num_modes 9 --seed $SEED \
           --out doxy_docked_${SEED}.pdbqt > doxy_log_${SEED}.txt 2>&1
      BEST=$(awk '$1=="1"{print $2; exit}' doxy_log_${SEED}.txt)
      echo "   seed $SEED : $BEST kcal/mol"
    done
    awk '$1=="1"{print $2}' doxy_log_*.txt \
      | awk '{s+=$1; ss+=$1*$1; c++} END{m=s/c;
             printf "   MEAN %.3f  SD %.3f  (n=%d)\n", m, sqrt(ss/c-m*m), c}'
  fi
fi

# ---------------------------------------------------------------
# PART 2 — sample the untrained prior
# ---------------------------------------------------------------
echo; echo "===== PART 2: sampling untrained prior ====="
cat > sample.toml <<'EOF'
run_type = "sampling"
device = "cuda:0"

[parameters]
model_file = "/home/donmigcage/Research/priors/reinvent.prior"
output_file = "prior_sample.csv"
num_smiles = 6400
unique_molecules = false
randomize_smiles = true
EOF
~/miniconda3/envs/reinvent4/bin/reinvent -l sample.log sample.toml 2>&1 | tail -20
if [ ! -f prior_sample.csv ]; then
  echo "!! sampling failed — see sample.log; aborting Part 3"; exit 1
fi
wc -l prior_sample.csv

# ---------------------------------------------------------------
# PART 3 — dock the prior sample
# ---------------------------------------------------------------
echo; echo "===== PART 3: docking prior sample ====="
$CONDA $R/dock_prior.py prior_sample.csv prior_docked.csv 2>&1 | tail -40

echo; echo "===== SUMMARY ====="
$CONDA - <<'PY'
import pandas as pd
d = pd.read_csv("prior_docked.csv").dropna(subset=["affinity"])
print("docked OK      :", len(d))
print("mean affinity  : %.3f kcal/mol" % d.affinity.mean())
print("median         : %.3f" % d.affinity.median())
print("best           : %.3f" % d.affinity.min())
for t in (-6.918, -7.5, -8.945):
    n = (d.affinity <= t).sum()
    print(f"<= {t:7.3f} : {n:5d}  ({100*n/len(d):.2f}%)")
PY
echo "########## JOB A COMPLETE $(date) ##########"
