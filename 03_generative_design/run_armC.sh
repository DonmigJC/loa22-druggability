#!/bin/bash
# =====================================================================
# ARM C — RL WITH THE BRIDGING OBJECTIVE
#
# Completes the three-arm generative ablation:
#
#   Arm A   untrained prior                              6,049 molecules
#   Arm B   RL, naive scoring at P_0                     6,319 molecules
#   Arm C   RL, corrected scoring at the functional site  this run
#
# Held constant from Arm B so the comparison is meaningful:
#   same prior file, sigma 120, learning rate 0.0001,
#   100 steps, batch size 64, geometric-mean aggregation.
#
# Changed:
#   docking box      P_0  ->  functional site
#   affinity term    linear ramp  ->  reverse sigmoid centred at -6.0
#   new component    bridging (simultaneous contact with both sites)
#   weights          0.6/0.2/0.2  ->  0.40/0.30/0.20/0.10
#   exhaustiveness   4 (unseeded)  ->  8 (seeded)
#
# THE COMPARISON THAT MATTERS
# ---------------------------
# Arm B's library has already been scored with this exact bridging
# metric by the dual-site screen: 22.53% of its molecules exceed
# bridge = 0.5.  Arm C succeeds if it produces a materially higher
# fraction.  That is a like-for-like measurement, not a circular one.
#
# Run:  bash ~/Research/run_armC.sh
# =====================================================================
set -u
R=~/Research
W=$R/arm_c
PY=~/miniconda3/envs/reinvent4/bin/python
mkdir -p "$W"; cd "$W" || exit 1

# Molecular dynamics may be running.  Keep worker count low if so.
if pgrep -f "gmx mdrun" > /dev/null; then
  export ARMC_WORKERS=${ARMC_WORKERS:-3}
  echo "NOTE: molecular dynamics is running."
  echo "      Using ARMC_WORKERS=$ARMC_WORKERS to limit contention."
  echo "      Both jobs will run slower than they would alone."
else
  export ARMC_WORKERS=${ARMC_WORKERS:-6}
  echo "No MD detected; using ARMC_WORKERS=$ARMC_WORKERS"
fi

echo
echo "########## ARM C  $(date) ##########"

# ---------------------------------------------------------------- config
cat > loa22_armC.toml <<'EOF'
run_type = "staged_learning"
device = "cuda:0"
tb_logdir = "tb_logs_armC"
json_out_config = "_armC_resolved.json"
seed = 42

[parameters]
prior_file = "/home/donmigcage/Research/priors/reinvent.prior"
agent_file = "/home/donmigcage/Research/priors/reinvent.prior"
summary_csv_prefix = "armC_generation"
use_checkpoint = false
batch_size = 64

[learning_strategy]
type = "dap"
sigma = 120
rate = 0.0001

[[stage]]
chkpt_file = "armC_stage1.chkpt"
max_steps = 100

[stage.scoring]
type = "geometric_mean"

[[stage.scoring.component]]
[stage.scoring.component.ExternalProcess]
[[stage.scoring.component.ExternalProcess.endpoint]]
name = "Bridge_Affinity"
weight = 0.7
[stage.scoring.component.ExternalProcess.endpoint.params]
executable = "/home/donmigcage/miniconda3/envs/reinvent4/bin/python"
args = "/home/donmigcage/Research/vina_scorer_armC.py"
property = "bridge_score"

[[stage.scoring.component]]
[stage.scoring.component.QED]
[[stage.scoring.component.QED.endpoint]]
name = "QED"
weight = 0.2

[[stage.scoring.component]]
[stage.scoring.component.MolecularWeight]
[[stage.scoring.component.MolecularWeight.endpoint]]
name = "MW"
weight = 0.1
[stage.scoring.component.MolecularWeight.endpoint.transform]
type = "double_sigmoid"
high = 500.0
low = 150.0
coef_div = 500.0
coef_si = 20.0
coef_se = 20.0
EOF

# ------------------------------------------------- PRE-FLIGHT VALIDATION
echo
echo "=== PRE-FLIGHT: does the scorer distinguish known molecules? ==="
echo
echo "  Expected ordering, from the dual-site screen:"
echo "    LOA22-B1/B2/B3   bridge 0.79-0.87, d_key 3.1-3.4 A   HIGH"
echo "    LOA22-1/2/3      bridge 0.01-0.03, d_key 6.7-7.4 A   LOW"
echo "    ethanol          trivially small                      ~0"
echo

printf '%s\n%s\n%s\n%s\n%s\n%s\n%s\n' \
  'O=C(Nc1ccc(O)c(C(=O)O)c1)c1cccc(C(F)(F)F)c1' \
  'CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1' \
  'COC(=O)c1sc(NC(=O)COc2ccc(Cl)cc2F)nc1C' \
  'CC(C)CNC(=O)Nc1ccc2nnc(-c3cc(C(F)(F)F)cc(C(F)(F)F)c3)n2n1' \
  'CC(CC(=O)N1CCc2ccc(C(F)(F)F)cc2C1)c1ccc(O)cc1' \
  'CC(=O)Nc1ccc(F)c(C(=O)NS(=O)(=O)c2ccccc2Cl)c1' \
  'CCO' \
  | $PY $R/vina_scorer_armC.py > preflight.json 2> preflight.err

if [ -s preflight.err ]; then
  echo "  stderr from the scorer:"
  head -10 preflight.err | sed 's/^/    /'
fi

$PY - <<'PY'
import json, sys
try:
    s = json.load(open("preflight.json"))["payload"]["bridge_score"]
except Exception as e:
    print(f"  PRE-FLIGHT FAILED: {e}")
    sys.exit(1)
names = ["LOA22-B1  (bridger)", "LOA22-B2  (bridger)", "LOA22-B3  (bridger)",
         "LOA22-1   (P_0 only)", "LOA22-2   (P_0 only)", "LOA22-3   (P_0 only)",
         "ethanol   (control)"]
for n, v in zip(names, s):
    bar = "#" * int(v * 40)
    print(f"  {n:<22} {v:.4f}  {bar}")
bridgers, old = s[:3], s[3:6]
if all(v == 0 for v in s):
    print("\n  ALL ZERO - the scorer is broken. Aborting."); sys.exit(1)
if min(bridgers) <= max(old):
    print("\n  ?? bridgers did not clearly outscore the P_0-only compounds.")
    print("     Review the transform parameters before trusting the run.")
    sys.exit(1)
print(f"\n  separation: bridgers {min(bridgers):.3f}-{max(bridgers):.3f}  "
      f"vs P_0-only {min(old):.3f}-{max(old):.3f}")
print("  PRE-FLIGHT PASSED")
PY
if [ $? -ne 0 ]; then echo; echo "Aborting Arm C."; exit 1; fi

# ---------------------------------------------------------------- run
echo
echo "=== LAUNCHING ARM C ==="
echo "    100 steps x 64 = 6,400 molecules"
echo "    exhaustiveness 8, seeded, $ARMC_WORKERS workers"
if pgrep -f "gmx mdrun" > /dev/null; then
  echo "    expect 10-16 h while MD is running"
else
  echo "    expect 5-7 h"
fi
echo

~/miniconda3/envs/reinvent4/bin/reinvent -l armC_run.log loa22_armC.toml

# ---------------------------------------------------------------- summary
echo
echo "=== SUMMARY ==="
$PY - <<'PY'
import glob, pandas as pd
f = sorted(glob.glob("armC_generation*.csv"))
if not f:
    print("no output CSV found"); raise SystemExit
d = pd.read_csv(f[-1])
v = d[d.SMILES_state == 1]
print(f"generated {len(d)}   valid {len(v)} ({100*len(v)/len(d):.1f}%)")
print(f"unique    {v.SMILES.nunique()}")
col = [c for c in v.columns if "Bridge_Affinity" in c and "raw" in c]
if col:
    print(f"mean combined score {v[col[0]].mean():.4f}")
print(f"mean QED            {v['QED (raw)'].mean():.4f}")
g = v.groupby("step")["Score"].mean()
print("\nmean composite score by step:")
for s in (1, 10, 25, 50, 75, 100):
    if s in g.index:
        print(f"  step {s:>3}: {g[s]:.4f}")
print("\nNext: score this library with the dual-site metric and compare")
print("      the bridge > 0.5 fraction against Arm B's 22.53%.")
PY
echo "########## ARM C COMPLETE $(date) ##########"
