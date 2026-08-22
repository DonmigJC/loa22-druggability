# Superseded scripts

These are earlier versions retained for transparency. Each was replaced
because it produced an incorrect result, and the manuscript discusses several
of these corrections directly. They are not required to reproduce the reported
results.

| Script | Replaced by | Reason |
|---|---|---|
| `transfer_poses.py` | `transfer_poses_v4.py` | Assumed matching atom order between the docked pose and the parameterised molecule. Format conversion permutes atom order, so chemically distinct atoms were superimposed, giving RMSD 4.0–5.8 Å. |
| `transfer_poses_v3.py` | `transfer_poses_v4.py` | Added an element-order guard that correctly refused to proceed, revealing the permutation. Matched via an Open Babel MOL2 intermediate that RDKit could not kekulise. |
| `prep_md_leads.py` | `param_relaxed.py` | Parameterised from the docked conformer. The strained geometry prevented convergence of the AM1-BCC charge calculation. |
| `md_analysis.sh` | `analyse_md_v2.sh` | Applied no periodic-boundary correction and fitted the ligand onto itself, which mathematically removes the translation being measured. Reported an apparent dissociation of 11.83 nm — the box width. |
| `vina_scorer_v2.py` | `vina_scorer_armC.py` | Intermediate scorer that docked into the P_0 volume, in which the Asp121 side chain lies outside the searched region. |
| `diag.py` | `diagnose_geometry.py` | Early diagnostic, superseded by the fuller version. |

## Note on the pose-transfer sequence

Four attempts were required. The first two failed silently or wrote output
despite a bad RMSD; the third added a check that refused to proceed and thereby
identified the cause; the fourth succeeded by matching atoms on molecular
connectivity rather than positional index, achieving 0.005 Å.

The verification threshold in `transfer_poses_v4.py` (reject above 0.05 Å) is
what made the difference. A check that cannot fail provides no information.
