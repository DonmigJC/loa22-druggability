# Loa22 druggability assessment: analysis code and data

[![DOI](https://zenodo.org/badge/1342212159.svg)](https://doi.org/10.5281/zenodo.22054886)

Code and data supporting:

> **Spatial dissociation between the druggable pocket and the functional
> peptidoglycan-binding site explains the intractability of *Leptospira* Loa22
> to conventional small-molecule design**
>
> [Authors]. [Journal], [year]. DOI: [pending]

---

## Contents

Every script used to produce the results reported in the manuscript, the
generated compound libraries with docking scores, the molecular dynamics
protocol files, and the structural inputs. Folder numbering follows the
pipeline order and maps onto the Methods subsections.

Molecular dynamics trajectory files (approximately 2 GB per 100 ns simulation,
1,250 ns simulated in total) are not included here and are available from the
corresponding author on request.

---

## Environments

Four Conda environments were used to avoid dependency conflicts between
pipeline stages. To recreate:

    conda env create -f environment/reinvent4_environment.yml
    conda env create -f environment/acpype_env_environment.yml
    conda env create -f environment/mmpbsa_environment.yml
    conda env create -f environment/pymol_env_environment.yml

**Hardware used:** AMD Ryzen 5 7600X (6 cores / 12 threads), NVIDIA RTX 4070
SUPER (12 GB), 15 GiB RAM, under WSL2 running Ubuntu 22.04.5. GROMACS 2021.6
was compiled from source with CUDA acceleration targeting compute capability
8.9.

---

## Pipeline

### 01_structure_validation - Methods 4.2, Structure retrieval, validation and residue numbering

`extract_plddt.py` reports AlphaFold per-residue confidence by structural
region. The OmpA domain (residues 78 to 186) has a mean pLDDT of 97.01 with a
minimum of 90.69. The lower whole-model mean of 83.16 reflects the disordered
N-terminal region, which also contains all twelve Ramachandran outliers.

### 02_cavity_detection - Methods 4.3, Cavity detection, search volumes and receptor preparation

`diagnose_geometry.py` measures the distance from a docking search volume
center to specified residues.

**This script identified the study's central problem.** The search volume used
in the initial generative campaign was centered on the highest-scoring cavity
(P_0) with a 25 Å edge, giving a 12.5 Å half-width. The Asp121 side-chain
centroid lies 16.27 Å from that center, outside the searched region entirely.

`find_target_pocket.py` compares all six detected cavities against the
pharmacophoric residues and the ten CDD-annotated peptidoglycan-binding
positions, identifying P_4 as the functional site (8/10 annotated residues) and
P_0 as containing none.

### 03_generative_design - Methods 4.4, Generative design, scoring and filtering

`loa22_config.toml` and `vina_scorer.py` specify the initial campaign.

The scorer's affinity transform is linear:

    reward = clip(affinity, -12, 0) / -12

This choice is consequential and is discussed in the manuscript: it yields a
reward increment of 0.25 for an affinity improvement from −6 to −9 kcal/mol,
while a QED improvement from 0.4 to 0.9 yields 0.50.

`loa22_armC.toml` and `vina_scorer_armC.py` specify the revised objective:
docking relocated to the functional site, a reverse-sigmoid affinity transform
centered at −6.0 kcal/mol, and an added geometric bridging term.

`dock_prior.py`, `run_prior_control.sh` and `compare_prior.py` implement the
untrained-prior control, which established that reinforcement learning improved
drug-likeness without improving predicted affinity (Cliff's delta 0.030).

### 04_dual_site_screening - Methods 4.5, Dual-site screening, lead selection and controls

`dual_site_dock.py` docks each library molecule at both search volumes and
computes the bridging score.

`separate_pharmacophore.py` repeats this with Asp121 and Arg142 measured
independently rather than as a combined group.

**This distinction matters.** The primary bridging metric reports distance to
whichever pharmacophoric residue is nearer. Measured separately, docking within
the functional-site volume places ligands against Asp121 (19.46% of molecules
within 4.5 Å) but essentially never against Arg142 (0.02%).

`evaluate_armC.py` re-scores the revised library with the same independent
metric applied to the initial library, permitting direct comparison.

### 05_lead_selection - Methods 4.5, Dual-site screening, lead selection and controls

`select_leads.py` applies seven simultaneous criteria and enforces Murcko
scaffold diversity. `redock_top20.py` performs validation re-docking at
exhaustiveness 32.

### 06_md_preparation - Methods 4.6, Molecular dynamics, free energy and statistics

`param_relaxed.py` derives GAFF2/AM1-BCC parameters from a relaxed conformer
rather than the docked pose. Force-field parameters are properties of molecular
topology rather than conformation, and the strained docked geometry prevented
convergence of the semi-empirical charge calculation.

`transfer_poses_v4.py` transfers docked coordinates onto the parameterised
molecule by **molecular graph matching** rather than positional index, because
format conversion permutes atom order. Transfer is verified by heavy-atom RMSD
against the docked pose, with a rejection threshold of 0.05 Å. Achieved values
were 0.005 Å.

`build_md_systems.sh` assembles, solvates and minimizes each system. The `.mdp`
files are the complete protocol specification and are included: CHARMM36
force-switched non-bonded treatment with dispersion correction disabled, TIP3P
water, 150 mM NaCl, position restraints active through both equilibration
stages, and the C-rescale barostat.

`apo_v2/` holds the ligand-free build that replaces the original apo run. The
original used SPC/E water with two chloride ions for neutrality and no added
salt, at 300 K, which does not match the protocol Methods states for every
other system. `build_apo_v2.sh` rebuilds it under the same protocol as the
ligand systems: CHARMM36 with TIP3P, 150 mM NaCl by `genion -neutral
-conc 0.15`, 500 ps NVT and 2 ns NPT restrained, 100 ns production at 310 K.
The only deliberate difference is the temperature-coupling group, `Protein
Non-Protein` rather than `Protein_LIG Water_and_ions`, because there is no
ligand to merge into the solute group, and these are default groups, so no custom
index file is needed. The build lands on 117,193 atoms with NA 109 and CL 112,
against 117,185, NA 109 and CL 112 for the ligand systems. The script verifies
the water model and the salt before it will run a replicate.

### 07_md_production - Methods 4.6, Molecular dynamics, free energy and statistics

`run_md_production.sh` runs equilibration and production with per-replicate
velocity seeds (42, 1337, 24601). Each replicate is independently equilibrated
from the minimized structure rather than branched from a shared equilibrated
state, so replicates differ in initial conditions and not only in trajectory
divergence.

`run_apo_md.sh` runs the ligand-free control simulation.

`md_docked_replicates.sh` adds the second and third docked-pose replicates,
reusing every build input from the first and changing only `gen_seed` in
`nvt.mdp`, so the three differ in initial velocities and nothing else.

`run_queue_resumable.sh` is the only thing that should be used to restart the
campaign after an interruption. Both driver scripts skip a replicate only when
`md_prod.gro` exists, and a run killed mid-production has none, so calling a
driver directly would restart that replicate from NVT and discard its
checkpoint. The wrapper resumes partial production first, then calls the
drivers in order. It takes an `flock` and refuses to start while any `gmx`
process is alive, so it cannot create a second `mdrun`.

`md_watch.sh` and `apo_watch.sh` report progress and an estimated finish time.
They only read files and are safe to run against a live campaign.

### 08_md_analysis - Methods 4.6, Molecular dynamics, free energy and statistics

`analyse_md_v2.sh` applies three sequential periodic-boundary corrections
before analysis (`-pbc whole`, `-pbc nojump`, then centring), and computes
ligand RMSD by fitting on protein backbone while measuring ligand deviation.

`bridging_by_rep.py` measures the same distances independently of the
`gmx mindist` series, by plain coordinate geometry on periodicity-corrected
static frames, with no structural fitting. It is a cross-check on that
series rather than a replacement for it, and the occupancies reported in the
manuscript come from the every-frame series. Residue identity is taken
from the PDB residue name and number jointly.

`where_is_ligand.py` was the first corrected analysis and is retained for
reference. `analyse_apo.sh` analyses the ligand-free control.

`md_replicate_analysis.sh` applies the same `gmx mindist` measurement to the
new replicates that was used for the first, so the series are comparable.

`summarize_docked_replicates.py` reports, per replicate, the trajectory length,
the first crossing of the 4.5 Å contact threshold, the final departure and the
fraction of frames in contact. It separates the first crossing from the final
departure deliberately, because a brief excursion is not a departure.

`analyse_conformational_stability_v2.sh` chooses its trajectory by frame count
rather than by file name. A replicate restarted with `-noappend` can hold more
than one concatenation and some of them are stale: LOA22-B2 replicate 1 keeps a
4,619-frame `md_prod_full.xtc` beside the 9,984-frame `md_prod_cat.xtc` that its
distance series used. Picking the longest is self-correcting.

`analyse_conformational_stability_v2.sh` and
`summarize_conformational_stability.py` rebuild the conformational-stability
table from the v2 trajectories. The summary script recomputes the four
originally submitted rows from the v1 files first, so the statistics code is
checked against known values before any new figure is trusted.

### 09_free_energy - Methods 4.6, Molecular dynamics, free energy and statistics

`run_mmpbsa.sh` and `mmpbsa.in` specify the end-state calculations.

Note `radiopt = 1` and `inp = 1`. Poisson-Boltzmann requires radii optimized
for that method. An initial calculation using generalized-Born-optimized radii
(mbondi2) returned a non-physical positive binding free energy and was
discarded. That superseded input file is included as
`mmpbsa_gbradii_superseded.in`.

### 10_validation_controls - Methods 4.5, Dual-site screening, lead selection and controls

`pgn_control.py` docks peptidoglycan fragments of increasing size (GlcNAc,
MurNAc, muramyl dipeptide), all retrieved from PubChem with molecular formulas
verified.

`decoy_validation.py` constructs property-matched but topologically dissimilar
decoys (Morgan fingerprint Tanimoto < 0.30) and tests whether the bridging
metric is confounded by molecular size.

`audit_geometry.py`, `audit_tables.py` and `audit_t3_mmpbsa.py` recompute every
claim in the manuscript that is reachable from the files in this repository,
including the cavity geometry, all of the docking tables and the end-state free
energies. They were written to check the paper rather than to produce it, and
are included so that a reader can run the same check.

### figures

Scripts generating the manuscript figures. `fig2.pml` requires the open-source
PyMOL build; the educational build's license excludes publication figures.

`fig6C_occupancy.py` replaces `fig6C.py` for panel C of the dynamics figure.
The earlier script hardcoded bridging counts from a 21-frame sample taken at
5 ns intervals, and the panel now plots occupancy evaluated at every frame, on the
same basis as the main-text table. Nothing is hardcoded: it reads the distance
series, computes the occupancies and writes no file unless all three compounds
reproduce the published values. It also refuses any pair of distance series
whose lengths disagree, which guards against the truncated analysis files noted
under Known limitations.

`figS6.py` draws the workflow figure. It replaces a version whose section
references came from a draft in which the Materials and Methods was section 2
and the Results was section 3, so every pointer on it was wrong, and which also
cited a table that does not exist. The script checks every number it prints
against the manuscript and writes no file if one is missing, and it refuses to
run if any draft-era reference reappears in its own content.

`frontiers_export.py` wraps a figure script and enforces the submission figure
rules before writing: minimum 8 pt text, 2 pt data lines, 180 mm width, RGB
TIFF with LZW compression. Export resolution is read from `EXPORT_DPI`. Every
submitted TIFF was produced through it, so it is required to regenerate them.

`fig1.py`, `figS1.py` and `fig2A_polish.pml` generate the workflow figure, the
sequence-alignment figure and the final version of panel A of the cavity
figure.

### superseded

Earlier script versions retained for transparency, each replaced because it
produced an incorrect result. See `superseded/README.md`.

`superseded/md_apo_replicates.sh` was written to add replicates to the original
apo run before that build was checked. It would have copied the SPC/E topology
and the missing salt into two new trajectories. It is kept because the attempt
is part of the record, and it now refuses to run against any build whose
topology includes `spce.itp` or carries no sodium.

---

## Data files

| File | Contents | Rows |
|---|---|---|
| `data/libraries/loa22_generation_1.csv` | Initial campaign | 6,400 |
| `data/libraries/armC_generation_1.csv` | Revised objective | 6,400 |
| `data/libraries/prior_docked.csv` | Untrained prior control | 6,049 |
| `data/docking/dual_site_results.csv` | Both sites, initial library | 6,319 |
| `data/docking/armC_dualsite_results.csv` | Both sites, revised library | 6,254 |
| `data/docking/armB_split.csv` | Asp121/Arg142 measured separately | 6,319 |
| `data/docking/armC_split.csv` | As above, revised library | 6,254 |
| `data/docking/bridging_candidates.csv` | Drug-like bridging candidates | 256 |
| `data/docking/final_leads.csv` | Selected compounds | 3 |
| `data/docking/decoy_results.csv` | Decoy validation | 60 |
| `data/docking/pgn_control.csv` | Peptidoglycan fragments | 3 |
| `tables/SupplementaryTableS15.csv` | Affinity against molecular weight and drug-likeness | 6,275 |
| `data/md_summary/` | MM-PBSA results and decomposition | n/a |
| `data/structures/` | AlphaFold model, prepared receptor, docked poses | n/a |

### Reading the generation CSV files

The column `Vina_Loa22 (raw)` stores the **transformed reward on the interval
[0, 1]**, not the raw binding affinity. Recover affinity in kcal/mol as:

    affinity = -12.0 * value

A comment in `05_lead_selection/extract_leads_correct.py` states otherwise and
is incorrect. This is noted because the deposited data would otherwise be
misinterpreted.

---

## Reproducing key results

Activate the generative environment first:

    conda activate reinvent4

**The druggability dichotomy** (Results 2.1):

    python 02_cavity_detection/find_target_pocket.py

**The search volume error** (Results 2.1, the docking box):

    python 02_cavity_detection/diagnose_geometry.py

**The prior control** (Results 2.2):

    python 03_generative_design/compare_prior.py

**Residue-resolved contact** (Results 2.3):

    python -c "
    import pandas as pd
    d = pd.read_csv('data/docking/armB_split.csv').dropna(subset=['d_asp','d_arg'])
    C = 4.5
    print(f'n = {len(d)}')
    print(f'contacting Asp121: {(d.d_asp<C).sum()} ({100*(d.d_asp<C).mean():.2f}%)')
    print(f'contacting Arg142: {(d.d_arg<C).sum()} ({100*(d.d_arg<C).mean():.2f}%)')
    print(f'contacting both  : {((d.d_asp<C)&(d.d_arg<C)).sum()}')
    "

**Library-wide site preference** (Results 2.3):

    python -c "
    import pandas as pd
    d = pd.read_csv('data/docking/dual_site_results.csv').dropna(subset=['aff_site','aff_p0'])
    print(f'n = {len(d)}')
    print(f'mean P_0            : {d.aff_p0.mean():.3f} kcal/mol')
    print(f'mean functional site: {d.aff_site.mean():.3f} kcal/mol')
    print(f'mean penalty        : {(d.aff_site-d.aff_p0).mean():.3f} kcal/mol')
    print(f'prefer functional   : {(d.aff_site<d.aff_p0).sum()} ({100*(d.aff_site<d.aff_p0).mean():.2f}%)')
    "

---

## Known limitations

- The initial generative run used no fixed random seed and is therefore not
  exactly reproducible. The library it produced is deposited here in full.
- The prior model file returns an invalid hash on load, with a `TL` annotation
  indicating possible transfer learning. It is deposited so that the exact
  model used can be recovered.
- One trajectory segment was truncated by a filesystem write conflict during
  a process collision. Analysis used the intact portions, which together span
  the full simulated interval.
- Superseded script versions are retained in `superseded/` with an explanation
  of why each was replaced.

---

## License

Code is released under the MIT License. Data are under CC-BY-4.0. See `LICENSE`.

## Citation

If you use this code or data, please cite the manuscript above.
