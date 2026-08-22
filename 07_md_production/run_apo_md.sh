#!/bin/bash
# =====================================================================
# JOB B (GPU) — Apo (ligand-free) Loa22 control simulation, 50 ns
# Uses CORRECTED CHARMM36 settings: force-switched LJ, DispCorr=no,
# position restraints during equilibration, fixed velocity seed.
# Run:  bash ~/Research/run_apo_md.sh
# =====================================================================
set -eu
MD=~/Research/04_Molecular_Dynamics
W=$MD/Apo/md_run
mkdir -p "$W"; cd "$W"
cp $MD/loa22_processed.gro $MD/loa22_topol.top $MD/loa22_posre.itp .
ln -sf $MD/charmm36-jul2022.ff .
mv loa22_topol.top topol.top

echo "########## JOB B START $(date) ##########"

NB="vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw            = 1.2
rcoulomb        = 1.2
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16
DispCorr        = no
cutoff-scheme   = Verlet"

CONS="constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4"

printf 'integrator = steep\nemtol = 1000.0\nemstep = 0.01\nnsteps = 50000\nnstlist = 10\n%s\npbc = xyz\n' "$NB" > em.mdp
printf 'integrator = steep\nemtol = 1000.0\nemstep = 0.01\nnsteps = 5000\nnstlist = 10\n%s\npbc = xyz\n' "$NB" > ions.mdp

printf 'define = -DPOSRES\nintegrator = md\nnsteps = 100000\ndt = 0.002\nnstxout-compressed = 5000\nnstenergy = 500\nnstlog = 500\ncontinuation = no\n%s\nnstlist = 10\n%s\ntcoupl = V-rescale\ntc-grps = Protein Non-Protein\ntau_t = 0.1 0.1\nref_t = 300 300\npcoupl = no\npbc = xyz\ngen_vel = yes\ngen_temp = 300\ngen_seed = 42\n' "$CONS" "$NB" > nvt.mdp

printf 'define = -DPOSRES\nintegrator = md\nnsteps = 500000\ndt = 0.002\nnstxout-compressed = 5000\nnstenergy = 500\nnstlog = 500\ncontinuation = yes\n%s\nnstlist = 10\n%s\ntcoupl = V-rescale\ntc-grps = Protein Non-Protein\ntau_t = 0.1 0.1\nref_t = 300 300\npcoupl = C-rescale\npcoupltype = isotropic\ntau_p = 2.0\nref_p = 1.0\ncompressibility = 4.5e-5\nrefcoord_scaling = com\npbc = xyz\ngen_vel = no\n' "$CONS" "$NB" > npt.mdp

printf 'integrator = md\nnsteps = 25000000\ndt = 0.002\nnstxout = 0\nnstxout-compressed = 5000\nnstenergy = 5000\nnstlog = 5000\ncontinuation = yes\n%s\nnstlist = 20\n%s\ntcoupl = V-rescale\ntc-grps = Protein Non-Protein\ntau_t = 0.1 0.1\nref_t = 300 300\npcoupl = C-rescale\npcoupltype = isotropic\ntau_p = 2.0\nref_p = 1.0\ncompressibility = 4.5e-5\npbc = xyz\ngen_vel = no\n' "$CONS" "$NB" > md_prod.mdp

echo "--- box, solvate, ionise ---"
gmx editconf -f loa22_processed.gro -o box.gro -c -d 1.2 -bt dodecahedron
gmx solvate  -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp   -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
printf 'SOL\n' | gmx genion -s ions.tpr -o ions.gro -p topol.top \
       -pname NA -nname CL -neutral

echo "--- energy minimisation ---"
gmx grompp -f em.mdp -c ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -ntmpi 1 -ntomp 6 -pin on
tail -12 em.log

echo "--- NVT 200 ps (restrained) ---"
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 6 -nb gpu -pin on

echo "--- NPT 1 ns (restrained) ---"
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 6 -nb gpu -pme gpu -pin on

echo "--- production 50 ns ---"
gmx grompp -f md_prod.mdp -c npt.gro -t npt.cpt -p topol.top -o md_prod.tpr
gmx mdrun -deffnm md_prod -ntmpi 1 -ntomp 6 -nb gpu -pme gpu -bonded gpu -pin on

echo "--- analysis ---"
printf 'q\n' | gmx make_ndx -f md_prod.tpr -o d.ndx >/dev/null 2>&1
NG=$(grep -c '^\[' d.ndx)
printf 'r 115-150 & a N CA C O\nname %d BindingSite\nq\n' "$NG" \
  | gmx make_ndx -f md_prod.tpr -o an.ndx >/dev/null 2>&1
gi() { awk -v g="$2" 'BEGIN{c=0} /^\[/{n=$0; gsub(/[][]/,"",n);
       gsub(/^ +| +$/,"",n); if(n==g){print c; exit} c++}' "$1"; }
P=$(gi an.ndx Protein); B=$(gi an.ndx Backbone)
S=$(gi an.ndx BindingSite); Y=$(gi an.ndx System)
printf '%s\n' "$Y" | gmx trjconv -s md_prod.tpr -f md_prod.xtc -n an.ndx -o w.xtc -pbc whole
printf '%s\n' "$Y" | gmx trjconv -s md_prod.tpr -f w.xtc -n an.ndx -o nj.xtc -pbc nojump
rm -f w.xtc
printf '%s\n%s\n' "$P" "$Y" | gmx trjconv -s md_prod.tpr -f nj.xtc -n an.ndx -o center.xtc -center -ur compact
rm -f nj.xtc
printf '%s\n%s\n' "$B" "$B" | gmx rms -s md_prod.tpr -f center.xtc -n an.ndx -o rmsd_backbone.xvg -tu ns
printf '%s\n%s\n' "$S" "$S" | gmx rms -s md_prod.tpr -f center.xtc -n an.ndx -o rmsd_bsite.xvg -tu ns
printf '%s\n' "$P" | gmx rmsf -s md_prod.tpr -f center.xtc -n an.ndx -o rmsf.xvg -res
awk '!/^[@#]/{s+=$2;ss+=$2*$2;c++} END{m=s/c;
     printf "APO binding-site RMSD: mean %.3f  SD %.3f nm\n", m, sqrt(ss/c-m*m)}' rmsd_bsite.xvg
echo "########## JOB B COMPLETE $(date) ##########"
