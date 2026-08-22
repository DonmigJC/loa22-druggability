load loa22_alphafold2_raw.pdb, loa22
hide everything
bg_color white
set ray_opaque_background, 1
set surface_quality, 1
set ray_shadows, 0

# P_0 residues
select p0, loa22 and resi 61+64+65+66+67+68+69+70+71+72+73+74+75+88+89+92+93+96+97+115+116+117+119+163+165+169+171+172+183+185+186+187
# functional site residues (CDD-annotated plus Arg142)
select fsite, loa22 and resi 80+81+120+121+124+128+135+138+178+182+142
select pharm, loa22 and resi 121+142

show surface, loa22
color grey85, loa22
color skyblue, p0
color firebrick, fsite
show sticks, pharm and not name N+C+O
color yellow, pharm and elem C

set transparency, 0.15
orient loa22 and resi 78-186
png fig2A.png, width=2400, height=1800, dpi=300, ray=1

# Panel B: close-up
hide surface
show cartoon, loa22
set cartoon_transparency, 0.7
color grey70, loa22
show sticks, p0
color skyblue, p0 and elem C
show sticks, pharm
color firebrick, pharm and elem C
orient pharm
zoom pharm, 8
png fig2B.png, width=2400, height=1800, dpi=300, ray=1
