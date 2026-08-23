# =====================================================================
# Figure S3 - all six detected cavities, three orientations.
#
# Run:  conda activate pymol_env
#       cd ~/Research
#       pymol -cq figS3.pml
#
# Replaces figS3.pml and figS3b.pml, which produced the three views in
# two separate runs at 0, 120 and 60 degrees. Fixed here:
#   1. all three views come from one run, so lighting and scale match.
#   2. the rotations are 0, 60 and 120 degrees, matching the caption.
#      The old view3 was at 240 degrees, not 120.
#   3. two_sided_lighting and backface_cull are set, which is what the
#      black artefact in the old view2 was caused by.
#   4. pockets are painted from P_5 up to P_0, so where two pockets
#      share a residue the higher-ranked one wins. Under the old order
#      P_5 overwrote P_0 at residues 61 and 187, P_1 at 54 and 58, and
#      P_3 at 112.
#   5. colours now match Figure 2 for the two pockets that appear in
#      both: P_0 blue, P_4 red. The six-colour set is CVD-checked.
# =====================================================================

load loa22_alphafold2_raw.pdb, loa22
hide everything
bg_color white
set ray_opaque_background, 1
set ray_shadows, 0
set two_sided_lighting, 1
set backface_cull, 0
set surface_quality, 1
set antialias, 2
set depth_cue, 0
set specular, 0.15
set label_size, 18
set label_color, black
set label_font_id, 7
set label_outline_color, white
set float_labels, 0

# structured domain only. Residues 1-47 are the disordered N-terminus.
create dom, loa22 and resi 48-195
delete loa22
show surface, dom
color grey90, dom

set_color pk0, [0.2000, 0.5020, 0.8000]
set_color pk1, [0.9020, 0.6235, 0.0000]
set_color pk2, [0.0000, 0.6196, 0.4510]
set_color pk3, [0.8000, 0.4745, 0.6549]
set_color pk4, [0.6980, 0.1333, 0.1333]
set_color pk5, [0.4157, 0.2392, 0.6039]

select p0, dom and resi 61+64+65+66+67+68+69+70+71+72+73+74+75+88+89+92+93+96+97+115+116+117+119+163+165+169+171+172+183+185+186+187
select p1, dom and resi 48+49+50+52+54+55+56+58+59+60+63+103+107
select p2, dom and resi 85+86+87+90+91+145+148+149+152+153
select p3, dom and resi 110+111+112+113+158+159+161+191+192+193+194
select p4, dom and resi 80+81+82+120+121+123+124+125+132+133+134+135+138+141+142+182
select p5, dom and resi 53+54+57+58+61+112+187+188+189+190
deselect

# painted lowest rank first so P_0 wins any shared residue
color pk5, p5
color pk4, p4
color pk3, p3
color pk2, p2
color pk1, p1
color pk0, p0

# direct labels, so identity never rests on colour alone
pseudoatom lab0, pos=[11.753, 13.665, 6.621], label=P_0
pseudoatom lab1, pos=[-23.795, 19.445, 15.928], label=P_1
pseudoatom lab2, pos=[12.945, -13.542, 16.821], label=P_2
pseudoatom lab3, pos=[-26.803, -10.017, 10.499], label=P_3
pseudoatom lab4, pos=[13.332, -8.329, -22.462], label=P_4
pseudoatom lab5, pos=[-26.917, 10.527, 6.846], label=P_5

orient dom
zoom dom, 16
ray 2000, 1600
png FigureS3_v1.png, dpi=300
turn y, 60
ray 2000, 1600
png FigureS3_v2.png, dpi=300
turn y, 60
ray 2000, 1600
png FigureS3_v3.png, dpi=300

print "=== FigureS3_v1/v2/v3 written at 0, 60 and 120 degrees ==="
