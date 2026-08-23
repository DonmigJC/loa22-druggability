# =====================================================================
# Figure 2 - version 4  (final)
#
# Run:  conda activate pymol_env
#       cd ~/Research
#       pymol -cq fig2.pml
#
# FIXES vs version 3
#   a. the panel A label read 12.90 because the pseudoatom coordinates
#      were mistyped. True P_0 centroid is 3.102657 4.779130 2.720982
#      and the true distance is 12.9144, which displays as 12.91.
#      Full-precision coordinates are now used.
#   b. orthoscopic projection is on, so a distance drawn in the image
#      plane is faithful to its true length rather than foreshortened.
#   c. panel D no longer attempts a cross-section. Clipping a PyMOL
#      surface reveals the inside of the far shell rather than a solid
#      cut face, and the enclosure contrast was not legible. Panel D is
#      now two face-on close-ups at IDENTICAL scale, one per site, which
#      shows a bounded groove against a flat exposed patch directly.
# =====================================================================

load loa22_alphafold2_raw.pdb, loa22
hide everything
remove hydrogens
remove loa22 and resi 1-47

bg_color white
set ray_opaque_background, 1
set ray_shadows, 0
set two_sided_lighting, 1
set backface_cull, 0
set ray_interior_color, grey50
set surface_quality, 1
set antialias, 2
set depth_cue, 0
set specular, 0.15
set orthoscopic, 1
set label_size, 36
set label_color, black
set label_font_id, 7
set label_distance_digits, 2
set float_labels, 1
set dash_color, black
set dash_gap, 0.35
set dash_width, 4.0
set dash_radius, 0.14

select p0,     loa22 and resi 61+64+65+66+67+68+69+70+71+72+73+74+75+88+89+92+93+96+97+115+116+117+119+163+165+169+171+172+183+185+186+187
select p0near, loa22 and resi 70+71+72+115+116+117+119
select cdd,    loa22 and resi 80+81+120+121+124+128+135+138+178+182
select fsite,  loa22 and resi 80+81+120+121+124+128+135+138+178+182+142
select pharm,  loa22 and resi 121+142 and not name N+C+O
deselect

pseudoatom cen_p0, pos=[3.102657, 4.779130, 2.720982]
pseudoatom cen_ph, pos=[7.108182, -1.803182, -7.642909]
distance d_centroid, cen_p0, cen_ph
distance d_nearest, /loa22//A/PRO`70/CD, /loa22//A/ARG`142/CD
distance d_pharm, /loa22//A/ASP`121/OD1, /loa22//A/ARG`142/NH1
hide everything, cen_p0
hide everything, cen_ph
hide dashes
hide labels
color black, d_centroid
color black, d_nearest
color black, d_pharm

python
SHARED = """
set_view (\
     0.310160,  -0.509688,  -0.802508,\
    -0.428285,   0.678713,  -0.596591,\
     0.848748,   0.528741,  -0.007782,\
     0.000000,   0.000000, -140.000000,\
     5.105400,   1.488000,  -2.461000,\
   100.000000, 180.000000,  -20.000000 )
"""
def shared_camera(zoom_sel="loa22", buf=4):
    cmd.do(SHARED)
    cmd.turn("y", 180)
    cmd.zoom(zoom_sel, buf)

def face_on(rows, origin, dist=55.0):
    # Look straight down the outward normal of one site. dist is fixed so
    # both panel D sub-panels are rendered at exactly the same scale.
    v = list(rows) + [0.0, 0.0, -dist] + list(origin) + [dist - 40.0, dist + 40.0, -20.0]
    cmd.set_view(v)

P0_ROWS = [-0.716568, 0.697517, 0.000000,
           -0.209243, -0.214958, 0.953945,
            0.665393, 0.683566, 0.299983]
FS_ROWS = [ 0.286178, 0.958176, 0.000000,
            0.682038, -0.203704, 0.702373,
            0.672998, -0.201004, -0.711809]
P0_CEN = [3.102657, 4.779130, 2.720982]
FS_CEN = [7.108182, -1.803182, -7.642909]
python end

# =====================  PANEL A  =====================================
hide everything
show surface, loa22
color grey85, loa22
color skyblue, p0
color firebrick, fsite
show dashes, d_centroid
show labels, d_centroid
shared_camera("loa22", 3)

set transparency, 0.0
ray 2400, 1800
png fig2A_opaque.png, dpi=300
set transparency, 0.25
ray 2400, 1800
png fig2A_t25.png, dpi=300
set transparency, 0.45
ray 2400, 1800
png fig2A_t45.png, dpi=300
set transparency, 0.0

# =====================  PANEL B  =====================================
hide everything
show cartoon, loa22
color grey80, loa22
set cartoon_transparency, 0.6
show sticks, p0near
color skyblue, p0near and elem C
set stick_radius, 0.16, p0near
show sticks, pharm
color yellow, pharm and elem C
set stick_radius, 0.30, pharm
show dashes, d_nearest
show dashes, d_pharm
show labels, d_nearest
show labels, d_pharm
shared_camera("pharm or p0near", 4)
ray 2400, 1800
png fig2B.png, dpi=300
set cartoon_transparency, 0.0
set stick_radius, 0.25

# =====================  PANEL C  =====================================
hide everything
show cartoon, loa22
color grey80, loa22
show surface, p0
set transparency, 0.55
color skyblue, p0
show spheres, cdd and name CA
set sphere_scale, 0.6, cdd and name CA
color orange, cdd
show sticks, pharm
color yellow, pharm and elem C
shared_camera("loa22", 3)
ray 2400, 1800
png fig2C.png, dpi=300
set transparency, 0.0

# =====================  PANEL D  =====================================
# Two face-on close-ups at identical scale. D1 is P_0, enclosure 0.24,
# volume 909.78 cubic angstrom. D2 is the functional site, enclosure
# 0.00, volume 145.69. Same camera distance, so the two are comparable.
hide everything
show surface, loa22
color grey85, loa22
color skyblue, p0
color firebrick, fsite

python
face_on(P0_ROWS, P0_CEN, 55.0)
python end
ray 1800, 1800
png fig2D1_P0.png, dpi=300

python
face_on(FS_ROWS, FS_CEN, 55.0)
python end
ray 1800, 1800
png fig2D2_fsite.png, dpi=300

# Slightly wider field, in case 55 crops either site
python
face_on(P0_ROWS, P0_CEN, 70.0)
python end
ray 1800, 1800
png fig2D1_P0_wide.png, dpi=300

python
face_on(FS_ROWS, FS_CEN, 70.0)
python end
ray 1800, 1800
png fig2D2_fsite_wide.png, dpi=300

print "=== done. A opaque/t25/t45, B, C, D1/D2 at two scales ==="
