# Panel A only, with a white outline behind the distance label so the
# 12.91 reads crisply through the 25 percent transparent surface.
load loa22_alphafold2_raw.pdb, loa22
hide everything
remove hydrogens
remove loa22 and resi 1-47
bg_color white
set ray_opaque_background, 1
set ray_shadows, 0
set two_sided_lighting, 1
set backface_cull, 0
set surface_quality, 1
set antialias, 2
set depth_cue, 0
set specular, 0.15
set orthoscopic, 1
set label_size, 36
set label_color, black
set label_font_id, 7
set label_outline_color, white
set label_distance_digits, 2
set float_labels, 1
set dash_color, black
set dash_gap, 0.35
set dash_width, 4.0
set dash_radius, 0.14

select p0,    loa22 and resi 61+64+65+66+67+68+69+70+71+72+73+74+75+88+89+92+93+96+97+115+116+117+119+163+165+169+171+172+183+185+186+187
select fsite, loa22 and resi 80+81+120+121+124+128+135+138+178+182+142
deselect
pseudoatom cen_p0, pos=[3.102657, 4.779130, 2.720982]
pseudoatom cen_ph, pos=[7.108182, -1.803182, -7.642909]
distance d_centroid, cen_p0, cen_ph
hide everything, cen_p0
hide everything, cen_ph
color black, d_centroid

show surface, loa22
color grey85, loa22
color skyblue, p0
color firebrick, fsite
show dashes, d_centroid
show labels, d_centroid
set_view (\
     0.310160,  -0.509688,  -0.802508,\
    -0.428285,   0.678713,  -0.596591,\
     0.848748,   0.528741,  -0.007782,\
     0.000000,   0.000000, -140.000000,\
     5.105400,   1.488000,  -2.461000,\
   100.000000, 180.000000,  -20.000000 )
turn y, 180
zoom loa22, 3
set transparency, 0.25
ray 2400, 1800
png fig2A_t25.png, dpi=300
print "=== fig2A_t25.png re-rendered with outlined label ==="
