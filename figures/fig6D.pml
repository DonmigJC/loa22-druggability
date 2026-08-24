# =====================================================================
# Figure 6 panel D - LOA22-B1 at 100 ns, relocated into P_0.
#
# Run:  conda activate pymol_env
#       cd ~/Research/05_MD_v2/LOA22-B1/rep1
#       pymol -cq ~/Research/fig6D.pml
#
# The script finds the production .gro itself. The run used -noappend after
# a restart, so the final frame may be md_prod.gro or md_prod.partNNNN.gro.
# If none is found it stops and says so rather than rendering an empty scene.
#
# Colours match Figure 2:
#   skyblue   P_0 lining residues
#   grey70    contact residues that are NOT part of P_0
#   yellow    LOA22-B1
# =====================================================================

python
import glob, os, sys
# prefer a periodic-boundary-corrected frame if one has been made
cands = sorted(glob.glob("md_prod.part*.gro")) + ["md_prod.gro", "frame_100ns_pbc.gro"]
found = None
for c in cands:
    if os.path.exists(c):
        found = c
if found is None:
    print("=== ERROR: no md_prod*.gro in %s ===" % os.getcwd())
    print("=== files present: %s ===" % ", ".join(sorted(os.listdir("."))[:40]))
    cmd.quit()
print("=== loading %s ===" % found)
cmd.load(found, "sys")
python end

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
set orthoscopic, 1
set label_size, 20
set label_color, black
set label_font_id, 7
set label_outline_color, white
set float_labels, 1
set transparency, 0.55

remove solvent
remove resn NA+CL+SOL+K
remove hydrogens

# P_0 lining residues, same set as Figures 2 and S3
select p0, polymer and resi 61+64+65+66+67+68+69+70+71+72+73+74+75+88+89+92+93+96+97+115+116+117+119+163+165+169+171+172+183+185+186+187

# residues within 5 A of the ligand in this frame
select contacts, byres (polymer within 5 of resn UNL)
# of those, the ones that belong to P_0
select contacts_p0, contacts and p0
select contacts_other, contacts and not p0

show cartoon, polymer
color grey85, polymer
set cartoon_transparency, 0.65

show surface, p0
color skyblue, p0

show sticks, contacts_p0
color skyblue, contacts_p0 and elem C
show sticks, contacts_other
color grey60, contacts_other and elem C
set stick_radius, 0.16, contacts

show sticks, resn UNL
color yellow, resn UNL and elem C
set stick_radius, 0.26, resn UNL

# report what was actually found, so the caption can be checked against it
python
mdl = cmd.get_model("contacts and name CA")
allc = sorted(set((a.resn, int(a.resi)) for a in mdl.atom), key=lambda t: t[1])
p0m = cmd.get_model("contacts_p0 and name CA")
inp0 = sorted(set(int(a.resi) for a in p0m.atom))
print("=== contact residues within 5 A of LOA22-B1 at 100 ns ===")
print("    " + ", ".join("%s%d" % (r.capitalize(), i) for r, i in allc))
print("=== of these, %d of %d are P_0 lining residues ===" % (len(inp0), len(allc)))
print("    " + ", ".join(str(i) for i in inp0))
python end

orient resn UNL or contacts
zoom resn UNL or contacts, 4
ray 2400, 1800
png Figure6D.png, dpi=300

turn y, 30
ray 2400, 1800
png Figure6D_alt.png, dpi=300

print "=== wrote Figure6D.png and Figure6D_alt.png ==="
