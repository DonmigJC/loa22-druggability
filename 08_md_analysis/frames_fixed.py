import os, math, subprocess, glob
R = os.path.expanduser("~/Research"); MD = os.path.join(R, "05_MD_v2")
KEY = {121:{"CB","CG","OD1","OD2"}, 142:{"CB","CG","CD","NE","CZ","NH1","NH2"}}
P0 = {61,64,65,66,67,68,69,70,71,72,73,74,75,88,89,92,93,96,97,115,116,117,
      119,163,165,169,171,172,183,185,186,187}
SKIP = {"SOL","HOH","WAT","NA","CL","TIP3"}
DOCKED = {"LOA22-B1":(3.43,3.04),"LOA22-B2":(3.33,3.13),"LOA22-B3":(3.12,2.88)}

def parse(path):
    key,p0,lig = [],[],[]
    for l in open(path):
        if not l.startswith(("ATOM","HETATM")): continue
        rs = l[17:20].strip()
        if rs in SKIP: continue
        nm = l[12:16].strip()
        if nm.startswith("H") or l[76:78].strip()=="H": continue
        try:
            rn = int(l[22:26]); c=(float(l[30:38]),float(l[38:46]),float(l[46:54]))
        except ValueError:
            continue                      # column overflow past 99999 atoms
        if rs in ("UNL","LIG"): lig.append(c); continue
        if rn in KEY and nm in KEY[rn]: key.append(c)
        if rn in P0: p0.append(c)
    return key,p0,lig

for lead in ("LOA22-B1","LOA22-B2","LOA22-B3"):
    run = os.path.join(MD, lead, "rep1"); out = os.path.join(run,"analysis")
    fr  = os.path.join(out,"frames5"); os.makedirs(fr, exist_ok=True)
    ndx = os.path.join(out,"an.ndx")
    if not os.path.exists(os.path.join(out,"center.xtc")): continue
    print(f"\n=== {lead} ===")
    dk0,dp0 = DOCKED[lead]
    print(f"{'time':>8}{'to key':>10}{'to P_0':>10}   verdict")
    print(f"{'docked':>8}{dk0:>10.2f}{dp0:>10.2f}   reference")
    nb=nk=0
    for t in range(0,100001,5000):
        f = os.path.join(fr, f"{t}.pdb")
        if not os.path.exists(f):
            subprocess.run(["gmx","trjconv","-s",os.path.join(run,"md_prod.tpr"),
                "-f",os.path.join(out,"center.xtc"),"-n",ndx,"-o",f,
                "-dump",str(t)], input="0\n", text=True,
                capture_output=True, cwd=run)
        if not os.path.exists(f): continue
        key,p0,lig = parse(f)
        if not lig or not key: print(f"{t//1000:>6} ns  (parse failed)"); continue
        dk = min(math.dist(a,b) for a in lig for b in key)
        dp = min(math.dist(a,b) for a in lig for b in p0) if p0 else 99
        v = "BRIDGING" if dk<4.5 and dp<4.5 else ("key only" if dk<4.5
            else ("P_0 only" if dp<4.5 else "neither"))
        if v=="BRIDGING": nb+=1
        if dk<4.5: nk+=1
        print(f"{t//1000:>6} ns{dk:>10.2f}{dp:>10.2f}   {v}")
    print(f"  --> bridging in {nb}/21 frames; pharmacophore contact in {nk}/21")
