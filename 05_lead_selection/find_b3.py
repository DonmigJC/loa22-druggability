import os, subprocess, tempfile, math, statistics as st
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, QED
from rdkit.Chem.Scaffolds import MurckoScaffold
from meeko import MoleculePreparation, PDBQTWriterLegacy
RDLogger.DisableLog("rdApp.*")

R = os.path.expanduser("~/Research")
REC = os.path.join(R, "loa22_receptor.pdbqt")
OUT = os.path.join(R, "final_leads")
BOX = dict(cx=7.381, cy=-2.952, cz=-10.291, s=27.6)
KEY = {121:{"CB","CG","OD1","OD2"},
       142:{"CB","CG","CD","NE","CZ","NH1","NH2"}}
key=[]
for l in open(REC):
    if l.startswith(("ATOM","HETATM")):
        try: rn=int(l[22:26])
        except ValueError: continue
        if rn in KEY and l[12:16].strip() in KEY[rn]:
            key.append((float(l[30:38]),float(l[38:46]),float(l[46:54])))

used = {MurckoScaffold.MurckoScaffoldSmiles(
        smi) for smi in [
        "O=C(Nc1ccc(O)c(C(=O)O)c1)c1cccc(C(F)(F)F)c1",
        "CC(=O)Nc1ccnc(C(=O)Nc2ccc(-c3nnco3)cc2)c1"]}

d = pd.read_csv(os.path.join(R,"dual_site","bridging_candidates.csv"))
q = d[(d.bridge>=0.75)&(d.d_key<=4.0)&(d.aff_site<=-5.5)&
      (d.QED>=0.55)&(d.MW<=500)&(d.logP<=5.0)&(d.LE>=0.20)]
q = q.sort_values("bridge", ascending=False)

def prep(smi,p):
    m=Chem.AddHs(Chem.MolFromSmiles(smi))
    if AllChem.EmbedMolecule(m,randomSeed=42)==-1: return False
    AllChem.MMFFOptimizeMolecule(m)
    s,ok,_=PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(m)[0])
    if not ok: return False
    open(p,"w").write(s); return True

def dock(lig,out,seed):
    r=subprocess.run(["vina","--receptor",REC,"--ligand",lig,
        "--center_x",str(BOX["cx"]),"--center_y",str(BOX["cy"]),
        "--center_z",str(BOX["cz"]),"--size_x",str(BOX["s"]),
        "--size_y",str(BOX["s"]),"--size_z",str(BOX["s"]),
        "--exhaustiveness","32","--num_modes","9","--seed",str(seed),
        "--out",out],capture_output=True,text=True,timeout=900)
    for line in r.stdout.splitlines():
        p=line.split()
        if p and p[0]=="1":
            try: return float(p[1])
            except (IndexError,ValueError): return None
    return None

def dk(pose):
    c=[(float(l[30:38]),float(l[38:46]),float(l[46:54]))
       for l in open(pose) if l.startswith(("ATOM","HETATM"))
       and l[77:79].strip() not in ("HD","H")]
    return min(math.dist(a,b) for a in c for b in key)

print(f"{'bridge':>8}{'aff32':>9}{'d_key32':>9}{'QED':>7}  verdict  SMILES")
print("-"*92)
for _, r in q.iterrows():
    m=Chem.MolFromSmiles(r.SMILES)
    sc=MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    if sc in used: continue
    with tempfile.TemporaryDirectory() as td:
        lig=os.path.join(td,"l.pdbqt")
        if not prep(r.SMILES,lig): continue
        vals,poses=[],[]
        for sd in (12345,67890,24680):
            po=os.path.join(OUT,f"cand_{sd}.pdbqt")
            v=dock(lig,po,sd)
            if v is not None: vals.append(v); poses.append(po)
        if not vals: continue
        dd=dk(poses[0])
        ok = dd<=4.5
        print(f"{r.bridge:>8.3f}{st.mean(vals):>9.3f}{dd:>9.2f}{r.QED:>7.3f}"
              f"  {'KEEP ' if ok else 'reject'}  {r.SMILES[:42]}")
        if ok:
            for sd,po in zip((12345,67890,24680),poses):
                os.rename(po, os.path.join(OUT,f"LOA22-B3_site_{sd}.pdbqt"))
            print(f"\n>>> LOA22-B3 = {r.SMILES}")
            print(f"    affinity {st.mean(vals):.3f} +/- "
                  f"{st.pstdev(vals):.3f}   d_key {dd:.2f} A")
            print(f"    QED {QED.qed(m):.3f}  MW {Descriptors.MolWt(m):.1f}  "
                  f"logP {Descriptors.MolLogP(m):.2f}")
            break
