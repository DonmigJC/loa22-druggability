import numpy as np, glob, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read(f):
    t,v = [],[]
    for l in open(f):
        if l.startswith(("@","#")): continue
        p=l.split()
        if len(p)>=2: t.append(float(p[0])); v.append(float(p[1])*10)
    return np.array(t), np.array(v)

COL = {"LOA22-B1":"#C0392B","LOA22-B2":"#2E86C1","LOA22-B3":"#27AE60"}
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)

for lead,c in COL.items():
    f=f"traj_docked_{lead}.xvg"
    if os.path.exists(f):
        t,v=read(f); ax[0].plot(t,v,color=c,lw=.8,label=lead)
    f=f"traj_alt_{lead}.xvg"
    if os.path.exists(f):
        t,v=read(f); ax[1].plot(t,v,color=c,lw=.8,label=lead)

for a,title in zip(ax, ["A  Docked pose (Asp121)","B  Alternative pose (Arg142)"]):
    a.axhline(4.5, ls="--", color="k", lw=1)
    a.set_xlabel("Time (ns)"); a.set_title(title, loc="left", fontweight="bold")
    a.spines[["top","right"]].set_visible(False)
    a.legend(frameon=False, fontsize=8)
ax[0].set_ylabel("Distance to pharmacophore (Å)")
plt.tight_layout(); plt.savefig("Figure6AB.png", dpi=300)
print("wrote Figure6AB.png")
