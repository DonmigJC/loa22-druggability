#!/usr/bin/env python3
"""
ACPYPE RETRY WITH DIAGNOSTICS
=============================

LOA22-B2 and LOA22-B3 failed during ligand parameterisation.  This script
captures the real error and then tries several alternative routes.

Why parameterisation fails, in rough order of frequency:

  1. sqm (the semi-empirical quantum program that computes AM1-BCC
     partial charges) fails to converge.  Common for molecules with
     unusual heterocycles or poor starting geometry.
  2. antechamber cannot assign GAFF atom types, typically for less
     common aromatic heterocycles.
  3. The MOL2 file written by Open Babel has atom or residue naming that
     antechamber rejects.
  4. A stale .acpype directory or temporary folder blocks the run.

Strategies attempted, in order:
  A. clean rerun from the Open Babel MOL2
  B. simplified atom and residue naming (LIG, sequential names)
  C. GAFF instead of GAFF2
  D. Gasteiger charges instead of AM1-BCC  (last resort - disclose it)

Run:  ~/miniconda3/envs/reinvent4/bin/python retry_acpype.py
"""
import os
import re
import shutil
import subprocess
import glob

R = os.path.expanduser("~/Research")
MD = os.path.join(R, "05_MD_v2")
TARGETS = ["LOA22-B2", "LOA22-B3"]


def sh(cmd, cwd=None, timeout=1800):
    return subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True,
                          text=True, timeout=timeout)


def clean(wd, name):
    """Remove artefacts from previous attempts."""
    for pat in (f"{name}.acpype", f".acpype_tmp_{name}", "ANTECHAMBER*",
                "ATOMTYPE.INF", "sqm.*", "NEWPDB.PDB", "PREP.INF",
                f"{name}_AC.*", "leap.log", "*.prmtop", "*.inpcrd"):
        for p in glob.glob(os.path.join(wd, pat)):
            shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) \
                else os.remove(p)


def show_errors(wd, res):
    """Surface the real cause from acpype output and sqm logs."""
    text = (res.stdout or "") + (res.stderr or "")
    keys = [l for l in text.splitlines()
            if re.search(r"error|fail|cannot|unable|Fatal|Traceback",
                         l, re.I)]
    for l in keys[-8:]:
        print(f"      {l.strip()[:110]}")
    for f in glob.glob(os.path.join(wd, "**", "sqm.out"), recursive=True):
        tail = open(f, errors="ignore").read().splitlines()[-15:]
        bad = [l for l in tail if re.search(r"error|conver|fail", l, re.I)]
        if bad:
            print(f"      sqm.out: {bad[-1].strip()[:100]}")


def check(wd, name):
    itp = os.path.join(wd, f"{name}.acpype", f"{name}_GMX.itp")
    gro = os.path.join(wd, f"{name}.acpype", f"{name}_GMX.gro")
    return os.path.exists(itp) and os.path.exists(gro)


def total_charge(wd, name):
    itp = os.path.join(wd, f"{name}.acpype", f"{name}_GMX.itp")
    tot, inb = 0.0, False
    for l in open(itp):
        if l.strip().startswith("[ atoms ]"):
            inb = True
            continue
        if inb:
            if l.strip().startswith("["):
                break
            p = l.split()
            if len(p) >= 7 and not l.strip().startswith(";"):
                try:
                    tot += float(p[6])
                except ValueError:
                    pass
    return tot


def simplify_mol2(src, dst):
    """Rewrite atom names as element+index and the residue as LIG."""
    lines = open(src).read().splitlines()
    out, section, counts = [], None, {}
    for l in lines:
        if l.startswith("@<TRIPOS>"):
            section = l.split(">")[1].strip()
            out.append(l)
            continue
        if section == "ATOM" and l.strip():
            p = l.split()
            if len(p) >= 9:
                el = re.match(r"[A-Za-z]{1,2}", p[5]).group(0)
                counts[el] = counts.get(el, 0) + 1
                p[1] = f"{el}{counts[el]}"
                p[7] = "1"
                p[8] = "LIG"
                out.append(f"{int(p[0]):>7} {p[1]:<8}{float(p[2]):>10.4f}"
                           f"{float(p[3]):>10.4f}{float(p[4]):>10.4f} "
                           f"{p[5]:<8}{p[7]:>4} {p[8]:<8}"
                           f"{float(p[9]) if len(p) > 9 else 0.0:>10.4f}")
                continue
        out.append(l)
    open(dst, "w").write("\n".join(out) + "\n")


print("=" * 78)
print("ACPYPE RETRY WITH DIAGNOSTICS")
print("=" * 78)

results = {}
for name in TARGETS:
    wd = os.path.join(MD, name)
    sdf = os.path.join(wd, f"{name}.sdf")
    if not os.path.exists(sdf):
        print(f"\n{name}: SDF missing — rerun prep_md_leads.py")
        continue

    print(f"\n{'=' * 78}\n{name}\n{'=' * 78}")
    done = False

    # ---------------- A: clean rerun ------------------------------
    print("\n  [A] clean rerun from the Open Babel MOL2")
    clean(wd, name)
    mol2 = os.path.join(wd, f"{name}.mol2")
    sh(f"obabel {sdf} -O {mol2} --partialcharge gasteiger", cwd=wd)
    if os.path.exists(mol2):
        res = sh(f"conda run -n acpype_env acpype -i {name}.mol2 "
                 f"-b {name} -n 0 -a gaff2 -c bcc -o gmx", cwd=wd)
        if check(wd, name):
            print("      SUCCESS")
            done = True
        else:
            show_errors(wd, res)
    else:
        print("      Open Babel produced no MOL2")

    # ---------------- B: simplified naming ------------------------
    if not done:
        print("\n  [B] simplified atom and residue naming")
        clean(wd, name)
        simple = os.path.join(wd, f"{name}_simple.mol2")
        try:
            simplify_mol2(mol2, simple)
            res = sh(f"conda run -n acpype_env acpype -i "
                     f"{os.path.basename(simple)} -b {name} -n 0 "
                     f"-a gaff2 -c bcc -o gmx", cwd=wd)
            if check(wd, name):
                print("      SUCCESS")
                done = True
            else:
                show_errors(wd, res)
        except Exception as e:
            print(f"      rewrite failed: {e}")

    # ---------------- C: GAFF instead of GAFF2 --------------------
    if not done:
        print("\n  [C] GAFF instead of GAFF2")
        clean(wd, name)
        res = sh(f"conda run -n acpype_env acpype -i {name}.mol2 "
                 f"-b {name} -n 0 -a gaff -c bcc -o gmx", cwd=wd)
        if check(wd, name):
            print("      SUCCESS  (note: GAFF, not GAFF2 — disclose this)")
            done = True
        else:
            show_errors(wd, res)

    # ---------------- D: Gasteiger charges ------------------------
    if not done:
        print("\n  [D] Gasteiger charges instead of AM1-BCC")
        clean(wd, name)
        res = sh(f"conda run -n acpype_env acpype -i {name}.mol2 "
                 f"-b {name} -n 0 -a gaff2 -c gas -o gmx", cwd=wd)
        if check(wd, name):
            print("      SUCCESS  (Gasteiger charges — LESS ACCURATE than")
            print("      AM1-BCC and inconsistent with LOA22-B1. Disclose,")
            print("      or preferably resolve the AM1-BCC failure.)")
            done = True
        else:
            show_errors(wd, res)

    results[name] = done
    if done:
        q = total_charge(wd, name)
        print(f"\n  topology total charge: {q:+.3f}  (target +0.000)")
        if abs(q) > 0.05:
            print("  !! charge mismatch")
            results[name] = False

print()
print("=" * 78)
print("SUMMARY")
print("=" * 78)
print(f"{'LOA22-B1':<12} already parameterised (charge -1.000)")
for n in TARGETS:
    print(f"{n:<12} {'READY' if results.get(n) else 'STILL FAILING'}")

if all(results.get(n) for n in TARGETS):
    print("\nAll three ligands parameterised. Next:")
    print("  bash ~/Research/build_md_systems.sh")
else:
    print("\nFor any still failing, send the [A] diagnostic lines above —")
    print("they identify whether the problem is sqm convergence, atom")
    print("typing, or file formatting.")
print("=" * 78)
