#!/bin/bash

set -e

ENV_NAME="dock_env"
OUTDIR="docking_results"

echo "Setting up environment..."

eval "$(conda shell.bash hook)"

if conda env list | grep -q "$ENV_NAME"; then
    conda remove -n $ENV_NAME --all -y
fi

conda create -n $ENV_NAME -c conda-forge \
    python=3.10 vina rdkit scipy numpy pandas gemmi scikit-learn pip -y

conda activate $ENV_NAME
pip install meeko

mkdir -p $OUTDIR/grids $OUTDIR/merged $OUTDIR/top_poses

# User Input

echo "Select docking mode:"
echo "1 = Blind docking"
echo "2 = Targeted docking"
read MODE

PROTEIN="protein.pdb"
LIGAND="ligand.sdf"

[ ! -f "$PROTEIN" ] && echo "❌ protein.pdb missing" && exit 1
[ ! -f "$LIGAND" ] && echo "❌ ligand.sdf missing" && exit 1

echo "Cleaning protein..."
awk '$1=="ATOM"' "$PROTEIN" > protein_clean.pdb


# Prep Receptor

echo "Preparing receptor..."
python - << EOF
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

mol = Chem.MolFromPDBFile("protein_clean.pdb", removeHs=False)
mol = Chem.AddHs(mol)

prep = MoleculePreparation()
setups = prep.prepare(mol)

writer = PDBQTWriterLegacy()

with open("protein.pdbqt","w") as f:
    for s in setups:
        r = writer.write_string(s)
        f.write(r[0] if isinstance(r, tuple) else r)
EOF

grep -v -E "ROOT|ENDROOT|BRANCH|ENDBRANCH|TORSDOF" protein.pdbqt > tmp && mv tmp protein.pdbqt


# Prep Ligand

echo "Preparing ligand..."
python - << EOF
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

mol = Chem.SDMolSupplier("ligand.sdf")[0]
mol = Chem.AddHs(mol)

AllChem.EmbedMolecule(mol)
AllChem.UFFOptimizeMolecule(mol)

prep = MoleculePreparation()
setups = prep.prepare(mol)

writer = PDBQTWriterLegacy()

with open("ligand.pdbqt","w") as f:
    for s in setups:
        r = writer.write_string(s)
        f.write(r[0] if isinstance(r, tuple) else r)
EOF


# Grid Generation (Blind + Safe Fpocket)

echo "Generating docking grids..."

if command -v fpocket &> /dev/null
then
    echo "fpocket detected"
    fpocket -f protein_clean.pdb
    USE_FPOCKET=1
else
    echo "fpocket not found, continuing without it"
    USE_FPOCKET=0
fi

python - << EOF
import numpy as np, json, os

USE_FPOCKET = int("$USE_FPOCKET")

coords = []

with open("protein_clean.pdb") as f:
    for line in f:
        if line.startswith("ATOM"):
            coords.append([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ])

coords = np.array(coords)

min_c = coords.min(axis=0)
max_c = coords.max(axis=0)

GRID = 30
STEP = GRID * 0.8

centers = []

if "$MODE" == "1":
    def has_atoms(c):
        d = np.linalg.norm(coords - c, axis=1)
        return np.sum(d < 15) > 20

    for x in np.arange(min_c[0], max_c[0], STEP):
        for y in np.arange(min_c[1], max_c[1], STEP):
            for z in np.arange(min_c[2], max_c[2], STEP):
                c = [x,y,z]
                if has_atoms(c):
                    centers.append({"center":c, "type":"blind"})
else:
    cx = float(input("center_x: "))
    cy = float(input("center_y: "))
    cz = float(input("center_z: "))
    sx = float(input("size_x: "))
    sy = float(input("size_y: "))
    sz = float(input("size_z: "))
    centers.append({"center":[cx,cy,cz], "type":"targeted", "size":[sx,sy,sz]})

def is_new(c, existing):
    for e in existing:
        if np.linalg.norm(np.array(c) - np.array(e["center"])) < 15:
            return False
    return True

if USE_FPOCKET:
    pdir = "protein_clean_out/pockets"
    added = 0

    if os.path.exists(pdir):
        for f in os.listdir(pdir):
            if f.endswith("_atm.pdb"):
                pts=[]
                for line in open(os.path.join(pdir,f)):
                    if line.startswith("ATOM"):
                        pts.append([
                            float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54])
                        ])
                if pts:
                    center = np.mean(pts, axis=0).tolist()
                    if is_new(center, centers):
                        centers.append({"center":center, "type":"pocket"})
                        added += 1

    print(f"Pocket grids added: {added}")

print("Total grids:", len(centers))

json.dump(centers, open("grid_centers.json","w"))
EOF


# Docking

echo "Starting docking..."

python - << EOF
import json, subprocess, pandas as pd, os, time

with open("grid_centers.json") as f:
    centers = json.load(f)

results = []
start = time.time()

for i, c in enumerate(centers):

    t0 = time.time()

    grid_dir = f"$OUTDIR/grids/grid_{i}"
    os.makedirs(grid_dir, exist_ok=True)

    cx,cy,cz = c["center"]

    if "size" in c:
        sx,sy,sz = c["size"]
    else:
        sx=sy=sz=30

    grid_type = c["type"]

    percent = round((i+1)/len(centers)*100,1)
    elapsed = round((time.time()-start)/60,2)

    print(f"▶ Grid {i+1}/{len(centers)} ({percent}%) | {grid_type} | elapsed {elapsed} min")

    out_file = f"{grid_dir}/docked.pdbqt"
    log_file = f"{grid_dir}/docking.log"

    cmd = f"""
    vina --receptor protein.pdbqt \
         --ligand ligand.pdbqt \
         --center_x {cx} --center_y {cy} --center_z {cz} \
         --size_x {sx} --size_y {sy} --size_z {sz} \
         --exhaustiveness 4 --cpu 4 \
         --out {out_file}
    """

    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    with open(log_file, "w") as f:
        f.write(res.stdout)

    if res.returncode != 0:
        print("⚠️ failed")
        continue

    print(f"   ✅ done in {round((time.time()-t0)/60,2)} min")

    for line in res.stdout.split("\\n"):
        if line.strip().startswith(tuple(str(i) for i in range(1,10))):
            parts = line.split()
            if len(parts) >= 2:
                results.append({
                    "grid": i,
                    "pose": int(parts[0]),
                    "affinity": float(parts[1]),
                    "type": grid_type,
                    "cx": cx, "cy": cy, "cz": cz,
                    "file": out_file
                })

df = pd.DataFrame(results).sort_values("affinity")
df.to_csv("$OUTDIR/all_results.csv", index=False)

print("\\n🏆 TOP 10 POSES:")
print(df.head(10))
EOF

conda activate dock_env

python << EOF
import pandas as pd
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdmolops
import tempfile, os, subprocess


# Config

OUTDIR = "docking_results"
TOPDIR = f"{OUTDIR}/top_poses"
PROTEIN_PDB = "protein.pdb"
TOP_N = 20

print("-" * 50)
print("Final conversion step")
print("-" * 50)

os.makedirs(TOPDIR, exist_ok=True)


# Clean Pdbqt

def clean_pdbqt(infile):
    tmp = tempfile.NamedTemporaryFile(delete=False, mode="w")

    for line in open(infile):
        if line.startswith(("ATOM","HETATM")):
            elem = line[77:79].strip()
            if len(elem) == 0 or not elem.isalpha() or len(elem) > 2:
                elem = "C"
            line = line[:77] + elem.rjust(2) + line[79:]
        tmp.write(line)

    tmp.close()
    return tmp.name


# Load Csv

df = pd.read_csv(f"{OUTDIR}/all_results.csv").sort_values("affinity").head(TOP_N)


# Load Protein (Critical)

protein_mol = Chem.MolFromPDBFile(PROTEIN_PDB, removeHs=False)

if protein_mol is None:
    raise ValueError("❌ Protein failed to load. Fix protein.pdb first.")

print(f"📊 Processing top {TOP_N} poses\n")

success = 0


# Loop

for i, row in df.iterrows():

    print(f"▶ Pose {i+1}/{TOP_N}")

    try:
        pdbqt_file = row["file"]

        if not os.path.exists(pdbqt_file):
            print("   ⚠️ Missing PDBQT")
            continue

        clean_file = clean_pdbqt(pdbqt_file)

        pdbqt_mol = PDBQTMolecule.from_file(clean_file)
        mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)

        if not mols:
            print("   ⚠️ No poses found")
            continue

        pose_idx = int(row["pose"]) - 1

        if pose_idx >= len(mols) or mols[pose_idx] is None:
            print("   ⚠️ Invalid pose → using BEST pose")
            pose_idx = 0

        lig = mols[pose_idx]

        aff = round(row["affinity"], 2)
        gid = int(row["grid"])

        base = f"{TOPDIR}/pose_{i+1}_aff_{aff}_grid_{gid}"

        pdb_path = base + ".pdb"
        sdf_path = base + ".sdf"
        mol2_path = base + ".mol2"
        complex_path = base.replace("pose_", "complex_") + ".pdb"


        # SAVE LIGAND

        Chem.MolToPDBFile(lig, pdb_path)
        Chem.MolToMolFile(lig, sdf_path)

        # optional MOL2 export
        result = subprocess.run(
            f'obabel "{sdf_path}" -O "{mol2_path}"',
            shell=True, capture_output=True, text=True
        )

        if result.returncode != 0 or not os.path.exists(mol2_path):
            print("   ⚠️ MOL2 failed (ignored)")
        else:
            print(f"   ✅ MOL2 → {mol2_path}")


        # combine protein and ligand

        complex_mol = rdmolops.CombineMols(protein_mol, lig)

        Chem.MolToPDBFile(complex_mol, complex_path)


        # rename ligand residue

        fixed_lines = []
        with open(complex_path) as f:
            for line in f:
                if line.startswith("HETATM"):
                    line = line[:17] + "LIG" + line[20:]
                fixed_lines.append(line)

        with open(complex_path, "w") as f:
            f.writelines(fixed_lines)

        print(f"   ✅ PDB → {pdb_path}")
        print(f"   ✅ SDF → {sdf_path}")
        print(f"   ✅ COMPLEX → {complex_path}")

        success += 1

    except Exception as e:
        print(f"   ❌ Failed: {e}")


# Summary

print("\n=========================================")
print("Done")
print(f"✅ Successful: {success}/{TOP_N}")
print(f"📁 Output: {TOPDIR}")
print("-" * 50)
EOF
# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------
print("\n=========================================")
print("Ready for CHARMM-GUI")
print(f"✅ Successful: {success}/{TOP_N}")
print(f"📁 Output: {TOPDIR}")
print("-" * 50)
EOF

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------
print("\n=========================================")
print("Conversion completed")
print(f"✅ Successful: {success}")
print(f"❌ Failed: {fail}")
print(f"📁 Output folder: {TOPDIR}")
print("-" * 50)
EOF

# ------------------------------------------------------------------
# Top Poses Export
# ------------------------------------------------------------------
echo "Exporting top poses..."

python - << EOF
import pandas as pd
from rdkit import Chem

df=pd.read_csv("$OUTDIR/all_results.csv").head(10)

def parse(file):
    poses=[];cur=[]
    for l in open(file):
        if l.startswith("MODEL"): cur=[]
        elif l.startswith("ENDMDL"): poses.append(cur)
        elif l.startswith(("ATOM","HETATM")): cur.append(l)
    return poses

def build(lines):
    m=Chem.RWMol()
    conf=Chem.Conformer(len(lines))
    for i,l in enumerate(lines):
        x,y,z=map(float,[l[30:38],l[38:46],l[46:54]])
        el=l[77:79].strip() or "C"
        idx=m.AddAtom(Chem.Atom(el))
        conf.SetAtomPosition(idx,(x,y,z))
    m.AddConformer(conf)
    return m

poses=parse("$OUTDIR/merged/all_docked.pdbqt")

for i,row in df.iterrows():
    try:
        m=build(poses[i])
        aff=round(row["affinity"],2)
        gid=int(row["grid"])
        Chem.MolToPDBFile(m,f"$OUTDIR/top_poses/pose_{i+1}_aff_{aff}_grid_{gid}.pdb")
    except: pass

print("Top poses exported")
EOF

echo "Docking pipeline finished."
