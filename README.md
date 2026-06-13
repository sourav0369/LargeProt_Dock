# LargeProt_Dock
Script for automating docking of large protein biomolecules through multi-grid docking protocol. The logic is at beginning stage. The code can be further upgraded for enzymes or other similar proteins that are also dependent on co-factors or metal ions for functioning.

## About the Workflow

Automated workflow for blind and targeted protein–ligand docking of large proteins using overlapping grid decomposition and optional fpocket-guided pocket detection.



## Features
- Automated receptor and ligand preparation
- Blind and targeted docking modes
- Multi-grid sampling strategy for large proteins
- Optional fpocket-assisted grid generation
- Automatic ranking of docking poses
- Export of PDB, SDF, MOL2, and protein–ligand complexes
- Compatible with downstream MD and CHARMM-GUI workflows



## Workflow
protein.pdb + ligand.sdf
            │
            ▼
   Receptor/Ligand Preparation
            │
            ▼
 Grid Generation (Blind + fpocket)
            │
            ▼
    Multi-grid AutoDock Vina
            │
            ▼
      Pose Ranking & Filtering
            │
            ▼
 PDB / SDF / MOL2 / Complex Export
            │
            ▼
    Downstream MD / CHARMM-GUI



## Requirements
- Conda
- AutoDock Vina
- RDKit
- Meeko
- Open Babel (optional)
- fpocket (optional)



## Usage

chmod +x pipeline.sh
./pipeline.sh



Required input files:
- protein.pdb
- ligand.sdf



## Output
The pipeline generates:
- Docking scores (.csv)
- Individual docked PDBQT files
- Top-ranked ligand structures (PDB/SDF/MOL2)
- Protein–ligand complex structures
