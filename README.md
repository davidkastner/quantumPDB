QuantumPDB
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/davidkastner/quantumpdb/workflows/CI/badge.svg)](https://github.com/davidkastner/quantumpdb/actions?query=workflow%3ACI)
[![Documentation Status](https://readthedocs.org/projects/quantumpdb/badge/?version=latest)](https://quantumpdb.readthedocs.io/en/latest/?badge=latest)

## Table of Contents
1. **Overview**
2. **Installation**
    * Download the package
    * Creating python environment
    * Command-line interface
3. **What is included?**
    * File structure
4. **Documentation**
    * Read the Docs
    * Examples
5. **Developer Guide**
    * GitHub refresher
6. **Areas of Active Development**
7. **QuantumPDB File Structure**  


## 1. Overview
The purpose of quantumPDB (qp) is to serve as a toolkit for working with our database of proteins to setup and facilitate DFT-calculation and cluster model creation.

![Software Diagram](https://raw.githubusercontent.com/davidkastner/quantumPDB/main/docs/_static/QuantumPDB.png)

## 2. Installation
Install the package by running the follow commands inside the repository. This will perform a developmental version install. It is good practice to do this inside of a virtual environment. A yaml environmental file has been created to facilitate the installation of dependencies.

### Setup developing environment
To begin working with quantumPDB, first clone the repo and then move into the top-level directory of the package.
Then perform a developer install.
Remember to update your GitHub [ssh keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).
```bash
git clone git@github.com:davidkastner/quantumPDB.git
```

### Creating python environment
All the dependencies can be loaded together using the prebuilt environment.yml file.
Compatibility is automatically tested for python versions 3.8 and higher.
If you are only going to be using the package run:
```bash
cd quantumPDB
conda env create -f environment.yml
source activate qp
python -m pip install -e .
```

### Command-line interface
All of the functionality of quantumPDB has been organized into a command-line interface (CLI).
After performing the developer install, the CLI can be called from anywhere using `qp`.


## 3. What is included?
### Package architecture
```
.
├── docs                           # Readthedocs documentation site
└── qp                             # quantumPDB subpackages and modules
    |── cli.py                     # Command-line interface entry point
    ├── checks                     # Perform quality and structural checks
    │   ├── fetch_pdb              # Get a PDB
    │   └── to_xyz                 # Save structure to a new XYZ
    ├── structure                  # Correct the PDB structure
    │   ├── missing_loops          # Use modeller to add missing loops
    │   └── add_hydrogens          # Get the structure with hydrogens
    ├── clusters                   # Generalizable plotting and vizualization
    │   └── coordination_spheres   # Select the first, second, etc. spheres
    └── manager  
        ├── failure_checkup
        ├── find_incomplete
        └── job_manager
```


## 4. Documentation
### Run the following commands to update the ReadTheDocs site
```bash
make clean
make html
```


## 5. Developer guide

### GitHub refresher
#### Push new changes
Use this Git sequence to make a quick push.

```
git status
git pull
git add -A .
git commit -m "Change a specific functionality"
git push -u origin main
```

#### Making a pull request
Use this Git sequence to make a branch and make a pull request.
Recommend for significant changes.

```
git checkout main
git pull

# Before you begin making changes, create a new branch
git checkout -b new-feature-branch
git add -A
git commit -m "Detailed commit message describing the changes"
git push -u origin new-feature-branch

# Visit github.com to add description, submit, merge the pull request

# Once finished on github.com, return to local
git checkout main
git pull

# Delete the remote branch
git branch -d new-feature-branch
```

#### Handle merge conflict

```
git stash push --include-untracked
git stash drop
git pull
```
## 6. Areas of active development
Currently working on handling all edge cases, including non-canonical amino acids. Additionally, support for mmCIFs will eventually needed to be added to work with newer and larger PDBs. Documentation is present for all functions in the code, but should be added with external examples for use. 

## 7. QuantumPDB generated file structure
An example file structure from a `qp` run given the parameters in `config.yaml` `input: qp_input.csv` and `output_dir: dataset/`.
If we run `qp run -c ./config.yaml`, then the file structure will be generated in `dataset` if `qp_input.csv` specfies `1a9s`.
```
.
├── config.yaml                           # The input yaml containing all `qp` job parameters 
├── qp_input.csv                          # List of PDB ID's to run `qp`
└── dataset                               # Specified with "output_dir: dataset/" in config.yaml
    └── 1a9s                              # One of these is generated for each entry in qp_input.csv
        ├── 1a9s_modeller.pdb             # Modeller optimized structure with added missing atoms and loops
        ├── 1a9s.ali                      # Alignment file needed by Modeller
        ├── 1a9s.pdb                      # Original 1a9s PDB downloaded directly from the Protein Data Bank
        ├── charge.csv                    # Generated file used to calculate the charge of the system
        ├── count.csv                     # Generated file that keeps track of the residue counts
        ├── Protoss                       # Directory storing all Protoss files
        │   ├── 1a9s_ligands.sdf          # Ligand structural files
        │   ├── 1a9s_log.txt              # Error log from Protoss server
        │   └── 1a9s_protoss.pdb          # Protonated result returned from the Protoss server
        └── A290                          # Generated for each cluster, named by chain (A) and the res number of center (290)
            ├── 0.pdb                     # Structure of the center residue
            ├── 1.pdb                     # Structure of the first sphere around the center
            ├── 2.pdb                     # Structure of the second sphere around the first
            ├── A290.pdb                  # Structure of the entire cluster in PDB format
            ├── A290.xyz                  # Structure of the entire cluster in XYZ format
            └── wpbeh                     # QM job file specified with "method: wpbeh" in "config.yaml"
                ├── A290.xyz              # Structure of the entire cluster in XYZ format
                ├── jobscript.sh          # SLURM/SGE submit script
                ├── ptchrges.xyz          # MM embedded point charges specified with "charge_embedding: true"
                ├── qmscript.in           # QM job input using TeraChem
                ├── qmscript.out          # QM job output details
                └── wpbeh                 # Results from the TeraChem QM calculation
                    ├── A290.basis
                    ├── A290.geometry
                    ├── A290.molden
                    ├── bond_order.list
                    ├── c0
                    ├── charge_mull.xls
                    ├── grad.xyz
                    ├── mullpop
                    ├── results.dat
                    └── xyz.xyz
```

### Copyright

Copyright (c) 2024, Kulik Group MIT

#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
