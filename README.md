QuantumPDB
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/davidkastner/quantumpdb/workflows/CI/badge.svg)](https://github.com/davidkastner/quantumpdb/actions?query=workflow%3ACI)
[![Documentation Status](https://readthedocs.org/projects/quantumpdb/badge/?version=latest)](https://quantumpdb.readthedocs.io/en/latest/?badge=latest)

## Table of Contents
1. **Overview**
    * Introduction
    * Purpose
2. **Installation**
    * Installing quantumPDB
    * Prerequisites
3. **What is included?**
    * File structure
4. **Documentation**
    * Read the Docs
    * Examples
5. **Developer Guide**
    * GitHub refresher


## 1. Overview
The purpose of quantumPDB (qp) is to serve as a toolkit for working with our database of DFT-calculated proteins.

File Structure overview: [Miro Board](https://miro.com/app/board/uXjVPSPcRKQ=/?share_link_id=972866970692)


## 2. Installation
Install the package by running the follow commands inside the repository. This will perform a developmental version install. It is good practice to do this inside of a virtual environment. A yaml environmental file has been created to automate the installation of dependencies.

### Creating python environment
All the dependencies can be loaded together using the prebuilt environment.yml file.
Compatibility is automatically tested for python versions 3.8 and higher.
If you are only going to be using the package run:
```bash
conda env create -f environment.yml
source activate qp
```

### Setup developing environment
To begin working with quantumPDB, first clone the repo and then move into the top-level directory of the package.
Then perform a developer install.
Remember to update your GitHub [ssh keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).
```bash
git clone git@github.com:davidkastner/quantumPDB.git
cd quantumPDB
python -m pip install -e .
```

### Command-line interface
All of the functionality of quantumPDB has been organized into a command-line interface (CLI).
With one additional step, the CLI can be called from anywhere.
We just have to setup a shortcut command in your BASHRC.
Add the following line to your BASHRC:
```bash
alias ml='python /path/to/the/quantumPDB/cli.py
```


## 3. What is included?
### File structure
```
.
|── cli.py                         # Command-line interface entry point
├── docs                           # Readthedocs documentation site
└── qp                             # quantumPDB subpackages and modules
    ├── setup                      # Get a PDB and perform necessary quality and structural checks
    │   ├── fetch_pdb              # Get a PDB
    │   ├── check_edia             # Check the quality of each chain
    │   ├── choose_conformer       # Choose the best conformer
    │   ├── to_pdb                 # Save structure to a new PDB
    │   └── to_xyz                 # Save structure to a new XYZ
    ├── structure                  # Correct the PDB structure
    │   ├── missing_hetatoms       # Add any missing heteroatoms
    │   ├── missing_loops          # Use modeller to add missing loops
    │   ├── protoss_upload         # Upload a PDB to Protoss to add hydrogens
    │   └── protoss_download       # Get the structure with hydrogens
    └── clusters                   # Generalizable plotting and vizualization
        ├── find_metal             # Identify the metal center
        └── coordination_shells    # Select the first, second, and tertiary coordination spheres

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

```
git status
git pull
git add -A .
git commit -m "Change a specific functionality"
git push -u origin main
```

#### Making a pull request
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


### Copyright

Copyright (c) 2023, David W. Kastner


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
