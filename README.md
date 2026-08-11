![Graphical Summary of README](docs/_static/header.jpg)
QuantumPDB
==============================

**QuantumPDB** (`qp`) automates generation of quantum mechanical (QM) cluster
models from protein structures: structure preparation, protonation, Voronoi
cluster extraction, TeraChem job setup, and post-processing analysis.

Full user documentation: [hjkgrpquantumpdb.readthedocs.io](https://hjkgrpquantumpdb.readthedocs.io/)

## Table of Contents
1. **Overview**
2. **Installation**
3. **Package layout**
4. **Documentation**
5. **Quick example**
6. **Developer guide**
7. **Areas of active development**
8. **Generated file structure**

## 1. Overview
QuantumPDB turns PDB / mmCIF inputs into ready-to-run QM cluster models for
metalloenzyme active sites. Unlike simple distance cutoffs, it builds
hierarchical interaction spheres with Voronoi tessellation.

![Software Diagram](https://raw.githubusercontent.com/davidkastner/quantumPDB/main/docs/_static/QuantumPDB.png)

## 2. Installation
Clone the repository and perform a developer install inside a conda environment:

```bash
git clone git@github.com:davidkastner/quantumPDB.git
cd quantumPDB
conda env create -f environment.yml
conda activate qp
python -m pip install -e .
```

Compatibility is tested for Python 3.10–3.12 (requires Python ≥ 3.8). After
install, the `qp` CLI is available: `qp run`, `qp submit`, and `qp analyze`.

Modeller needs a free academic license key:

```bash
export KEY_MODELLER="XXXX"
```

## 3. Package layout
```
.
├── docs/                  # Sphinx / Read the Docs sources
├── example/               # Runnable examples (see example/basics/)
├── config.yaml            # Annotated config template (valid keys)
└── qp/
    ├── cli.py             # CLI entry point
    ├── structure/         # Fetch, mmCIF conversion, Modeller, NHIE-oxo
    ├── protonate/         # Protoss API and active-site fixes
    ├── cluster/           # Voronoi spheres and cluster I/O
    ├── manager/           # TeraChem job creation and submission
    ├── analyze/           # Job checkup and Multiwfn post-processing
    ├── resources/         # Bundled assets and helper scripts
    └── tests/             # Pytest suite and golden samples
```

## 4. Documentation
User guides and API docs are hosted on Read the Docs. To build locally:

```bash
cd docs
make clean
make html
# open _build/html/index.html
```

See also `docs/README.md` for the Sphinx / conda docs environment.

## 5. Quick example
```bash
qp run -c example/basics/cluster_only.yaml
```

This downloads 1OS7, runs Modeller + Protoss, and writes cluster models under
`example/basics/output/`. See [Quickstart](https://hjkgrpquantumpdb.readthedocs.io/en/latest/quickstart.html)
and `example/basics/README.md`.

## 6. Developer guide

### Push new changes
```
git status
git pull
git add -A .
git commit -m "Change a specific functionality"
git push -u origin main
```

### Making a pull request
```
git checkout main
git pull
git checkout -b new-feature-branch
git add -A
git commit -m "Detailed commit message describing the changes"
git push -u origin new-feature-branch
# Open the PR on GitHub, then:
git checkout main
git pull
git branch -d new-feature-branch
```

### Handle merge conflict
```
git stash push --include-untracked
git stash drop
git pull
```

## 7. Areas of active development
Temporary mmCIF support is available via `qp.structure.mmcif_to_pdb`: local
`.cif` / `.mmcif` inputs and RCSB entries without a classic PDB file are
converted to `{id}.pdb` before Modeller / Protoss / clustering. Multi-character
chain IDs and >3-character residue names are remapped (see
`{id}_mmcif_remap.json`); structures that exceed classic PDB limits
(>99999 atoms or >62 chains) are skipped with a warning in batch runs.
Center-residue selection accepts either the original or remapped residue names
when a remap sidecar is present. Native mmCIF handling for larger entries is
still planned.

## 8. Generated file structure
Example layout after `qp run` (and optionally `qp submit`) for PDB `1a9s` with
`output_dir: dataset/`:

```
.
├── config.yaml
├── proteins.csv
└── dataset
    └── 1a9s
        ├── 1a9s_modeller.pdb
        ├── 1a9s.ali
        ├── 1a9s.pdb
        ├── charge.csv
        ├── count.csv
        ├── Protoss
        │   ├── 1a9s_ligands.sdf
        │   ├── 1a9s_log.txt
        │   └── 1a9s_protoss.pdb
        └── A290                      # metal_id directory
            ├── 0.pdb                 # center
            ├── 1.pdb                 # first sphere
            ├── 2.pdb                 # second sphere
            ├── A290.pdb              # combined cluster
            ├── A290.xyz
            └── wpbeh                 # method directory from qp submit
                ├── A290.xyz
                ├── jobscript.sh
                ├── ptchrges.xyz      # if charge_embedding: true
                ├── qmscript.in       # TeraChem input
                └── .submit_record
```

See `docs/output.rst` and `qp/tests/samples/` for authoritative formats.

### Copyright

Copyright (c) 2024, Kulik Group MIT

#### Acknowledgements

Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
