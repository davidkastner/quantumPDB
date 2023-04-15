Getting Started
===============

*Welcome to quantumPDB!*

1 Overview
----------
The purpose of quantumPDB (qp) is to serve as a toolkit for working with our database of DFT-calculated proteins.


2 Installation
------------
Install the package by running the follow commands inside the repository. This will perform a developmental version install. It is good practice to do this inside of a virtual environment. A yaml environmental file has been created to automate the installation of dependencies.

**Creating python environment**
All the dependencies can be loaded together using the prebuilt environment.yml or environment_dev.yml files.
We provide two YAML files. The dev version contains additional packages for code maintenance.
Compatibility is automatically tested for python versions 3.8 and higher.
If you are only going to be using the package run:

::

    bash
    conda env create -f environment.yml
    source activate ml


**Setup developing environment**
To begin working with quantumPDB, first clone the repo and then move into the top-level directory of the package.
The perform a developer install.
Remember to update your GitHub [ssh keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).

::

    bash
    git clone git@github.com:davidkastner/quantumPDB.git
    cd quantumPDB
    python -m pip install -e .


3 What is included?
-------------------
**File structure**


::

    .
    |── cli.py          # Command-line interface entry point
    ├── docs            # Readthedocs documentation site
    ├── ml              # Directory containing the quantumAllostery modules
    │   ├── process     # Processes raw dynamics data
    │   ├── predict     # Machine learning analysis
    │   ├── manage      # File management functionality and routines
    │   ├── analyze     # Data analysis to combine process and plot routines
    │   └── plot        # Automated plotting and vizualization 
    └── ...



4 Documentation
---------------
Run the following commands to update the ReadTheDocs site:

::

    bash
    make clean
    make html



5 Developer guide
-----------------

GitHub refresher for those who would like to contribute
**Push new changes**

::
    
    bash
    git status
    git pull
    git add .
    git commit -m "Change a specific functionality"
    git push -u origin main



Copyright
---------

Copyright (c) 2023, David W. Kastner