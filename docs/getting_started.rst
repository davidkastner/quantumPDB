Getting Started
===============

This page covers how to install QuantumPDB and its dependencies.

Prerequisites
-------------

- **Python 3.8+**
- **Conda** (recommended for managing the environment and installing Modeller)
- **Internet access** (required for Protoss protonation via the ProteinsPlus web server)

Installation
------------

1. **Clone the repository:**

   .. code-block:: bash

      git clone https://github.com/davidkastner/quantumPDB.git
      cd quantumPDB

2. **Create the conda environment:**

   The provided ``environment.yml`` installs all required dependencies, including
   Modeller from the ``salilab`` channel:

   .. code-block:: bash

      conda env create -f environment.yml
      conda activate qp

3. **Install the package:**

   .. code-block:: bash

      pip install -e .

   This performs a development install so that the ``qp`` command is available
   on your path.

Modeller Setup
--------------

`Modeller <https://salilab.org/modeller/>`_ is used for modeling missing atoms,
residues, and loops (Stage 1). It is installed automatically via the conda
environment, but requires a **free academic license key**.

1. Register at https://salilab.org/modeller/ to obtain a license key.
2. Set the key as an environment variable:

   .. code-block:: bash

      export KEY_MODELLER="XXXX"

   Replace ``XXXX`` with your key. Add this line to your ``~/.bashrc`` or
   ``~/.zshrc`` so it persists across sessions.

Protoss
-------

`Protoss <https://proteins.plus/>`_ assigns protonation states and resolves
alternate conformations (Stage 2). QuantumPDB accesses Protoss through the
ProteinsPlus web API, so no local installation is needed --- just an internet
connection.

.. note::

   The Protoss web server may throttle users who submit too many requests in a
   short period. If you are processing a large batch, expect some requests to be
   rate-limited. QuantumPDB handles retries automatically, but very large batches
   may require patience.

Optional Dependencies
---------------------

These are only needed for specific pipeline stages and are **not** required to
generate cluster models.

**TeraChem** (Stage 4 --- QM calculations)
   TeraChem is the primary GPU-accelerated quantum chemistry code supported by
   QuantumPDB. A license is required. See https://www.petachem.com/.

**ORCA** (Stage 4 --- QM calculations)
   ORCA is supported as a CPU-based alternative. See https://orcaforum.kofo.mpg.de/.

**Multiwfn** (Stage 5 --- post-processing)
   `Multiwfn <http://sobereva.com/multiwfn/>`_ is required only for charge
   scheme analysis and dipole moment calculations in the analyze stage. Install
   it separately and ensure the ``Multiwfn`` executable is on your ``PATH``, or
   set the ``multiwfn_path`` parameter in your config file.

Verifying the Installation
--------------------------

After installation, verify that the ``qp`` command is available:

.. code-block:: bash

   qp --help

This should display the three available commands: ``run``, ``submit``, and ``analyze``.
