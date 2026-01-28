Quickstart
==========

This guide walks through a minimal example: generating a QM cluster model for
**taurine dioxygenase (TauD)**, a mononuclear non-heme iron enzyme (PDB: `1OS7
<https://www.rcsb.org/structure/1OS7>`_).

1. Create the Configuration File
---------------------------------

Create a file called ``config.yaml``:

.. code-block:: yaml

   # RUN
   input: 1OS7
   output_dir: output
   modeller: true
   protoss: true
   coordination: true
   skip: all
   center_residues: [FE]
   number_of_spheres: 2
   smoothing_method: 2
   capping_method: 1

This tells QuantumPDB to:

- Download 1OS7 from the RCSB PDB
- Model any missing atoms or residues with Modeller
- Assign protonation states with Protoss
- Extract cluster models centered on all iron (``FE``) atoms
- Build 2 interaction spheres using Voronoi tessellation with dummy atom smoothing
- Cap chain breaks with hydrogens

2. Run the Pipeline
--------------------

.. code-block:: bash

   qp run -c config.yaml

QuantumPDB processes the structure through three stages:

1. **Structure preparation** --- Modeller rebuilds missing atoms, residues, and
   loops, producing ``1os7_modeller.pdb``
2. **Protonation** --- Protoss assigns hydrogen positions and resolves alternate
   conformations, producing ``Protoss/1os7_protoss.pdb``
3. **Cluster extraction** --- Voronoi-based spheres are constructed around each
   ``FE`` center, producing PDB and XYZ cluster files

3. Expected Output
-------------------

After a successful run, the output directory will contain:

.. code-block:: text

   output/
   └── 1os7/
       ├── 1os7.pdb                  # Original PDB
       ├── 1os7.ali                  # Modeller alignment
       ├── 1os7_modeller.pdb         # Rebuilt structure
       ├── Protoss/
       │   ├── 1os7_protoss.pdb      # Protonated structure
       │   ├── 1os7_protoss_orig.pdb # Pre-active-site-fix copy
       │   ├── 1os7_ligands.sdf      # Ligand structures (SDF)
       │   └── 1os7_log.txt          # Protoss log
       ├── charge.csv                # Per-residue charges
       ├── count.csv                 # Residue counts per sphere
       └── FE_A501/                  # Cluster directory (one per center)
           ├── 1/                    # Sphere 1
           │   ├── cluster.pdb
           │   └── cluster.xyz
           └── 2/                    # Sphere 2 (includes sphere 1)
               ├── cluster.pdb
               └── cluster.xyz

Each cluster directory is named by the center residue, chain, and residue
number (e.g., ``FE_A501`` for iron at position 501 on chain A). Numbered
subdirectories contain progressively larger clusters.

4. Inspecting the Results
--------------------------

Open the cluster PDB files in a molecular viewer such as PyMOL or VMD to
verify that the active site is correctly captured. The ``charge.csv`` file
records the computed net charge for each residue, and ``count.csv`` lists
the number of residues in each interaction sphere.

Next Steps
----------

- **Multiple structures:** Use a CSV input file to process many PDBs at once
  (see :doc:`input_formats`)
- **QM calculations:** Set up TeraChem or ORCA job files with ``qp submit``
  (see :doc:`cli`)
- **Analysis:** Run charge scheme analysis with ``qp analyze``
  (see :doc:`cli`)
- **Configuration:** Fine-tune parameters for your system
  (see :doc:`configuration`)
