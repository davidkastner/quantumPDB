Quickstart
==========

This guide walks through a minimal example: generating a QM cluster model for
**taurine dioxygenase (TauD)**, a mononuclear non-heme iron enzyme (PDB: `1OS7
<https://www.rcsb.org/structure/1OS7>`_).

The same configuration also lives under ``example/basics/`` so you can run it
directly from a clone of the repository.

1. Create the Configuration File
---------------------------------

Create a file called ``config.yaml`` (or use ``example/basics/cluster_only.yaml``):

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

Or, from the repository root:

.. code-block:: bash

   qp run -c example/basics/cluster_only.yaml

QuantumPDB processes the structure through three stages:

1. **Structure preparation** --- Modeller rebuilds missing atoms, residues, and
   loops, producing ``1os7_modeller.pdb``
2. **Protonation** --- Protoss assigns hydrogen positions and resolves alternate
   conformations, producing ``Protoss/1os7_protoss.pdb``
3. **Cluster extraction** --- Voronoi-based spheres are constructed around each
   ``FE`` center, producing PDB and XYZ cluster files

3. Expected Output
-------------------

After a successful run, the output directory will contain something like:

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
       │   └── 1os7_log.txt          # Protoss clash log (empty under API v2)
       ├── charge.csv                # Per-cluster / per-sphere charges
       ├── count.csv                 # Residue counts per sphere
       └── A501/                     # Cluster directory named by metal_id
           ├── 0.pdb                 # Center residue(s)
           ├── 1.pdb                 # First interaction sphere
           ├── 2.pdb                 # Second sphere (includes sphere 1)
           ├── A501.pdb              # Combined cluster PDB
           └── A501.xyz              # Combined cluster XYZ

Each cluster directory is named by the matched center residue ID
(``metal_id``), for example ``A501`` for the iron at position 501 on chain A.
Numbered ``{i}.pdb`` files are cumulative sphere models; the combined
``{metal_id}.pdb`` / ``.xyz`` files are what ``qp submit`` uses.

For ground-truth examples of this layout, see the regression samples under
``qp/tests/samples/`` (for example ``1lm6/A204/``).

4. Inspecting the Results
--------------------------

Open the cluster PDB files in a molecular viewer such as PyMOL or VMD to
verify that the active site is correctly captured. The ``charge.csv`` file
records net charge per cluster and sphere, and ``count.csv`` lists residue
composition in each interaction sphere. See :doc:`output` for exact formats.

Next Steps
----------

- **Multiple structures:** Use a CSV input file to process many PDBs at once
  (see :doc:`input_formats`)
- **Cluster options:** Tune spheres, smoothing, capping, and naming
  (see :doc:`cluster_models`)
- **QM calculations:** Create TeraChem job files with ``qp submit``
  (see :doc:`qm_jobs`)
- **Analysis:** Run checkup and Multiwfn post-processing with ``qp analyze``
  (see :doc:`analysis`)
- **Configuration:** Fine-tune parameters for your system
  (see :doc:`configuration`)
