CLI Reference
=============

QuantumPDB provides three commands, each corresponding to a stage of the
pipeline. All commands require a YAML configuration file passed with the
``-c`` / ``--config`` flag.

qp run
------

Generate QM cluster models from protein structures.

.. code-block:: bash

   qp run -c config.yaml

This command executes up to three stages depending on the configuration:

1. **Structure preparation** (``modeller: true``) --- Models missing atoms,
   residues, and loops using Modeller.
2. **Protonation** (``protoss: true``) --- Assigns hydrogen positions and
   resolves alternate conformations via the Protoss web server.
3. **Cluster extraction** (``coordination: true``) --- Constructs Voronoi-based
   interaction spheres around the specified center residues.

Each stage can be enabled or disabled independently. The ``skip`` parameter
controls whether existing output files are reused (``'all'``) or specific
stages are re-run.

**Example --- single PDB:**

.. code-block:: yaml

   input: 1OS7
   output_dir: output
   modeller: true
   protoss: true
   coordination: true
   center_residues: [FE]

**Example --- batch from CSV:**

.. code-block:: yaml

   input: proteins.csv
   output_dir: dataset
   modeller: true
   protoss: true
   coordination: true
   center_residues: [FE]

See :doc:`input_formats` for CSV column requirements.

qp submit
---------

Create and submit QM calculation job files.

.. code-block:: bash

   qp submit -c config.yaml

This command has two sub-steps controlled by separate flags:

1. **Create job files** (``create_jobs: true``) --- Generates TeraChem or ORCA
   input scripts for each cluster, including scheduler submission scripts
   (SLURM or SGE).
2. **Submit jobs** (``submit_jobs: true``) --- Submits the generated jobs to the
   scheduler, up to ``job_count`` concurrent jobs.

**Example --- single-point calculations with TeraChem:**

.. code-block:: yaml

   input: proteins.csv
   output_dir: dataset
   optimization: false
   method: wpbeh
   basis: lacvps_ecp
   gpus: 1
   memory: 8G
   dielectric: 10
   scheduler: slurm
   create_jobs: true
   submit_jobs: true
   job_count: 80

**Example --- geometry optimization with charge embedding:**

.. code-block:: yaml

   input: proteins.csv
   output_dir: dataset
   optimization: true
   method: wpbeh
   basis: lacvps_ecp
   charge_embedding: true
   charge_embedding_cutoff: 20
   capping_method: 2
   gpus: 1
   memory: 8G
   dielectric: 10
   scheduler: slurm
   create_jobs: true

.. note::

   The ``submit`` command requires a CSV input file with ``oxidation`` and
   ``multiplicity`` columns so that the correct electronic state can be set in
   the QM input files.

qp analyze
----------

Check job status and run post-processing analysis.

.. code-block:: bash

   qp analyze -c config.yaml

This command supports two independent functions:

1. **Job checkup** (``job_checkup: true``) --- Classifies all QM jobs by
   status (done, running, queued, backlog, error) and generates summary plots
   and CSVs in a ``checkup/`` directory.
2. **Charge analysis** (``calc_charge_schemes: true``) --- Runs Multiwfn to
   compute partial atomic charges using a specified scheme.
3. **Dipole calculation** (``calc_dipole: true``) --- Computes substrate dipole
   moments using the center of mass as the reference point.

**Example --- check job status:**

.. code-block:: yaml

   input: proteins.csv
   output_dir: dataset
   method: wpbeh
   job_checkup: true

**Example --- compute Hirshfeld charges:**

.. code-block:: yaml

   input: proteins.csv
   output_dir: dataset
   method: wpbeh
   charge_scheme: Hirshfeld
   calc_charge_schemes: true

.. note::

   Charge analysis and dipole calculations require `Multiwfn
   <http://sobereva.com/multiwfn/>`_ to be installed. Set the
   ``multiwfn_path`` parameter if the executable is not on your ``PATH``.
