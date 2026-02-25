Configuration Reference
=======================

All QuantumPDB parameters are set in a YAML configuration file passed to
the CLI with ``-c config.yaml``. Parameters are organized into three
sections corresponding to the three CLI commands.

Only include the parameters you need --- all parameters have sensible defaults.

RUN Parameters
--------------

These parameters control ``qp run`` (structure preparation, protonation, and
cluster extraction).

**Pipeline control:**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``input``
     - ``[]``
     - Path to input CSV, single PDB code, local PDB path, or list of PDB
       codes (e.g., ``[1OS7, 1FYZ]``).
   * - ``output_dir``
     - ``''``
     - Directory where all output files are saved.
   * - ``modeller``
     - ``false``
     - Enable Modeller to rebuild missing atoms, residues, and loops.
   * - ``protoss``
     - ``false``
     - Enable Protoss to assign protonation states via the web server.
   * - ``coordination``
     - ``false``
     - Enable QM cluster model generation.
   * - ``skip``
     - ``'all'``
     - Skip existing results. ``'modeller'`` skips only Modeller,
       ``'protoss'`` skips only Protoss, ``'all'`` skips both if their
       output files already exist. To force re-processing, delete the
       output files.
   * - ``optimize_select_residues``
     - ``1``
     - Modeller refinement level: ``0`` = none, ``1`` = missing residues
       only, ``2`` = all residues.
   * - ``max_clash_refinement_iter``
     - ``5``
     - Maximum iterations of the Modeller--Protoss feedback loop for
       resolving steric clashes.

**Cluster model parameters:**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``center_residues``
     - ``[]``
     - Three-letter residue code(s) for the cluster center. Use YAML list
       syntax: ``[FE]``, ``[FE, FE2]``. See :doc:`input_formats` for the
       full center definition syntax.
   * - ``additional_ligands``
     - ``[]``
     - List of three-letter residue codes to force into the first interaction
       sphere (e.g., ``[AKG, SIN]``). These residues are also protected from
       pruning.
   * - ``number_of_spheres``
     - ``2``
     - Number of interaction spheres to construct. More spheres capture more
       of the protein environment at increased computational cost.
   * - ``radius_of_first_sphere``
     - ``4.0``
     - Distance cutoff in angstroms for the first (distance-based) sphere.
       The default of 4.0 Å works well for most systems. Increase for
       large ligands or weak interactions.
   * - ``include_ligands``
     - ``2``
     - Ligand inclusion level: ``0`` = first sphere only, ``1`` = non-water
       ligands, ``2`` = all ligands and waters.
   * - ``max_atom_count``
     - ``None``
     - Maximum number of atoms in the QM cluster. When set, the most distant
       residues are pruned until the count is below the threshold. Capping
       atoms are not counted toward this limit.
   * - ``merge_distance_cutoff``
     - ``0.0``
     - Distance in angstroms for merging nearby center residues into a single
       center. Set slightly larger than the metal--metal distance. Examples:
       4.0 for MMO diiron (~3.4 Å), 12.0 for PHM dicopper (~11 Å).
   * - ``capping_method``
     - ``1``
     - How to cap chain breaks: ``0`` = none, ``1`` = hydrogen caps,
       ``2`` = ACE/NME caps. Use hydrogen caps for single-point energy
       calculations and ACE/NME for geometry optimizations.
   * - ``smoothing_method``
     - ``2``
     - Voronoi regularization method: ``0`` = box plot, ``1`` = DBSCAN,
       ``2`` = dummy atom (recommended), ``3`` = none.

**Output control:**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``compute_charges``
     - ``true``
     - Compute per-residue charges from protonation states and ligand
       properties.
   * - ``count_residues``
     - ``true``
     - Count the number of residues in each interaction sphere and write
       ``count.csv``.
   * - ``write_xyz``
     - ``true``
     - Write XYZ coordinate files alongside PDB cluster files.


SUBMIT Parameters
-----------------

These parameters control ``qp submit`` (QM job creation and submission).

**QM calculation settings:**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``optimization``
     - ``false``
     - ``true`` for geometry optimization, ``false`` for single-point energy.
   * - ``method``
     - ``'wpbeh'``
     - DFT functional. TeraChem supports: wpbeh, b3lyp, ub3lyp, pbe0, etc.
   * - ``basis``
     - ``'lacvps_ecp'``
     - Basis set. Common choices: ``lacvps_ecp``, ``6-31g*``, ``def2-svp``.
   * - ``guess``
     - ``'generate'``
     - Initial wavefunction guess method.
   * - ``gpus``
     - ``1``
     - Number of GPUs for TeraChem calculations.
   * - ``memory``
     - ``'8G'``
     - Memory allocation per job.
   * - ``dielectric``
     - ``10``
     - Dielectric constant for implicit solvent (PCM). Use ``10`` for
       protein interior, ``78.4`` for aqueous solvent.
   * - ``use_implicit_solvent``
     - ``not charge_embedding``
     - Enable PCM implicit solvent (COSMO). Defaults to ``true`` when
       ``charge_embedding`` is ``false``, and ``false`` when
       ``charge_embedding`` is ``true``, preserving backward-compatible
       behavior. Set explicitly to ``true`` alongside
       ``charge_embedding: true`` to use both MM point charges and
       implicit solvent simultaneously.
   * - ``charge_embedding``
     - ``false``
     - Include MM point charges around the QM cluster. By default uses
       ff14SB partial charges (also TIP3P for water, standard charges for
       common ions). The charge source can be customized with
       ``charge_embedding_charges``.
   * - ``charge_embedding_cutoff``
     - ``20``
     - Distance cutoff in angstroms for MM point charges, measured from the
       QM cluster centroid. If any atom of a residue falls within this
       distance, all atoms of that residue are included to preserve
       integer per-residue charge contributions.
   * - ``charge_embedding_charges``
     - ``null``
     - Path to a JSON file with custom partial charges, keyed by residue
       name then atom name (e.g., ``{"ALA": {"N": -0.4157, ...}}``).
       When ``null``, uses the built-in AMBER ff14SB charges.
   * - ``scheduler``
     - ``'slurm'``
     - Job scheduler: ``'slurm'`` or ``'sge'``.
   * - ``pcm_radii_file``
     - ``'pcm_radii'``
     - Path to a custom PCM radii file (TeraChem-specific).

**Job management:**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``create_jobs``
     - ``false``
     - Generate QM input files for each cluster.
   * - ``submit_jobs``
     - ``false``
     - Submit generated jobs to the scheduler.
   * - ``job_count``
     - ``80``
     - Maximum number of concurrent jobs on the scheduler.


ANALYZE Parameters
------------------

These parameters control ``qp analyze`` (job monitoring and post-processing).

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``job_checkup``
     - ``false``
     - Check the status of all QM jobs and generate summary reports.
   * - ``delete_queued``
     - ``false``
     - Delete ``.submit_record`` files for unfinished jobs, allowing them
       to be resubmitted.
   * - ``multiwfn_path``
     - ``'Multiwfn'``
     - Path to the Multiwfn executable. If Multiwfn is on your ``PATH``,
       the default works.
   * - ``charge_scheme``
     - ``'Hirshfeld'``
     - Charge scheme for partial charge analysis. Options: ``Hirshfeld``,
       ``Voronoi``, ``Mulliken``, ``ADCH``, ``Hirshfeld-I``, ``CM5``.
   * - ``calc_charge_schemes``
     - ``false``
     - Run Multiwfn to compute partial atomic charges.
   * - ``calc_dipole``
     - ``false``
     - Run Multiwfn to compute substrate dipole moments using the center
       of mass as the reference point.

Minimal Configuration Examples
------------------------------

**Cluster generation only** (most common starting point):

.. code-block:: yaml

   input: 1OS7
   output_dir: output
   modeller: true
   protoss: true
   coordination: true
   center_residues: [FE]

**Full pipeline** (generation + QM + analysis):

.. code-block:: yaml

   # RUN
   input: proteins.csv
   output_dir: dataset
   modeller: true
   protoss: true
   coordination: true
   center_residues: [FE]
   number_of_spheres: 2

   # SUBMIT
   method: wpbeh
   basis: lacvps_ecp
   dielectric: 10
   create_jobs: true
   submit_jobs: true

   # ANALYZE
   charge_scheme: Hirshfeld
   calc_charge_schemes: true
