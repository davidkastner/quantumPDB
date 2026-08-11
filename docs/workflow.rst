Pipeline Overview
=================

QuantumPDB processes protein structures through five sequential stages. The
first three stages are handled by ``qp run``, job management by ``qp submit``,
and post-processing by ``qp analyze``.

.. code-block:: text

   PDB file ──► Structure Prep ──► Protonation ──► Cluster Extraction ──► Job Management ──► Analysis
                 (qp run)          (qp run)         (qp run)              (qp submit)        (qp analyze)

Stage 1: Structure Preparation
------------------------------

**Module:** ``qp.structure``

Raw PDB files from the Protein Data Bank often have missing atoms, residues,
or entire loops due to crystallographic disorder. QuantumPDB uses
`Modeller <https://salilab.org/modeller/>`_ to rebuild these regions.

**What happens:**

1. The structure is fetched from the RCSB or read from a local path. Classic
   PDB is preferred; if only mmCIF is available (or a local ``.cif`` /
   ``.mmcif`` is supplied), QuantumPDB converts it to classic PDB and may
   write a ``*_mmcif_remap.json`` sidecar. Oversized structures (>99999 atoms
   or >62 chains) are skipped in batch mode. See :doc:`input_formats`.
2. REMARK 465 (missing residues) and REMARK 470 (missing atoms) records are
   parsed to identify gaps.
3. Unresolved terminal residues are trimmed (they cannot be reliably modeled).
4. Modeller builds the missing regions using homology modeling against the
   known portions of the structure.

**Output:** ``{pdb}_modeller.pdb``

**Key parameter:** ``optimize_select_residues`` controls whether Modeller
refines only missing residues (1, default), all residues (2), or none (0).

Stage 2: Protonation
--------------------

**Module:** ``qp.protonate``

Crystal structures lack hydrogen atoms. QuantumPDB submits the prepared
structure to the `Protoss <https://proteins.plus/>`_ web server, which assigns
hydrogen positions and resolves alternate conformations.

**What happens:**

1. Partial occupancy is resolved by selecting a self-consistent coordinate set,
   prioritizing center residues and canonical amino acids.
2. The structure is uploaded to Protoss for protonation.
3. Historically, the Protoss clash log was checked for steric clashes so that
   problematic residues could be removed and rebuilt by Modeller (up to
   ``max_clash_refinement_iter`` times). ProteinsPlus API v2 no longer
   exposes this log, so that feedback loop is currently inactive; the
   ``*_log.txt`` file is written empty as a placeholder. If residues appear
   to be missing after protonation, inspect the structure manually.
4. Metalloenzyme-specific corrections are applied:

   - Histidine ring flips for metal coordination
   - Backbone nitrogen deprotonation near metals (e.g., nitrile hydratase)
   - Removal of hydrogen atoms too close to metal centers
   - Correction of Protoss hydroxylamine artifacts

5. Ligand charges and spins are computed from Protoss SDF output.

**Output:** ``Protoss/{pdb}_protoss.pdb``, ``Protoss/{pdb}_ligands.sdf``,
``charge.csv``

Stage 3: Cluster Extraction
----------------------------

**Module:** ``qp.cluster``

This is the core algorithmic stage. QuantumPDB constructs hierarchical
interaction spheres around the specified center residue(s) to carve out a QM
cluster model.

**What happens:**

1. Center residues are identified by name, chain, and residue number
   (see :doc:`input_formats` for syntax).
2. The **first sphere** is built using a distance cutoff (default 4.0 Å),
   with chemistry-aware filtering that excludes non-coordinating atoms.
3. **Second and higher spheres** are built using Voronoi tessellation ---
   residues are included based on topological contact rather than distance,
   which captures non-spherical active-site environments.
4. Voronoi cells in sparse regions are regularized using dummy atoms (or
   other smoothing methods).
5. If ``max_atom_count`` is set, the most distant residues are pruned.
6. Chain breaks are capped (hydrogens or ACE/NME groups).
7. Net charges are computed from protonation states, ionizable side chains,
   and user-provided metal oxidation states.

**Output:** Per-center directories named by ``metal_id`` containing cumulative
``0.pdb``, ``1.pdb``, ... sphere models plus combined ``{metal_id}.pdb`` /
``.xyz``, along with ``charge.csv`` and ``count.csv``. See :doc:`output` and
:doc:`cluster_models`.

Cluster extraction is performed **per matched center set**, not strictly per
chain. A multi-chain protein produces a separate cluster for each matching
center (or merged center group). Matching centers on different chains remain
separate unless you merge them with dash syntax or ``merge_distance_cutoff``.

Stage 4: Job Management
------------------------

**Module:** ``qp.manager``

QuantumPDB generates ready-to-run QM input files and manages job submission.
See :doc:`qm_jobs` for a full walkthrough.

**What happens:**

1. For each cluster, a TeraChem input file (``qmscript.in``) is generated with
   the specified method, basis set, and electronic state. The current job
   generator is TeraChem-only.
2. The electronic state (charge and spin multiplicity) is computed from:

   - The cluster's ``charge.csv`` and ``spin.csv`` (from Stages 2--3)
   - The user-provided ``oxidation`` and ``multiplicity`` columns in the
     input CSV

3. If charge embedding is enabled, MM point charges are placed around the
   QM cluster within the cutoff distance. Selection is residue-based: if
   any atom of a residue is within the cutoff, all atoms of that residue
   are included to ensure integer per-residue charge sums. By default
   AMBER ff14SB charges are used, but users can supply a custom charge
   file via ``charge_embedding_charges``.
4. Scheduler scripts (SLURM or SGE) are generated alongside the QM input.
5. Jobs are submitted up to the ``job_count`` limit, with a ``.submit_record``
   file preventing duplicate submissions.

**Output:** QM input files, scheduler scripts, and ``.submit_record`` markers.

Stage 5: Analysis
-----------------

**Module:** ``qp.analyze``

Post-processing extracts electronic properties from completed QM calculations.
See :doc:`analysis` for details.

**What happens:**

1. **Job checkup** classifies all jobs by status (done, running, queued,
   backlog, error) and generates summary CSVs and plots.
2. **Charge analysis** uses `Multiwfn <http://sobereva.com/multiwfn/>`_ to
   compute partial atomic charges. Available schemes:

   - Hirshfeld, Voronoi (deformation density), Mulliken, ADCH, Hirshfeld-I, CM5

3. **Dipole calculation** computes substrate dipole moments using the center
   of mass as the reference point.

.. note::

   The "Voronoi" charge scheme refers to Voronoi deformation density charges
   computed by Multiwfn. This is unrelated to the Voronoi tessellation used for
   cluster construction in Stage 3.

**Output:** Charge data files, dipole results, and job status summaries in
``checkup/``.
