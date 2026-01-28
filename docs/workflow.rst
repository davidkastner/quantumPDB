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

1. The PDB file is fetched from the RCSB (or read from a local path).
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
3. The Protoss log is checked for steric clashes. If clashes are found, the
   problematic residues are removed and rebuilt by Modeller. This feedback
   loop repeats up to ``max_clash_refinement_iter`` times (default: 5).
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

**Output:** Per-center directories containing ``cluster.pdb`` and
``cluster.xyz`` for each sphere level, plus ``charge.csv`` and ``count.csv``.

Cluster extraction is performed **per chain** --- a multi-chain protein produces
separate clusters for each chain that contains a matching center residue.

Stage 4: Job Management
------------------------

**Module:** ``qp.manager``

QuantumPDB generates ready-to-run QM input files and manages job submission.

**What happens:**

1. For each cluster, a QM input file is generated for TeraChem (primary) or
   ORCA (alternative) with the specified method, basis set, and electronic
   state.
2. The electronic state (charge and spin multiplicity) is computed from:

   - The cluster's ``charge.csv`` and ``spin.csv`` (from Stages 2--3)
   - The user-provided ``oxidation`` and ``multiplicity`` columns in the
     input CSV

3. If charge embedding is enabled, MM point charges from an ff14SB force
   field are placed around the QM cluster within the cutoff distance.
4. Scheduler scripts (SLURM or SGE) are generated alongside the QM input.
5. Jobs are submitted up to the ``job_count`` limit, with a ``.submit_record``
   file preventing duplicate submissions.

**Output:** QM input files, scheduler scripts, and ``.submit_record`` markers.

Stage 5: Analysis
-----------------

**Module:** ``qp.analyze``

Post-processing extracts electronic properties from completed QM calculations.

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
