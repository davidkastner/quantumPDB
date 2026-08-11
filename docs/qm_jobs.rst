QM Jobs
=======

``qp submit`` turns finished cluster models into TeraChem inputs and optional
scheduler submissions. The implementation is in ``qp.manager``.

Prerequisites
-------------

1. Clusters already exist under ``output_dir`` from ``qp run``.
2. ``input`` is a path to an **existing CSV file** (not a PDB ID or YAML list).
3. The CSV includes ``oxidation`` and ``multiplicity`` columns.
4. Each cluster directory contains exactly one ``.xyz`` file (usually
   ``{metal_id}.xyz``).

Electronic State
----------------

For each cluster QuantumPDB computes:

- **Cluster charge** from ``charge.csv`` (amino-acid / sphere contribution plus
  ligand formal charges present in the spheres)
- **Total charge** as cluster charge + CSV ``oxidation``
- **Multiplicity** as CSV ``multiplicity`` plus any radical ligand spin from
  ``spin.csv`` (for example NO or O2)

Incorrect oxidation / multiplicity values are a common source of failed or
meaningless QM calculations.

Generated Files
---------------

With ``create_jobs: true``, each cluster gets a method subdirectory:

.. code-block:: text

   {metal_id}/{method}/
   ├── {metal_id}.xyz
   ├── qmscript.in
   ├── jobscript.sh
   └── ptchrges.xyz          # only if charge_embedding: true

- ``qmscript.in`` --- TeraChem input written by ``job_scripts.write_qm``
- ``jobscript.sh`` --- SLURM or SGE wrapper from ``scheduler``
- ``ptchrges.xyz`` --- MM embedding charges

The current generator is **TeraChem-only**. ORCA ``.inp`` writing is not
implemented.

Solvent and Embedding
---------------------

- ``use_implicit_solvent`` defaults to the opposite of ``charge_embedding`` for
  backward compatibility, but both may be enabled together.
- ``dielectric`` and ``pcm_radii_file`` configure PCM / COSMO.
- ``charge_embedding_cutoff`` is residue-based: if any atom of a residue is
  inside the cutoff from the QM centroid, the whole residue is included.
- Default charges come from AMBER ff14SB (plus common ions / TIP3P water).
  Supply ``charge_embedding_charges`` JSON of the form
  ``{"RES": {"ATOM": q, ...}, ...}`` to override.

Submission Control
------------------

- ``submit_jobs: true`` submits via ``sbatch`` or ``qsub``.
- ``job_count`` caps how many jobs are submitted in one manage pass.
- A hidden ``.submit_record`` marker prevents duplicate submission.
- To clear unfinished markers and allow resubmission, run ``qp analyze`` with
  ``delete_queued: true``, then submit again.

.. note::

   Generated SLURM / SGE scripts include site-specific partitions, modules, and
   sleep / queue settings. Edit ``qp.manager.job_scripts`` or the produced
   ``jobscript.sh`` files for your local cluster.

Related Pages
-------------

- :doc:`cli` --- command examples
- :doc:`configuration` --- submit keys
- :doc:`analysis` --- checkup and resubmission helpers
- :doc:`faq` --- common submission failures
