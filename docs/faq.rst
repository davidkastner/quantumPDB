FAQ & Troubleshooting
=====================

Installation
------------

**Modeller fails to import or run.**
   Modeller requires a license key. Register at https://salilab.org/modeller/
   to obtain a free academic key, then set it as an environment variable:

   .. code-block:: bash

      export KEY_MODELLER="XXXX"

   Ensure you installed Modeller from the ``salilab`` conda channel
   (``conda install -c salilab modeller``). The ``environment.yml`` file
   handles this automatically.

**Multiwfn is not found.**
   Multiwfn must be installed separately. Download it from
   http://sobereva.com/multiwfn/ and ensure the executable is on your
   ``PATH``, or set the ``multiwfn_path`` parameter in your config file:

   .. code-block:: yaml

      multiwfn_path: /path/to/Multiwfn

   Multiwfn is only needed for the ``qp analyze`` stage.

Protoss
-------

**Protoss returns an error or times out.**
   Protoss is accessed via the ProteinsPlus web API and requires an internet
   connection. Common issues:

   - **Rate limiting:** The Protoss server may throttle users who submit too
     many requests in a short period. If you are processing a large batch,
     the server may reject some requests. QuantumPDB handles retries, but
     large batches may require running over multiple sessions.
   - **File size:** PDB files larger than 4 MB cannot be uploaded to Protoss.
   - **Server downtime:** The ProteinsPlus server may occasionally be
     unavailable. Check https://proteins.plus/ to verify the service is
     running.

**Protoss removes residues due to steric clashes.**
   Protoss may still drop clashing residues during protonation. Under the
   older ProteinsPlus API, QuantumPDB read the Protoss clash log and fed
   those residues back to Modeller for rebuilding (up to
   ``max_clash_refinement_iter`` times). ProteinsPlus API v2 no longer
   provides that log, so automatic clash remodeling is currently inactive
   and ``Protoss/{pdb}_log.txt`` is empty. Temporary workaround: compare
   the Modeller and Protoss PDBs for missing residues and inspect unusual
   conformations near the active site. A future release may restore this
   behavior via a local input/output residue diff.

Cluster Extraction
------------------

**The cluster is too large for my QM calculation.**
   Set the ``max_atom_count`` parameter to cap the cluster size:

   .. code-block:: yaml

      max_atom_count: 500

   The most distant residues are pruned until the atom count is below the
   threshold. You can also reduce ``number_of_spheres`` to limit the number
   of interaction layers.

**Important residues are missing from the cluster.**
   Several options:

   - Increase ``number_of_spheres`` to capture more of the environment.
   - Increase ``radius_of_first_sphere`` beyond the default 4.0 Å.
   - Use ``additional_ligands`` to force specific residues into the cluster:

     .. code-block:: yaml

        additional_ligands: [AKG, SIN]

     This accepts a list of three-letter residue codes. Residues matching
     these codes are included in the first sphere and protected from pruning.

**Multi-metallic centers are generating separate clusters instead of one.**
   Use merged center syntax in your CSV or increase ``merge_distance_cutoff``:

   .. code-block:: yaml

      # Automatic merging by distance
      center_residues: [FE2]
      merge_distance_cutoff: 4.0

   Or specify explicitly in a CSV ``center`` column::

      FE2_A5001-FE2_A5002

   See :doc:`input_formats` for the full center syntax reference.

Job Submission
--------------

**Jobs are not being submitted.**
   Check that both flags are set:

   .. code-block:: yaml

      create_jobs: true
      submit_jobs: true

   The ``submit`` command requires a CSV input file with ``oxidation`` and
   ``multiplicity`` columns so that the electronic state can be set correctly.

**Jobs are being resubmitted.**
   QuantumPDB creates a hidden ``.submit_record`` file in each job directory
   after submission. If this file exists, the job is skipped. To resubmit a
   job, set ``delete_queued: true`` in the analyze config to clear the records,
   then run ``qp submit`` again.

General
-------

**Can I use structures from MD simulations or generative models?**
   Yes. QuantumPDB works with any PDB-format file, including MD snapshots,
   cryo-EM structures, NMR ensembles, and generative model outputs. However,
   there is no trajectory parsing --- you must provide a single PDB frame.

**Which QM code should I use?**
   TeraChem (GPU-accelerated) is the only QM code currently supported by the
   job generator. The ``method`` and ``basis`` parameters follow TeraChem
   conventions. ORCA input writing is not implemented.

**No clusters were generated / the run exits asking for a center.**
   Provide either YAML ``center_residues`` or a CSV ``center`` column. If both
   are present, the CSV value wins for each PDB. Batch YAML lists such as
   ``[FE, FE2]`` are consumed **one item per PDB**, not as two centers inside
   one structure. For multiple centers in one structure, use fuzzy syntax
   (``FE_FE2``), CSV center definitions, or ``merge_distance_cutoff``.

**An mmCIF structure was skipped with an oversized warning.**
   Classic PDB conversion supports at most 99,999 atoms and 62 chains.
   Larger entries raise ``OversizedStructureError`` and are skipped in batch
   mode. See :doc:`input_formats`.

**``qp submit`` fails even though ``qp run`` worked.**
   ``submit`` requires ``input`` to be an existing CSV path (not a PDB ID or
   YAML list). Each cluster directory must contain exactly one ``.xyz`` file.
   Missing ``oxidation`` / ``multiplicity`` columns will produce incorrect
   electronic states.

**Multiwfn cannot find molden files.**
   Charge and dipole analysis look under ``{cluster}/{method}/scr/`` for
   ``*.molden`` outputs from completed TeraChem jobs. Confirm ``method``
   matches the submit config and that the calculation finished successfully.
   ``charge_scheme`` must be a single scheme name (not a comma-separated list).
