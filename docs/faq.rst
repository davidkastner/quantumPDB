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
   This is expected behavior. Protoss identifies steric clashes in the
   structure, and QuantumPDB automatically feeds these back to Modeller for
   rebuilding. This loop repeats up to ``max_clash_refinement_iter`` times
   (default: 5). If clashes persist, check the structure for unusual
   conformations near the active site.

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
   TeraChem (GPU-accelerated) is the primary target and best tested. ORCA
   (CPU-based) is supported as an alternative. The ``method`` and ``basis``
   parameters follow TeraChem conventions by default.
