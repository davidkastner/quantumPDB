Input Formats
=============

QuantumPDB accepts several input formats for specifying which protein structures
to process. All inputs are specified in the YAML configuration file.

Single PDB Code
----------------

Process a single structure by its PDB accession code:

.. code-block:: yaml

   input: 1OS7

The structure is automatically downloaded from the RCSB PDB.

Local PDB File
--------------

Process a local PDB file by providing the path:

.. code-block:: yaml

   input: /path/to/structure.pdb

List of PDB Codes
------------------

Process multiple structures by listing PDB codes:

.. code-block:: yaml

   input: [1OS7, 1FYZ, 1PHM]

CSV Input File
--------------

For high-throughput processing, use a CSV file. This is the recommended format
when processing many structures or when per-PDB parameters are needed.

.. code-block:: yaml

   input: proteins.csv

**Required column:**

- ``pdb_id`` --- PDB accession code or path to a local PDB file

**Optional columns:**

- ``center`` --- Center residue definition for this specific PDB (overrides
  the ``center_residues`` parameter in the YAML config)
- ``oxidation`` --- Metal oxidation state (required for ``qp submit``)
- ``multiplicity`` --- Spin multiplicity (required for ``qp submit``)

**Example CSV:**

.. code-block:: text

   pdb_id,center,oxidation,multiplicity
   1OS7,FE,3,6
   1FYZ,FE2_A5001-FE2_A5002,6,11
   1PHM,CU_A357-CU_A358,4,3

Center Residue Syntax
---------------------

The center residue defines which atoms are at the core of the QM cluster model.
The syntax differs between YAML and CSV formats.

YAML Format (``center_residues``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use brackets in the YAML config for **general mode**, which matches all HETATM
records with the given residue name(s):

.. code-block:: yaml

   # Single residue type
   center_residues: [FE]

   # Multiple residue types (generates separate clusters for each)
   center_residues: [FE, FE2]

General mode only matches HETATM records, not protein ATOM records. This
prevents generating a cluster for every instance of a common amino acid.

CSV Format (``center`` column)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CSV ``center`` column supports **specific mode** with exact residue
identification:

.. code-block:: text

   # Single residue (name only, general mode)
   FE

   # Specific residue (name_chain+number)
   FE_A199

   # Merged centers (dash-separated)
   CU_A357-CU_A358

   # Complex merged center
   FE_A155-OXY_A555-HEM_A155-HIS_A93

The format for specific residues is ``RESNAME_CHAINID``, where ``RESNAME`` is
the three-letter residue code and ``CHAINID`` is the chain letter followed by
the residue number (e.g., ``FE_A199`` for iron at position 199 on chain A).

**Merged centers** treat multiple residues as a single center. This is used for:

- Multi-metallic active sites (e.g., ``CU_A357-CU_A358`` for dicopper)
- Heme groups (e.g., ``FE_A155-HEM_A155-OXY_A555-HIS_A93``)
- Oligomeric substrates where you want the cluster centered on a subset

.. important::

   When both ``center_residues`` (YAML) and a ``center`` column (CSV) are
   provided, the **CSV value takes priority** for each PDB.

Merging Nearby Centers
^^^^^^^^^^^^^^^^^^^^^^

For multi-metallic sites where centers are close together, use the
``merge_distance_cutoff`` parameter to automatically merge centers within a
given distance:

.. code-block:: yaml

   center_residues: [FE2]
   merge_distance_cutoff: 4.0

This merges any ``FE2`` atoms within 4.0 Å of each other into a single center.
Set the cutoff slightly larger than the metal--metal distance. Examples:

- **MMO diiron** (Fe--Fe ~3.4 Å): ``merge_distance_cutoff: 4.0``
- **PHM dicopper** (Cu--Cu ~11 Å): ``merge_distance_cutoff: 12.0``

Configuration File Structure
-----------------------------

The YAML configuration file has three sections corresponding to the three CLI
commands. You only need to include the parameters relevant to the command you
are running.

.. code-block:: yaml

   #######
   # RUN #
   #######
   input: proteins.csv
   output_dir: output
   modeller: true
   protoss: true
   coordination: true
   center_residues: [FE]
   number_of_spheres: 2

   ##########
   # SUBMIT #
   ##########
   method: wpbeh
   basis: lacvps_ecp
   create_jobs: true

   ###########
   # ANALYZE #
   ###########
   charge_scheme: Hirshfeld
   calc_charge_schemes: true

See :doc:`configuration` for the complete parameter reference.
