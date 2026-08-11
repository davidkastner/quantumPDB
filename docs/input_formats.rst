Input Formats
=============

QuantumPDB accepts several input formats for specifying which protein structures
to process. All inputs are specified in the YAML configuration file via the
``input`` key.

Single PDB Code
---------------

Process a single structure by its PDB accession code:

.. code-block:: yaml

   input: 1OS7

The structure is automatically downloaded from the RCSB PDB. If a classic PDB
file is unavailable, QuantumPDB falls back to mmCIF and converts it (see
:ref:`mmcif-support`).

Local PDB File
--------------

Process a local PDB file by providing the path:

.. code-block:: yaml

   input: /path/to/structure.pdb

Local mmCIF File
----------------

Local ``.cif`` / ``.mmcif`` paths are also accepted:

.. code-block:: yaml

   input: /path/to/structure.cif

QuantumPDB converts the mmCIF to a classic PDB under ``output_dir`` before
Modeller, Protoss, and clustering. See :ref:`mmcif-support`.

List of PDB Codes
-----------------

Process multiple structures by listing PDB codes:

.. code-block:: yaml

   input: [1OS7, 1FYZ, 1PHM]

Plain-Text PDB List
-------------------

A text file with one PDB ID or path per line is also accepted:

.. code-block:: yaml

   input: pdb_list.txt

CSV Input File
--------------

For high-throughput processing, use a CSV file. This is the recommended format
when processing many structures or when per-PDB parameters are needed.

.. code-block:: yaml

   input: proteins.csv

**Required column:**

- ``pdb_id`` --- PDB accession code or path to a local PDB / mmCIF file

**Optional columns:**

- ``center`` --- Center residue definition for this specific PDB (overrides
  the ``center_residues`` parameter in the YAML config)
- ``oxidation`` --- Metal oxidation state (**required for** ``qp submit``)
- ``multiplicity`` --- Spin multiplicity (**required for** ``qp submit``)

**Example CSV:**

.. code-block:: text

   pdb_id,center,oxidation,multiplicity
   1OS7,FE,3,6
   1FYZ,FE2_A5001-FE2_A5002,6,11
   1PHM,CU_A357-CU_A358,4,3

.. important::

   ``qp run`` accepts PDB IDs, local paths, lists, text files, and CSVs.
   ``qp submit`` and electronic-state setup require ``input`` to be a path to
   an **existing CSV file**.

.. _mmcif-support:

mmCIF Support
-------------

QuantumPDB currently uses classic PDB as the working format for Modeller,
Protoss, and clustering. mmCIF is supported through conversion:

1. **RCSB fallback** --- If fetching a classic PDB returns 404, QuantumPDB
   downloads the mmCIF and converts it.
2. **Local mmCIF** --- Paths ending in ``.cif`` or ``.mmcif`` are converted
   via ``qp.structure.mmcif_to_pdb``.
3. **Remapping** --- Multi-character chain IDs and residue names longer than
   three characters are remapped to classic PDB-compatible values. The mapping
   is written as ``{stem}_mmcif_remap.json`` beside the converted PDB.
4. **Center matching** --- When a remap sidecar is present, center residue
   specifications may use either the original mmCIF names or the remapped
   names.
5. **Size limits** --- Structures exceeding classic PDB limits (**>99999
   atoms** or **>62 chains**) raise ``OversizedStructureError`` and are
   skipped with a warning during batch ``qp run``.

Native large-mmCIF handling without classic PDB conversion is still planned.

Center Residue Syntax
---------------------

The center residue defines which atoms are at the core of the QM cluster model.
The syntax differs between YAML and CSV formats, but matching always follows
two underlying modes implemented by ``CenterResidue``.

Matching Modes (HETATM vs ATOM)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**1. Name-based / fuzzy mode** (no ``exact:`` prefix and no ``-``)
   Tokens are residue names, optionally joined by underscores
   (``FE``, ``FE_FE2``, ``BGC_GAL_NGA``). QuantumPDB matches **all residues
   with those names that are recorded as HETATM** (BioPython
   ``res.id[0] != ' '``). Protein-chain ATOM residues are ignored, even if
   the three-letter name matches.

   This is intentional: ``center_residues: [HIS]`` must not build a cluster
   for every histidine in the fold. The same chemistry can still be selected
   when it appears as a hetero ligand:

   - ``HIS`` written as a **HETATM ligand** → matched in fuzzy mode
   - ``HIS`` in the **protein backbone/side chain as ATOM** → not matched
     in fuzzy mode

   Nearby matched HETATM centers can be combined afterward with
   ``merge_distance_cutoff`` (see below). With the default cutoff ``0.0``,
   each matched HETATM typically becomes its own cluster.

**2. Exact-list / strict mode**
   Tokens are exact ``RESNAME_CHAINID`` keys. Strict mode does **not**
   require HETATM: coordinating amino acids stored as ATOM records can be
   included.

   - **Single residue** must use the ``exact:`` prefix:
     ``exact:FE_A199``. Bare ``FE_A199`` stays fuzzy (names ``FE`` and
     ``A199``), because ``A199`` could be a real CCD ligand code.
   - **Multiple residues** may omit the prefix when joined by dashes:
     ``CU_A357-CU_A358`` or ``FE_A155-OXY_A555-HEM_A155-HIS_A93``.
   - The prefix is also allowed on multi-residue lists:
     ``exact:FE_A155-HIS_A93``.

.. important::

   Use ``exact:FE_A199`` to pin one residue. Do not rely on the shape of
   ``FE_A199`` alone --- without ``exact:`` or a dash-separated list it is
   fuzzy HETATM name matching.

YAML Format (``center_residues``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use brackets in the YAML config for **general / fuzzy** name matching:

.. code-block:: yaml

   # Single residue type (all HETATM FE)
   center_residues: [FE]

   # Batch mode: one center type per PDB in processing order
   center_residues: [FE, FE2]

.. warning::

   A YAML list such as ``[FE, FE2]`` is **not** interpreted as "build FE and
   FE2 centers inside every structure". In batch processing, QuantumPDB pops
   one center definition per PDB. To place multiple center **names** in one
   structure, use fuzzy underscore syntax (``FE_FE2``), an explicit CSV
   ``center`` value, or ``merge_distance_cutoff``.

CSV Format (``center`` column)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CSV ``center`` column accepts fuzzy names, ``exact:`` singles, and
dash-separated strict lists:

.. code-block:: text

   # Fuzzy: all HETATM residues named FE
   FE

   # Fuzzy multi-type: all HETATM FE and FE2
   FE_FE2

   # Still fuzzy: names FE and A199 (A199 may be a CCD code)
   FE_A199

   # Strict single residue (required exact: prefix)
   exact:FE_A199

   # Strict merged center (dash-separated; prefix optional)
   CU_A357-CU_A358

   # Strict complex center including a protein HIS (ATOM)
   FE_A155-OXY_A555-HEM_A155-HIS_A93

Each strict token uses ``RESNAME_CHAINID``, where ``RESNAME`` is the
residue code and ``CHAINID`` is the chain letter followed by the residue
number (e.g., ``FE_A199`` for iron at position 199 on chain A).

**Dash-merged centers** treat the listed residues as a single center. Typical
uses:

- Multi-metallic active sites (e.g., ``CU_A357-CU_A358`` for dicopper)
- Heme / cofactor assemblies that must include a specific protein residue
  (e.g., ``FE_A155-HEM_A155-OXY_A555-HIS_A93``)
- Any case where the center must include ATOM amino acids that fuzzy mode
  would skip
- A single exact residue via ``exact:HIS_A93``

.. important::

   When both ``center_residues`` (YAML) and a ``center`` column (CSV) are
   provided, the **CSV value takes priority** for each PDB.

Merging Nearby Centers
^^^^^^^^^^^^^^^^^^^^^^

For multi-metallic sites selected in **fuzzy** mode, use
``merge_distance_cutoff`` to automatically merge matched HETATM centers that
lie within a given distance:

.. code-block:: yaml

   center_residues: [FE2]
   merge_distance_cutoff: 4.0

This merges any matched ``FE2`` centers within 4.0 Å of each other into a
single center. Set the cutoff slightly larger than the metal--metal distance.
Examples:

- **MMO diiron** (Fe--Fe ~3.4 Å): ``merge_distance_cutoff: 4.0``
- **PHM dicopper** (Cu--Cu ~11 Å): ``merge_distance_cutoff: 12.0``

Strict dash lists already define one combined center; they do not depend on
``merge_distance_cutoff`` to stay together.

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
