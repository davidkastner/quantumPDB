Cluster Models
==============

Stage 3 of ``qp run`` builds hierarchical QM cluster models around catalytic
centers. The implementation lives in ``qp.cluster.spheres`` and is exercised by
the golden samples under ``qp/tests/samples/``.

Center Matching
---------------

Centers are parsed by ``CenterResidue``. Keep the user-facing syntax from
:doc:`input_formats`, but remember the HETATM rule:

- **General / fuzzy names** --- ``FE`` or underscore lists such as ``FE_FE2``
  match **all HETATM** residues with those names. ATOM records in the protein
  chain are skipped. Example: a free ``HIS`` ligand stored as HETATM can be a
  center; the same ``HIS`` as an ATOM residue in the fold cannot.
- **Fuzzy lookalikes** --- bare ``FE_A199`` is still fuzzy (names ``FE`` and
  ``A199``), because ``A199`` may be a CCD code.
- **Strict single residue** --- ``exact:FE_A199`` pins one
  ``RESNAME_CHAINID`` (ATOM or HETATM). The ``exact:`` prefix is required.
- **Strict dash lists** --- ``CU_A357-CU_A358`` or
  ``FE_A155-OXY_A555-HEM_A155-HIS_A93`` select exact keys and may omit
  ``exact:``. Those residues form one ``metal_id``.
- **Optional distance merge** --- after fuzzy matching, ``merge_distance_cutoff``
  can combine nearby HETATM centers.
- **mmCIF aliases** --- when ``*_mmcif_remap.json`` is present, original or
  remapped residue names both work.


Sphere Construction
-------------------

1. **Sphere 0** contains the matched center residue(s).
2. **Sphere 1** uses a distance cutoff (``radius_of_first_sphere``, default
   4.0 Å) plus chemistry-aware filtering that prefers coordinating N/C/O/S
   environments around metals.
3. **Higher spheres** grow by Voronoi adjacency: residues that share a face
   with the previous shell are included, which captures non-spherical active
   sites better than a pure radial cutoff.
4. ``include_ligands`` controls ligand / water inclusion:

   - ``0`` --- ligands/waters only if already selected for the first sphere
   - ``1`` --- include non-water ligands in outer spheres
   - ``2`` --- include ligands and waters (default)

5. ``additional_ligands`` forces named residues into the first sphere and
   protects them from pruning.

Smoothing Methods
-----------------

Voronoi cells in sparse regions can produce long, chemically unreasonable
contacts. ``smoothing_method`` maps to:

======= ================= ==========================================
Value   Method            CLI parameters
======= ================= ==========================================
``0``   ``box_plot``      none
``1``   ``dbscan``        ``eps=6``, ``min_samples=3``
``2``   ``dummy_atom``    ``mean_distance=3`` (recommended default)
``3``   disabled          no smoothing
======= ================= ==========================================

Dummy-atom smoothing inserts auxiliary points to regularize Voronoi cells and
is the default for most metalloenzyme sites.

Merging, Pruning, and Oligomers
-------------------------------

- ``merge_distance_cutoff`` recursively merges centers closer than the cutoff
  into one ``metal_id`` (for example ``C601_C602``).
- ``max_atom_count`` prunes the most distant residues until the atom count is
  under the limit. Cap atoms are not counted.
- Oligomeric ligands are completed when needed so formal charges stay
  consistent; partial oligomer inclusion may emit charge warnings.

Capping
-------

Chain breaks created by carving the cluster are capped with
``capping_method``:

- ``0`` --- no caps
- ``1`` --- hydrogen caps (default; good for single points)
- ``2`` --- ACE / NME caps (preferred for geometry optimizations)

Enabling capping (or ``compute_charges``) forces Protoss on so protonation
states are available.

Outputs and Naming
------------------

Each matched center writes a directory named by default ``metal_id``:

.. code-block:: text

   {metal_id}/
   ├── 0.pdb
   ├── 1.pdb
   ├── 2.pdb
   ├── {metal_id}.pdb
   └── {metal_id}.xyz

Optional controls:

- ``write_xyz`` --- write the combined XYZ (needed for ``qp submit``)
- ``write_hetero_pdb`` --- keep HETATM records in the combined PDB
- ``cluster_name_template`` --- custom directory/file names; also writes
  ``cluster_name_map.csv``
- ``compute_charges`` / ``count_residues`` --- write ``charge.csv`` and
  ``count.csv`` (see :doc:`output`)

Charge and Spin Checks
----------------------

After extraction, QuantumPDB can sanity-check that the electron count implied
by XYZ elements is consistent with the charge / spin bookkeeping used later by
``qp submit``. Inspect ``charge.csv``, ``spin.csv``, and CLI warnings if a
cluster looks electronically inconsistent.

Related Pages
-------------

- :doc:`configuration` --- parameter defaults
- :doc:`output` --- file formats
- :doc:`qm_jobs` --- how clusters become TeraChem jobs
- :doc:`special_workflows` --- NHIE-oxo conversion before clustering
