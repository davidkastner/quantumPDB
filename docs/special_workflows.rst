Special Workflows
=================

Most users only need the standard Modeller → Protoss → cluster pipeline. This
page documents optional paths that are easy to miss in the code.

NHIE Oxo / Succinate Conversion
-------------------------------

Set ``convert_to_nhie_oxo: true`` (requires ``modeller`` and ``protoss``) to run
the specialized AKG → oxo / succinate NHIE pathway implemented in
``qp.structure.convert_nhie_oxo``.

**Intended use**

- Non-heme iron / α-ketoglutarate systems where you want the post-reaction
  succinate + oxo state represented before clustering

**What happens**

1. Standard Modeller + Protoss preparation runs.
2. QuantumPDB rewrites relevant ligand coordinates / identities toward oxo and
   succinate (SIN) representations.
3. Protoss may be invoked again for the updated ligand set / SDF bookkeeping.
4. Clustering proceeds on the converted structures.

**Limitations**

- This is **not** a general ligand mutation tool.
- It is unused by default and has no dedicated regression test in
  ``qp/tests``.
- Prefer leaving it ``false`` unless you know the target chemistry matches.

Preset Protoss Cache
--------------------

For a small set of PDB IDs, QuantumPDB may copy precomputed Protoss outputs
from ``qp/resources/prepared/{pdb}/Protoss/`` instead of calling the web API.

Treat this as a performance / reproducibility optimization for known entries,
**not** as a general offline Protoss mode. Most structures still require network
access to ProteinsPlus.

Helper Scripts
--------------

``qp/resources/scripts/`` contains offline utilities (failure walking, PDB
metadata parsing, residue statistics, etc.). They are not exposed as ``qp``
subcommands and are aimed at maintainers / dataset curation rather than the
standard user pipeline.

Related Pages
-------------

- :doc:`workflow` --- default five-stage pipeline
- :doc:`configuration` --- ``convert_to_nhie_oxo``
- :doc:`cluster_models` --- clustering after conversion
