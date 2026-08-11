Output Structure
================

QuantumPDB organizes output into a directory hierarchy under the path
specified by ``output_dir``. Each PDB gets its own subdirectory containing
all pipeline outputs.

Directory Layout
----------------

After running all three stages (``qp run``, ``qp submit``, ``qp analyze``),
a typical output directory looks like:

.. code-block:: text

   output_dir/
   └── {pdb}/
       ├── {pdb}.pdb                    # Original fetched / converted PDB
       ├── {pdb}.cif                    # Present when input was mmCIF
       ├── {pdb}_mmcif_remap.json       # Present when mmCIF remapping occurred
       ├── {pdb}.ali                    # Modeller alignment file
       ├── {pdb}_modeller.pdb           # Rebuilt structure (Stage 1)
       ├── Protoss/
       │   ├── {pdb}_protoss.pdb        # Protonated structure (Stage 2)
       │   ├── {pdb}_protoss_orig.pdb   # Copy before active-site fixes
       │   ├── {pdb}_ligands.sdf        # Ligand structures (SDF format)
       │   └── {pdb}_log.txt            # Protoss clash log (empty under API v2)
       ├── charge.csv                   # Per-cluster / per-sphere charges
       ├── count.csv                    # Residue counts per sphere
       ├── spin.csv                     # Radical species spins (if present)
       ├── cluster_name_map.csv         # Written when cluster_name_template is set
       └── {metal_id}/                  # One directory per matched center
           ├── 0.pdb                    # Center residue(s) only
           ├── 1.pdb                    # First interaction sphere
           ├── 2.pdb                    # Second sphere (includes sphere 1)
           ├── {metal_id}.pdb           # Combined cluster PDB
           ├── {metal_id}.xyz           # Combined cluster XYZ
           └── {method}/                # After qp submit (e.g. wpbeh)
               ├── {metal_id}.xyz
               ├── qmscript.in          # TeraChem input
               ├── jobscript.sh         # SLURM / SGE submit script
               ├── ptchrges.xyz         # MM point charges (if enabled)
               └── .submit_record       # Marker after submission

Cluster Directory Naming
------------------------

By default, cluster directories are named by the matched center residue ID(s)
(``metal_id``), **not** ``{RESNAME}_{CHAIN}{RESID}``. Examples from the
regression samples:

- ``A204`` --- iron center on chain A, residue 204
- ``A446`` / ``B446`` --- separate centers on two chains
- ``C601_C602`` --- merged dicopper center
- ``B501_B502_B503`` --- merged multi-metal center

Set ``cluster_name_template`` to override naming (for example
``"A_{radius}"``). When a template is used, QuantumPDB also writes
``cluster_name_map.csv`` mapping the generated names back to ``metal_id``.

Each numbered ``{i}.pdb`` file is a cumulative sphere model:

- ``0.pdb`` --- the center residue(s)
- ``1.pdb`` --- first coordination / interaction sphere
- ``2.pdb`` --- second sphere (superset of sphere 1)
- higher numbers --- additional Voronoi shells

The combined ``{metal_id}.pdb`` / ``{metal_id}.xyz`` files contain the full
outermost cluster used for QM job setup.

File Descriptions
-----------------

**Structure files:**

- ``{pdb}.pdb`` --- The original or converted classic PDB file.
- ``{pdb}.cif`` / ``{pdb}_mmcif_remap.json`` --- Retained when the structure
  arrived as mmCIF; the sidecar maps original chain IDs and residue names to
  classic PDB-compatible values.
- ``{pdb}_modeller.pdb`` --- Structure after Modeller rebuilds missing atoms,
  residues, and loops.
- ``{pdb}_protoss.pdb`` --- Protonated structure with hydrogens, alternate
  conformations resolved, and active-site corrections applied.
- ``{pdb}_protoss_orig.pdb`` --- Protoss output before metalloenzyme-specific
  active-site corrections.

**Ligand data:**

- ``{pdb}_ligands.sdf`` --- SDF file containing ligand structures from Protoss.
  Used to compute ligand charges and spins.

**Charge and count data:**

- ``charge.csv`` --- Matrix of net charges per cluster and sphere. The header
  is ``Name,0,1,2,...`` where numeric columns are sphere indices. Cluster rows
  use the ``metal_id`` (for example ``A204``); ligand rows follow as
  ``RESNAME_CHAINRESID`` with their formal charges. Example:

  .. code-block:: text

     Name,0,1,2
     A204,0,-1,-2

     GOL_A301,0

- ``count.csv`` --- Residue composition of each interaction sphere. Header is
  ``Name,1,2,...``; each cell lists residue counts such as
  ``"GLN 1, HIS 2, HOH 3"``.
- ``spin.csv`` --- Spin contributions from radical species (for example NO as
  doublet, O2 as triplet). Only generated when radical ligands are present.
- ``cluster_name_map.csv`` --- Optional mapping written when
  ``cluster_name_template`` is set.

**Cluster models:**

- ``{i}.pdb`` --- Sphere-level PDB models inside the cluster directory.
- ``{metal_id}.pdb`` --- Combined QM cluster model in PDB format with caps.
- ``{metal_id}.xyz`` --- Combined QM cluster model in XYZ format (used by
  TeraChem job creation).

**QM job files** (after ``qp submit``):

- ``qmscript.in`` --- TeraChem input (current job generator is TeraChem-only)
- ``jobscript.sh`` --- Scheduler submission script (SLURM or SGE)
- ``.submit_record`` --- Hidden marker created after submission to prevent
  duplicate submissions
- ``ptchrges.xyz`` --- MM point charges for charge embedding (if enabled)

**Analysis output** (after ``qp analyze``):

Job checkup reports are written under ``checkup/`` in the **current working
directory** (not necessarily ``output_dir``):

- ``checkup/failure_modes.csv`` --- Failure mode classification for all jobs
- ``checkup/job_status.csv`` --- Status of each job
- ``checkup/failure_modes.png`` --- Plot of failure modes
- ``checkup/job_status.png`` --- Plot of job statuses

Multiwfn charge and dipole results are written beside the QM outputs under
``{cluster}/{method}/``.
