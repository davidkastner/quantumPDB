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
       ├── {pdb}.pdb                    # Original fetched PDB
       ├── {pdb}.ali                    # Modeller alignment file
       ├── {pdb}_modeller.pdb           # Rebuilt structure (Stage 1)
       ├── Protoss/
       │   ├── {pdb}_protoss.pdb        # Protonated structure (Stage 2)
       │   ├── {pdb}_protoss_orig.pdb   # Copy before active-site fixes
       │   ├── {pdb}_ligands.sdf        # Ligand structures (SDF format)
       │   └── {pdb}_log.txt            # Protoss processing log
       ├── charge.csv                   # Per-residue charges
       ├── count.csv                    # Residue counts per sphere
       ├── spin.csv                     # Radical species spins (if present)
       └── {CENTER}_{CHAIN}{ID}/        # One directory per center residue
           ├── 1/                       # Sphere 1 (innermost)
           │   ├── cluster.pdb
           │   └── cluster.xyz
           ├── 2/                       # Sphere 2 (includes sphere 1)
           │   ├── cluster.pdb
           │   └── cluster.xyz
           └── ...                      # Additional spheres

Cluster Directory Naming
------------------------

Cluster directories are named after the center residue using the format
``{RESNAME}_{CHAIN}{RESID}``. For example:

- ``FE_A501`` --- iron at position 501 on chain A
- ``CU_A357-CU_A358`` --- merged dicopper center
- ``HEM_A155`` --- heme group

Each numbered subdirectory contains a progressively larger cluster model.
Sphere 1 is the innermost (first coordination shell), sphere 2 adds the
next layer of interacting residues, and so on. Higher-numbered spheres are
supersets --- sphere 2 contains all atoms from sphere 1 plus the second
shell.

File Descriptions
-----------------

**Structure files:**

- ``{pdb}.pdb`` --- The original PDB file, either downloaded from the RCSB
  or copied from a local path.
- ``{pdb}_modeller.pdb`` --- The structure after Modeller has rebuilt missing
  atoms, residues, and loops.
- ``{pdb}_protoss.pdb`` --- The protonated structure with hydrogens added,
  alternate conformations resolved, and active-site corrections applied.
- ``{pdb}_protoss_orig.pdb`` --- A copy of the Protoss output before
  metalloenzyme-specific active-site corrections.

**Ligand data:**

- ``{pdb}_ligands.sdf`` --- SDF file containing ligand structures as output
  by Protoss. Used to compute ligand charges and spins.

**Charge and count data:**

- ``charge.csv`` --- Net charge for each residue in the cluster. Columns:
  chain ID, residue name, residue ID, net charge. Includes both amino acid
  charges (from protonation state analysis) and ligand charges (from Protoss
  SDF).
- ``count.csv`` --- Number of residues in each interaction sphere.
- ``spin.csv`` --- Spin contributions from radical species (e.g., NO as
  doublet, O2 as triplet). Only generated when radical ligands are present.

**Cluster models:**

- ``cluster.pdb`` --- QM cluster model in PDB format with chain break caps.
- ``cluster.xyz`` --- QM cluster model in XYZ format (suitable for direct
  use with QM codes).

**QM job files** (after ``qp submit``):

- QM input scripts (TeraChem ``.in`` or ORCA ``.inp``)
- Scheduler submission scripts (``.sh``)
- ``.submit_record`` --- Hidden marker file created after job submission to
  prevent duplicate submissions
- ``ptchrges.xyz`` --- MM point charges for charge embedding (if enabled)

**Analysis output** (after ``qp analyze``):

- ``checkup/failure_modes.csv`` --- Failure mode classification for all jobs
- ``checkup/job_status.csv`` --- Status of each job
- ``checkup/failure_modes.png`` --- Plot of failure modes
- ``checkup/job_status.png`` --- Plot of job statuses
