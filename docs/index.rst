.. image:: _static/logo-white.svg

QuantumPDB
==========

**QuantumPDB** (``qp``) is a fully integrated, open-source Python package that automates
the generation of quantum mechanical (QM) cluster models from protein structures. It takes
raw PDB files and produces ready-to-run QM calculation inputs through a five-stage pipeline:
structure preparation, protonation, cluster extraction, job management, and post-processing
analysis.

QM cluster models capture the electronic structure of enzyme active sites by carving out a
region of the protein surrounding the catalytic center. Unlike simple distance-based cutoffs,
QuantumPDB uses Voronoi tessellation to construct hierarchical interaction spheres that
faithfully capture the chemical environment, even for non-spherical active sites.

Key Features
------------

- **Automated structure preparation** --- models missing atoms, residues, and loops via Modeller
- **Protonation state assignment** --- assigns hydrogen positions and resolves alternate conformations via Protoss
- **Voronoi-based cluster extraction** --- contact-based interaction spheres instead of distance cutoffs
- **Flexible center definitions** --- single metals, multi-metallic sites, heme groups, oligomeric substrates
- **QM job management** --- generates input files for TeraChem (GPU) and ORCA (CPU) with optional MM charge embedding
- **Post-processing analysis** --- partial charge schemes and dipole moments via Multiwfn
- **High-throughput** --- processes hundreds of enzymes from a single CSV input file

QuantumPDB works with any PDB-format structure, including those from X-ray crystallography,
cryo-EM, NMR, MD simulation snapshots, and generative models.

.. container:: .buttons

   `GitHub <https://github.com/davidkastner/quantumpdb>`_
   `Zenodo Dataset <https://doi.org/10.5281/zenodo.17239737>`_

.. toctree::
   :maxdepth: 2
   :caption: Getting Started
   :hidden:

   getting_started
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   :hidden:

   workflow
   input_formats
   configuration
   cli
   output

.. toctree::
   :maxdepth: 2
   :caption: Reference
   :hidden:

   api
   faq
   citation
