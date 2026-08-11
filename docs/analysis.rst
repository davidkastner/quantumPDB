Analysis
========

``qp analyze`` covers job monitoring and Multiwfn post-processing. The code
lives in ``qp.analyze``.

Job Checkup
-----------

Set ``job_checkup: true`` to classify jobs under each
``{cluster}/{method}/`` directory. Status categories include done, running,
queued, backlog, and error / failure modes parsed from scheduler markers and
``qmscript.out``.

Reports are written to ``checkup/`` in the **current working directory**
(process CWD), not necessarily ``output_dir``:

- ``failure_modes.csv`` / ``failure_modes.png``
- ``job_status.csv`` / ``job_status.png``

``delete_queued: true`` removes ``.submit_record`` files for unfinished jobs so
``qp submit`` can enqueue them again.

Multiwfn Charge Schemes
-----------------------

Set ``calc_charge_schemes: true`` and choose a single ``charge_scheme``:

- ``Hirshfeld`` (default)
- ``Voronoi`` (Voronoi deformation density; unrelated to cluster Voronoi spheres)
- ``Mulliken``
- ``ADCH``
- ``Hirshfeld-I``
- ``CM5``

QuantumPDB looks for completed TeraChem molden files under
``{cluster}/{method}/scr/*.molden``. ``multiwfn_path`` must point to a Multiwfn
executable. Atomic radii / settings come from bundled ``qp.resources`` assets.

.. warning::

   ``charge_scheme`` must be **one** scheme name. Comma-separated values such as
   ``Hirshfeld,CM5`` raise ``ValueError``.

Dipole Moments
--------------

``calc_dipole: true`` runs Multiwfn dipole analysis using the center of mass as
the reference. There is no separate ``calc_dipole_coc`` config key in the
current CLI.

Common Failure Modes
--------------------

- Wrong ``method`` subdirectory name relative to submit config
- Missing molden / incomplete QM output
- Multiwfn not on ``PATH``
- Running checkup from a directory other than where you expect ``checkup/``
- Analyzing before jobs finish

Related Pages
-------------

- :doc:`qm_jobs` --- job creation and markers
- :doc:`output` --- where analysis files land
- :doc:`faq` --- troubleshooting
- :doc:`configuration` --- analyze parameters
