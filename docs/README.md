# Compiling QuantumPDB documentation

Docs are built with [Sphinx](https://www.sphinx-doc.org/) and the
`revitron_sphinx_theme`. Dependencies are pinned in `docs/requirements.yaml`
(also used by Read the Docs via `readthedocs.yml`).

## Local build

From a conda / micromamba environment that already has the package importable:

```bash
# optional dedicated docs env
conda env create -f docs/requirements.yaml
conda activate qp-docs   # or whatever name you chose / already use

cd docs
make clean
make html
```

Open `_build/html/index.html`.

If you are developing QuantumPDB itself, install the package first so autodoc
can import `qp`:

```bash
pip install -e ..
```

## Read the Docs

`readthedocs.yml` at the repository root builds with Python 3.11, installs from
`docs/requirements.yaml`, then `pip install .`. Hosted docs:
https://hjkgrpquantumpdb.readthedocs.io/

## Notes

- Autosummary stubs under `docs/autosummary/` are checked in and regenerated
  when `autosummary_generate` runs during the Sphinx build.
- Legacy module pages such as `qp.cluster.spheres_bad` may still exist on disk
  but are excluded from the user-facing `qp.cluster` index.
