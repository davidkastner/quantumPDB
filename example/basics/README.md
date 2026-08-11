# Basics

Minimal configs shared with the Sphinx quickstart / user guides. Run them from
the repository root after installing `qp` (`pip install -e .`) and setting
`KEY_MODELLER` if you enable Modeller.

## Files

| File | Purpose |
|------|---------|
| `cluster_only.yaml` | Download 1OS7 and build FE-centered clusters |
| `proteins.csv` | Batch CSV with oxidation / multiplicity for submit |
| `batch_submit.yaml` | Create TeraChem jobs from existing clusters |
| `charge_embedding.yaml` | Submit with MM point-charge embedding + PCM |
| `analyze.yaml` | Job checkup + Multiwfn charge / dipole flags |
| `custom_charges.json` | Tiny custom-charge JSON shape for embedding |

## Cluster-only quickstart

```bash
qp run -c example/basics/cluster_only.yaml
```

Expected layout under `example/basics/output/{pdb}/` matches
[Output Structure](https://hjkgrpquantumpdb.readthedocs.io/en/latest/output.html):
`Protoss/`, `charge.csv`, `count.csv`, and `{metal_id}/0.pdb`, `1.pdb`, …
plus `{metal_id}.xyz`.

## Batch submit

1. Produce clusters with `qp run` using a CSV `input`.
2. Create jobs (does not require a live scheduler if `submit_jobs: false`):

```bash
qp submit -c example/basics/batch_submit.yaml
```

`qp submit` requires `input` to be an existing CSV path with `oxidation` and
`multiplicity` columns.

## Charge embedding

```bash
qp submit -c example/basics/charge_embedding.yaml
```

## Analyze

```bash
qp analyze -c example/basics/analyze.yaml
```

`job_checkup` writes `checkup/` in the **current working directory**. Multiwfn
steps need completed TeraChem `scr/*.molden` files and a Multiwfn executable.

## Ground-truth samples

Regression outputs used by the test suite live in `qp/tests/samples/` and are
useful references for `charge.csv` / cluster naming.
