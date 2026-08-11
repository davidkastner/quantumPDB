# QuantumPDB examples

Contributor examples live under this directory. Each subfolder is a self-contained
mini-workflow.

| Directory | Maintainer focus |
|-----------|------------------|
| [`basics/`](basics/) | Minimal starter configs (cluster → submit → analyze) |

## Adding your own

1. Create a new subdirectory, e.g. `example/my_workflow/`.
2. Include a short `README.md` with purpose, commands, and expected outputs.
3. Prefer paths relative to the repository root in YAML/CSV so commands work as
   `qp run -c example/my_workflow/config.yaml`.
