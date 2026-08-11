## Installing testing dependencies
```bash
pip install -e ".[test]"
```

## Sample PDBs
Outputs generated with default parameters unless otherwise stated. For more details, see [samples.yaml](samples.yaml).

- `1lm6` (boxplot smoothing)
- `1sp9`
- `2chb` (merge cutoff 4.0)
- `2fd8` (merge cutoff 2.0, max atom count 102)
- `2q4a`
- `2r6s` (DBSCAN smoothing)
- `3a8g`
- `3x20`
- `4ilv` (ACE/NME capping)
- `4z42` (merge cutoff 4.0)
- `6f2a`

When making changes to core methods, use [samples.yaml](samples.yaml) to update the ground truth files as necessary.

## Coordination integration cases
End-to-end coordination-step checks live under [coordination_cases/](coordination_cases/) with metadata in [coordination_cases/cases.yaml](coordination_cases/cases.yaml). Each case ships precomputed modeller + Protoss outputs so runtime skips those stages (`skip: all`) and only regenerates clusters/charges.

```bash
pytest -v qp/tests/test_coordination_cases.py
# or: pytest -m integration
```

Cases cover protonation quirks (ASP/GLU/HIS), Cys/CSS ligands, cluster smoothing, nitrile hydratase variants, O₂ spin/charge, and partial-occupancy selection.

## Running tests
To generate coverage report
```bash
pytest --cov=qp --cov-report=html
open htmlcov/index.html
```
