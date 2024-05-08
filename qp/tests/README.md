## Installing testing dependencies
```bash
pip install -e ".[test]"
```

## Sample PDBs
Outputs generated with default parameters unless otherwise stated. For more details, see [samples.yaml](samples.yaml).

- `1lm6` (boxplot smoothing)
- `1sp9`
- `2chb` (merge cutoff 2.0)
- `2fd8` (merge cutoff 2.0, max atom count 102)
- `2q4a`
- `2r6s` (DBSCAN smoothing)
- `3a8g`
- `3x20`
- `4ilv` (ACE/NME capping)
- `4z42` (merge cutoff 2.0)
- `6f2a`

When making changes to core methods, use [samples.yaml](samples.yaml) to update the ground truth files as necessary. 

## Running tests
To generate coverage report
```bash
pytest --cov=qp --cov-report=html
open htmlcov/index.html
```