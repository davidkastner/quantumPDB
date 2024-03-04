## Installing testing dependencies
```bash
pip install -e ".[test]"
```

## Sample PDBs
Outputs generated with default parameters unless otherwise stated

- `1lm6` (boxplot smoothing)
- `1sp9`
- `2q4a`
- `2r6s` (DBSCAN smoothing)
- `3a8g`
- `3x20`
- `4ilv` (ACE/NME capping)
- `6f2a`

## Running tests
To generate coverage report
```bash
pytest --cov=qp --cov-report=html
open htmlcov/index.html
```