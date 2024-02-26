## Installing testing dependencies
```bash
pip install -e ".[test]"
```

## Sample PDBs
Outputs generated with default parameters unless otherwise stated

- `1sp9`
- `2q4a`
- `3a8g`
- `4qm8` (ACE/NME capping)
- `1lm6` (boxplot smoothing)
- `2r6s` (DBSCAN smoothing)

## Running tests
To generate coverage report
```bash
pytest --cov=qp --cov-report=html
open htmlcov/index.html
```