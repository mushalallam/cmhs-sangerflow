# Contributing

Bug reports and pull requests are welcome. Do not attach patient chromatograms or identifiers to
issues. Use synthetic or explicitly redistributable test data.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
ruff check .
pytest
```

Changes to trimming, alignment, consensus, or variant logic must include regression tests and a
clear explanation of the biological assumption being changed. Keep generated data and `.ab1`
files out of commits.

