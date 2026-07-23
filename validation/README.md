# Validation Helpers

This package contains repeatable physics-preservation checks for modernization work.

Use:

```bash
python scripts/write_physics_baseline.py
python scripts/check_physics_baseline.py
```

Generate the baseline before refactoring. Re-run the check after each phase. The scripts are intentionally independent of the GUI so they can run on a headless cluster.

