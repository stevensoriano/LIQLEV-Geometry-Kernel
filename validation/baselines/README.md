# Baseline Output Directory

`write_physics_baseline.py` writes JSON baselines here.

Recommended first baseline:

```bash
python scripts/write_physics_baseline.py --output validation/baselines/physics_baseline.json
```

The generated JSON captures the current solver output for canonical cases. Keep a copy available during modernization so future changes can prove result parity.

