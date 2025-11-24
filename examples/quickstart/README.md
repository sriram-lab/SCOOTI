Quickstart (End-to-End)
=======================

This quickstart runs a minimal, reproducible pipeline that:
- generates demo flux predictions (MATLAB), and
- infers metabolic objectives (Python),
with all inputs/configs kept inside this folder. Outputs are written under `out/`.

Prereqs: MATLAB + COBRA Toolbox configured; Python env installed per top-level README.

Run everything
--------------
```
bash examples/quickstart/run_all.sh
```

Run steps manually
------------------
- Flux (unconstrained + constrained):
```
bash scooti/run_flux.sh examples/quickstart/configs/unconstrained.json
bash scooti/run_flux.sh examples/quickstart/configs/constrained.json
```

- Inference:
```
bash scooti/run_trainer.sh examples/quickstart/configs/inference.json
```

Outputs
-------
- Flux predictions: `examples/quickstart/out/unconstrained_models/` and `.../constrained_models/`
- Inference results: `examples/quickstart/out/regression_models/`

Cleanup
-------
```
rm -rf examples/quickstart/out
```

