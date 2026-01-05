# Inference Reproduce

Run meta-regression using your own unconstrained models while keeping the constrained models fixed.

- Unconstrained models: `./examples/quickstart/unconstrained_models/`
- Constrained models: `./examples/constrained_demo/out/constrained_models/`
- Outputs: `./examples/quickstart/inference_reproduce/out/regression_models/`

## Usage

From the SCOOTI repo root:

```
bash examples/quickstart/inference_reproduce/run_reproduce_inference.sh
# or with a custom config JSON
bash examples/quickstart/inference_reproduce/run_reproduce_inference.sh /abs/path/to/config.json
```

Ensure the unconstrained directory contains pairs like `*_metadata.json` and `*_fluxes.csv.gz` (or `*_fluxes.csv`).
