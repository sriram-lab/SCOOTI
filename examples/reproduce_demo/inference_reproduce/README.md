# Inference Reproduce

Run meta-regression using your own unconstrained models while keeping the constrained models fixed.

- Unconstrained models: `/home/daweilin/sriram-lab/SCOOTI/examples/reproduce_demo/unconstrained_models/`
- Constrained models: `./examples/constrained_demo/out/constrained_models/`
- Outputs: `./examples/reproduce_demo/inference_reproduce/out/regression_models/`

## Usage

From the SCOOTI repo root:

```
bash examples/reproduce_demo/inference_reproduce/run_reproduce_inference.sh
# or with a custom config JSON
bash examples/reproduce_demo/inference_reproduce/run_reproduce_inference.sh /abs/path/to/config.json
```

Ensure the unconstrained directory contains pairs like `*_metadata.json` and `*_fluxes.csv.gz` (or `*_fluxes.csv`).
