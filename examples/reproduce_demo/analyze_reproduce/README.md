# Analyze Reproduce (SCOOTI)

Mirror of analyze_demo that reads coefficients from your inference reproduce output under the SCOOTI repo.

- Coefficients: `./examples/reproduce_demo/inference_reproduce/out/regression_models/`
- Constrained fluxes: `./examples/constrained_demo/out/constrained_models/`
- Outputs: `./examples/reproduce_demo/analyze_reproduce/out/`

## Usage

Full (legacy) engine:

```
bash examples/reproduce_demo/analyze_reproduce/run_analyze_reproduce.sh
```

Minimal engine:

```
bash examples/reproduce_demo/analyze_reproduce/run_analyze_reproduce_minimal.sh
```

Pass a custom config JSON as an argument to override paths or labeling.
