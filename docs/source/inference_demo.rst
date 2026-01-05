Inference of Metabolic Objectives (Demo)
========================================

Overview
--------
Given unconstrained flux predictions (ideal objectives) and constrained fluxes (omics-informed),
SCOOTI infers objective coefficients via meta-regression: the constrained flux is modeled as a
linear combination of unconstrained flux distributions. The coefficient for each objective
captures its contribution to the observed metabolic state.

Inputs and Settings
-------------------
Prepare a config with:

::

  {
    "unconModel": "./examples/unconstrained_demo/out/unconstrained_models/",
    "conModel":   "./examples/constrained_demo/out/constrained_models/",
    "savePath":   "./examples/inference_demo/out/regression_models/",
    "kappaArr":   "10,1,0.1",
    "rhoArr":     "10,1,0.1",
    "dkappaArr":  "10,1,0.1",
    "expName":    "inference_demo",
    "unconNorm":  "T",
    "conNorm":    "F",
    "medium":     "DMEMF12",
    "method":     "cfr",
    "model":      "recon1",
    "inputType":  "flux",
    "clusterPath": "",
    "objListPath": "",
    "rank":      "F",
    "stackModel": "F",
    "sampling":  "F",
    "learner":   "L",
    "geneKO":    "F",
    "geneListPath": "",
    "learningRate": 0.001,
    "epo": 5000,
    "fileSuffix": "_fluxes.csv.gz"
  }

We provide this at ``examples/inference_demo/demo_inference_config.json``.

Run the Pipeline
----------------

Demo run (one‑liner)
^^^^^^^^^^^^^^^^^^^^
From the repository root:

::

  # Demo wrapper
  bash examples/inference_demo/run_inference.sh

  # Or run the full quickstart (generates fluxes then runs inference)
  bash examples/quickstart/run_all.sh

Real use (one‑liner)
^^^^^^^^^^^^^^^^^^^^
Call the generic trainer with your JSON config:

::

  bash scooti/run_trainer.sh examples/inference_demo/demo_inference_config.json

Notes
-----
- Ensure that the medium tag (e.g., ``DMEMF12``) matches between your flux outputs and the
  inference config; otherwise, no files are found for regression.
- Ensure the suffix (``fileSuffix``) matches your flux outputs; quickstart/demo uses ``_fluxes.csv.gz``.
- The run prints a device summary (CPU/GPU) at start.

Outputs and Interpretation
--------------------------
- The ``savePath`` folder contains CSV files of objective coefficients for each sample/condition.
- Each sample obtains a vector of coefficients (default 52), reflecting the contribution of each
  objective to the observed fluxes.
- Coefficients can be profiled across samples or compared between conditions to infer metabolic
  priorities (e.g., growth, redox balance, antioxidant capacity).

Timing
------
- Training meta‑regression (demo scale): 10–45 minutes CPU; faster with GPU.
- Larger cohorts: up to a few hours depending on sample count.
