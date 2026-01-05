Prepare Unconstrained Models Optimizing Ideal Objectives
=======================================================

This section shows how to generate “unconstrained” flux models, i.e., models optimized
for each ideal objective without applying omics-based constraints. These fluxes serve as
predictors for meta-regression analyses.

What “Unconstrained” Means
--------------------------
- No gene/protein expression constraints are applied.
- Only the default bounds and model structure are used.
- SCOOTI optimizes the GEM independently for each objective (e.g., ATP, NADH, glutathione).

Timing
------
Depending on the solver and hardware, this can take from a couple of hours up to a day.

Configuration
-------------
Create or edit a config JSON with these key fields:

::

  {
    "COBRA_path": "~/cobratoolbox",
    "GEM_path": "./scooti/metabolicModel/GEMs/Shen2019.mat",
    "model_name": "Recon1",
    "obj_candidate_list_file": "./scooti/metabolicModel/GEMs/obj52_metabolites_shen2019.csv",
    "pfba": 1,
    "data_dir": "",
    "prefix_name": "model",
    "medium": "DMEMF12",
    "simulation": "CFR",
    "constraint": 0,
    "save_root_path": "./examples/unconstrained_demo/out/unconstrained_models/"
  }

We include an example config:

- ``examples/unconstrained_demo/unconstrained_demo_config.json``

Run The Pipeline
----------------

Demo run (one‑liner)
^^^^^^^^^^^^^^^^^^^^
From the repository root:

::

  # Demo wrapper
  bash examples/unconstrained_demo/run_unconstrained.sh

  # Or run the full quickstart (generates both unconstrained and constrained)
  bash examples/quickstart/run_all.sh

Real use (one‑liner)
^^^^^^^^^^^^^^^^^^^^
Call the generic runner with your JSON config:

::

  bash scooti/run_flux.sh examples/unconstrained_demo/unconstrained_demo_config.json

Which script or command is run to obtain this output?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use either the demo wrapper or the direct runner (both are valid):

::

  bash examples/unconstrained_demo/run_unconstrained.sh
  # or
  bash scooti/run_flux.sh examples/unconstrained_demo/unconstrained_demo_config.json

Example Output (truncated)
--------------------------

::

  [CFR] Flux results saved to ./examples/unconstrained_demo/out/unconstrained_models//[NovxxYYYYHHMM]CFR_ct1_obj1_data1_fluxes.csv
  [CFR] ... one file per objective

Notes
-----
- The output folder contains one flux distribution per objective; filenames include objective IDs.
- If an optimization returns infeasible, check the GEM and solver setup (e.g., Gurobi)
  and ensure ``COBRA_path`` and ``GEM_path`` are correct.
- Medium can be set to ``DMEMF12`` or ``KSOM`` (or left empty) as appropriate for your system; keep
  it consistent with how you plan to normalize/compare later.
