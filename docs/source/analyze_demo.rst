Analyze Metabolic Objectives (Demo)
===================================

Overview
--------
Use the analyzer to examine inferred objective coefficients, generate visualizations,
and compute comparative metrics.

Config (JSON)
-------------
Minimal example at ``examples/analyze_demo/analyze_config.minimal.json``:

::

  {
    "flux_paths": {"exp": "./examples/constrained_demo/out/constrained_models/"},
    "coef_paths": {"exp": "./examples/inference_demo/out/regression_models/"},
    "save_root_path": "./examples/analyze_demo/out/",
    "engine": "minimal",
    "reduction": "auto",
    "GEM_path": "./scooti/metabolicModel/GEMs/Shen2019.mat",
    "uncon_model_path": "./examples/unconstrained_demo/out/unconstrained_models/",
    "col_map": {},
    "prefix": "analyze_demo",
    "medium": "DMEMF12",
    "get_flux": {"kappa": [10, 1, 0.1], "rho": [10, 1, 0.1], "rank": false},
    "flux_analysis": {"kappaArr": [10, 1, 0.1], "rhoArr": [10, 1, 0.1]},
    "get_coef": {"metType_cluster": false},
    "coef_analysis": {
      "clustering": false, "entropy": false, "distance": true, "compare": true,
      "umap_para": [5, 50], "method": "average", "ref_col": null
    }
  }

Run (Demo)
----------
From the repository root:

::

  bash examples/analyze_demo/run_analyze_minimal.sh

Real Use (One‑liner)
--------------------
Call the generic analyzer with your config:

::

  bash scooti/run_analyzer.sh examples/analyze_demo/analyze_config.json

Notes
-----
- Ensure the ``coef_paths`` point to CSVs written by the inference step.
- The analyzer writes figures and tables into ``save_root_path``.
- UMAP safety: both engines automatically clamp ``n_components`` and ``n_neighbors`` to ``N-1``
  samples to avoid ARPACK errors (you’ll see a short adjustment log). To override the UMAP
  parameters, set ``umap_para: [n_neighbors, n_components]``.

Label configuration
-------------------

By default, both engines label samples by the ``coef_paths`` key for each column (e.g., ``exp``).
You can override labeling in the JSON without code using the ``labels`` section:

::

  {
    "labels": {
      "mode": "regex",              # one of: group_key (default), column, regex, split
      "regex": "^(.*?)_",           # only for mode=regex
      "regex_group": 1,              # which capture group to use (default 1)
      "split_delim": "_",           # only for mode=split
      "split_index": 0,              # only for mode=split
      "group_map": {"exp": "Experiment"}  # optional renaming of outputs
    }
  }

Modes
^^^^^
- ``group_key``: use ``coef_paths`` keys to label points (default).
- ``column``: use exact column names.
- ``regex``: labels derived by applying the regex to each column name.
- ``split``: split each column name by a delimiter and take a specific index.

These options apply to both minimal and legacy engines (legacy applies overrides after loading coefficients).

Timing
------
- Minimal engine: 5–20 minutes depending on number of samples and UMAP neighbors.
- Legacy engine (Pareto overlays, trait analysis): 15–60 minutes depending on feature set.
