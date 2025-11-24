Trade-off Analysis (Demo)
=========================

Overview
--------
Evaluate pairwise trade-offs between metabolic objectives. Trade-offs arise when two
objectives cannot be simultaneously maximized (cells strike a balance along a Pareto front).

Two paths are supported:

- Minimal: correlation-based trade-off CSVs and pairwise scatter plots.
- Legacy: full `metObjAnalyzer.tradeoff_analysis` with Pareto overlays.

Minimal Demo
------------

::

  conda activate scooti
  bash examples/tradeoff_demo/run_tradeoff_minimal.sh

Config (`examples/tradeoff_demo/tradeoff_config.minimal.json`):

::

  {
    "engine": "minimal",
    "coef_paths": {"exp": "./examples/inference_demo/out/regression_models/"},
    "save_root_path": "./examples/tradeoff_demo/out/",
    "tradeoff": true,
    "pareto_mets": ["glutathione", "biomass", "cholesterol"],
    "tradeoff_corr_threshold": -0.4,
    "umap_para": [10, 50]
  }

Outputs
^^^^^^^
- ``*_tradeoffs.csv``: strongest negative correlations across objectives.
- ``*_tradeoff_A_vs_B.png``: pairwise scatter plots among the listed objectives.
- Embedding plots/CSV from UMAP/PCA (safe UMAP automatically clamps components to N-1).

Label configuration
-------------------
Trade-off plots and embeddings use the same labeling controls as the analyze demo. Add a ``labels``
section to your JSON to switch between ``group_key`` (default), ``column``, ``regex``, or ``split``
strategies. See :doc:`analyze_demo` for details and examples.

Legacy Demo
-----------

::

  bash examples/tradeoff_demo/run_tradeoff_legacy.sh

Config (`examples/tradeoff_demo/tradeoff_config.legacy.json`):

::

  {
    "engine": "legacy",
    "tradeoff_analysis": {
      "pareto_mets": ["glutathione", "biomass", "cholesterol"],
      "pareto2DAnalysis": true,
      "traitAnalysis": true,
      "corr_heatmap": true
    }
  }

Outputs
^^^^^^^
- Correlation matrices (Pearson heatmap); objective pair scatter plots.
- Pareto front overlays for listed objectives; optional 3D Pareto analysis.
- If ``traitAnalysis`` is enabled: archetype-level summaries of metabolic trade-offs.

Timing
------
- Minimal trade‑off (correlations + pairwise plots): 10–30 minutes.
- Legacy trade‑off (Pareto overlays/3D): 20–90 minutes depending on grid density.
