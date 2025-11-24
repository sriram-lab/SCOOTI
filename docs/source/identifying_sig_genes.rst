Identifying Significant Genes or Proteins
========================================

This short demo shows how to generate up-/down-regulated gene (or protein) lists from an example
omics table using two-sample t-tests (optionally with FDR correction) and a fold-change
criterion. The output CSVs can be used directly by the constrained flux pipeline.

Inputs
------

- Bulk expression table (genes × samples):
  ``examples/example_omics/johnson_18_GSE117444_prolif_qui_count.csv``
  with columns ``P_rep1..3`` (proliferative) and ``7dCI_rep1..3`` (quiescent)

How to Run (Demo)
-----------------

1. Activate the environment:

   ``conda activate scooti``

2. Run the demo script with its config:

   ``bash examples/siggenes_demo/run_siggenes.sh examples/siggenes_demo/siggenes_config.json``

Real use (Flux/Inference)
-------------------------
Use the demo to prepare gene lists, then run the modeling/inference pipelines with one‑liners:

::

  # Generate unconstrained models (demo wrapper or direct)
  bash examples/unconstrained_demo/run_unconstrained.sh
  # or
  bash scooti/run_flux.sh examples/unconstrained_demo/unconstrained_demo_config.json

  # Run inference (demo wrapper or direct)
  bash examples/inference_demo/run_inference.sh
  # or
  bash scooti/run_trainer.sh examples/inference_demo/demo_inference_config.json

   This writes four CSVs into ``examples/example_sigGenes``:

   - ``johnson_18_GSE117444_upgenes.csv``
   - ``johnson_18_GSE117444_dwgenes.csv``
   - ``johnson_18_GSE117444_P_upgenes.csv``
   - ``johnson_18_GSE117444_P_dwgenes.csv``

Configuration
-------------

The JSON file controls group identification and thresholds::

  {
    "tablePath": "examples/example_omics/johnson_18_GSE117444_prolif_qui_count.csv",
    "refPrefix": "7dCI_rep",
    "expPrefix": "P_rep",
    "alpha": 0.05,
    "fdr": false,
    "equal_var": true,
    "outDir": "examples/example_sigGenes",
    "baseStem": "johnson_18_GSE117444",
    "emitReverseP": true
  }

- Set ``fdr`` to ``true`` to apply Benjamini–Hochberg correction before thresholding.
- Adjust ``alpha`` for the significance threshold.

Use in Constrained Flux
-----------------------

Point your constrained config to the output folder and file suffixes, for example::

  {
    "data_dir": "./examples/example_sigGenes/",
    "uplabel": "upgenes",
    "dwlabel": "dwgenes",
    ...
  }

Then run::

  bash scooti/run_flux.sh examples/quickstart/configs/constrained.json

Single‑cell embryogenesis (findSigGenes)
----------------------------------------

Use the new demo to identify DEGs for embryogenesis transitions with `findSigGenes`:

::

  bash examples/identifySigGenes_demo/run_identify_siggenes.sh

Config keys include:
- ``sc_path``: 10x folder (``matrix.mtx``, ``genes.tsv``, ``barcodes.tsv``)
- ``sc_transitions``: list of {label, ref_regex, exp_regex} (e.g., Zygote→2cell)
- ``sc_method``: ``AVGSTD`` (default) or ``CI``; ``sc_std_num``: 2

Outputs are written under ``examples/identifySigGenes_demo/out/``.

Timing
------
- Bulk CSVs: 5–20 minutes
- Single‑cell transitions: 10–60 minutes (depends on cells/genes)
