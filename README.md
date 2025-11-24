![ICON](./icon.png)
# SCOOTI
Welcome to `SCOOTI: Single-Cell Optimization Objective and Tradeoff Inference`.
This guide walks you through setting up the environment and running SCOOTI's MATLAB and Python modules smoothly.

# Python Environment Setup

## Option 1: Install via conda (Recommended)

We provide an environment.yml file.

```
conda env create -f environment.yml
conda activate scooti
```

Note: Different versions of Gurobi and MATLAB (and solver settings) can lead
to different feasible flux distributions that reach the same objective
value. Consequently, inferred objective coefficients may differ numerically
across environments, while qualitative trends should remain similar.

This will install all required packages, including `numpy`, `pandas`, `scanpy`, `cobra`, and others. PyTorch is optional and only needed for GPU‑based learners (e.g., LassoTorch/MLP).

Note (Ubuntu 22.04): If `conda env create -f environment.yml` fails due to phate/adjustText, use pip instead (Option 2): `pip install -r requirements.txt` or `pip install .`.

### Optional: Install PyTorch (GPU or CPU)

If you plan to use the GPU‑based learners, install the official PyTorch wheel after creating the env. We recommend the cu118 build (bundles CUDA/cuDNN):

```
# Install PyTorch 2.0.1 with CUDA 11.8 wheel (bundles CUDA/cuDNN)
pip install --force-reinstall --no-deps --index-url https://download.pytorch.org/whl/cu118 torch==2.0.1+cu118

# Verify
python - <<'PY'
import torch
print('torch', torch.__version__, 'cuda', torch.version.cuda, 'avail', torch.cuda.is_available())
import torch.backends.cudnn as cudnn
print('cudnn', cudnn.version())
PY
```


## Option 2: Install via pip

If you prefer pip, you can install manually (PyTorch remains optional):

```
pip install -r requirements.txt
# Then (optionally) install PyTorch (choose ONE of the following):
# CPU-only wheel
# pip install --index-url https://download.pytorch.org/whl/cpu torch==2.0.1
# GPU wheel (CUDA 11.8)
# pip install --index-url https://download.pytorch.org/whl/cu118 torch==2.0.1+cu118
```
Or if installing SCOOTI as a local Python package:
```
pip install .
```
Make sure your Python version is 3.10.

## Import SCOOTI
``
python -c "from scooti._version import version; print(version)"
``


# MATLAB Environment Setup (optional if you don't need metabolic modeling)

You need MATLAB installed with:
1. Optimization Toolbox (`gurobi` is required)
2. Parallel Computing Toolbox (optional, for speedup)

Inside MATLAB:
```
addpath(genpath('scooti/metabolicModel/'));
addpath(genpath('scooti/metabolicModel/utils/'));
```
This loads all necessary `.m` functions (e.g., `CFRinterface`, `run_simulation`).


# Python Package Development Setup

If you are developing SCOOTI (not just using):

- Install in editable mode:
```
pip install -e .
```
Now any code changes will take effect immediately without reinstalling.


# Running the Modules
## Run Flux Modeling (MATLAB)

Before you start, please revise the `.json` files to access your cobratoolbox:
```
"COBRA_path": ".YOUR_COBRATOOLBOX_PATH", ...
```


Quickstart (one command):
```
bash examples/quickstart/run_all.sh
```

Configurations (JSON)

Unconstrained (examples/quickstart/configs/unconstrained.json)
```json
{
  "jj": 1,
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
  "save_root_path": "./examples/quickstart/out/unconstrained_models/"
}
```

Constrained (examples/quickstart/configs/constrained.json)
```json
{
  "jj": 1,
  "COBRA_path": "~/cobratoolbox",
  "GEM_path": "./scooti/metabolicModel/GEMs/Shen2019.mat",
  "model_name": "Recon1",
  "obj_candidate_list_file": "./scooti/metabolicModel/GEMs/obj52_metabolites_shen2019.csv",
  "data_dir": "./examples/run_flux/example_sigGenes/",
  "prefix_name": "ProQui",
  "medium": "DMEMF12",
  "uplabel": "upgenes",
  "dwlabel": "dwgenes",
  "simulation": "CFR",
  "constraint": 1,
  "save_root_path": "./examples/quickstart/out/constrained_models/",
  "kappa": 0.1,
  "rho": 10,
  "DFA_kappa": -1
}
```

Inference (examples/quickstart/configs/inference.json)
```json
{
  "unconModel": "./examples/quickstart/out/unconstrained_models/",
  "conModel": "./examples/quickstart/out/constrained_models/",
  "savePath": "./examples/quickstart/out/regression_models/",
  "kappaArr": "10,1,0.1",
  "rhoArr": "10,1,0.1",
  "dkappaArr": "10,1,0.1",
  "expName": "quickstart",
  "unconNorm": "T",
  "conNorm": "F",
  "medium": "DMEMF12",
  "method": "cfr",
  "model": "recon1",
  "inputType": "flux",
  "fileSuffix": "_fluxes.csv.gz"
}
```

Example for unconstrained models:

```
bash scooti/run_flux.sh examples/run_flux/unconstrained_demo_config.json
```

Example for constrained models:

```
bash scooti/run_flux.sh examples/run_flux/constrained_demo_config.json
```

This will:
- Load a config .json
- Launch CFRinterface(config) to infer fluxes.

# Run Regression Training (Python)

Example:
```
bash scooti/run_trainer.sh examples/run_inference/demo_inference_config.json
```

This will:
- Train regressors with unconstrained and constrained fluxes.

# Prepare Unconstrained Models (Ideal Objectives)

We include a minimal demo to generate “unconstrained” flux models that optimize
each objective independently (no omics-based constraints).

Demo run:

```
bash examples/unconstrained_demo/run_unconstrained.sh
#!/usr/bin/env bash
# Or run the full quickstart (generates flux then inference)
bash examples/quickstart/run_all.sh
```

Outputs are written under:

```
examples/unconstrained_demo/out/unconstrained_models/
```

Each objective produces a ``*_fluxes.csv(.gz)`` file. These serve as predictors for
the downstream meta-regression.

Real use (one‑liner):

```
bash scooti/run_flux.sh examples/unconstrained_demo/unconstrained_demo_config.json
```


## Troubleshooting Tips

| Issue | Solution |
|:---|:---|
| ImportError: numpy.dtype size changed | Make sure numpy and numba versions match (use our `environment.yml`). |
| Undefined function (e.g. `validate_config`) | Make sure you added `./SCOOTI/metabolicModels/` to MATLAB path. |
| Python package conflicts | Please make sure you installed packages listed in `requirements.txt` |

# Identifying Significant Genes or Proteins

We include a small demo to generate up-/down-regulated gene (or protein) lists from an example
bulk RNA-seq table using a two-sample t-test (optionally with FDR) and a
fold-change criterion. The output CSVs can be used directly by the constrained
flux pipeline.

Demo run (uses the Johnson dataset under `examples/example_omics`):

```
conda activate scooti
bash examples/siggenes_demo/run_siggenes.sh examples/siggenes_demo/siggenes_config.json
```

This writes four CSVs into `examples/example_sigGenes/`:

- `johnson_18_GSE117444_upgenes.csv`
- `johnson_18_GSE117444_dwgenes.csv`
- `johnson_18_GSE117444_P_upgenes.csv`
- `johnson_18_GSE117444_P_dwgenes.csv`

Point your constrained config to the folder and suffixes:

```
"data_dir": "./examples/example_sigGenes/",
"uplabel": "upgenes",
"dwlabel": "dwgenes"
```

Then run constrained flux as shown above.

Real use (one‑liners):

```
# Unconstrained (demo wrapper or direct)
bash examples/unconstrained_demo/run_unconstrained.sh
# or
bash scooti/run_flux.sh examples/unconstrained_demo/unconstrained_demo_config.json

# Inference (demo wrapper or direct)
bash examples/inference_demo/run_inference.sh
# or
bash scooti/run_trainer.sh examples/inference_demo/demo_inference_config.json
```

# Project Structure Overview
```
SCOOTI/
├── scooti/                     # Python package and entry scripts
│   ├── __init__.py
│   ├── metabolicModel/         # MATLAB code (flux modeling)
│   │   └── utils/
│   ├── regressors/             # Python regressors (objective inference)
│   ├── GeneralMethods/         # Shared utilities
│   ├── utils/                  # Python helpers and plotting
│   ├── SCOOTI_trainer.py       # Objective inference entry (Python)
│   ├── SCOOTI_sampler.py
│   ├── run_flux.sh             # Flux pipeline runner (MATLAB)
│   ├── run_trainer.sh          # Trainer runner (Python)
│   └── examples -> ../examples # Back-compat symlink to top-level examples
│
├── examples/                   # Example configs and tiny data
│   ├── run_flux/
│   │   ├── constrained_demo_config.json
│   │   └── unconstrained_demo_config.json
│   └── run_inference/
│       └── demo_inference_config.json
│
├── docs/                       # Sphinx documentation
├── environment.yml             # Conda environment (see README note)
├── requirements.txt            # pip-based dependency list
├── setup.py                    # Packaging
├── pyproject.toml
└── README.md
```
Please cite [this paper](https://doi.org/10.1016/j.cels.2024.12.005) if you use SCOOTI in your work.


# Analyze Inferred Objectives

SCOOTI provides two analysis engines that read the coefficient tables from the inference step and produce figures/tables:

- minimal: lightweight, minimal dependencies; safe UMAP with automatic fallback; optional HDBSCAN clustering; entropy, distance to biomass, pairwise comparisons, metabolite-type grouping, and trade-off CSVs/plots.
- legacy: full feature set via `metObjAnalyzer`/`regressionAnalyzer` (e.g., Pareto front overlays), requires more optional dependencies.

Recommended: start with the minimal engine.

Run the demo (minimal):

```
conda activate scooti
bash examples/analyze_demo/run_analyze_minimal.sh
```

Or the legacy analyzer:

```
bash examples/analyze_demo/run_analyze.sh
```

Internals and configuration

- Both demos call the generic runner: `bash scooti/run_analyzer.sh <config.json>`.
- The JSON has an `engine` key: `"minimal"` or `"legacy"`.
- UMAP safety (both engines): when projecting to 2D, the code inspects the number of samples `N` and clamps `n_components` and `n_neighbors` to `N-1` to avoid SciPy ARPACK errors (`k >= N`). You’ll see a log message like: `[UMAP] Adjusted n_components from 50 to 28 (N=29).`
- Minimal engine options (subset):
  - `clustering`, `entropy`, `distance`: booleans controlling analysis outputs.
  - `reduction`: `"umap"`, `"pca"`, or `"auto"` (try UMAP, else PCA).
  - `umap_para`: `[n_neighbors, n_components]` (auto-clamped by N-1).
  - `cluster_umap`: boolean to run HDBSCAN on the 2D embedding and save labels.
  - `metType_cluster` + `metType_map`: aggregate objectives by metabolite category and cluster heatmap.
  - `tradeoff` + `pareto_mets`: compute negative-correlation trade-offs and plot pairwise scatters for listed objectives.


# Advanced: Trade-off Analysis

SCOOTI supports exploring pairwise and multi-objective trade-offs between metabolic objectives.

Two demo paths are provided under `examples/tradeoff_demo/`:

- Minimal trade-off demo (correlations + pairwise scatters):

```
bash examples/tradeoff_demo/run_tradeoff_minimal.sh
```

Config highlight (`examples/tradeoff_demo/tradeoff_config.minimal.json`):

```
{
  "engine": "minimal",
  "coef_paths": {"exp": "./examples/inference_demo/out/regression_models/"},
  "save_root_path": "./examples/tradeoff_demo/out/",
  "tradeoff": true,
  "pareto_mets": ["glutathione", "biomass", "cholesterol"],
  "tradeoff_corr_threshold": -0.4,
  "umap_para": [10, 50]
}
```

Outputs include:
- Correlation-based trade-off CSV (`*_tradeoffs.csv`)
- Pairwise scatter plots for the listed objectives (`*_tradeoff_<A>_vs_<B>.png`)
- UMAP/PCA scatter and embedding CSV; optional clustering labels (`cluster_umap: true`).

- Legacy trade-off demo (Pareto analysis via `metObjAnalyzer.tradeoff_analysis`):

```
bash examples/tradeoff_demo/run_tradeoff_legacy.sh
```

Config highlight (`examples/tradeoff_demo/tradeoff_config.legacy.json`):

```
{
  "engine": "legacy",
  "tradeoff_analysis": {
    "pareto_mets": ["glutathione", "biomass", "cholesterol"],
    "pareto2DAnalysis": true,
    "traitAnalysis": true
  }
}
```

This path can generate correlation matrices, objective pair scatter plots, and Pareto front overlays; 3D Pareto surfaces are supported when enabled.


## Labeling (config.json)

Scatter plots (UMAP/PCA) and related outputs use labels to color points. By default, labels come
from the `coef_paths` keys (e.g., `"exp"`). You can override labeling without writing code by adding
the `labels` section to your config JSON. Supported modes:

- mode: `"group_key"` (default): label by the `coef_paths` key for each column.
- mode: `"column"`: use the exact column names (original sample/condition IDs).
- mode: `"regex"`: extract labels from column names using a regex pattern.
  - `regex`: pattern string (e.g., `"^(.*?)_"`)
  - `regex_group`: which capture group to use (default: 1)
- mode: `"split"`: split column names by a delimiter and pick a field.
  - `split_delim`: delimiter string (default: `_`)
  - `split_index`: integer index (default: 0)
- Optional remap: `group_map` is a dict to rename resulting labels.

Example snippet:

```
{
  "engine": "minimal",
  "coef_paths": {"exp": "./examples/inference_demo/out/regression_models/"},
  "labels": {
    "mode": "regex",
    "regex": "^(.*?)_",
    "regex_group": 1,
    "group_map": {"exp": "Experiment"}
  }
}
```

Notes
- Works for both engines. For the legacy engine, these overrides are applied after coefficients are loaded.
- If `labels` is omitted, both engines default to `group_key` labeling.

# Identify Significant Genes (bulk + single‑cell)

We provide a unified demo to generate up-/down-regulated gene lists from bulk CSVs (Sharma’21, Min’18, Johnson’18) and single‑cell embryogenesis transitions using `findSigGenes`.

Quick run:
```
bash examples/identifySigGenes_demo/run_identify_siggenes.sh
```
Config (bulk + single‑cell) at `examples/identifySigGenes_demo/identify_siggenes_config.json`:
- Bulk: `sharma_csv`, `min_csv`, `qp_csv` with groups/regex definitions
- Single‑cell: `sc_path` (10x folder), `sc_transitions` (e.g., Zygote→2cell, 2cell→32cell), `sc_method: AVGSTD`

Outputs (CSV):
- `.../<prefix>_upgenes.csv`, `.../<prefix>_dwgenes.csv` and their reverse “_P_” counterparts
- For constrained flux, point `data_dir` to this output folder and set `uplabel/dwlabel`

Timing: 10–60 minutes depending on dataset size and CPU.
