# Inference of Metabolic Objectives and Trade‑offs using SCOOTI

## Summary
A step‑by‑step protocol to infer metabolic objectives and analyze trade‑offs from omics data using SCOOTI. The workflow covers DEG identification (bulk and single‑cell), flux modeling (unconstrained and constrained), objective inference, and downstream analysis including dimensionality reduction and trade‑off exploration.

## Before You Begin
- Innovation: SCOOTI couples omics‑constrained flux prediction with meta‑regression over idealized single‑objective fluxes to infer condition‑specific metabolic priorities, enabling interpretable trade‑off analyses (correlation/Pareto) across samples.
- Approvals: Not applicable for computational workflows. If using human/animal datasets, ensure institutional permissions for data use.

## Key Resources
- Software: MATLAB (Gurobi), Python 3.10, `scanpy`, `torch`, `cobra` (see `environment.yml`).
- GEM: `Shen2019.mat` (includes metabolite “accoa” in the objective list) under `./scooti/metabolicModel/GEMs/`.
- Demo data and configs under `examples/*_demo/`.

## Folder Structure (essentials)
```
SCOOTI/
├── environment.yml               # Core Python deps (torch optional)
├── requirements.txt              # pip deps (no torch)
├── README.md                     # Setup + quickstart + config examples
├── protocol.md                   # This protocol
├── scooti/
│   ├── run_flux.sh               # MATLAB flux pipeline launcher
│   ├── run_trainer.sh            # Objective inference launcher (Python)
│   ├── run_analyzer.sh           # Objective analysis launcher (Python)
│   ├── metabolicModel/
│   │   └── GEMs/
│   │       └── Shen2019.mat      # GEM used in demos (includes 'accoa')
│   ├── regressors/
│   │   ├── regressorMetaLearner.py  # sklearn-based learners
│   │   ├── LassoTorch.py            # optional (requires PyTorch)
│   │   └── MLPRegressor.py          # optional (requires PyTorch)
│   ├── metObjAnalyzer.py         # Analysis utilities
│   └── utils/                    # Helpers
└── examples/
    ├── quickstart/
    │   └── configs/              # unconstrained.json, constrained.json, inference.json
    ├── identifySigGenes_demo/    # bulk + single-cell DEG demo
    ├── unconstrained_demo/       # unconstrained flux demo (MATLAB)
    ├── constrained_demo/         # constrained flux demo (MATLAB)
    ├── inference_demo/           # objective inference demo (Python)
    ├── analyze_demo/             # legacy analyzer demo
    └── tradeoff_demo/            # legacy trade-off demo (Pareto/3D)
```

## Environment Setup
- Python (conda recommended):
```
conda env create -f environment.yml
conda activate scooti
```
- Optional: install PyTorch (only if using GPU learners like Lasso/MLP)
  - GPU (CUDA 11.8; bundles CUDA/cuDNN):
    - `pip install --force-reinstall --no-deps --index-url https://download.pytorch.org/whl/cu118 torch==2.0.1+cu118`
  - CPU only:
    - `pip install --index-url https://download.pytorch.org/whl/cpu torch==2.0.1`
- MATLAB: Install and license Gurobi; set COBRA toolbox path in JSON configs.

Timing: 30–90 minutes (once) depending on download/solver installs.

## 1. Identify Significant Genes/Proteins (DEGs)
Use the unified demo for bulk CSVs and single‑cell embryogenesis transitions.

- Run:
```
bash examples/identifySigGenes_demo/run_identify_siggenes.sh
```
- Config: `examples/identifySigGenes_demo/identify_siggenes_config.json`
  - Bulk: Sharma’21 (`sharma_csv`), Min’18 (`min_csv`), Johnson’18 (`qp_csv`)
  - Single‑cell: `sc_path` (10x folder), `sc_transitions` (e.g., Zygote→2cell, 2cell→32cell), `sc_method: AVGSTD`
- Outputs: up/down gene CSVs under `examples/identifySigGenes_demo/out/` (and reverse “_P_” files). Point constrained flux configs to this folder (`data_dir`) with `uplabel/dwlabel`.

Timing: 10–60 minutes (bulk per dataset); 10–60 minutes (single‑cell transitions).

Config (bulk: identifySigGenes_demo/identify_siggenes_config.json)
```json
{
  "out_dir": "./examples/identifySigGenes_demo/out/",
  "fdr_alpha": 0.05,
  "fc_threshold": 1.0,
  "sharma_csv": "/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sharma_21_prolif_pmid_34274479.csv",
  "min_csv": "/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/min_18_GSE122927_ReadCount.csv",
  "qp_csv": "/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/johnson_18_GSE117444_prolif_qui_count.csv",
  "sharma_groups": {
    "EGF":  {"exp": [2,6,10],  "ctrl": [1,5,9,15]},
    "TPA":  {"exp": [8,12],    "ctrl": [1,5,9,15]},
    "H89-EGF": {"exp": [13,16],  "ctrl": [7,11,17]},
    "H89-TPA": {"exp": [14,18],  "ctrl": [7,11,17]}
  },
  "min_pairs": [
    {"label": "p21_2E2",            "regex": "(?=.*p21_high)(?=.*2E2)|(?=.*p21_low)(?=.*2E2)"},
    {"label": "p21_3B6",            "regex": "(?=.*p21)(?=.*3B6)"},
    {"label": "SerumStarvation_2E2","regex": "(?=.*SerumStarvation)(?=.*2E2)|(?=.*2E2)(?=.*Control)"},
    {"label": "SerumStarvation_3B6","regex": "(?=.*SerumStarvation)(?=.*3B6)|(?=.*3B6)(?=.*Control)"},
    {"label": "Meki_2E2",           "regex": "(?=.*Meki)(?=.*2E2)|(?=.*2E2)(?=.*Control)"},
    {"label": "Meki_3B6",           "regex": "(?=.*Meki)(?=.*3B6)|(?=.*3B6)(?=.*Control)"},
    {"label": "CDK46i_2E2",         "regex": "(?=.*CDK46i)(?=.*2E2)|(?=.*2E2)(?=.*Control)"},
    {"label": "CDK46i_3B6",         "regex": "(?=.*CDK46i)(?=.*3B6)|(?=.*3B6)(?=.*Control)"},
    {"label": "ContactIn_2E2",      "regex": "(?=.*ContactIn)(?=.*2E2)|(?=.*2E2)(?=.*Control)"},
    {"label": "ContactIn_3B6",      "regex": "(?=.*ContactIn)(?=.*3B6)|(?=.*3B6)(?=.*Control)"}
  ],
  "qp_split": {"ctrl_n": 3, "exp_n": 3}
}
```

Config (single‑cell embryo: identifySigGenes_demo/identify_siggenes_scembryo_config.json)
```json
{
  "out_dir": "./examples/identifySigGenes_demo/out_sc/",
  "sc_path": "/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/single_cell/",
  "sc_label_split": "_",
  "sc_label_index": 0,
  "sc_method": "AVGSTD",
  "sc_std_num": 2,
  "sc_transitions": [
    {"label": "scEmbryo_sc1C2C", "ref_regex": "^Zygote", "exp_regex": "^2cell"},
    {"label": "scEmbryo_sc2CBC", "ref_regex": "^2cell",  "exp_regex": "^32cell"}
  ]
}
```

Runs
- Bulk/quiescence:
  - `bash examples/identifySigGenes_demo/run_identify_siggenes.sh examples/identifySigGenes_demo/identify_siggenes_config.json`
- Single‑cell embryo:
  - `bash examples/identifySigGenes_demo/run_identify_siggenes.sh examples/identifySigGenes_demo/identify_siggenes_scembryo_config.json`

## 2. Generate Unconstrained Flux Models (ideal objectives)
- Command:
```
bash examples/unconstrained_demo/run_unconstrained.sh
# or
bash scooti/run_flux.sh examples/unconstrained_demo/unconstrained_demo_config.json
```
- Notes: Uses `Shen2019.mat`, objective list `obj52_metabolites_shen2019.csv`, `medium: DMEMF12`.

Timing: 2–12 hours depending on solver and number of objectives.

Config (unconstrained_demo/unconstrained_demo_config.json)
```json
{
  "jj": 1,
  "COBRA_path": "~/cobratoolbox",
  "GEM_path": "./scooti/metabolicModel/GEMs/Shen2019.mat",
  "model_name": "Recon1",
  "obj_candidate_list_file": "./scooti/metabolicModel/GEMs/obj52_metabolites_shen2019.csv",
  "pfba": 1,
  "medium_perturbation": 0,
  "data_dir": "",
  "prefix_name": "model",
  "medium": "DMEMF12",
  "simulation": "CFR",
  "constraint": 0,
  "save_root_path": "./examples/unconstrained_demo/out/unconstrained_models/"
}
```

## 3. Generate Constrained Flux Models (omics‑informed)
- Command:
```
bash examples/constrained_demo/run_constrained.sh
# or
bash scooti/run_flux.sh examples/constrained_demo/constrained_demo_config.json
```
- Required fields: `data_dir` points to DEG CSVs, `uplabel/dwlabel`, `CFR_kappa: 0.1`, `CFR_rho: 10`, `medium: DMEMF12`.

Timing: 2–8 hours depending on samples and solver.

Config (constrained_demo/constrained_demo_config.json)
```json
{
  "jj": 1,
  "COBRA_path": "~/cobratoolbox",
  "GEM_path": "./scooti/metabolicModel/GEMs/Shen2019.mat",
  "model_name": "Recon1",
  "obj_candidate_list_file": "./scooti/metabolicModel/GEMs/obj52_metabolites_shen2019.csv",
  "input_obj_tb": "",
  "paraLen": 1,
  "random_para": 0,
  "init_objective": 1,
  "genekoflag": 0,
  "rxnkoflag": 0,
  "FSflag": 0,
  "pfba": 1,
  "medium_perturbation": 0,
  "data_dir": "./examples/example_sigGenes/",
  "prefix_name": "example",
  "medium": "DMEMF12",
  "uplabel": "upgenes",
  "dwlabel": "dwgenes",
  "simulation": "CFR",
  "constraint": 1,
  "save_root_path": "./examples/constrained_demo/out/constrained_models/",
  "CFR_model_path": "",
  "pairwise_CFR_model": 0,
  "extraWeight": 0,
  "algorithm": "cfr",
  "data_series": "",
  "prefix_series": "",
  "medium_series": "",
  "CFR_kappa": 0.1,
  "CFR_rho": 10,
  "DFA_kappa": -1
}
```

## 4. Infer Metabolic Objectives (meta‑regression)
- Command:
```
bash examples/inference_demo/run_inference.sh
# or
bash scooti/run_trainer.sh examples/inference_demo/demo_inference_config.json
```
- Notes: Ensure `medium` and `fileSuffix` match your flux outputs. Outputs written to `examples/inference_demo/out/regression_models/`.

Timing: 10–45 minutes (demo scale); up to a few hours for large cohorts.

Config (inference_demo/demo_inference_config.json)
```json
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
```

## 5. Analyze Objectives
- Legacy (Pareto overlays, trait analysis):
```
bash examples/analyze_demo/run_analyze.sh
```
- Notes: Label configuration supported; UMAP safety handled in code.

Timing: 15–60 minutes depending on data size.

Config (analyze_demo/analyze_config.json)
```json
{
  "flux_paths": {"exp": "./examples/constrained_demo/out/constrained_models/"},
  "coef_paths": {"exp": "./examples/inference_demo/out/regression_models/"},
  "save_root_path": "./examples/analyze_demo/out/",
  "engine": "legacy",
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
```

## 6. Trade‑off Analysis
- Legacy (Pareto overlays/3D):
```
bash examples/tradeoff_demo/run_tradeoff_legacy.sh
```
- Notes: Configure pareto_mets, correlation thresholds, UMAP parameters as needed.

Timing: 20–90 minutes depending on feature set.

Config (tradeoff_demo/tradeoff_config.legacy.json)
```json
{
  "engine": "legacy",
  "coef_paths": {"exp": "./examples/inference_demo/out/regression_models/"},
  "save_root_path": "./examples/tradeoff_demo/out/",
  "tradeoff_analysis": {
    "pareto_mets": ["glutathione", "biomass", "cholesterol"],
    "pareto2DAnalysis": true,
    "traitAnalysis": true
  }
}
```
## Troubleshooting
- If you select a GPU learner (`--learner Lasso` or `--learner MLP`) without PyTorch installed,
  the trainer will print an instruction to install the torch wheel via pip. Otherwise, use the
  default sklearn learner (`--learner L`).
- Torch import error on some HPC nodes: install the PyTorch cu118 wheel (see README) or use CPU wheel.
- No files found by trainer/analyzer: ensure `medium` and `fileSuffix` match the flux outputs.
- UMAP ARPACK errors (`k >= N`): use lower `umap_para`; built‑in auto‑clamping is enabled.

## Resource Availability
- Lead contact: Sriram Chandrasekaran (csriram@umich.edu)
- Technical contact: Da‑Wei Lin (daweilin@umich.edu)
- Materials availability: Not applicable.
- Data/code availability: SCOOTI code and demo configs are included in this repository; example data paths are specified in configs.

## Figures and Files
- Provide figures (300 dpi) as separate files upon submission; the analyzer and trade‑off demos produce publication‑quality plots.

## Highlights (examples)
- Legacy analyzer provides Pareto overlays and trait analysis.
- End‑to‑end pipeline from omics to objectives and trade‑offs.
- Single‑cell embryogenesis transitions supported via `findSigGenes`.
