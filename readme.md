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

This will install all required packages, including `numpy`, `pandas`, `torch`, `scanpy`, `cobra`, and others.

## Option 2: Install via pip

If you prefer pip, you can install manually:

```
pip install -r requirements.txt
```
Or if installing SCOOTI as a local Python package:
```
pip install .
```
Make sure your Python version is 3.10.


# MATLAB Environment Setup (optional if you don't need metabolic modeling)

You need MATLAB installed with:
1. Optimization Toolbox (`gurobi` is required)
2. Parallel Computing Toolbox (optional, for speedup)

Inside MATLAB:
```
addpath(genpath('SCOOTI/metabolicModels/'));
addpath(genpath('SCOOTI/metabolicModels/utils/'));
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

Example for unconstrained models:

```
bash ./SCOOTI/run_flux.sh ./SCOOTI/example/unconstrained_demo_config.json
```

Example for constrained models:

```
bash ./SCOOTI/run_flux.sh ./SCOOTI/example/constrained_demo_config.json
```

This will:
- Load a config .json
- Launch CFRinterface(config) to infer fluxes.

# Run Regression Training (Python)

Example:
```
bash example/run_trainer.sh example/demo_trainer_config.json
```

This will:
- Train regressors with unconstrained and constrained fluxes.


## Troubleshooting Tips

| Issue | Solution |
|:---|:---|
| ImportError: numpy.dtype size changed | Make sure numpy and numba versions match (use our `environment.yml`). |
| Undefined function (e.g. `validate_config`) | Make sure you added `./SCOOTI/metabolicModels/` to MATLAB path. |
| Python package conflicts | Please make sure you installed packages listed in `requirements.txt` |

# Project Structure Overview
SCOOTI/
├── setup.py               # Python packaging file
├── environment.yml        # Conda environment setup
├── requirements.txt       # pip-based dependency list
├── README.md               # Project overview and usage
├── setup_instructions.md   # Detailed installation and usage guide
│
├── metabolicModels/        # MATLAB modeling code
│   ├── utils/              # Helper MATLAB functions
│   └── (core .m files like CFRinterface.m, multiObj_CBM.m)
│
├── SCOOTI/                 # Python modules
│   ├── __init__.py
│   ├── trainer/            # Training scripts (SCOOTI_trainer.py, regressorTraining.py)
│   └── (other Python modules)
│
├── example/                # Example configs and scripts
│   ├── run_flux.sh         # Bash script to run MATLAB flux pipeline
│   ├── run_trainer.sh      # Bash script to run Python trainer
│   ├── demo_flux_config.json
│   └── demo_trainer_config.json

Please cite [this paper](https://doi.org/10.1016/j.cels.2024.12.005) if you use SCOOTI in your work.



