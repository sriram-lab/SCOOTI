Repository Structure
====================

Overview of the SCOOTI repository layout after recent cleanups.

Top-level
---------
::

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
  ├── examples/                   # Example configs (runnable)
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

Notes
-----
- MATLAB code is kept inside the Python package under ``scooti/metabolicModel`` so that packaged distributions may include it.
- Examples are moved to a top-level ``examples/`` directory and symlinked back to ``scooti/examples`` for backward compatibility.
- Two entry scripts are provided for minimal workflows:
  - ``scooti/run_flux.sh`` — runs the MATLAB-based flux modeling pipeline.
  - ``scooti/run_trainer.sh`` — runs the Python-based objective inference pipeline.

Quick Commands
--------------
Run flux modeling with the unconstrained example::

  bash scooti/run_flux.sh examples/run_flux/unconstrained_demo_config.json

Run objective inference example::

  bash scooti/run_trainer.sh examples/run_inference/demo_inference_config.json

