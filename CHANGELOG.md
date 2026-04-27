# Development Log

## [2026-04-27]

### Fixed
- Fixed all-metabolite objective candidate scans in `scooti/metabolicModel/CFRinterface.m` by skipping demand-objective candidates that do not create any valid objective reaction in the selected GEM. This allows broad candidate lists such as all-metabolite lists to continue through valid metabolites instead of failing on absent metabolites.
- Updated `scooti/run_flux.sh` so MATLAB exceptions return a non-zero shell exit code and print the MATLAB error report, making SCOOTI CLI failures visible to scripts and batch jobs.

### Notes
- All-metabolite flux scans still require a working MATLAB solver setup. On Great Lakes, load MATLAB and Gurobi before running flux modeling, for example: `module load matlab/R2024b gurobi/10.0.2`.

## [2025-05-12]

### Fixed
- Fixed the bugs in `./SCOOTI/regressors/LassoTorch.py`

### Changed
- Moved supporting functions from `./SCOOTI/metabolicModel/utils/`
- Moved `regressorMetaLearner.py` to `./SCOOTI/regressors/`
- Modularized MATLAB-based flux modeling in `./SCOOTI/metabolicModel/utils/`
- Changed `./SCOOTI/metabolicModel/CFRinterface.m` to input `config`
- Changed `./SCOOTI/metabolicModel/multiObj_CBM.m` to input `config`

### Added
- Added CLI tools `./SCOOTI/run_flux.sh` and `./SCOOTI/run_trainer.sh`
- Added example `.json` file for unconstrained models `./SCOOTI/examples/run_flux/unconstrained_demo_config.json`
- Added example `.json` file for constrained models `./SCOOTI/examples/run_flux/constrained_demo_config.json`
- Added example `.json` file for objective inference `./SCOOTI/examples/run_inference/demo_inference_config.json`
- Added example lists of significant genes `./SCOOTI/examples/example_sigGenes/`
- Added example flux predictions `./SCOOTI/examples/example_fluxPreduction/`

### Planned
- Modularize `./SCOOTI/metabolicModel/DFAinterface.m`
- Update the uses of `DFAinterface.m` in `multiObj_CBM.m`

## [2025-04-28]

### Fixed
- Fixed the bugs in `./SCOOTI/metabolicModel/multiObj_CBM.m`
- Fixed the bugs in `./SCOOTI/metabolicModel/CFRinterface.m`

### Changed
- Copied all files from `./SCOOTI/GeneralMethods/` to `./SCOOTI/utils/`

### Added
- Added directory `./SCOOTI/utils/`
    - Added `plot_func.py` (optimized Python classes for visualization)
    - Added `fluxModels.py` (optimized Python class for loading flux data)
    - Added `dataPreprocessing.py` (optimized Python class for normalization and imputation)
- Added directory `./SCOOTI/regressors/`
    - Added `LassoTorch.py` (Lasso regression using GPU deployed by PyTorch)
    - Added `MLPRegressor.py` (MLP-regression using GPU deployed by PyTorch)
