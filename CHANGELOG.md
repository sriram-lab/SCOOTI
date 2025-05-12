# Development Log

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
