# Development Log

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
