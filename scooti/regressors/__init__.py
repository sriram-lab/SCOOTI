#__init__.py
# Subpackage exports (PyTorch learners are optional)

# Set a fixed random seed for reproducibility across regressors.
# Default seed is 8, can be overridden via env var SCOOTI_SEED.
import os
import random
import numpy as np

_SCOOTI_SEED = int(os.environ.get("SCOOTI_SEED", "8"))
try:
    random.seed(_SCOOTI_SEED)
    np.random.seed(_SCOOTI_SEED)
except Exception:
    pass

try:  # Torch is optional
    import torch
    torch.manual_seed(_SCOOTI_SEED)
    try:
        torch.cuda.manual_seed_all(_SCOOTI_SEED)
    except Exception:
        pass
except Exception:
    pass

from .regressorMetaLearner import *

# Optional exports (require PyTorch). Import softly so torch is not mandatory.
try:
    from .LassoTorch import *  # noqa: F401,F403
except Exception:
    pass

try:
    from .MLPRegressor import *  # noqa: F401,F403
except Exception:
    pass
