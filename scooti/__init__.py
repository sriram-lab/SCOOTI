"""SCOOTI package init."""

from ._version import version  # noqa: F401


def load():
    """Import the full SCOOTI API on demand."""
    import pandas as pd  # noqa: F401
    import json  # noqa: F401
    import numpy as np  # noqa: F401
    import matplotlib.pyplot as plt  # noqa: F401
    import seaborn as sns  # noqa: F401
    import cobra  # noqa: F401
    import sys  # noqa: F401
    import warnings
    import scipy.stats as ss  # noqa: F401
    import os  # noqa: F401
    from tqdm.notebook import tqdm, trange  # noqa: F401
    from adjustText import adjust_text  # noqa: F401
    from scipy.stats import mannwhitneyu  # noqa: F401

    warnings.simplefilter("ignore")

    # Set cobra solver to glpk in order to avoid err msg
    config = cobra.Configuration()
    config.solver = "glpk"

    import importlib

    def _export_public(module):
        public = getattr(module, "__all__", None)
        if public is None:
            public = [name for name in module.__dict__ if not name.startswith("_")]
        for name in public:
            globals().setdefault(name, module.__dict__[name])

    for module_name in [
        "GeneralMethods",
        "regressionAnalyzer",
        "metObjAnalyzer",
        "regressors",
        "utils",
    ]:
        module = importlib.import_module(f"{__name__}.{module_name}")
        globals()[module_name] = module
        _export_public(module)
