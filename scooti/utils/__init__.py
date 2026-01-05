#__init__.py
# import essentail modules for shortening the code for other modules
from .stat_tests import *
from .MatplotProp import CanvasStyle, PltProps, Significance
from .findSigGenes import *
from .dataPreprocessing import *
try:
    from .plot_func import *
except Exception as e:
    # Optional plotting dependencies (umap, hdbscan, phate, etc.) may be missing or incompatible.
    # Defer failures to when plotting is actually invoked.
    print(f"[scooti.utils] Optional plotting imports unavailable: {e}")
