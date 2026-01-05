# Shim: re-export findSigGenes from utils to unify API and ensure get_transition_genes is available
try:
    from scooti.utils.findSigGenes import *  # re-export everything
except ImportError as e:
    raise
