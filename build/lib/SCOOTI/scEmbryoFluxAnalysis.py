"""
scEmbryoFluxAnalysis.py
==============================================================
Constrained model analysis for single-cell mouse embryogenesis
- The model predicts the fluxes of ideal optimality
- We aim to model fluxes that can separate cell types
- Clustering performance was applied to the clusters of UMAP of flux predictions
- Best parameter combination is firstly determined by silhouette score
- If the score is all negative, mutual info and rand index were considered
"""
#packages
import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
import warnings; warnings.simplefilter('ignore')
import scipy.stats as ss
import os
from tqdm.notebook import tqdm, trange
from adjustText import adjust_text
import sys, os
from SCOOTI.metObjAnalyzer import metObjAnalyzer

# function to change labels
def label_func(df):
    return pd.Series(df.columns).apply(
        lambda x: x.split('_')[0]
        ).replace({'1C':'sc1C2C', '2C':'sc2CBC'})

# initiate the object
method = 'sl' # superLearner
moa = metObjAnalyzer(
        flux_paths={
        'sc1C2C':'/nfs/turbo/umms-csriram/daweilin/fluxPrediction/sc1C2C/paraScan/',
        'sc2CBC':'/nfs/turbo/umms-csriram/daweilin/fluxPrediction/sc2CBC/paraScan/',
            },
        coef_paths={},
        save_root_path='/home/daweilin/StemCell/Project_mESC_JinZhang/test_output/',
        GEM_path='/nfs/turbo/umms-csriram/daweilin/data/models/Shen2019.mat',
        uncon_model_path='/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/pfba/KSOM/',
        col_map={'sc1C2C':'sc1C2C', 'sc2CBC':'sc2CBC'},
        label_func=label_func,
        samplingFlux_path='/nfs/turbo/umms-csriram/daweilin/fluxPrediction/RandomObjCoef/embryo_ksom_50percent/',
        sel_para='',
        prefix='scEmbryo',
        medium='KSOM',
      )

# flux analysis
moa.fluxAnalysis(kappa_arr=[10, 1, 0.1, 0.01, 0.001], rho_arr=[10, 1, 0.1, 0.01, 0.001])
