"""Analysis sinlge-cell embryogenesis

This script is to analyze the metabolic objectives of single-cell embryogenesis.
Note that the flux prediction was based on the algorithm COMPASS instead of the linear version of iMAT.
"""


#packages
import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
import sys
sys.path.append('/home/daweilin/StemCell/GeneralMethods/')
from regressionAnalyzer import *
from stat_tests import *
from AnalysisKits import *
from MatplotProp import CanvasStyle, PltProps
PltProps()
import warnings; warnings.simplefilter('ignore')
import scipy.stats as ss
import os
from tqdm.notebook import tqdm, trange
from adjustText import adjust_text
# Regression models
from regressorCollection import *

# Set cobra solver to glpk in order to avoid err msg
config = cobra.Configuration()
config.solver = "glpk"


from scipy.stats import mannwhitneyu

# +++++++++++++++++++
# + Data processing +
# +++++++++++++++++++
# access to the data
coef_paths = {
        'scEmbryo':'./regression_models/scEmbryo/flux_sl_scZygote2CBC_input_nonorm_outcome_nonorm_compass_recon1.csv',
        }


# get fluxes of multi-objective models with constraints
coefs = {}
for k in coef_paths.keys():
    path = coef_paths[k]
    coefs[k] = pd.read_csv(path, index_col=0)
    #coefs[k].index = pd.Series(coefs[k].index).replace(mnames)

# Merge
coefs_df = pd.concat((coefs[k] for k in coefs.keys()), axis=1)
# Separate dataframes by parameters
coef_sel = coefs_df.copy()
coef_sel = coef_sel[
            coef_sel.any(axis=1)
            ]
para = 'k0.1_r0.01'
# get pca results and labels
labels = pd.Series(coef_sel.columns).apply(
        lambda x: x.split('_')[0]
        )
labels = ['Zygote' if 'Zygote' in ele else ele for ele in labels]
print(labels)


# +++++++++++++++++++++++++++++++++++
# + Distances to biomass objectives +
# +++++++++++++++++++++++++++++++++++
# distance to biomass
biomass_model_coef = getRxnCoefficients(
        model_path='./models/Shen2019.mat',
        model_name='Recon1',
        )
# Get labels
coef_dist = coef_distance_to_biomassObj(
        coef_sel, labels, biomass_model_coef, 'scEmbryo',
        norm=True, rank=False, func='euclidean',
        save_root_path='./result_figures/',
        histplot=True, boxplot_cols=[],
        boxplot_order=['Zygote', '2cell', '32cell']
        )
