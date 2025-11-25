#!/usr/bin/env python3
import pandas as pd
from scooti.metObjAnalyzer import metObjAnalyzer

# function to change labels
def label_func(df):
    return pd.Series(df.columns).apply(
        lambda x: x.split('_')[0]
    ).replace({'1C':'sc1C2C', '2C':'sc2CBC'})

# initiate the object (paths assume running from repo root)
moa = metObjAnalyzer(
    flux_paths={
        'sc1C2C':'./examples/constrained_demo/out/constrained_models/',
        'sc2CBC':'./examples/constrained_demo/out/constrained_models/',
    }, # optional if you only want to analyze objectives
    coef_paths={
        # Update these to your produced coefficients if different
        'sc1C2C':'./examples/endToend_demo/out/regression_models/flux_sc1C2C_input_norm_outcome_nonorm_k0.1_r0.01.csv',
        'sc2CBC':'./examples/endToend_demo/out/regression_models/flux_sc2CBC_input_norm_outcome_nonorm_k0.1_r0.01.csv',
    },
    save_root_path='./examples/analyze_demo/out_scEmbryo/',
    GEM_path='./scooti/metabolicModel/GEMs/Shen2019.mat',
    uncon_model_path='./examples/unconstrained_demo/out/unconstrained_models/',
    col_map={'sc1C2C':'sc1C2C', 'sc2CBC':'sc2CBC'},
    label_func=label_func,
    sel_para='k0.1_r0.01',
    prefix='scEmbryo',
    medium='KSOM',
)

coefs = moa.get_coef()
try:
    zeros = coefs.columns[(coefs == 0).all(0)]
    print('Columns with zero coefficients:', zeros)
except Exception:
    pass
