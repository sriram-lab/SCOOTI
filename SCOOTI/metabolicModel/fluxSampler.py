"""
regressorDiseaseAnalysis_stackModels.py
=======================================================
Analysis of metabolic objectives and fluxes in diseases
"""

# for ipython
%load_ext autoreload
%autoreload 2

#packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
import warnings; warnings.simplefilter('ignore')
import scipy.stats as ss
from tqdm.notebook import tqdm, trange
import sys, os


# Set cobra solver to glpk in order to avoid err msg
config = cobra.Configuration()
config.solver = "glpk"

# set up the number of samples
sample_num = 3


def processing_metadata():
    # Save settings in metadata
    metadata['obj'] = obj
    metadata['obj_type'] = obj_type
    metadata['obj_c'] = obj_c
    metadata['output_path'] = root_path
    metadata['input_path'] = ''
    metadata['file_name'] = out_name
    metadata['with_constraint'] = ctrl
    metadata['CFR_kappa'] = np.nan
    metadata['CFR_rho'] = np.nan
    metadata['medium'] = medium
    metadata['genekoflag'] = 0
    metadata['rxnkoflag'] = 0
    metadata['media_perturbation'] = 0
    obj_check = model.c!=0 # find objective functions
    if length(obj_check)~=0,
      metadata.objWeights = model.c(obj_check);
      metadata.objRxns = model.rxns{obj_check};
    else
      metadata.objWeights = 0;
      metadata.objRxns = '';
    end
    metadata.model_path = model_path;
    metadata.upStage = upsheet;
    metadata.dwStage = dwsheet;

    % Load reconstructed models
    if length(CFR_model)>0,
      recon_model = load(CFR_model);
      fn = fieldnames(recon_model);
      recon_model = getfield(recon_model, fn{1});
    else,
      recon_model = '';
    end
    metadata.CFRModel = CFR_model;
    metadata.extraWeight = extra_weight;
    metadata.algorithm = algorithm;
    
%%   convert structure to json files
    encodedJSON = jsonencode(metadata);
    JSONFILE_name= sprintf('%s_metadata.json', excelname);
    fid = fopen(JSONFILE_name,'w');
    fprintf(fid, encodedJSON);
    fclose('all')





from cobra import Model, Reaction, Metabolite
# path to access the metabolic model
GEM_path='./GEMs/Shen2019.mat'
gem = cobra.io.load_matlab_model(GEM_path)


# change culture medium
media = pd.read_excel('FINAL_MEDIUM_MAP_RECON1.xlsx',sheet_name='DMEMF12')
EX_lbs = media.iloc[:,2]
EX_rxns = media.iloc[:,5]
with gem:
    medium = gem.medium
    for EX_rxn, new_lb in zip(EX_rxns, EX_lbs):
        try:
            medium[EX_rxn] = new_lb
            gem.medium = medium
        except Exception as e:
            print(e)


# get objective candidates
objectives = pd.read_csv('./obj52_metabolites_shen2019.csv', index_col=0)
compartments = ['c', 'm', 'n', 'x', 'r', 'g', 'l']
# create new folder for each sampled flux models
rootpath = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/fva/'
for i in range(sample_num):
    if not os.path.exists(rootpath+f'fva_{i}'):
        os.makedirs(rootpath+f'fva_{i}')

# sampled fluxes
fluxSamples = {}
for obj in objectives['metabolites']:
    if obj=='gh':
        obj = 'biomass_objective'
        gem.objective = obj
    else:
        demand_rxn_mets = {}
        for c in compartments:
            met = f'{obj}[{c}]'
            try:
                gem.metabolites.get_by_id(met)
                demand_rxn_mets[met] = -1
            except Exception as e:
                print('skip')
                print(e)
        
        reaction = Reaction(f'{obj}_demand')
        reaction.name = f'Objective candidate {obj}'
        gem.add_reactions([reaction])
        reaction.add_metabolites(demand_rxn_mets)
        gem.objective = f'{obj}_demand'

    #print(gem.objective.expression)
    optgp = cobra.sampling.OptGPSampler(gem, seed=1)
    optgp_samples = optgp.sample(sample_num)
    optgp_samples = pd.DataFrame(optgp_samples, columns=[r.id for r in gem.reactions])
    optgp_samples = optgp_samples.sample(frac=1)
    optgp_samples.index = np.arange(len(optgp_samples))
    optgp_samples = optgp_samples.T
    for col in optgp_samples.columns:
        df = pd.DataFrame(optgp_samples[col])
        df.to_csv(
                '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/fva/fva_{col}/'
                )


  
