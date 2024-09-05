"""
regressorDiseaseAnalysis_stackModels.py
=======================================================
Analysis of metabolic objectives and fluxes in diseases
"""

# for ipython
#%load_ext autoreload
#%autoreload 2

#packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
import warnings; warnings.simplefilter('ignore')
import scipy.stats as ss
from tqdm.notebook import tqdm, trange
import sys, os, json
from datetime import datetime
import time
import cobra.flux_analysis.parsimonious as pFBA
from cobra.sampling import sample
# Set cobra solver to glpk in order to avoid err msg
config = cobra.Configuration()
config.solver = "glpk"

# set up the number of samples
sample_num = 10


def processing_metadata(
        candidate,#model,
        objectives,
        obj_c,
        root_path,
        model_path,
        data_path='',
        sample_name='',
        out_name='',
        upsheet='',
        dwsheet='',
        ctrl=0,
        kappa=1,
        rho=1,
        medium='DMEMF12',
        genekoflag=False,
        rxnkoflag=False,
        media_perturbation=False,
        ):
    # Save settings in metadata
    metadata = {}
    #objname = str(model.objective.expression).split('*')[1].split(' ')[0]
    metadata['obj'] = list(objectives)
    metadata['obj_type'] = 'demand'
    metadata['obj_c'] = list(obj_c)
    metadata['output_path'] = root_path
    if ctrl==1:
      metadata['input_path'] = data_path
    else:
      metadata['input_path'] = sample_name
    metadata['file_name'] = out_name
    metadata['with_constraint'] = ctrl
    metadata['CFR_kappa'] = kappa
    metadata['CFR_rho'] = rho
    metadata['medium'] = medium
    metadata['genekoflag'] = genekoflag
    metadata['rxnkoflag'] = rxnkoflag
    metadata['media_perturbation'] = media_perturbation
    metadata['objWeights'] = 1
    metadata['objRxns'] = candidate#str(objname)
    metadata['model_path'] = model_path
    metadata['upStage'] = upsheet
    metadata['dwStage'] = dwsheet
    

    # manage file names
    # File name with datetime as prefix
    current_time = datetime.now()
    file_prefix = time.strftime("%b%d%Y%H%M%S")
    # Output file
    filename = f'{root_path}/[{file_prefix}]{out_name}'
    print('The CFR result has been saved in', filename);
    # output file
    excelname = filename
    # convert structure to json files
    with open(f'{excelname}_metadata.json', 'w') as f:
        json.dump(metadata, f)

    return excelname
    





def objective_setting_functionPy(
        root_path='/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/fluxSampling/',
        GEM_path='/home/daweilin/StemCell/Project_mESC_JinZhang/SCOOTI/SCOOTI/metabolicModel/GEMs/Shen2019.mat'
        ):

    from cobra import Model, Reaction, Metabolite
    # path to access the metabolic model
    gem = cobra.io.load_matlab_model(GEM_path)
    # change culture medium
    media = pd.read_excel('/home/daweilin/StemCell/Project_mESC_JinZhang/SCOOTI/SCOOTI/metabolicModel/FINAL_MEDIUM_MAP_RECON1.xlsx',sheet_name='DMEMF12')
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
    objectives = pd.read_csv('/home/daweilin/StemCell/Project_mESC_JinZhang/SCOOTI/SCOOTI/metabolicModel/obj52_metabolites_shen2019.csv', index_col=0)
    objectives = objectives.iloc[1:,:] # testing
    compartments = ['c', 'm', 'n', 'x', 'r', 'g', 'l']


    # sampled fluxes
    gem_tmp = gem.copy()
    objective_candidates = []
    for i, obj in enumerate(objectives['metabolites']):
        # attributes for the metadata
        #obj_c = np.zeros(len(objectives['metabolites']))
        #obj_c[i] = 1
        # create a copy of the model
        if obj=='gh':
            obj = 'biomass_objective'
            objective_candidates.append(obj)
            #gem_tmp.objective = obj
        else:
            demand_rxn_mets = {}
            for c in compartments:
                met = f'{obj}[{c}]'
                try:
                    gem_tmp.metabolites.get_by_id(met)
                    demand_rxn_mets[met] = -1
                except Exception as e:
                    print('skip')
                    print(e)
            
            reaction = Reaction(f'{obj}_demand')
            reaction.name = f'Objective candidate {obj}'
            gem_tmp.add_reactions([reaction])
            reaction.add_metabolites(demand_rxn_mets)
            objective_candidates.append(f'{obj}_demand')
            print(len(gem_tmp.reactions))
            #print(gem_tmp.reactions.get_by_id(f'{obj}_demand'))
    
        #print(gem.objective.expression)
    for i, candidate in enumerate(objective_candidates):
        gem_tmp2 = gem_tmp.copy()
        pFBA.add_pfba(gem_tmp2, candidate)
        #gem_tmp.objective = candidate
        obj_c = np.zeros(len(objectives['metabolites']))
        obj_c[i] = 1
        #optgp = cobra.sampling.OptGPSampler(gem_tmp2, seed=1)
        #optgp.n_samples = 100
        #optgp_samples = optgp.sample(sample_num)
        optgp_samples = sample(gem_tmp2, 1000, processes=30)
        rxns = [r.id for r in gem_tmp2.reactions]
        optgp_samples = pd.DataFrame(
                optgp_samples,
                columns=rxns
                )
        optgp_samples = optgp_samples.sample(frac=1)
        optgp_samples['Obj'] = optgp_samples[candidate].to_numpy()
        optgp_samples.index = np.arange(len(optgp_samples))
        optgp_samples = optgp_samples.T
        print(optgp_samples)
        for col in optgp_samples.columns:
            # create new folder for each sampled flux models
            rootpath = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/fluxSampling1000/'
            if not os.path.exists(rootpath+f'fs_{col}'):
                os.makedirs(rootpath+f'fs_{col}')
            # get metadata
            out_name = f'model_ct1_obj{i+1}_data1'
            excelname = processing_metadata(
                    candidate,#gem_tmp2,
                    objectives['metabolites'].to_numpy(),
                    obj_c,
                    f'/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/fluxSampling1000/fs_{col}/',
                    GEM_path,
                    data_path='',
                    sample_name='',
                    out_name=out_name,
                    upsheet='',
                    dwsheet='',
                    ctrl=0,
                    kappa=1,
                    rho=1,
                    medium='DMEMF12',
                    genekoflag=False,
                    rxnkoflag=False,
                    media_perturbation=False,
                    )
            # save flux data
            df = pd.DataFrame(optgp_samples[col])
            df.columns = ['upgene']
            print(df)
            df.to_csv(
                    f'{excelname}_fluxes.csv.gz', compression='gzip'
                    )

    
if __name__=="__main__":
    objective_setting_functionPy()
