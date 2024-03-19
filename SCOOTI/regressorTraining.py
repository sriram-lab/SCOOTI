"""
regressorTraining.py
=================================================================
A module enables regression models to learn from metabolic fluxes
"""

# packages
import pandas as pd
import pickle
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
import sys
from SCOOTI.regressionAnalyzer import *
from SCOOTI.stat_tests import *
from SCOOTI.MatplotProp import CanvasStyle, PltProps, Significance
PltProps()
import warnings; warnings.simplefilter('ignore')
import scipy.stats as ss
import os
from tqdm.notebook import tqdm, trange
# Regression models
#from SCOOTI.regressorCollection import *
from SCOOTI.regressorMetaLearner import *

# Set cobra solver to glpk in order to avoid err msg
config = cobra.Configuration()
config.solver = "glpk"


class regressorTraining:


    def __init__(
            self,
            unconstrained_models_path,
            constrained_models_path,
            save_root_path,
            kappa_arr=[10, 1, 0.1, 0.01, 0.001],
            rho_arr=[10, 1, 0.1, 0.01, 0.001],
            expName='regression',
            uncon_norm=True,
            con_norm=False,
            medium='DMEMF12',
            method='cfr',
            model='recon1',
            input_type='flux',
            cluster_path='',
            rank=False,
            objList_path=''
            ):
        self.save_root_path = save_root_path
        self.input_type = input_type
        self.unconstrained_models_path = unconstrained_models_path
        self.constrained_models_path = constrained_models_path
        self.kappa_arr = kappa_arr
        self.rho_arr = rho_arr
        self.unconstrained_models = pd.DataFrame()
        self.constrained_models = pd.DataFrame()
        self.uncon_norm = uncon_norm
        self.con_norm = con_norm
        self.medium = medium
        self.expName = expName
        self.method = method
        self.cluster_path = cluster_path
        self.objList_path = objList_path
        # parameters for regressions
        #suffix = 'cfr_recon3d_dmemf12_paraScan_k10_r0.001'
        kappa_arr_s = f'{kappa_arr}'.replace(' ', '').replace(',', '_')
        rho_arr_s = f'{rho_arr}'.replace(' ', '').replace(',', '_')
        if rank:
            self.suffix = f'rank_{method}_{uncon_norm}_{con_norm}_{model}_{medium}_k{kappa_arr_s}_rho{rho_arr_s}'
        else:
            self.suffix = f'{method}_{uncon_norm}_{con_norm}_{model}_{medium}_k{kappa_arr_s}_rho{rho_arr_s}'

        # run regression
        self.metaLearner = self.run_training(rank=rank)

    @staticmethod
    def get_unconstrained_models(root_path, norm):
        print('Start processing the unconstrained models...')
        # unconstrained models
        # root_path = '/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/pfba/Recon3D/DMEMF12/ori_models/'
        uncon_res = unconstrained_models(
                root_path,
                norm=norm,
                medium='DMEMF12'
                )

        return uncon_res
   
    @staticmethod
    def get_constrained_models(root_path, norm, kappa_arr, rho_arr):
        print('Start processing the constrained models...')
        # GSE159929 scRNAseq datasets
        # root_paths = f'/nfs/turbo/umms-csriram/daweilin/fluxPrediction/Enrichr/Disease_Perturbations_Recon3D/'
        con_res = constrained_models(
                root_path+'/', 
                CFR_paraScan=True,
                norm=norm,
                CFR_k=kappa_arr,
                CFR_r=rho_arr,# input_path_pattern='NCI60'
                )
        return con_res

    
    def run_training(self, rank=True):
        # get unconstrained models
        self.unconstrained_models = self.get_unconstrained_models(
                self.unconstrained_models_path, self.uncon_norm
                )
        if len(self.objList_path)>0:
            objList = pd.read_csv(self.objList_path, index_col=0)
            self.unconstrained_models = self.unconstrained_models[objList.iloc[:, 0]]
        # get constrained models
        if self.method=='compass':
            # load fluxes
            self.constrained_models = pd.read_csv(self.constrained_models_path, index_col=0)
            # check if the dimension of flux data are the same
            overlaps = np.intersect1d(
                    self.constrained_models.index.to_numpy(),
                    self.unconstrained_models.index.to_numpy()
                    )
            print(len(self.unconstrained_models), len(overlaps))
            # just in case
            #con_res = con_res[con_res.index.isin(overlaps)].reindex(overlaps)
            #uncon_res = uncon_res[uncon_res.index.isin(overlaps)].reindex(overlaps)
        else:
            self.constrained_models = self.get_constrained_models(
                    self.constrained_models_path,
                    self.con_norm,
                    self.kappa_arr,
                    self.rho_arr
                    )

        # convert data into ranks
        if rank:
            self.unconstrained_models = self.unconstrained_models.rank()
            self.constrained_models = self.constrained_models.rank()

        # get data
        uncon_res = self.unconstrained_models
        con_res = self.constrained_models
        # parameters
        input_type = self.input_type
        save_root_path = self.save_root_path
        suffix = self.suffix
        cluster_path = self.cluster_path
        # execute regression
        from multiprocessing import Pool, cpu_count, set_start_method
        
        # parameter settings
        #suffix = 'cfr_recon3d_dmemf12_paraScan_k10_r0.001'
        #input_type = 'flux' # flux/pageRank/mixture
        suffix = self.suffix
        expName = self.expName
         
        if input_type=='pageRank':
            prefix = 'pg'
            print('Start loading models...')
            # loading datasets
            uncon_models, con_models = remove_all_zeros(uncon_res, con_res)
            pg_uncon_models, pg_con_models = load_pageRankData(suffix, expName)
            pg_uncon_models, pg_con_models = remove_all_zeros(pg_uncon_models, pg_con_models)
        
        elif input_type=='flux':
            prefix = 'flux'
            print('Start loading models...')
            # loading datasets
            uncon_models, con_models = remove_all_zeros(uncon_res, con_res)
            pg_uncon_models, pg_con_models = uncon_models, con_models
        else:
            prefix = 'st'
            print('Start loading models...')
            # loading datasets
            uncon_models, con_models = remove_all_zeros(uncon_res, con_res)
            pg_uncon_models, pg_con_models = load_pageRankData(suffix, expName)
            pg_uncon_models, pg_con_models = remove_all_zeros(pg_uncon_models, pg_con_models)
        
        
        print(f'{save_root_path}/{prefix}_rfelm_{self.expName}_{self.suffix}.csv')
        print('Start training regression models...')
        # convert fba results into mass flux graph
        #pg_res_df = process_Pandas_data(topology, ccle1, num_processes=36)
        import time   
        start = time.time()
        
        num_processes = min(con_models.shape[1], cpu_count())
        
        # 'with' context manager takes care of pool.close() and pool.join() for us
        with Pool(num_processes) as pool:
            
            # we need a sequence to pass pool.map; this line creates a generator (lazy iterator) of columns
            seq = [(
                uncon_models,
                pd.DataFrame(con_models[col_name]),
                expName,
                suffix,
                pg_uncon_models,
                pd.DataFrame(pg_con_models[col_name]),
                input_type,
                cluster_path
                ) for col_name in con_models.columns]
        
            # +-------------+
            # +Testing block+
            # +-------------+
            # model_training(uncon_models, con_models, expName, suffix, pg_uncon_models, pg_con_models, input_type)
            print('Mapping...')
            # pool.map returns results as a list
            results_list = pool.starmap(model_training, seq)
            if cluster_path=='':
                # return list of processed columns, concatenated together as a new dataframe
                RFElm_df = pd.concat((resl[0] for resl in results_list), axis=1)
                lasso_df = pd.concat((resl[1] for resl in results_list), axis=1)
                EN_df = pd.concat((resl[2] for resl in results_list), axis=1)
                LL_df = pd.concat((resl[3] for resl in results_list), axis=1)
                integrate_df = pd.concat((resl[4] for resl in results_list), axis=1)
            else:
                # return list of processed columns, concatenated together as a new dataframe
                cluster_dfs = []
                for jj in range(len(results_list[0])):
                   cluster_dfs.append(pd.concat((resl[jj] for resl in results_list), axis=1))

        
        end = time.time()

        if cluster_path=='':
            # output the models
            RFElm_df.columns = con_models.columns
            RFElm_df.to_csv(f'{save_root_path}/{prefix}_rfelm_{self.expName}_{self.suffix}.csv')
            
            lasso_df.columns = con_models.columns
            lasso_df.to_csv(f'{save_root_path}/{prefix}_lasso_{self.expName}_{self.suffix}.csv')
            
            EN_df.columns = con_models.columns
            EN_df.to_csv(f'{save_root_path}/{prefix}_EN_{self.expName}_{self.suffix}.csv')
            
            LL_df.columns = con_models.columns
            LL_df.to_csv(f'{save_root_path}/{prefix}_LL_{self.expName}_{self.suffix}.csv')
            
            integrate_df.columns = con_models.columns
            integrate_df.to_csv(f'{save_root_path}/{prefix}_sl_{self.expName}_{self.suffix}.csv')
        else:
            for jj in range(len(cluster_dfs)):
                cluster_dfs[jj].columns = con_models.columns
                if jj==len(cluster_dfs)-1:
                    cluster_dfs[jj].to_csv(f'{save_root_path}/{prefix}_cluster_weight_{self.expName}_{self.suffix}.csv')
                elif jj==len(cluster_dfs)-2:
                    cluster_dfs[jj].to_csv(f'{save_root_path}/{prefix}_sl_{self.expName}_{self.suffix}.csv')
                    integrate_df = cluster_dfs[jj]
                else:
                    cluster_dfs[jj].to_csv(f'{save_root_path}/{prefix}_cluster{jj}_{self.expName}_{self.suffix}.csv')


        return integrate_df
    

#unconModel='/nfs/turbo/umms-csriram/daweilin/fluxPrediction/unconstrained_models/pfba/DMEMF12/'
#conModel='/nfs/turbo/umms-csriram/daweilin/fluxPrediction/BulkRNAseq/AVGSTD_1/NCI60/'
#savePath='/nfs/turbo/umms-csriram/daweilin/regression_models/BulkRNAseq_paraScan/'
#kappaArr=[0.001]
#rhoArr=[10]
#expName='NCI60'
#medium='DMEMF12'
#method='cfr'
#model='recon1'
#inputType='flux'
#
#
#regression_exe = regressorTraining(
#    unconModel,
#    conModel,
#    savePath,
#    kappa_arr=kappaArr,
#    rho_arr=rhoArr,
#    expName=expName,
#    uncon_norm=True,
#    con_norm=False,
#    medium=medium,
#    method=method,
#    model=model,
#    input_type=inputType
#)







