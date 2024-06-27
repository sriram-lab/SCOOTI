"""
metObjAnalyzer.py
==================================================================================
Object-oriented classes that serve as a user interface to run analysis/simulations
"""

#packages
import pandas as pd
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
from adjustText import adjust_text
# Regression models
#from SCOOTI.regressorCollection import *
from SCOOTI.regressorMetaLearner import *

# Set cobra solver to glpk in order to avoid err msg
config = cobra.Configuration()
config.solver = "glpk"


from scipy.stats import mannwhitneyu


class metObjAnalyzer:
    """Analysis of optimized fluxes and metabolic objectives.

    The object consists of four parts of analysis. 
    1) Optimized flux: 
        suggests the best parameters and visualization.
    2) Metabolic objective: 
        compare to different samples and allocations.
    3) Tradeoffs and traits: 
        identify anticorrelated metabolites and metabolic traits.
    4) Flux reconstruction:
        reconstruct flux prediction with inferred coefficients.
    
    
    Parameters
    ----------
    flux_paths : {dictionary},
        Paths to access the directories with flux prediction data.
        Keys should be the labels (sample names) associated with the directories.
    
    coef_paths : {dictionary},
        Paths to access the files (.csv) with inferred metabolic objectives.
        Keys should be the labels (sample names) associated with the files.
    
    save_root_path : {str},
        Path to save the outputs (including files and figures).
        Absolute path is recommended.

    GEM_path : {str},
        Paths to access the genome-scale metabolic models (.mat files).
        Only Shen et al., Recon2.2, and Recon3D are applicable for this version.

    uncon_model_path: {str},
        Paths to access the directories with flux predictions of single-objective models.

    col_map : {dictionary},
        Serialize unique column names (samples) and corresponding preferred names.

    label_func : {function},
        A function to convert column names into unique sample names.

    samplingFlux_path : {str},
        Paths to access the directories with flux predictions 
        models with randomly sampled objective coefficients.

    sel_para : {str},
        A parameter combination selected for CFR (e.g. k0.1_r0.01)

    prefix : {str},
        Project name (e.g. scEmbryo)
   
    medium : {string}, options=['KSOM', 'DMEMF12']
        medium used to model the fluxes with FBA
    
    Attributes
    ----------
    flux_df : {pandas.DataFrame},
        Flux data predicted with constrained models without objectives

    coef_df : {pandas.DataFrame},
        Coefficients of metabolic objectives

    labels : {pandas.DataFrame},
        label for each column of flux or coefficient data
    
    
    Examples
    --------
    >>> from SCOOTI.metObjAnalyzer import metObjAnalyzer
    >>> # set up a function to replace column names
    >>> def label_func(df):
    >>>     return pd.Series(df.columns).apply(
    >>>             lambda x: x.split('_')[0]
    >>>         ).replace({'1C':'sc1C2C', '2C':'sc2CBC'})
    >>> # initiate the object
    >>> moa = metObjAnalyzer(
    >>>     flux_paths={
    >>>         'sc1C2C':'./fluxPrediction/sc1C2C/paraScan/',
    >>>         'sc2CBC':'./fluxPrediction/sc2CBC/paraScan/',
    >>>                 },
    >>>     coef_paths={
    >>>         'sc1C2C':f'./regression_models/scEmbryo_paraScan/flux_sl_sc1C2C_input_norm_outcome_nonorm_k0.1_r0.01.csv',
    >>>         'sc2CBC':f'./regression_models/scEmbryo_paraScan/flux_sl_sc2CBC_input_norm_outcome_nonorm_k0.1_r0.01.csv',
    >>>         },
    >>>     save_root_path='/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
    >>>     GEM_path='./models/Recon1.mat',
    >>>     uncon_model_path='./fluxPrediction/unconstrained_models/pfba/KSOM/',
    >>>     col_map={'sc1C2C':'sc1C2C', 'sc2CBC':'sc2CBC'},
    >>>     label_func=label_func,
    >>>     samplingFlux_path='./fluxPrediction/RandomObjCoef/embryo_ksom_50percent/',
    >>>     sel_para='k0.1_r0.01',
    >>>     prefix='scEmbryo',
    >>>     medium='KSOM',
    >>>     )
    >>> # Analyze fluxes and look for the best parameter combinations
    >>> moa.fluxAnalysis(
    >>>         kappa_arr=[10, 1, 0.1],
    >>>         rho_arr=[10, 1, 0.1],
    >>>     )
    >>> # read coefficients
    >>> moa.get_coef()
    >>> # analyze coefficients
    >>> coefAnalysis(
    >>>     clustering=True,
    >>>     entropy=True,
    >>>     distance=True,
    >>>     umap_para=[5, 50]
    >>>     )
    """

    def __init__(
            self,
            flux_paths,
            coef_paths,
            save_root_path,
            GEM_path,
            uncon_model_path,
            col_map={},
            label_func=None,
            samplingFlux_path='',
            sel_para='',
            prefix='',
            medium='KSOM',
            GEM_name=''
            ):
        # save input parameters into class attributes
        self.flux_paths = flux_paths
        self.coef_paths = coef_paths
        self.save_root_path = save_root_path
        self.label_func = label_func
        self.samplingFlux_path = samplingFlux_path
        self.sel_para = sel_para
        self.col_map = col_map
        self.prefix = prefix
        self.medium = medium
        self.GEM_path = GEM_path
        self.uncon_model_path = uncon_model_path
        self.GEM_name = GEM_name

        # class attributes that will be updated after calling functions
        self.flux_scores = {'Coverage':[], 'MI':[], 'RI':[], 'SI':[], 'kappa':[], 'rho':[], 'label':[]}
        self.clustering_scores = pd.DataFrame()
        self.flux_df = pd.DataFrame()
        self.coef_df = pd.DataFrame()
        self.labels = pd.Series()

        # reconstructed fluxes
        self.wrxns = pd.DataFrame()
        self.wobjs = pd.DataFrame()

    def get_unconstrained_model(self, sel_mets=[], flux_plot=False, corr_plot=False, obj_sanity_check=False):
        """Sanity checks and analysis of unconstrained models

        Three goals of the function:
        1) Output flux plots of the unconstrained models for publications and sanity check
        2) Compute correlations of unconstrained models of metabolites
        3) Sanity checks of objective values which should align with the corresponding maximized values
        
        Parameters
        ----------
        sel_mets : {array-like}, default=[],
            selected metabolites for publishing tables.
            It will output the entire tables if leaving the parameter default.
        
        flux_plots : {bool}, default=False,
            Plot the table of fluxes for selected metabolites with clustermaps.
        
        corr_plots : {bool}, default=False,
            Plot the correlation of fluxes between unconstrained models with clustermaps.

        obj_sanity_check : {bool}, default=False,
            Make plots of objective values of each unconstrained model with clustermaps.
            It separates the plots based on the category of metabolites.


        Returns
        -------
        self : object
            returns self

        uncon_res : {pandas.DataFrame} of shape (n_reactions, n_metabolites),
            The data table of flux predictions generated from either pFBA or FBA without constraints.
            
        uncon_res_tmp : {pandas.DataFrame} of shape (n_innercellular_reactions, n_metabolites),
            The data table of flux predictions is actually the `uncon_res` without exchange reactions.

        """
        medium = self.medium
        # import unconstrained models (without sampling)
        uncon_res = unconstrained_models(self.uncon_model_path, norm=True, medium=f'{medium}')
        uncon_res = uncon_res[uncon_res.columns[uncon_res.any(axis=0)]]
        uncon_res = uncon_res[uncon_res.any(axis=1)]
        #uncon_res.to_csv(f'{self.save_root_path}/unconstrained_models_pFBA_{self.medium}_{self.GEM_name}.csv')
        
        # output the flux predictions of selected metabolites
        sel_mets = uncon_res.columns.to_numpy() if len(sel_mets)==0 else sel_mets
        uncon_res_tmp = uncon_res[sel_mets][uncon_res[sel_mets].any(axis=1)]
        uncon_res_tmp = uncon_res_tmp.sort_values(by=list(sel_mets))
        # remove exchange fluxes
        uncon_res_tmp[(pd.Series(uncon_res_tmp.index).str.contains('EX_')==0).to_numpy()]
        # get demand fluxes (objective candidate fluxes)
        uncon_dm = unconstrained_models(self.uncon_model_path, return_variables=False, norm=True, medium=f'{medium}')
        uncon_dm = uncon_dm[uncon_dm.columns[uncon_dm.any(axis=0)]]
        uncon_dm.index = pd.Series(uncon_dm.index).apply(lambda x: x.split('[')[0])
        uncon_dm = uncon_dm.groupby(uncon_dm.index).sum()
        
        # add a column of metabolite types
        unconT = uncon_dm.copy().T
        unconT['metTypes'] = unconT.index.map(mets_category())
        unconT = unconT.sort_values(by=['metTypes'])

        # clustermap of flux prediction
        if flux_plot==True:
            sns.set(font_scale=2.5)
            cmap = sns.diverging_palette(220, 20, as_cmap=True)
            sns.clustermap(uncon_res_tmp, cmap=cmap, figsize=(12,40))
            plt.savefig(f'{self.save_root_path}/unconstrained_models_rxns_{medium}.png')
        
        if corr_plot==True:
            # clustermap of flux correlation
            sns.set(font_scale=2.5)
            cmap = sns.diverging_palette(220, 20, as_cmap=True)
            sns.clustermap(uncon_res.corr(), cmap=cmap, figsize=(40,40))
            plt.savefig(f'{self.save_root_path}/unconstrained_models_corrmap_{medium}.png')

        if obj_sanity_check==True:   
            # make plots of objective fluxes by categories of metabolite types
            for metType in unconT['metTypes'].unique():
                # clustermap of flux correlation
                metdf = unconT[unconT['metTypes']==metType]
                metdf = metdf[metdf.columns[pd.Series(metdf.columns).str.contains('h2o|Obj|metType')==0]]
                metdf = metdf[metdf.columns[metdf.any(axis=0)]]
                print(metdf)
                kws = dict(cbar_kws=dict(ticks=[int(0), np.round(metdf.max().max(), 1)], orientation='horizontal'))
                sns.set(font_scale=2.5)
                g = sns.clustermap(metdf, cmap='viridis', figsize=(3+metdf.shape[1],10), **kws)
                x0, _y0, _w, _h = g.cbar_pos
                g.ax_cbar.set_position([x0, 0.9, g.ax_row_dendrogram.get_position().width, 0.02])
                g.ax_cbar.tick_params(axis='x')#, length=10)
                plt.savefig(f'{self.save_root_path}/uncon_demand_flux_norm_{metType}.png')
        
        return uncon_res, uncon_res_tmp, unconT


    def output_table(self):
        """Export dataframe into a .csv file for publications or other usages

        Save .csv files of fluxes or coefficients
        
        Parameters
        ----------
        self : {object}
            the function will gather the flux data and/or coefficient data for outputs

        Returns
        -------
        self : object
            returns self 

        """
        # output flux models
        if len(self.flux_df)>0:
            self.flux_df.to_csv(f'{self.save_root_path}/{self.prefix}_fluxPrediction_published.csv')
        # output coefficient models
        if len(self.coef_df)>0:
            self.coef_df.to_csv(f'{self.save_root_path}/{self.prefix}_objectiveCoef_published.csv')
        
        

    def getSpecificSamples(self, labels, sel_label_arr, update_attr=False):
        """Slice table with columns associated with correspondingly selected labels
        
        Users select specific columns via given labels to chop the data table.
        Since the column names are not the same as the label names,
        selected labels will compare with the entire array of labels to index the columns.

        Parameters
        ----------
        labels : {array-like} of shape (n_columns,)
            name array that should has the same length of columns

        sel_label_arr : {array-like}
            names chosen from labels that should be unique

        update_attr : {bool}, default=False
            update class attributes of flux/coefficients dataframe and labels
            if update_attr is set true


        Returns
        -------
        self : object
            returns self 

        """
        # get specific cells
        remain_inds = [(labels==label) for label in sel_label_arr]
        remain_inds = sum(remain_inds)
        
        remain_labels = np.array(labels)[remain_inds]
        try:
            remain_flux = self.flux_df[self.flux_df.columns[remain_inds]]
            if update_attr:
                self.labels = remain_labels
                self.flux_df = remain_flux
        except:
            print('Cannot slice flux data.')
        try:
            remain_coef = self.coef_df[self.coef_df.columns[remain_inds]]
            if update_attr:
                self.labels = remain_labels
                self.coef_df = remain_coef
        except:
            print('Cannot slice coefficient data.')
    


    @staticmethod
    def label_setup(df, label_func):
        # get labels
        return label_func(df)


    def get_flux(self, kappa=1, rho=1, rank=False, stack_model=False):
        """Load inferred metabolic objectives.

        The coefficients of metabolic objectives were obtained by
        the unconstrained models with single objectives regressiing on
        the condition/cell-specific constrained models

        Parameters
        ----------
        kappa_arr : {array-like} of shape (n_parameters,), default=[10, 1, 0.1, 0.01, 0.001]
            parameter of interest.

        rho_arr : {array-like} of shape (n_parameters,), default=[10, 1, 0.1, 0.01, 0.001]
            parameter of interest.


        Returns
        -------
        self : object
            Returns self.
        """
        # collect models
        fdfs = {}
        for k in self.flux_paths.keys():
            print(k)
            res = load_multiObj_models(
                    self.flux_paths[k], medium=self.medium,
                    return_variables=True, norm=False, stack_model=stack_model,
                    CFR_paraScan=True, CFR_k=[kappa], CFR_r=[rho],
                    file_suffix='_fluxes.csv.gz'
                    )
            fdfs[k] = res
        
        # get fluxes of constrained models without objectives
        flux_df = pd.concat((fdfs[k] for k in fdfs.keys()), axis=1)
        flux_df = flux_df[flux_df.columns[flux_df.any(axis=0)]].replace(
                [np.inf, -np.inf], [0, 0]
                )
        if rank:
            flux_df = flux_df.rank()

        # get labels
        self.labels = self.label_setup(flux_df, self.label_func)
        print(self.labels)

        self.flux_df = flux_df



    def fluxAnalysis(
            self,
            kappa,
            rho,
            umap_para=[5, 10],
            label=''
            ):
        """Cluster fluxes of UMAP with HDBSCAN and output clustering perfomance.

        Cell-specific constrained models without objective functions
        1) The model predicts the fluxes of ideal optimality
        2) We aim to model fluxes that can separate cell types

        Parameters
        ----------
        kappa_arr : {array-like} of shape (n_parameters,), default=[10, 1, 0.1, 0.01, 0.001]
            parameter of interest.

        rho_arr : {array-like} of shape (n_parameters,), default=[10, 1, 0.1, 0.01, 0.001]
            parameter of interest.

        umap_para : {array-like}, default=[5, 10]
            the first number of the array stands for the number of neighbors and
            the second number of the array is the total dimension of UMAP


        Returns
        -------
        self : object
            Returns self.

        flux_df : {pandas.DataFrame} of shape (n_reactions, n_samples), optional
            The data table of flux predictions generated from either CFR or DFA.
            Only provided if `self.sel_para` is given in the class attributes.

        scores : {pandas.DataFrame} of shape (n_parameters, 4), optional
            The table of clustering evaluations which comes with coverage rate of HDBSCAN,
            rand index, mutual information, and silhouette scores.
            It is only provided if `self.sel_para` is not given in the class attributes.
            A .csv file with the evaluations is output at the same time.
        """
        # Evaluation of constrained models w/o objectives
        # clustering methods
        cf = clustering_func(
                    self.flux_df,
                    self.save_root_path,
                    f'{self.prefix}_Flux_{kappa}_{rho}',
                    mets_category(),
                )
        
        # initialize the clustering objective
        cf.corr_clustermap(self.labels, show_cluster=False)
        # call umap function to reduce the dimension of data
        umaps = cf.reduction_scatter(
                self.labels, continuous=False, func='UMAP', para=umap_para
                )
        
        # cluster the umap results
        recluster = cf.reclustering(umaps, ['UMAP1', 'UMAP2'], min_size=10)
        clustered_res = recluster[recluster['cluster']>=0] # remove outliers
        if len(clustered_res):
            # evaluation 
            scoreSI = cf.clustering_evaluation(
                    clustered_res, 'Cell type', 'cluster', method='SI'
                    )
            scoreRI = cf.clustering_evaluation(
                    clustered_res, 'Cell type', 'cluster', method='RI'
                    )
            scoreMI = cf.clustering_evaluation(
                    clustered_res, 'Cell type', 'cluster', method='MI'
                    )
            self.flux_scores['Coverage'].append(len(clustered_res)/len(recluster))
            self.flux_scores['SI'].append(scoreSI)
            self.flux_scores['RI'].append(scoreRI)
            self.flux_scores['MI'].append(scoreMI)
            self.flux_scores['kappa'].append(kappa)
            self.flux_scores['rho'].append(rho)
            self.flux_scores['label'].append(label)

        # clustering report
        clustering_scores = pd.DataFrame(self.flux_scores)
        self.clustering_scores = clustering_scores

        return clustering_scores

    def paraScanFlux(
            self,
            kappa_arr=[10, 1, 0.1, 0.01, 0.001],
            rho_arr=[10, 1, 0.1, 0.01, 0.001],
            umap_para=[5, 10],
            extra_suffix=''
            ):
        """Scan parameters for the evaluation of clustering/separations of fluxes

        This function is only applicable when the flux data `self.flux_df` is loaded.
        It only performs a for-loop to scan parameters and save the evaluation results.


        Parameters
        ----------
        kappa_arr : {array-like} of shape (n_parameters,), default=[10, 1, 0.1, 0.01, 0.001]
            parameter of interest.

        rho_arr : {array-like} of shape (n_parameters,), default=[10, 1, 0.1, 0.01, 0.001]
            parameter of interest.

        
        Returns
        -------
        self : object
            Returns self.
        clustering_score_df : {pandas.DataFrame}

        Notes
        -----
        The function will update the evaluation score for clustering in class attribute

        """
        # iterate thru all the parameters
        for kappa in kappa_arr:
            for rho in rho_arr:
                # execute the function of flux analysis
                clustering_res = self.fluxAnalysis(
                        kappa, rho, umap_para=umap_para
                        )
        # save the evaluation results in the end
        self.save_clustering_result(extra_suffix)

    def save_clustering_result(self, extra_suffix=''):
        # save results
        print('Save the clustering evaluation...')
        self.flux_scores.to_csv(
                f'{self.save_root_path}/{self.prefix}_clustering_scores{extra_suffix}.csv'
                )
        

    def get_coef(self, metType_cluster=False):
        """Load inferred metabolic objectives.

        The coefficients of metabolic objectives were obtained by
        the unconstrained models with single objectives regressiing on
        the condition/cell-specific constrained models

        Parameters
        ----------
        kappa_arr : {array-like} of shape (n_parameters,), default=[10, 1, 0.1, 0.01, 0.001]
            parameter of interest.

        rho_arr : {array-like} of shape (n_parameters,), default=[10, 1, 0.1, 0.01, 0.001]
            parameter of interest.


        Returns
        -------
        self : object
            Returns self.
        """
        # get fluxes of multi-objective models with constraints
        coefs = {}
        for k in self.coef_paths.keys():
            path = self.coef_paths[k]
            coefs[k] = pd.read_csv(path, index_col=0)
        
        # Merge
        coefs_df = pd.concat((coefs[k] for k in coefs.keys()), axis=1)
        # remove metabolites that were not selected
        coef_sel = coefs_df.copy()
        coef_sel = coef_sel[
                    coef_sel.any(axis=1)
                    ]
        # remove columns with 0 coefficients
        keep_cols = coef_sel.columns[coef_sel.any(axis=0)]
        print(
                'Columns with zero coefficients:', coef_sel.columns[coef_sel.any(axis=0)==0]
                )
        coef_sel = coef_sel[keep_cols]
        coef_sel = coef_sel[coef_sel.abs()>1E-16].fillna(0) # remove metabolites with extremely low coef
        # get labels
        labels = self.label_setup(coef_sel, self.label_func)
        print(labels)

        self.coef_df = coef_sel
        self.labels = labels
        
        if metType_cluster:
            # cluster the metabolites by their types
            metDict = mets_category()
            self.coef_df['metTypes'] = metDict[self.coef_df.index]
            self.coef_df = self.coef_df.groupby('metTypes').sum()
        

    def get_randomObj_coef(self, th=50, sample_size=100000):
        """get random coefficients of designated objectives

        Enable randomly sampling coefficients of selected metabolites based on how many portion
        of cells choose a metabolite as a feature. The coefficients will be saved into a .csv file

        Parameters
        ----------
        th : {float}, default=50
            the threshold to select metabolites as features. the number represents the percentage of
            cells that choose the metabolites; thus, it ranges from [0, 100]
        sample_size : {integer}, default=100000
            number of random samples used to generate coefficients.

        Returns
        -------
        None
        """
        # Separate dataframes by parameters
        coef_sel = self.coef_df.copy()
        coef_sel = coef_sel[
                    coef_sel.any(axis=1)
                    ]
        # get selected metabolites
        coef_ratio = coef_sel.div(coef_sel.sum(axis=0), axis=1)
        single_obj = np.array([]) # initialize
        for unique_label in np.unique(self.labels):
            coef_unique_label = coef_ratio[coef_ratio.columns[self.labels==unique_label]]
            single_obj = np.unique(np.append(
                    single_obj,
                    coef_unique_label.index[(coef_unique_label>0).mean(axis=1)>=th],
                    ))
        
        # sampling and save
        df = pd.DataFrame(np.random.rand(len(single_obj), sample_size))
        df.index = single_obj
        df.columns = [f'Sample_{ind}' for ind in np.arange(len(df.columns))]
        df.to_csv(f'{self.save_root_path}/samplingObjCoef_{self.prefix}_{th}percent.csv')
        coef = pd.read_csv(f'{self.save_root_path}/samplingObjCoef_{self.prefix}_{th}percent.csv')

        return coef


        
    def coefAnalysis(
            self,
            norm=True,
            unknown_clustering=True,
            clustering=True,
            entropy=True,
            distance=True,
            compare=True,
            umap_para=[5, 50],
            recluster_min_size=10,
            method='average',
            ):
        """Analysis of inferred metabolic objectives.

        The coefficients of metabolic objectives were analyzed by
        1) comparing to different samples,
        2) dimension reduction methods for metabolic objectives,
        3) objective entropy (degree of order),
        4) distances from the biomass objective to the inferred objectives  

        Parameters
        ----------
        norm : {bool}, default=True
            normalize the coefficients by the sum of all coefficients of a sample.

        clustering : {bool}, default=True
            enable dimension reduction methods or clustering methods for metabolic objectives.
            UMAP is the default method in the class, but it is allowed to be swapped by
            PCA, tSNE, and PHATE.

        entropy : {bool}, default=True
            enable the calculations of entropy based on the coefficients of metabolic objectives
            which follows the equation derived from `p(log(p))`. `p` is converted from 

        distance : {bool}, default=True
            enable dimension reduction methods or clustering methods for metabolic objectives

        umap_para : {array-like}, default=[5, 10]
            the first number of the array stands for the number of neighbors and
            the second number of the array is the total dimension of UMAP

        recluster_min_size : {int}, default=10
            the parameter for HDBSCAN to cluster the UMAP results

        method : {str}, default='average'
            method used to cluster samples with scipy.cluster.hierarchy.linkage. 
            Options include 'average', 'single', 'complete', 'weighted', 'centroid'.
            For more details, one can visit the official document.
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html

        Returns
        -------
        self : object
            Returns self.

        Notes
        -----
        The figures of the analysis will be saved in the path given by `self.save_root_path`.
        """
        # get unique column names
        cols = np.unique(self.labels)
        # +++++++++++++++++++++++++++++++++
        # + clustering for unknown labels +
        # +++++++++++++++++++++++++++++++++
        if unknown_clustering:
            # copy the coef data
            plot_df = self.coef_df.copy() 
            if norm:
                plot_df = plot_df.div(plot_df.sum(axis=0), axis=1)
                cf = clustering_func(
                            plot_df,
                            self.save_root_path,
                            f'{self.prefix}_coef_norm_{self.sel_para}',
                            mets_category(),
                        )
            else:
                cf = clustering_func(
                            plot_df,
                            self.save_root_path,
                            f'{self.prefix}_coef_{self.sel_para}',
                            mets_category(),
                        )
            
            # clustering
            dn = cf.corr_clustermap(self.labels, show_cluster=True, method=method)
            
            # initialize the clustering objective
            # cf.corr_clustermap(self.labels)
            # call umap function to reduce the dimension of data
            umaps = cf.reduction_scatter(
                    self.labels, continuous=False, func='UMAP', para=umap_para
                    )
        
            # cluster the umap results
            recluster = cf.reclustering(umaps, ['UMAP1', 'UMAP2'], min_size=recluster_min_size)

            return recluster, dn

        # ++++++++++++++
        # + Clustering +
        # ++++++++++++++
        if clustering:
            # initiate the objective for clustering and analysis
            plot_df = self.coef_df.copy()#[coef_sel.index!='h2o']
            if norm:
                plot_df = plot_df.div(plot_df.sum(axis=0), axis=1)
                cf = clustering_func(
                            plot_df,
                            self.save_root_path,
                            f'{self.prefix}_coef_norm_{self.sel_para}',
                            mets_category(),
                        )
            else:
                cf = clustering_func(
                            plot_df,
                            self.save_root_path,
                            f'{self.prefix}_coef_{self.sel_para}',
                            mets_category(),
                        )
            
            # clustering
            cf.corr_clustermap(self.labels, show_cluster=False)
            #cf.reduction_scatter(labels, continuous=False, func='PCA')
            cf.reduction_scatter(
                    self.labels,
                    continuous=False,
                    func='UMAP',
                    para=umap_para
                    )


        # +----------------------+
        # + Entropy calculations +
        # +----------------------+
        if entropy:
            plot_df = self.coef_df.copy() 
            Enz = logVariance_allocation(
                    plot_df,
                    self.labels,
                    self.prefix,
                    self.col_map,
                    dataType='Objective'
                    )

            # Overlay the entropy with the UMAP plot
            # initiate the objective for clustering and analysis
            #plot_df = plot_df[Enz.index]
            cf = clustering_func(
                        plot_df,#.div(plot_df.sum(axis=0), axis=1),
                        self.save_root_path,
                        f'{self.prefix}_Enz_{self.sel_para}',
                        mets_category(),
                    )
            
            # clustering
            cf.reduction_scatter(
                    Enz['Entropy of Objective'].to_numpy(),
                    continuous=True, func='UMAP', para=umap_para
                    )

        # +++++++++++++++++++++++++++++++++++
        # + Distances to biomass objectives +
        # +++++++++++++++++++++++++++++++++++
        if distance:
            # distance to biomass
            # get biomass coefficients
            biomass_coef_df = getRxnCoefficients(
                    model_path='/home/daweilin/StemCell/cancer_model.mat',
                    model_name='Recon1'
                    )
            # Get labels
            coef_dist = coef_distance_to_biomassObj(
                    self.coef_df, 
                    self.labels, 
                    biomass_coef_df,
                    self.prefix,
                    norm=True, 
                    rank=False, 
                    func='euclidean',
                    save_root_path=self.save_root_path,
                    histplot=True, 
                    boxplot_cols=[],
                    boxplot_order=[col for col in cols]
                    )

        # Comparisons of metabolic objectives in 2 or 3 different samples
        if len(cols)==2:
            # get labels
            col2 = cols[1] 
            col1 = cols[0]


            # +++++++++++++++++++++++++++++
            # + Allocaiton of metabolites +
            # +++++++++++++++++++++++++++++
            # comparison between 1C2C and 2CBC
            ref_col = col2
            overlap_allocation_plot(
                    self.coef_df,
                    self.labels,
                    ref_col,
                    prefix=self.prefix+f'_{col1}_{col2}',
                    norm=True,
                    cutoff=0.0#.0001
                    )
            
            # rank of coefficients
            for col in cols:
                allocation_plot(
                        self.coef_df,
                        self.labels,
                        col,
                        prefix=self.prefix,
                        norm=False,
                        cutoff=0.001
                        )

            if compare==True:
                # +----------------------+
                # + Significant Features +
                # +----------------------+
                # boxplot
                boxplot_fig(
                        self.coef_df,
                        mets_category(),
                        self.labels, col2, col1, col2,
                        f'{self.prefix}_coef_{col1}vs{col2}',
                        value_ordering=True, fc=1, portion=0.1,
                        norm=True, plottype='stripplot',
                        save_root_path=self.save_root_path
                        )
                
                
                # +------------------------+
                # + Proportion of Features +
                # +------------------------+
                # count how many cells choose a metabolite as a feature
                portion_df = lollipop_fig(
                    self.coef_df,
                    mets_category(),
                    self.labels,
                    f'{self.prefix}_coef_{col1}vs{col2}',
                    [col1, col2],
                    value_ordering=True,
                    save_root_path=self.save_root_path,
                    cutoff=0.1,
                    )

        elif len(cols)==3:

            # get labels
            col1, col2, col3 = cols[0], cols[1], cols[2]
            
            if compare==True:

                plot_df = self.coef_df.copy() 
                # +----------------------+
                # + Significant Features +
                # +----------------------+
                # boxplot
                boxplot_fig(
                    plot_df,
                    mets_category(),
                    np.array(self.labels),
                    col2, col1, col2,
                    f'{self.prefix}_coef_{col1}vs{col2}vs{col3}',
                    col3=col3, norm=True, fc=1, portion=0.1,
                    plottype='stripplot', value_ordering=True,
                    xlabel='Normalized coefficient',
                    save_root_path=self.save_root_path
                        )

                # +------------------------+
                # + Proportion of Features +
                # +------------------------+
                # count how many cells choose a metabolite as a feature
                lollipop_fig(
                    plot_df,
                    mets_category(),
                    np.array(self.labels),
                    f'{self.prefix}_coef_{col1}vs{col2}vs{col3}',
                    cols, cutoff=0.1,
                    save_root_path=self.save_root_path
                    )

        else:
            raise Exception(
                    'The length of unique column names should be either 2 or 3.'
                    )


    def tradeoff_analysis(
            self,
            input_type='coef',
            corr_th=0.5,
            tri_mets=['gh', 'chsterol', 'glutathione'],
            pareto_mets=['gh'],
            pairplot=False,
            pairScatter=False,
            triangle=False,
            pareto3D=False,
            pareto2DAnalysis=False,
            archetypeAnalysis=False,
            sample_trait_func=None,
            control_trait_func=None,
            norm=True,
            ):
        """Analysis of inferred metabolic objectives.

        The coefficients of metabolic objectives were analyzed by
        1) correlation of coefficients: a pair of metabolites that compete for allocations,  
        2) triangle plot: tradeoffs and relationships among three metabolites,
        3) pareto analysis: the realistic allocations compared to the optimal solutions in 2D or 3D,
        4) trait analysis: PCA-reduced allocations used to identify metabolic traits

        Parameters
        ----------
        input_type : {string}, default='coef'
            the data used to perform trade-off analysis. Options include 'coef' and 'fluxRecon'.
            The first type of option will get `coef_df` from the class attribute.
            The second type of option will get `wobjs` from the class attribute.

        corr_th : {float}, default=0.5
            Threshold value used to filter out pairs of metabolites with low magnitude of correlations.

        tri_mets : {string array}, default=['gh', 'chsterol', 'glutathione']
            Input of metabolites for the triangle plot. The length of the input array has to be 3.
            `glutathione`, the label, is only applicable when using Recon1 or Shen model which merges
            the coefficients of oxidized and reduced glutathione (`gthox` and `gthrd`). 
            Only useful when the parameter `triangle` is set true.

        pareto_mets : {string array}, default=['gh']
            Input of metabolites for the plots of i) scatters of relationships or ii) Pareto fronts between pairs of metabolites.
            The length of the input array has to be at least 1.

        triangle : {bool}, default=False
            enable triangle plots of allocations among the three metabolites defined in `tri_mets`.

        pareto3D : {bool}, default=False
            enable Pareto analysis among three metabolites provided in the array of `pareto_mets`
            which includes the calculation of Euclidean distance from the data points to the surface.

        traitAnalysis : {bool}, default=False
            enable PCA-based analysis of allocations in high dimensional objective coefficients.
            The vertices of traits represent the extreme traits of the metabolic objectives
            based on the mechanistic limitation.

        sample_trait_func : {bool}, default=False
            additional function that gets omics-based data points closer to traits on vertices

        control_trait_func : {bool}, default=False
            additional function that gets mechanistic data points closet to traits on vertices

        Returns
        -------
        self : object
            Returns self.
        """
        # get input data
        input_df = self.coef_df if input_type=='coef' else self.wobjs
        # get unique column names
        cols = np.unique(self.labels)
        # +-----------------------+
        # + Load and process data +
        # +-----------------------+
        # Trade-off of objective coefficients
        tradeoff_df = input_df.copy().T
        tradeoff_df['cellType'] = self.labels.to_numpy()
        # normalize coefficients in each cell
        if norm:
            tradeoff_df.iloc[:,:-1] = tradeoff_df.iloc[:, :-1].div(
                    tradeoff_df.iloc[:, :-1].sum(axis=1), axis=0
                    ).fillna(0)
        # merge rows for glutathione
        if 'glutathione' in tri_mets:
            # merge oxidized and reduced glutathione
            tradeoff_df['glutathione'] = tradeoff_df['gthox'].add(
                    tradeoff_df['gthrd']
                    ).to_numpy()
            tradeoff_df = tradeoff_df.drop(
                    columns=['gthox', 'gthrd']
                    )
            tradeoff_df['cellType'] = labels.to_numpy()
        
        
        # +--------------------------------+
        # + Relationship among metabolites +
        # +--------------------------------+
        # Trade off plots without reference points
        ObjectiveTradeOffs(
                tradeoff_df,
                f'{self.prefix}_{self.sel_para}',
                'cellType',
                sampling_objFlux=[],
                corrplot=True,
                th=corr_th,
                save_root_path=self.save_root_path,
                compare_mets=[met for met in pareto_mets],
                pairplot=False,
                single_scatter=False,
                theory_ref=False,
                cellType_curvefit=True,
                pareto_line='connected',
                hue_order=[col for col in cols]
                )

        # +------------------------------+
        # + Triangle plot of coefficient +
        # +------------------------------+
        if triangle:
            # ternary plot
            triangle_plot(
                    tradeoff_df_tmp,
                    [],
                    tri_mets[0],
                    tri_mets[1],
                    tri_mets[2],
                    self.prefix
                    )

        # +--------------------------------+
        # + 3D Pareto surface and distance +
        # +--------------------------------+
        if pareto3D:
            # point to Pareto surface distance
            # distance to the 3d surface
            dist3d = pareto_surface_to_points(
                    tradeoff_df_norm,
                    sampling_objFlux,
                    [met for met in pareto3D_mets],
                    save_root_path=self.save_root_path,
                    prefix=self.prefix,
                    )
            # pvalues
            dist_p = []
            for phase in cols:
                dist_p.append(ss.ttest_ind(
                       dist3d[dist3d['cellTypes']=='Random']['Distances'],
                        dist3d[dist3d['cellTypes']==phase]['Distances']
                        )[1])
            
            # 3d pareto surface plots
            pareto_surface_3d(
                    tradeoff_df_norm,
                    sampling_objFlux,
                    [met for met in pareto3D_mets],
                    self.prefix,
                    rot=[10, 80]
                    )
        
        # +-----------------------------------------+
        # + 2D pareto analysis and metabolic traits +
        # +-----------------------------------------+
        if pareto2DAnalysis:
            # get obj fluxes modeled by random sample of objective coefficients
            print('Loading sampling flux data...')
            sampling_objFlux = read_sampling_objFlux(
                    path=self.samplingFlux_path,
                    medium=self.medium
                    )
            sampling_objFlux.columns = pd.Series(
                    sampling_objFlux.columns
                    ).replace({'biomass_objective':'gh'})
            # analysis of 2D pareto fronts
            print('Start running Pareto analysis...')
            ObjectiveTradeOffs(
                    tradeoff_df,
                    f'{self.prefix}_{self.sel_para}',
                    'cellType',
                    sampling_objFlux=sampling_objFlux,
                    save_root_path=self.save_root_path,
                    corrplot=False,
                    pairplot=True,
                    single_scatter=True,
                    compare_mets=[met for met in pareto_mets],
                    theory_ref=True,
                    cellType_curvefit=False,
                    pareto_line='connected',
                    hue_order=[col for col in cols]
                    )


        # +---------------------------------+
        # + Archetype analysis with control +
        # +---------------------------------+
        if archetypeAnalysis_w_control:
            try:
                print(sampling_objFlux)
            except:
                # get obj fluxes modeled by random sample of objective coefficients
                print('Loading sampling flux data...')
                sampling_objFlux = read_sampling_objFlux(
                        path=self.samplingFlux_path,
                        medium=self.medium
                        )
                sampling_objFlux.columns = pd.Series(
                        sampling_objFlux.columns
                        ).replace({'biomass_objective':'gh'})
            # normalize coefficients in each cell
            tradeoff_df_norm = tradeoff_df.copy()
            if norm:
                tradeoff_df_norm.iloc[:,:-1] = tradeoff_df.iloc[:, :-1].div(
                        tradeoff_df.iloc[:, :-1].sum(axis=1), axis=0
                        ).fillna(0)
            
                   
            # +---------------------+
            # + PCA Metabolic trait +
            # +---------------------+
            # dimension reduction of randomly simulated fluxes
            archetype_df = sampling_objFlux.copy().T
            archetype_df = archetype_df.iloc[:-1, :]
            
            #archetype_df = archetype_df[archetype_df.index!='h2o']
            archetype_df = archetype_df.T[archetype_df.any(axis=0)].T
            # merge the ratio of fluxes and ratio of inferred coefficients
            # [updates: outer join for archetypes]
            sel_df = input_df.copy()
            if norm:
                sel_tmp_df = sel_df.div(sel_df.sum(axis=0), axis=1)
                arch_tmp_df = archetype_df.div(archetype_df.sum(axis=0), axis=1)
            merge_tmp = pd.concat((
                sel_tmp_df, arch_tmp_df
                ), axis=1, join='outer')
            labels_tmp = self.labels.to_list()+['Control']*arch_tmp_df.shape[1]
            # visualize with PCA plot
            cf = clustering_func(
                    merge_tmp.fillna(0),
                        self.save_root_path,
                        f'{self.prefix}_50percent',#_woH2O',
                        mets_category(),
                    )
            
            # 3 traits
            pc_df = cf.reduction_scatter(
                    labels_tmp, continuous=False, func='PCA', para=[2,2]
                    )
            
            # 4 traits
            plot_order = ['Control']+[col for col in cols]
            pc_df = cf.reduction_scatter3D(
                    labels_tmp, continuous=False, func='PCA',
                    save_animate=False, projection=False,
                    para=[3,3], plot_order=plot_order,
                    rot=[30, 120+90+30], alpha=[0.5, 0.5, 0.5, 0.5]
                    )

            if sample_trait_func!=None and control_trait_func!=None:
                # ++++++++++++++++++++++++++++
                # + Metabolic trait analysis +
                # ++++++++++++++++++++++++++++
                # trait analysis
                merge_tmp, labels_tmp, arch_df = trait_analysis(
                        sampling_objFlux,
                        input_df,
                        self.labels,
                        sample_trait_func,
                        control_trait_func,
                        n_pc=2,
                        sample_name=self.prefix,
                        plot=True
                        )
                
                # analyze the subsystems 
                subsys_report = subsystems_of_traits(
                        merge_tmp,
                        arch_df,
                        uncon_model_path=self.uncon_model_path,
                        medium=self.medium,
                        gem_path=self.GEM_path,
                        sample_name=self.prefix,
                        plot=True,
                        )
                
                
                # Overlay the coefficients of each metabolite to the pca plot
                for met in sampling_objFlux.columns[:-1]:
                    # visualize with PCA plot
                    cf = clustering_func(
                            merge_tmp,
                            self.save_root_path,
                            f'{self.prefix}_50percent_{met}',#_woH2O',
                            mets_category(),
                            )
                    coef = merge_tmp[merge_tmp.index==met].to_numpy()[0].astype(float)
                    # calculate values only based on embryo data
                    switch_arr = pd.Series(merge_tmp.columns).apply(
                            lambda x: 0 if x.split('_')[0]=='data1' else 1
                            )
                    #coef = coef*switch_arr
                    min_v = np.min(coef[coef>0])
                    met_sum = merge_tmp[merge_tmp.index==met].to_numpy()[0].astype(float)
                    #np.log10(
                    #        merge_tmp[merge_tmp.index==met].to_numpy()[0].astype(float)+min_v/10
                    #        )
                    
                    # create an alpha array
                    alpha_arr = pd.Series(merge_tmp.columns).apply(
                            lambda x: 0 if x.split('_')[0]=='data1' else 0.5
                            )
                    # 3 traits
                    pc_df = cf.reduction_scatter(
                            met_sum, continuous=True, func='PCA', para=[2,2],
                            alpha=alpha_arr
                            )

                return subsys_report
        

        # +---------------------+
        # + PCA Metabolic trait +
        # +---------------------+
        if archetypeAnalysis_RSS:
            # archetypes computing
            AA_input = tradeoff_df.iloc[:,:-1].T
            AA_input.columns = tradeoff_df['cellType']
            AA = archetypal_computing(
                            AA_input,
                            n_arch=6,
                            select_labels=[],
                            show_other=0,
                            prefix='',
                            n_arch_scan=False,
                            simplex_plot=True,
                            archetype_profile_plot=False,
                            norm=norm,
                            save_root_path=self.save_root_path
                            )

        # +---------------------+
        # + PCA Metabolic trait +
        # +---------------------+
        if archetypeAnalysis_PCA:

            AA_input = tradeoff_df.iloc[:,:-1].T
            AA_input.columns = tradeoff_df['cellType']
            # get data
            cf = clustering_func(
                    AA_input,
                    self.save_root_path,
                        f'{self.prefix}_{input_type}_norm_{norm}',
                        mets_category(),
                    )
            # 3 traits
            pc_df = cf.reduction_scatter(
                    tradeoff_df['cellType'],
                    continuous=False, func='PCA', para=[2,2]
                    )
            # 4 traits
            #plot_order = [disease]+['gtex']
            pc_df = cf.reduction_scatter3D(
                    tradeoff_df['cellType'],
                    continuous=False, func='PCA',
                    high_dimension=True,
                    save_animate=False, projection=True,
                    para=[3,3], #plot_order=plot_order,
                    rot=[30, 120+90+30], alpha=[0.5, 0.5, 0.5, 0.5]
                    )

    def fluxRecon(
            self,
            objNorm=True,
            clustering=False,
            validation=False,
            ):
        """Flux reconstruction
        
        Predict the fluxes by  coefficients*unconstrained models

        Parameters
        ----------
        objNorm : {bool}, default=True
            normalize the flux values with objective values

        clustering : {bool}, default=True
            enable dimension reduction methods or clustering methods for metabolic objectives.
            UMAP is the default method in the class, but it is allowed to be swapped by
            PCA, tSNE, and PHATE.

        validation : {bool}, default=False
            enable the evaluation of regression models with R2, Pearson and Spearman correlations.

        Returns
        -------
        self : object
            Returns self.
        res : {pandas.DataFrame}, optional
            a table of evaluation results. Only provided if `validation` is set true.

        """
        # get flux reconstruction
        wrxns, wobjs = flux_reconstruction(
                self.coef_df,
                root_path=self.uncon_model_path,
                medium=self.medium,
                )
        
        self.wrxns = wrxns
        wobjs.index = pd.Series(wobjs.index).apply(lambda x: x.split('[')[0])
        wobjs = wobjs.groupby(wobjs.index).sum()
        if objNorm:
            wobjs = wobjs.div(wobjs[wobjs.index=='Obj'].values, axis=1)
        wobjs = wobjs[wobjs.index!='Obj']
        wobjs = wobjs[wobjs.any(axis=1)]
        self.wobjs = wobjs
        
        if clustering:
            # initiate the objective for clustering and analysis
            cf = clustering_func(
                    wrxns,#.div(wobjs.iloc[-1,:], axis=1),
                    self.save_root_path,
                    f'{self.prefix}_fluxRecon',
                    mets_category(),
                    )
            
            # clustering
            cf.corr_clustermap(self.labels)
            cf.reduction_scatter(self.labels, continuous=False, func='UMAP')

        if validation:
            # report of evaluation
            cols = np.intersect1d(
                    self.coef_df.columns,
                    self.flux_df.columns
                    )
            res = evaluation_score(
                    self.coef_df[cols],
                    self.labels,
                    self.flux_df[cols],
                    plot=False,
                    uncon_path=self.uncon_model_path,
                    medium=self.medium
                    )


            return res

    def functionalAnalysis(
            self,
            cellTypes,
            sheetnames,
            expression_tb_path,
            gesa_name,
            pluri_gene_path,
            sigGenes_root_path,
            ):
        """Functional analysis overlaying with single-cell data
        
        Refer to functional genes from MSigDB for hypothesis validation.
        
        Parameters
        ----------
        cellTypes : {string array},
            names of the labels used for analysis. The array should contain unique elements.

        sheetnames : {string array},
            sheetnames to access the data of significant genes corresponding to `cellTypes`.

        expression_tb_path : {string},
            path to access the file of gene expression data.

        gesa_name : {string},
            directory to access all the functional genes.

        pluri_gene_path : {string},
            path to access the directories where save different types of functional genes.

        sigGenes_root_path : {string},
            path to access the files of significant genes.
        
        Returns
        -------
        self : object
            Returns self.

        Notes
        -----
        This method is only applicable for specific structure of directories and data in Version 0.0.1.

        Example
        -------
        >>> # umap clustering for cluster separation
        >>> cf = clustering_func(
        >>>             coef_sel.div(coef_sel.sum(axis=0), axis=1),
        >>>             './result_figs/',
        >>>             f'scEmbryo_coef_for_gsea',
        >>>             mets_category(),
        >>>         )
        >>> # get umap components
        >>> umaps = cf.reduction_scatter(
        >>>         coef_sel.columns, continuous=False, func='UMAP'
        >>>         )
        >>> # cell types 
        >>> cells = umaps['Cell type']
        >>> exptb = pd.read_excel('/home/daweilin/StemCell/Project_mESC_JinZhang/validation_data/GSE136714/GSE136714_raw.xlsx')
        >>> # cluster the umap results
        >>> recluster = cf.reclustering(
        >>>         umaps, ['UMAP1', 'UMAP2'], min_size=5
        >>>         )
        >>> # functional analysis
        >>> GESA_analyais(
        >>>         umaps, recluster, 'Cell type', exptb,
        >>>         gesa_name='pluripotency',
        >>>         cellTypes=['1C2C', '2CBC'],
        >>>         sheetnames=['2C', 'BC'], reduction_func='UMAP',
        >>>         pluri_gene_path='./mcdb/',
        >>>         sigGenes_root_path='./data/scEmbryo/GSE136714/sigGenes/',
        >>>         save_root_path='./result_figs/'
        >>>         )
        """
        # umap clustering for cluster separation
        cf = clustering_func(
                    self.coef_df.div(self.coef_df.sum(axis=0), axis=1),
                    self.save_root_path,
                    f'{self.prefix}_coef_for_gsea',
                    mets_category(),
                )
        
        # get umap components
        umaps = cf.reduction_scatter(
                coef_sel.columns, continuous=False, func='UMAP'
                )
        # cell types 
        cells = umaps['Cell type']#.map({'sc1C2C':'2cell', ''})
        exptb = pd.read_excel()
        
        
        # cluster the umap results
        recluster = cf.reclustering(
                umaps, ['UMAP1', 'UMAP2'], min_size=5
                )
        
        # functional analysis
        GESA_analyais(
                umaps, recluster, 'Cell type', exptb,
                gesa_name=gesa_name,
                cellTypes=cellTypes,
                sheetnames=sheetnames,
                reduction_func='UMAP',
                pluri_gene_path=pluri_gene_path,
                sigGenes_root_path=sigGenes_root_path,
                save_root_path=self.save_root_path
                )


