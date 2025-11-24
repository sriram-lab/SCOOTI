"""
regressionAnalyzer.py
=================================================================================
A collection of methods for analysis of metabolic objectives and flux predictions
"""
# packages
import pandas as pd
import json
import numpy as np
import os
from tqdm import tqdm
import cobra
import scipy.io

config = cobra.Configuration()
config.solver = "glpk"

class FluxModels:
    """
    A class to handle different types of metabolic models including unconstrained, constrained, multi-objective, and gene knockout models.

    Methods
    -------
    load_files_gen(filenames):
        Generate pandas DataFrames from a list of file names.

    read_gene_list(geneList_path='/home/daweilin/StemCell/unique_gene_list.mat', targetGene_path=''):
        Read a gene list from a .mat file and optionally from a target gene CSV file.

    load_KOmat_gen(files, colnames, geneList, targetGeneList=[]):
        Generate pandas DataFrames from gene knockout .mat files.

    load_models(root_path, scenario='unconstrained', norm=False, return_variables=True, medium='DMEMF12', biomass_reaction_name='biomass_objective', CFR_paraScan=False, DFA_paraScan=False, randomScan=False, stack_model=False, file_suffix='_fluxes.csv.gz', input_path_pattern='', ind_labels=False, geneList=None, targetGeneList=None, CFR_k=None, CFR_r=None, DFA_k=None, unconstrained_models=False):
        Load various types of metabolic models based on the specified scenario.

    Examples
    --------

    # Usage example
    if __name__ == "__main__":
        flux_models = FluxModels()
    
        # Example for loading unconstrained models
        unconstrained_df = flux_models.load_models(
            root_path='/path/to/unconstrained/models',
            scenario='unconstrained',
            norm=True,
            return_variables=True
        )
        print(unconstrained_df)
    
        # Example for loading constrained models
        constrained_df = flux_models.load_models(
            root_path='/path/to/constrained/models',
            scenario='constrained',
            norm=True,
            return_variables=True,
            CFR_paraScan=True,
            CFR_k=[10, 1, 0.1]
        )
        print(constrained_df)
    
        # Example for loading multi-objective models
        multiObj_df = flux_models.load_models(
            root_path='/path/to/multiObj/models',
            scenario='multiObj',
            norm=True,
            return_variables=True,
            medium='KSOM'
        )
        print(multiObj_df)
    
        # Example for loading gene knockout models
        geneKO_df = flux_models.load_models(
            root_path='/path/to/geneKO/models',
            scenario='geneKO',
            norm=True,
            return_variables=True,
            geneList=['gene1', 'gene2'],
            targetGeneList=['targetGene1', 'targetGene2']
        )
        print(geneKO_df)
    """

    def __init__(self):
        """
        Initialize the FluxModels class by setting the COBRA solver to 'glpk'.
        """
        config = cobra.Configuration()
        config.solver = "glpk"


    @staticmethod
    def load_files_gen(filenames):
        """
        Generate pandas DataFrames from a list of file names.

        Parameters
        ----------
        filenames : list of str
            List of file paths to be read.

        Yields
        ------
        DataFrame
            DataFrame generated from each file.
        """
        for filename in tqdm(filenames):
            yield pd.read_csv(filename)

    

    @staticmethod
    def read_gene_list(geneList_path='/home/daweilin/StemCell/unique_gene_list.mat', targetGene_path=''):
        """
        Read a gene list from a .mat file and optionally from a target gene CSV file.

        Parameters
        ----------
        geneList_path : str, optional
            Path to the .mat file containing the gene list (default is '/home/daweilin/StemCell/unique_gene_list.mat').
        targetGene_path : str, optional
            Path to the target gene CSV file (default is '').

        Returns
        -------
        tuple
            Tuple containing the gene columns and target gene list.
        """
        geneList = scipy.io.loadmat(geneList_path)
        key = 'geneList' if 'recon' in geneList_path else 'unique_gene_list'
        genes_columns = [g[0][0] if len(g[0]) > 0 else 'unknown' for g in geneList[key]]
        targetGeneList = pd.read_csv(targetGene_path).iloc[:, 0].to_numpy() if targetGene_path else []
        return genes_columns, targetGeneList

        

    @staticmethod
    def load_KOmat_gen(files, colnames, geneList, targetGeneList=[], process_mode='', biomass_reaction_name='biomass_objective'):
        """
        Generate pandas DataFrames from gene knockout .mat files.

        Parameters
        ----------
        files : list of str
            List of .mat file paths to be read.
        colnames : list of str
            List of column names for the DataFrames.
        geneList : list of str
            List of gene names.
        targetGeneList : list of str, optional
            List of target gene names (default is []).

        Yields
        ------
        DataFrame
            DataFrame generated from each file.
        """
        geneList = np.append('WT', geneList)
        gene_to_id = {gene: idx for idx, gene in enumerate(geneList)}


        def WTRepeat(res, gene_to_id, biomass_reaction_name):
            """
            Calculate flux differences between WT and KO predictions.
        
            Parameters:
                res (pd.DataFrame): Flux prediction DataFrame with columns named like "sample_gene".
        
            Returns:
                out_df (pd.DataFrame): KO fluxes reshaped to long format.
                geneIDs (np.ndarray): Array of gene names (excluding WT).
            """

            res = res[~res.index.str.contains(
                        f'Obj|_demand|{biomass_reaction_name}'
                        )]
            # get sample names and gene names
            samples = pd.Series(res.columns).apply(lambda x: '_'.join(x.split('_')[:-1]))
            genes = pd.Series(res.columns).apply(lambda x: x.split('_')[-1])
            # Split column names into samples and genes
            res.columns = genes
            WT = res.iloc[:,0].to_numpy() # put WT back
            res = res.mul(0, axis=0).add(WT, axis=0) # cat[WT, WT, ...]
            res.iloc[:,0] = WT # put WT back
            res_long = res.melt()
            res_long.columns = ['gene_ids', np.unique(samples)[0]]
            res_long['gene_ids'] = res_long['gene_ids'].replace(gene_to_id)
        
            return res_long

        def fluxDifference(res, gene_to_id, biomass_reaction_name):
            """
            Calculate flux differences between WT and KO predictions.
        
            Parameters:
                res (pd.DataFrame): Flux prediction DataFrame with columns named like "sample_gene".
        
            Returns:
                out_df (pd.DataFrame): KO fluxes reshaped to long format.
                geneIDs (np.ndarray): Array of gene names (excluding WT).
            """

            res = res[~res.index.str.contains(
                        f'Obj|_demand|{biomass_reaction_name}'
                        )]
            # get sample names and gene names
            samples = pd.Series(res.columns).apply(lambda x: '_'.join(x.split('_')[:-1]))
            genes = pd.Series(res.columns).apply(lambda x: x.split('_')[-1])
            # Split column names into samples and genes
            res.columns = genes
            WT = res.iloc[:,0].to_numpy() # put WT back
            res = -res.subtract(res.iloc[:,0], axis=0) # WT, WT-KO1, WT-KO2...
            res.iloc[:,0] = WT # put WT back
            res_long = res.melt()
            res_long.columns = ['gene_ids', np.unique(samples)[0]]
            res_long['gene_ids'] = res_long['gene_ids'].replace(gene_to_id)
        
            return res_long

        def mat_to_df_process(path, genes_columns, colnames, targetGeneList, process_mode='', biomass_reaction_name='biomass_objective'):
            print('Process..')
            mat = scipy.io.loadmat(path)
            rxns_index = [r[0][0] if len(r[0]) > 0 else 'unknown' for r in mat['rxn_labels']]
            geneko_flux = pd.DataFrame(mat['geneko_flux'], index=rxns_index, columns=genes_columns)
            #print(geneko_flux)
            geneko_flux.columns = [f'{colname}_{x}' for x in geneko_flux.columns]
            if len(targetGeneList)>0:
                targetGeneList = np.append(targetGeneList, 'WT')
                cols = geneko_flux.columns[
                        [x.split('_')[-1] in targetGeneList for x in geneko_flux.columns]
                        ]
                geneko_flux = geneko_flux[cols]
                
            #print(geneko_flux.sum(axis=0))
                if process_mode=='diff':
                    geneko_flux = fluxDifference(geneko_flux, gene_to_id, biomass_reaction_name='biomass_objective')
                elif process_mode=='repeat':
                    geneko_flux = WTRepeat(geneko_flux, gene_to_id, biomass_reaction_name='biomass_objective')

                return geneko_flux
            

        for file, colname in tqdm(zip(files, colnames)):
            yield mat_to_df_process(file, geneList, colname, targetGeneList, process_mode, biomass_reaction_name)

    def load_models(self,
                    root_path,
                    scenario='unconstrained',
                    norm=False,
                    return_variables=True,
                    medium='DMEMF12',
                    biomass_reaction_name='biomass_objective',
                    CFR_paraScan=False,
                    DFA_paraScan=False,
                    randomScan=False,
                    stack_model=False,
                    file_suffix='_fluxes.csv.gz',
                    input_path_pattern='',
                    ind_labels=False,
                    geneList=None,
                    targetGeneList=None,
                    CFR_k=None,
                    CFR_r=None,
                    DFA_k=None,
                    unconstrained_models=False,
                    process_mode='',
                    ):
        """
        Load various types of metabolic models based on the specified scenario.

        Parameters
        ----------
        root_path : str
            Path to the root directory containing the model files.
        scenario : str, optional
            The type of model to load: 
            'unconstrained', 'constrained', 'multiObj', 'geneKO', geneKO_unconstrained 
            (default is 'unconstrained').
        norm : bool, optional
            Whether to normalize the data (default is False).
        return_variables : bool, optional
            Whether to return flux variables or objective values (default is True).
        medium : str, optional
            Medium used for input metabolites (default is 'DMEMF12').
        biomass_reaction_name : str, optional
            Name of the biomass reaction (default is 'biomass_objective').
        CFR_paraScan : bool, optional
            Whether to scan for CFR parameters (default is False).
        DFA_paraScan : bool, optional
            Whether to scan for DFA parameters (default is False).
        randomScan : bool, optional
            Whether to read all models without considering parameters (default is False).
        stack_model : bool, optional
            Whether to stack models with constraints (default is False).
        file_suffix : str, optional
            Suffix of the filenames to read (default is '_fluxes.csv.gz').
        input_path_pattern : str, optional
            Pattern to match in the input path (default is '').
        ind_labels : bool, optional
            Whether to add index labels (default is False).
        geneList : list of str, optional
            List of gene names (default is None).
        targetGeneList : list of str, optional
            List of target gene names (default is None).
        CFR_k : list of float, optional
            List of CFR kappa values of interest (default is None).
        CFR_r : list of float, optional
            List of CFR rho values of interest (default is None).
        DFA_k : list of float, optional
            List of DFA kappa values of interest (default is None).
        unconstrained_models : bool, optional
            Whether to load unconstrained models (default is False).

        Returns
        -------
        DataFrame
            DataFrame containing the loaded model data.
        """
        if geneList is None:
            geneList = []
        if targetGeneList is None:
            targetGeneList = []
        if CFR_k is None:
            CFR_k = [10, 1, 0.1, 0.01, 0.001]
        if CFR_r is None:
            CFR_r = [10, 1, 0.1, 0.01, 0.001]
        if DFA_k is None:
            DFA_k = [10, 1, 0.1, 0.01, 0.001]

        print(f'Start processing the data of {scenario} models...')
        flux_files = [f for f in os.listdir(root_path) if os.path.isfile(os.path.join(root_path, f))]

        df_collect = []
        colnames = []

        for i, file in tqdm(enumerate(flux_files)):
            if '_metadata' in file:
                with open(os.path.join(root_path, file), 'rb') as J:
                    JSON = json.load(J)
                    switch = (randomScan or
                              (CFR_paraScan and JSON['CFR_kappa'] in CFR_k and JSON['CFR_rho'] in CFR_r) or
                              (DFA_paraScan and JSON['DFA_kappa'] in DFA_k))
                    if input_path_pattern:
                        switch = switch and input_path_pattern in JSON['input_path']
                    if switch and (scenario != 'unconstrained' or sum(np.array(JSON['obj_c'])) > 0):
                        objMet = next((obj for obj, c in zip(JSON['obj'], JSON['obj_c']) if c != 0), '')
                        cellType = JSON['input_path'].split('/')[-1].split('.xlsx')[0].split('.csv')[0]
                        if scenario == 'unconstrained':
                            medium_from_model = JSON['medium']
                            if medium_from_model != medium or not objMet:
                                continue
                        if scenario == 'geneKO' or scenario == 'geneKO_unconstrained':
                            fn = f"{file.split('_metadata')[0]}_CFR-geneDel.mat"
                            if JSON['medium'] != medium or fn not in flux_files:
                                continue
                        else:
                            fn = f"{file.split('_metadata')[0]}{file_suffix}"
                        if fn in flux_files:
                            if CFR_paraScan:
                                cellType += f"_k{JSON['CFR_kappa']}_r{JSON['CFR_rho']}"
                            if DFA_paraScan:
                                cellType += f"_dk{JSON['DFA_kappa']}"
                            if stack_model:
                                with open(JSON['CFRModel'].replace(
                                    '_model_CFR.mat', '_metadata.json'
                                    ), 'rb') as J2:
                                    JSON2 = json.load(J2)
                                    stack_model_name = JSON2[
                                            'input_path'
                                            ].split('/')[-1].split('.xlsx')[0].split('.csv')[0]
                                    cellType = f"[{stack_model_name}]{cellType}"
                            df_collect.append(os.path.join(root_path, fn))

                        if scenario=='geneKO_unconstrained':
                            objMet = np.array(JSON['obj'])[np.array(JSON['obj_c'])!=0][0]
                            colnames.append(objMet)
                        elif scenario=='unconstrained':
                            objMet = np.array(JSON['obj'])[np.array(JSON['obj_c'])!=0][0]
                            colnames.append('rxns')
                            colnames.append(objMet)
                        elif scenario=='geneKO':
                            colnames.append(cellType if cellType else file.split('_')[3])
                        else:
                            colnames.append('rxns')
                            colnames.append(cellType if cellType else file.split('_')[3])

        print('Start merging tables')
        if process_mode=='diff' or process_mode=='repeat':
            res_df = pd.concat(self.load_KOmat_gen(
                df_collect,
                colnames,
                geneList=geneList,
                targetGeneList=targetGeneList,
                process_mode=process_mode,
                biomass_reaction_name=biomass_reaction_name
                ), axis=1
                               )
            return res_df


        if scenario == 'geneKO' or scenario=='geneKO_unconstrained':
            print('error length', len(df_collect))
            res_df = pd.concat(self.load_KOmat_gen(
                df_collect,
                colnames,
                geneList=geneList,
                targetGeneList=targetGeneList
                ), axis=1
                               )
        else:
            res_df = pd.concat(self.load_files_gen(df_collect), axis=1)
            res_df.columns = colnames
        res_df = res_df.loc[:, ~res_df.columns.duplicated(keep='first')]
        if scenario != 'geneKO' and scenario != 'geneKO_unconstrained':
            res_df.index = res_df['rxns']
            res_df = res_df.iloc[:, 1:]

        if norm:
            res_df = res_df.div(res_df[~res_df.index.str.contains(f'Obj|_demand|{biomass_reaction_name}|WT')].abs().sum(axis=0), axis=1)
            res_df = res_df[~res_df.index.str.contains(f'Obj|_demand|{biomass_reaction_name}')]
        else:
            res_df = res_df[~res_df.index.str.contains(f'Obj|_demand|{biomass_reaction_name}')]

        if return_variables:
            return res_df
        else:
            dm_df = res_df[res_df.index.str.contains(f'Obj|_demand|{biomass_reaction_name}')]
            return dm_df




#def fluxDifference(res):
#    """
#    Calculate flux differences between WT and KO predictions.
#
#    Parameters:
#        res (pd.DataFrame): Flux prediction DataFrame with columns named like "sample_gene".
#
#    Returns:
#        out_df (pd.DataFrame): KO fluxes reshaped to long format.
#        geneIDs (np.ndarray): Array of gene names (excluding WT).
#    """
#    # Split column names into samples and genes
#    col_split = res.columns.str.rsplit('_', n=1, expand=True)
#    samples = col_split[0].values
#    genes = col_split[1].values
#
#    df_collect = []
#    geneIDs = []
#
#    for sample in np.unique(samples):
#        # Boolean index for columns of current sample
#        idx = (samples == sample)
#        sample_cols = res.columns[idx]
#        sample_res = res[sample_cols]
#
#        sample_genes = genes[idx]
#
#        # Separate WT and KO columns
#        wt_mask = (sample_genes == 'WT')
#        ko_mask = ~wt_mask
#
#        WT = sample_res.loc[:, sample_cols[wt_mask]]
#        KO = sample_res.loc[:, sample_cols[ko_mask]]
#
#        # Align and subtract WT from each KO column
#        diff = KO.subtract(WT.iloc[:, 0], axis=0).sum(axis=0)
#
#        # Collect KO fluxes reshaped
#        longKO = KO.copy()
#        longKO.columns = [sample] * KO.shape[1]
#        df_collect.append(longKO)
#
#        geneIDs.extend(sample_genes[ko_mask])
#
#    out_df = pd.concat(df_collect, axis=0)
#    geneIDs = np.array(geneIDs)
#
#    return out_df, geneIDs



# remove the zeros from the datatables
def remove_all_zeros(unconstrained_models, constrained_models):
    
    # rearrange datasets
    ideal_model = unconstrained_models[unconstrained_models.index.isin(constrained_models.index)]
    con1 = (ideal_model.sum(axis=1)==0)
    con2 = (constrained_models.sum(axis=1)==0)
    unconstrained_clean = ideal_model[(con1 & con2)==0].replace([np.inf, -np.inf, np.nan], [0,0,0])
    constrained_clean = constrained_models[(con1 & con2)==0].replace([np.inf, -np.inf, np.nan], [0,0,0])
    
    print(
            'Check the size of the models:',
            unconstrained_clean.shape,
            constrained_clean.shape
            )
    
    return unconstrained_clean, constrained_clean



