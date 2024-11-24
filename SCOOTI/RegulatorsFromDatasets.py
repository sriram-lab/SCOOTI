"""
RegulatorsFromDatasets.py
==============================================================
records the process to output the regulators from each dataset
"""

# import essential packages
import numpy as np
import pandas as pd
import anndata as ad
import os
import scanpy as sc
from scipy.io import mmwrite, mmread, loadmat, savemat
import gzip
from tqdm import tqdm
# essential packages for reading datasets and identifying regulators
import sys
from SCOOTI.GeneralMethods.findSigGenes import findSigGenes, posterior_probability, transition_expression
from SCOOTI.regressionAnalyzer import *
from statsmodels.stats.multitest import fdrcorrection
# datasets path
datasets_repo_path = '/nfs/turbo/umms-csriram/daweilin/data/'

def update_medium_compassRecon1():
    with open('/home/daweilin/Compass/compass/Resources/Metabolic Models/RECON1_mat/model/model.lb.json', 'r') as f:
        lb = json.load(f)
        lb = pd.Series(lb)

    with open('/home/daweilin/Compass/compass/Resources/Metabolic Models/RECON1_mat/model/model.rxns.json', 'r') as f:
        rxns = json.load(f)
        rxns = pd.Series(rxns)

    media = pd.read_excel(
            '/home/daweilin/StemCell/Project_mESC_JinZhang/SCOOTI/SCOOTI/metabolicModel/FINAL_MEDIUM_MAP_RECON1.xlsx',
            sheet_name='KSOM', index_col=0
            )
    EX_mets = media.iloc[:,1]
    EX_rxns = media.iloc[:,4]
    
    for i in range(len(EX_rxns)):
        exrxn = EX_rxns.iloc[i]
        if '_L(e)' in exrxn:
            exrxn = exrxn.replace(
                    '_L(e)', '_DASH_L_LPAREN(e)RPAREN_'
                    )
        ex_rxn_ind = rxns==exrxn
        if sum(ex_rxn_ind)>0:
            print(exrxn)
            print(lb[ex_rxn_ind])
            print(EX_mets.iloc[i])
            lb[ex_rxn_ind] = EX_mets.iloc[i]

    lb = lb.to_list()
    with open('/home/daweilin/StemCell/Project_mESC_JinZhang/model.lb.json', 'w') as f:
        json.dump(lb, f)



def scHeartFailure():
    
    # single cell embryo path
    path = '/nfs/turbo/umms-csriram/daweilin/data/GSE183852_humanHeartFailure/GSE183852_Integrated_Counts.csv.gz'
    usecols = ['gene', 'TWCM-11-103', 'TWCM-13-285']
    with pd.read_csv(path, chunksize=100) as reader:
        for chunk in reader:
            print(chunk)
            cols = pd.Series(
                    chunk.columns
                    ).apply(
                            lambda x: '_'.join(x.split('_')[:-1])
                            )
            cols[0] = 'gene'
            chunk = chunk[chunk.columns[cols.isin(usecols)]]
            chunk.head()
            break
    
    fr = findSigGenes(path)
    fr.table_to_10xformat(
        gene_cols=[],
        barcode_cols=[],
        suffix='',
        sep='\t',
        transpose=False,
        chunksize=1000
    )

    path = '/nfs/turbo/umms-csriram/daweilin/data/GSE183852_humanHeartFailure/GSE183852_Integrated_Counts/barcodes.tsv'
    df = pd.read_csv(path)


def scEmbryoTransitions_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):


    embryoSC_paths = 'scEmbryo/GSE136714/single_cell_copy/'

    # single cell embryo path
    embryoSC_paths = 'scEmbryo/GSE136714/single_cell/'
    
    # find regulators
    fr = findSigGenes(datasets_repo_path+embryoSC_paths)
    adata = fr.read_scRNAseq()
    genedf = fr.get_genedf(transpose=False)

    savedf = genedf.T.copy()
    genedf.T.to_csv('/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/expression.tsv', )

    savedf = pd.read_csv('/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/expression.tsv', sep='\t')
    savedf.to_csv('/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/expression.csv', index=None)

    # transitions for 1C2C
    ref_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[0]).str.contains('Zygote').to_numpy()

    exp_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[0]).to_numpy()=='2cell'

    up1, dw1 = fr.get_transition_genes(
            ref_cells, exp_cells, split_str='_', method='AVGSTD', save_files=False
            )
    # transitions for 2CBC
    ref_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[0]).to_numpy()=='2cell'

    exp_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[0]).to_numpy()=='32cell'

    up2, dw2 = fr.get_transition_genes(
            ref_cells, exp_cells, split_str='_', method='AVGSTD', save_files=False
            )

    # GEM
    glist = pd.read_csv('/home/daweilin/StemCell/glist.csv')
    cols = dw1.columns#[pd.Series(dw1.columns).isin(glist['Var1'])]
    # boolean gene states
    plot_dw = []
    for dw in [dw1, dw2]:
        plot_dw.append(dw[cols].copy())
    plot_up = []
    for up in [up1, up2]:
        plot_up.append(up[cols].copy())
    ups = pd.concat(plot_up, axis=0)
    dws = pd.concat(plot_dw, axis=0)

    labels = pd.Series(ups.index).apply(lambda x: x.split('_')[0])

    cf = clustering_func(
            ups.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scEmbryo_upGene',
            mets_category(),
        )

    # show correlations
    cf.corr_clustermap(labels)

    cf = clustering_func(
            dws.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scEmbryo_dwGene',
            mets_category(),
        )

    # show correlations
    cf.corr_clustermap(labels)

    # +1, 0, -1 for significantly up-regulated, neutral, and down-regulated genes
    merge_df = ups.astype(int)-dws.astype(int)
    cf = clustering_func(
            merge_df.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scEmbryo_sigGene',
            mets_category(),
        )
    # show correlations
    labels = labels.replace(['2cell', '32cell'], ['sc1C2C', 'sc2CBC'])
    cf.corr_clustermap(labels.to_numpy())
    cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[5,50])
#
def scEmbryo_sigGenes(datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'):

    # single cell embryo path
    embryoSC_paths = '/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/single_cell/'

    # find regulators
    fr = findSigGenes(embryoSC_paths)
    adata = fr.read_scRNAseq()
    genedf = fr.get_genedf(transpose=False)
    fr.get_transition_genes(split_str='~')
    

    up1, dw1 = fr.get_top_last_genes(
            split_str='_',
            ratio=0.4,
            prefix_define=f'scEmbryo',
            save_files=False,
            zscores=False,
            th=1
            )
#
def scEmbryo_transition_expression():

    
    fr = findSigGenes('/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/single_cell/')
    fr.read_scRNAseq(
            folder='',
            rename_cellNames=False
            )
    gdf = fr.get_genedf(transpose=True)

    # get columns of interest
    cols = gdf.columns[
            pd.Series(gdf.columns).str.contains(
                'Zygote|2cell|32cell'
                )
            ]
    dfsel = gdf[cols]
    # find regulators
    fr = findSigGenes('/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/')
    # dataframe
    fr.df_to_10xformat(
            gdf,
            prefix='stageEmbryoKSOM',
            )
    zygote = dfsel.columns[
            pd.Series(dfsel.columns).str.contains(
                'Zygote'
                )
            ]
    twocell = dfsel.columns[
            pd.Series(dfsel.columns).apply(
                lambda x: x.split('_')[0]=='2cell'
                )
            ]

    bc = dfsel.columns[
            pd.Series(dfsel.columns).apply(
                lambda x: x.split('_')[0]=='32cell'
                )
            ]
    # get transition of expression between two states
    sc1C2C = transition_expression(dfsel[zygote], dfsel[twocell])
    sc2CBC = transition_expression(dfsel[twocell], dfsel[bc])
    # merge df
    mergedf = pd.concat((sc1C2C, sc2CBC), axis=1)
    # find regulators
    fr = findSigGenes('/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/')
    # dataframe
    fr.df_to_10xformat(
            mergedf,
            prefix='transitionEmbryoKSOM',
            )


    # single cell embryo path
    embryoSC_paths = '/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/GSE136714_raw.xlsx'
    df = pd.read_excel(embryoSC_paths, index_col=0)
    # get columns of interest
    cols = df.columns[
            pd.Series(df.columns).str.contains(
                'Zygote|2cell|32cell'
                )
            ]
    dfsel = df[cols]
    # find regulators
    fr = findSigGenes('/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/')
    # dataframe
    fr.df_to_10xformat(
            dfsel,
            prefix='stageEmbryoKSOM',
            )
    zygote = dfsel.columns[
            pd.Series(dfsel.columns).str.contains(
                'Zygote'
                )
            ]
    twocell = dfsel.columns[
            pd.Series(dfsel.columns).apply(
                lambda x: x.split('_')[0]=='2cell'
                )
            ]

    bc = dfsel.columns[
            pd.Series(dfsel.columns).apply(
                lambda x: x.split('_')[0]=='32cell'
                )
            ]
    # get transition of expression between two states
    sc1C2C = transition_expression(dfsel[zygote], dfsel[twocell])
    sc2CBC = transition_expression(dfsel[twocell], dfsel[bc])
    # merge df
    mergedf = pd.concat((sc1C2C, sc2CBC), axis=1)
    # find regulators
    fr = findSigGenes('/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/')
    # dataframe
    fr.df_to_10xformat(
            mergedf,
            prefix='transitionEmbryoKSOM',
            )


# function to process the labels of the columns
def embryo_labels(df):
    # split the sample names
    strsplits = pd.Series(df.columns).apply(
            lambda x: x.split('_')[0]
            )
    # replace sample names
    strsplits = strsplits.apply(
            lambda x: 'Zygote' if 'Zygote' in x else x
            ).replace({
                '32cell':'Blastocyst'
                })
    return strsplits


def transition_expression(df_early, df_late):
    # get fold changes of gene expression
    expDiff = df_late.div
            df_early.mean(axis=1), axis=0
            ).fillna(0)
    ss.f_oneway(
            dfsel[zygote].iloc[0,:].to_numpy(), [0.62]
            )
    rv.sf(0.62)
    #infexpDiff = expDiff==np.inf
    #expDiff[infexpDiff] = expDiff.replace(
    #        np.inf, 0
    #        ).max().max()
    from sklearn.preprocessing import quantile_transform
    normExpDiff = quantile_transform(
            expDiff, n_quantiles=1000
            )
    normExpDiff = pd.DataFrame(normExpDiff)
    normExpDiff.index = expDiff.index
    normExpDiff.columns = expDiff.columns
    return normExpDiff


def COMPASS_flux_to_SCOOTI_format():

    # KSOM
    # transitions
    # organize fluxes
    transition = pd.read_csv('/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/transitionEmbryoKSOM_10x/reactions.tsv', sep='\t', index_col=0)
    transition[transition.index.str.contains('_neg')] = -transition[transition.index.str.contains('_neg')]
    index = pd.Series(transition.index).apply(
            lambda x: x.split('_pos')[0].split('_neg')[0]
            )
    transition.index = index
    transition = transition.groupby(transition.index).sum()
    transition.to_csv('/nfs/turbo/umms-csriram/daweilin/fluxPrediction/scEmbryo_compass/transition_KSOM.csv')

    # clustering fluxes
    cf = clustering_func(
                    transition,
                    '/nfs/turbo/umms-csriram/daweilin/result_files/embryoProject/',
                    f'scEmbryo_COMPASS_transition_fluxKSOM',
                    mets_category(),
                )
    remain_labels = embryo_labels(transition)
    umaps = cf.reduction_scatter(
            remain_labels, continuous=False, func='UMAP', para=[5,100]
            )
    # organize fluxes
    stage = pd.read_csv('/nfs/turbo/umms-csriram/daweilin/data/scEmbryo/GSE136714/stageEmbryoKSOM_10x/reactions.tsv', sep='\t', index_col=0)
    stage[stage.index.str.contains('_neg')] = -stage[stage.index.str.contains('_neg')]
    index = pd.Series(stage.index).apply(
            lambda x: x.split('_pos')[0].split('_neg')[0]
            )
    stage.index = index
    stage = stage.groupby(stage.index).sum()
    stage.to_csv('/nfs/turbo/umms-csriram/daweilin/fluxPrediction/scEmbryo_compass/scEmbryo_KSOM.csv')

    # clustering fluxes
    cf = clustering_func(
                    stage,
                    '/nfs/turbo/umms-csriram/daweilin/result_files/embryoProject/',
                    f'scEmbryo_COMPASSfluxKSOM',
                    mets_category(),
                )
    


    remain_labels = embryo_labels(stage)
    umaps = cf.reduction_scatter(
            remain_labels, continuous=False, func='UMAP', para=[5,100]
            )

def scCellCycleHela_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):

    # single cell cycle path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/cellCycleMetabolism/scMatrix/'
    
    # find regulators
    fr = findSigGenes(SC_paths)
    fr.merge_adata(rename_cellNames=True)
    genedf = fr.genedf

    # get top and last genes for each cell
    up1, dw1 = fr.get_top_last_genes(
            split_str='_', ratio=0.4, prefix_define=f'scCellCycle', save_files=False, zscores=False, th=1
            )

    # significant genes
    glist = pd.read_csv('/home/daweilin/StemCell/glist.csv')
    cols = dw1.columns[pd.Series(dw1.columns).isin(glist['Var1'])]
    
    plot_dw = []
    for dw in [dw1,]:
        plot_dw.append(dw[cols].copy())
    plot_up = []
    for up in [up1,]:
        plot_up.append(up[cols].copy())
    ups = pd.concat(plot_up, axis=0)
    dws = pd.concat(plot_dw, axis=0)

    labels = pd.Series(dws.index).apply(lambda x: x.split('_')[-2])
    

    # +1, 0, -1 for significantly up-regulated, neutral, and down-regulated genes
    merge_df = ups.astype(int)-dws.astype(int)
    cf = clustering_func(
            merge_df.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scCellCycle_sigGene',
            mets_category(),
        )
    # show correlations
    cf.corr_clustermap(labels.to_numpy())
    cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[5,500])

    # up-regulated genes 
    cf = clustering_func(
            ups.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scCellCycle_upGene',
            mets_category(),
        )

    # show correlations
    cf.corr_clustermap(labels)
    cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[50,50])

    # down-regulated genes
    cf = clustering_func(
            dws.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scCellCycle_dwGene',
            mets_category(),
        )

    # show correlations
    cf.corr_clustermap(labels)
    cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[50,50])


#####################

    # transitions for G1
    ref_cells = (pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[2]).str.contains('S')==0).to_numpy()

    exp_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[2]).str.contains('S').to_numpy()

    up1, dw1 = fr.get_transition_genes(
            ref_cells, exp_cells, split_str='_', method='CI', alpha=0.68,
            std_num=0.3, correction=False, save_files=False
            )

    # S
    ref_cells = (pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[2]).str.contains('G2')==0).to_numpy()

    exp_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[2]).str.contains('G2').to_numpy()

    up2, dw2= fr.get_transition_genes(
            ref_cells, exp_cells, split_str='_', method='CI', alpha=0.68,
            std_num=0.3, correction=False, save_files=False
            )

    # G2
    ref_cells = (pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[2]).str.contains('G1')==0).to_numpy()

    exp_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[2]).str.contains('G1').to_numpy()

    up3, dw3 = fr.get_transition_genes(
            ref_cells, exp_cells, split_str='_', method='CI', alpha=0.68,
            std_num=0.3, correction=False, save_files=False
            )

    #up4, dw4 = fr.get_compared_genes(
    #        split_str='_', method='CI', alpha=0.68,
    #        std_num=1, correction=False
    #        )
    glist = pd.read_csv('/home/daweilin/StemCell/glist.csv')
    cols = dw1.columns#[pd.Series(dw1.columns).isin(glist['Var1'])]
    
    plot_dw = []
    for dw in [dw1, dw2, dw3]:
        plot_dw.append(dw[cols].copy())
    plot_up = []
    for up in [up1, up2, up3]:
        plot_up.append(up[cols].copy())
    ups = pd.concat(plot_up, axis=0)
    dws = pd.concat(plot_dw, axis=0)
    #fr.regulator_clustering(ups, dws, 'CI', True, 0.95)

    labels = pd.Series(dws.index).apply(lambda x: x.split('_')[2])

    cf = clustering_func(
            ups.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scCellCycle_CI55_upGene',
            mets_category(),
        )

    # show correlations
    cf.corr_clustermap(labels)
    cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[50,50])

    cf = clustering_func(
            dws.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scCellCycle_CI55_dwGene',
            mets_category(),
        )

    # show correlations
    cf.corr_clustermap(labels)
    cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[50,50])


def sharma_21_cells_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/', save_to_mat=False
        ):

    # single cell cycle path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sharma_21_prolif_pmid_34274479.csv'
    # read data
    expData = pd.read_csv(SC_paths, index_col=0)
    expData.columns = pd.Series(expData.columns).apply(lambda x: int(x.split('c')[1]))
    expData = expData[sorted(expData.columns)]
    #expData = expData.div(expData.sum(axis=0), axis=1)
    # get samples
    ctrl1 = [1,5,9,15]
    ctrl2 = [7,11,17]
    exp1 = [2,6,10]
    exp2 = [13,16]
    exp3 = [8,12]
    exp4 = [14,18]
    # get sig genes
    df_dict = {
            'EGF':[exp1, ctrl1],
            'TPA':[exp2, ctrl1],
            'H89-EGF':[exp3, ctrl2],
            'H89-TPA':[exp4, ctrl2]
            }
 
    # calculate pvalues
    genes = expData.index
    for k in df_dict.keys():
        genedf = expData[df_dict[k][0]+df_dict[k][1]].T
        print(genedf.shape)
        # calculate pvalues
        pvalues = []
        foldChanges = []
        for i in range(genedf.shape[1]):
            _, p = ss.ttest_ind(
                    genedf.iloc[:len(df_dict[k][0]), i],
                    genedf.iloc[len(df_dict[k][0]):, i]
                    )
            #print(genedf.iloc[:len(df_dict[k][0]), i].shape)
            #print(genedf.iloc[len(df_dict[k][0]):, i].shape)
            p = p if np.isnan(p)==0 else 1
            fc = genedf.iloc[:len(df_dict[k][0]), i].mean()/genedf.iloc[len(df_dict[k][0]):, i].mean()
            pvalues.append(p) 
            foldChanges.append(fc)
        # FDR
        pvalues = fdrcorrection(pvalues)
        # fix issues of inf and nan
        adj_foldChanges = []
        for ele in foldChanges:
            if np.isinf(ele): 
                adj_foldChanges.append(max(np.array(foldChanges)[np.isinf(foldChanges)==0]))
            elif np.isnan(ele):
                adj_foldChanges.append(1)
            else:
                adj_foldChanges.append(ele)
        # replace the original foldchanges
        foldChanges = adj_foldChanges
        # save to .mat files
        if save_to_mat:
            genes = np.array(genes)#.reshape(len(genes), 1)
            prior = sum(np.array(pvalues)<0.05)/len(pvalues)
            ppde = [
                    posterior_probability(
                        p,prior=prior
                        ) for p in pvalues
                    ]
            # save results into a .mat file
            mdic = {
                    'GeneID':genes,
                    'PPDE':ppde,
                    'FC':foldChanges
                    }
            cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigMat/'
            isExist = os.path.exists(cell_dir)
            if not isExist:
                # create a new dir if not existing
                os.makedirs(cell_dir)
                print(f'Create a folder for {cell_dir}')
            
            # save up-/down-regulated genes for each cells
            savemat(cell_dir+'Sharma_21'+f'_{k}_.mat', mdic)


            # save results into a .mat file
            mdic_P = {
                    'GeneID':genes,
                    'PPDE':ppde,
                    'FC':1/np.array(foldChanges)
                    }
            cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigMat/'
            # save up-/down-regulated genes for each cells
            savemat(
                    cell_dir+'Sharma_21_P'+f'_{k}_.mat',
                    mdic_P
                    )


        # get up- or down-regulated genes for each cells
        dwgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)<1)
        upgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)>1)
        
        print(upgenedf.sum())
        cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigGenes/prolif_qui/'
        isExist = os.path.exists(cell_dir)
        if not isExist:
            # create a new dir if not existing
            os.makedirs(cell_dir)
            print(f'Create a folder for {cell_dir}')
        
        # save up-/down-regulated genes for each cells
        pd.DataFrame({
            'upgenes':genedf.columns[upgenedf].to_numpy()
            }).to_csv(cell_dir+'Sharma_21'+f'_{k}_upgenes.csv')

        pd.DataFrame({
            'dwgenes':genedf.columns[dwgenedf].to_numpy()
            }).to_csv(cell_dir+'Sharma_21'+f'_{k}_dwgenes.csv')

        # save regulators for reverse conditions
        pd.DataFrame({
            'dwgenes':genedf.columns[upgenedf].to_numpy()
            }).to_csv(cell_dir+'Sharma_21_P'+f'_{k}_dwgenes.csv')

        pd.DataFrame({
            'upgenes':genedf.columns[dwgenedf].to_numpy()
            }).to_csv(cell_dir+'Sharma_21_P'+f'_{k}_upgenes.csv')


def min_18_cells_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/', save_to_mat=False
        ):

    # single cell cycle path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/min_18_GSE122927_ReadCount.csv'
    # read data
    expData = pd.read_csv(SC_paths)
    expData.index = expData.gene
    # get sample names
    samples = pd.Series(expData.columns[2:]).apply(
            lambda x: '_'.join(x.split('_')[:-1])
            )
    df_dict = {}
    # p21_2E
    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*p21_high)(?=.*2E2)|(?=.*p21_low)(?=.*2E2)', regex=True)]
    cols = sorted(cols)
    genedf = expData[cols].T
    df_dict['p21_2E2'] = genedf
    # p21_3B6
    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*p21)(?=.*3B6)', regex=True)]
    cols = sorted(cols)
    genedf = expData[cols].T
    df_dict['p21_3B6'] = genedf
    # Serum Starvation 2E2
    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*SerumStarvation)(?=.*2E2)|(?=.*2E2)(?=.*Control)', regex=True)]
    genedf = expData[cols].T
    df_dict['SerumStarvation_2E2'] = genedf
    # Serum Starvation 3B6
    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*SerumStarvation)(?=.*3B6)|(?=.*3B6)(?=.*Control)', regex=True)]
    genedf = expData[cols].T
    df_dict['SerumStarvation_3B6'] = genedf
    # Meki 2E2
    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*Meki)(?=.*2E2)|(?=.*2E2)(?=.*Control)', regex=True)]
    genedf = expData[cols].T
    df_dict['Meki_2E2'] = genedf
    # Meki 3B6
    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*Meki)(?=.*3B6)|(?=.*3B6)(?=.*Control)', regex=True)]
    genedf = expData[cols].T
    df_dict['Meki_3B6'] = genedf
    # CDK46i 2E2
    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*CDK46i)(?=.*2E2)|(?=.*2E2)(?=.*Control)', regex=True)]
    genedf = expData[cols].T
    df_dict['CDK46i_2E2'] = genedf
    # CDK46i 3B6
    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*CDK46i)(?=.*3B6)|(?=.*3B6)(?=.*Control)', regex=True)]
    genedf = expData[cols].T
    df_dict['CDK46i_3B6'] = genedf
    # ContactIn 2E2
    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*ContactIn)(?=.*2E2)|(?=.*2E2)(?=.*Control)', regex=True)]
    genedf = expData[cols].T
    df_dict['ContactIn_2E2'] = genedf
    # ContactIn 3B6
    cols = expData.columns[pd.Series(expData.columns).str.contains(r'(?=.*ContactIn)(?=.*3B6)|(?=.*3B6)(?=.*Control)', regex=True)]
    genedf = expData[cols].T
    df_dict['ContactIn_3B6'] = genedf
    # get genes
    genes = expData.index
    for k in df_dict.keys():
        #k='ContactIn_2E2'
        genedf = df_dict[k]
        # calculate pvalues
        pvalues = []
        foldChanges = []
        for i in tqdm(range(genedf.shape[1])):
            _, p = ss.ttest_ind(genedf.iloc[len(genedf)//2:, i], genedf.iloc[:len(genedf)//2, i])
            p = p if np.isnan(p)==0 else 1

            pvalues.append(p) 
            fc = genedf.iloc[len(genedf)//2:, i].mean()/genedf.iloc[:len(genedf)//2, i].mean()
            foldChanges.append(fc)

        # fix issues of inf and nan
        adj_foldChanges = []
        for ele in foldChanges:
            if np.isinf(ele): 
                adj_foldChanges.append(max(np.array(foldChanges)[np.isinf(foldChanges)==0]))
            elif np.isnan(ele):
                adj_foldChanges.append(1)
            else:
                adj_foldChanges.append(ele)

        # FDR
        pvalues = fdrcorrection(pvalues)
        # replace the original foldchanges
        foldChanges = adj_foldChanges
        # save to .mat files
        if save_to_mat:
            genes = np.array(genes)#.reshape(len(genes), 1)
            prior = sum(np.array(pvalues)<0.05)/len(pvalues)
            ppde = [
                    posterior_probability(
                        p,prior=prior
                        ) for p in pvalues
                    ]
            # save results into a .mat file
            mdic = {
                    'GeneID':genes,
                    'PPDE':ppde,
                    'FC':foldChanges
                    }
            cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigMat/'
            isExist = os.path.exists(cell_dir)
            if not isExist:
                # create a new dir if not existing
                os.makedirs(cell_dir)
                print(f'Create a folder for {cell_dir}')
            
            # save up-/down-regulated genes for each cells
            savemat(
                    cell_dir+'min_18_GSE122927'+f'_{k}_.mat',
                    mdic
                    )

            # save results into a .mat file
            mdic_P = {
                    'GeneID':genes,
                    'PPDE':ppde,
                    'FC':1/np.array(foldChanges)
                    }
            cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigMat/'
            # save up-/down-regulated genes for each cells
            savemat(
                    cell_dir+'min_18_GSE122927_P'+f'_{k}_.mat',
                    mdic_P
                    )

        # get up- or down-regulated genes for each cells
        dwgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)<1)
        upgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)>1)
        
        print(upgenedf.sum())
        cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigGenes/prolif_qui/'
        isExist = os.path.exists(cell_dir)
        if not isExist:
            # create a new dir if not existing
            os.makedirs(cell_dir)
            print(f'Create a folder for {cell_dir}')
        
        # save up-/down-regulated genes for each cells
        pd.DataFrame({
            'upgenes':genedf.columns[upgenedf].to_numpy()
            }).to_csv(cell_dir+'min_18_GSE122927'+f'_{k}_upgenes.csv')

        pd.DataFrame({
            'dwgenes':genedf.columns[dwgenedf].to_numpy()
            }).to_csv(cell_dir+'min_18_GSE122927'+f'_{k}_dwgenes.csv')

        # save regulators for reverse conditions
        pd.DataFrame({
            'dwgenes':genedf.columns[upgenedf].to_numpy()
            }).to_csv(cell_dir+'min_18_GSE122927_P'+f'_{k}_dwgenes.csv')

        pd.DataFrame({
            'upgenes':genedf.columns[dwgenedf].to_numpy()
            }).to_csv(cell_dir+'min_18_GSE122927_P'+f'_{k}_upgenes.csv')



def scCellCycleHelaTransitions_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):

    # single cell cycle path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/cellCycleMetabolism/scMatrix/'
    
    # find regulators
    fr = findSigGenes(SC_paths)
    fr.merge_adata(rename_cellNames=True)
    genedf = fr.genedf

    # get top and last genes for each cell
    up1, dw1 = fr.get_top_last_genes(
            split_str='_', ratio=0.2, prefix_define=f'A549-TGFB1', save_files=True, zscores=False, th=1
            )
    # transitions for G1S
    ref_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[2]).str.contains('G1').to_numpy()

    exp_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[2]).str.contains('S').to_numpy()

    fr.get_transition_genes(ref_cells, exp_cells, split_str='_', method='AVGSTD', std_num=2)

    # transitions for SG2
    ref_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[2]).str.contains('S').to_numpy()

    exp_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[2]).str.contains('G2').to_numpy()

    fr.get_transition_genes(ref_cells, exp_cells, split_str='_', method='AVGSTD', std_num=2)


    # transitions for G2G1
    ref_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[2]).str.contains('G2').to_numpy()

    exp_cells = pd.Series(
            genedf.index
            ).apply(lambda x: x.split('_')[2]).str.contains('G1').to_numpy()

    fr.get_transition_genes(ref_cells, exp_cells, split_str='_', method='AVGSTD', std_num=2)

def QPcells_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/', save_to_mat=False
        ):

    # single cell cycle path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/johnson_18_GSE117444_prolif_qui_count.csv'


    expData = pd.read_csv(SC_paths)
    expData.index = expData.gene
    expData = expData.iloc[:, 1:7]
    

    # assign the expression table to the objective
    genedf = expData
    # get genes
    genes = expData.index#.columns
    # calculate pvalues
    pvalues = []
    foldChanges = []
    for i in tqdm(range(genedf.shape[0])):
        _, p = ss.ttest_ind(genedf.iloc[i, :3], genedf.iloc[i, 3:])
        pvalues.append(p)
        fc = genedf.iloc[i, 3:].mean()/genedf.iloc[i, :3].mean()
        foldChanges.append(fc)

    # FDR
    pvalues = fdrcorrection(pvalues)
    # fix issues of inf and nan
    adj_foldChanges = []
    for ele in foldChanges:
        if np.isinf(ele): 
            adj_foldChanges.append(max(np.array(foldChanges)[np.isinf(foldChanges)==0]))
        elif np.isnan(ele):
            adj_foldChanges.append(1)
        else:
            adj_foldChanges.append(ele)

    # replace the original foldchanges
    foldChanges = adj_foldChanges
    # save to .mat files
    if save_to_mat:
        genes = np.array(genes)#.reshape(len(genes), 1)
        prior = sum(np.array(pvalues)<0.05)/len(pvalues)
        ppde = [
                posterior_probability(
                    p,prior=prior
                    ) for p in pvalues
                ]
        # save results into a .mat file
        mdic = {
                'GeneID':genes,
                'PPDE':ppde,
                'FC':foldChanges
                }
        cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigMat/'
        isExist = os.path.exists(cell_dir)
        if not isExist:
            # create a new dir if not existing
            os.makedirs(cell_dir)
            print(f'Create a folder for {cell_dir}')
        
        # save up-/down-regulated genes for each cells
        savemat(
                cell_dir+'johnson_18_GSE117444'+f'_{k}_.mat',
                mdic
                )

        # save results into a .mat file
        mdic_P = {
                'GeneID':genes,
                'PPDE':ppde,
                'FC':1/np.array(foldChanges)
                }
        cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigMat/'
        # save up-/down-regulated genes for each cells
        savemat(
                cell_dir+'johnson_18_GSE117444_P'+f'_{k}_.mat',
                mdic_P
                )



    # get up- or down-regulated genes for each cells
    dwgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)<1)
    upgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)>1)
    cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigGenes/prolif_qui/'
    isExist = os.path.exists(cell_dir)
    if not isExist:
        # create a new dir if not existing
        os.makedirs(cell_dir)
        print(f'Create a folder for {cell_dir}')
    
    # save up-/down-regulated genes for each cells
    pd.DataFrame({
        'upgenes':genedf.columns[upgenedf].to_numpy()
        }).to_csv(cell_dir+'johnson_18_GSE117444'+'_upgenes.csv')

    pd.DataFrame({
        'dwgenes':genedf.columns[dwgenedf].to_numpy()
        }).to_csv(cell_dir+'johnson_18_GSE117444'+'_dwgenes.csv')

    # save reverse conditions
    pd.DataFrame({
        'dwgenes':genedf.columns[upgenedf].to_numpy()
        }).to_csv(cell_dir+'johnson_18_GSE117444_P'+'_dwgenes.csv')

    pd.DataFrame({
        'upgenes':genedf.columns[dwgenedf].to_numpy()
        }).to_csv(cell_dir+'johnson_18_GSE117444_P'+'_upgenes.csv')


def QPcells_to_json():
    # single cell cycle path
    path = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/DESeq_res/'
    files = os.listdir(path)
    for file in files:
        # read data
        expData = pd.read_csv(path+file)
        # save results into a .mat file
        mdic = {
                'GeneID':expData.Gene.to_numpy(),
                'PPDE':expData.PPDE.to_numpy(),
                'FC':expData.RealFC.to_numpy()
                }
        cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigMat/'
        isExist = os.path.exists(cell_dir)
        if not isExist:
            # create a new dir if not existing
            os.makedirs(cell_dir)
            print(f'Create a folder for {cell_dir}')
                
        # save up-/down-regulated genes for each cells
        k = file.split('_EBSeq')[0]
        if "_P_" in file:
            k = '_'.join(k.split('_')[:2]+['P']+k.split('_')[2:-1])
        savemat(cell_dir+f'{k}_.mat', mdic)
        print('Save file', cell_dir+f'{k}_.mat')

