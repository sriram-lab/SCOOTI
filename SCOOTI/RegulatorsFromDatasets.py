"""

  ______ _           _   _____                  _       _                 
 |  ____(_)         | | |  __ \                | |     | |                
 | |__   _ _ __   __| | | |__) |___  __ _ _   _| | __ _| |_ ___  _ __ ___ 
 |  __| | | '_ \ / _` | |  _  // _ \/ _` | | | | |/ _` | __/ _ \| '__/ __|
 | |    | | | | | (_| | | | \ \  __/ (_| | |_| | | (_| | || (_) | |  \__ \
 |_|    |_|_| |_|\__,_| |_|  \_\___|\__, |\__,_|_|\__,_|\__\___/|_|  |___/
   __                        ____    __/ |   _                            
  / _|                      / __ \  |___/   (_)                           
 | |_ _ __ ___  _ __ ___   | |  | |_ __ ___  _  ___ ___                   
 |  _| '__/ _ \| '_ ` _ \  | |  | | '_ ` _ \| |/ __/ __|                  
 | | | | | (_) | | | | | | | |__| | | | | | | | (__\__ \                  
 |_| |_|  \___/|_| |_| |_|  \____/|_| |_| |_|_|\___|___/                  
                                                                          
                                                                          

==========
The script
==========
1) realizes the instance of `findRegulators`.
2) records the processing methods of each dataset.
3) outputs the significant regulators.

"""


#%load_ext autoreload
#%autoreload 2

# import essential packages
import numpy as np
import pandas as pd
import anndata as ad
import os
import scanpy as sc
from scipy.io import mmwrite, mmread, loadmat
import gzip
# essential packages for reading datasets and identifying regulators
import sys
sys.path.append('/home/daweilin/StemCell/GeneralMethods/')
from findRegulators import findRegulators
from regressionAnalyzer import *
# datasets path
datasets_repo_path = '/nfs/turbo/umms-csriram/daweilin/data/'





def re_replace(s):
    import re
    s = re.sub('[^0-9a-zA-Z]+', '-', s)
    return s

#xl = pd.ExcelFile(nci60)
#xl.sheet_names
"""
Bulk RNAseq datasets including HumanHPA, CCLE, and NCI60 datasets
"""

def HumanHPA_sigGenes(datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'):

    hpa = '/nfs/turbo/umms-csriram/daweilin/data/HumanProteinAtlas/rna_tissue_hpa.tsv.zip'
    
    hpa_df = pd.read_csv(hpa, sep='\t')
    hpaW_df = hpa_df.pivot(index='Gene', columns='Tissue', values='nTPM')
    hpaW_df.columns = pd.Series(hpaW_df.columns).apply(lambda x: re_replace(x))
    hpaW_df = hpaW_df.T
    hpaW_df.columns = pd.Series(hpaW_df.columns).replace(dict(zip(hpa_df['Gene'], hpa_df['Gene name'])))
    
    # single cell embryo path
    bulkRNAseq_paths = 'BulkRNAseq/'
    
    # find regulators
    fr = findRegulators(datasets_repo_path+bulkRNAseq_paths)
    fr.genedf = hpaW_df
    #fr.get_top_last_genes(split_str='_', prefix='HumanHPA')
    fr.get_top_last_genes(split_str='_', prefix='HumanHPANeg', flip=True)

def ccle_sigGenes(datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'):

    ccle = '/nfs/turbo/umms-csriram/daweilin/data/CCLE/CCLE_expression.csv.gz'
    
    ccle_df = pd.read_csv(ccle, index_col=0)
    ccle_df.columns = pd.Series(ccle_df.columns).apply(lambda x: x.split(' ')[0])
    
    # single cell embryo path
    bulkRNAseq_paths = 'BulkRNAseq/'
    
    # find regulators
    fr = findRegulators(datasets_repo_path+bulkRNAseq_paths)
    fr.genedf = ccle_df
    fr.get_top_last_genes(split_str='_', prefix='ccleNeg', flip=True)
    #fr.get_top_last_genes(split_str='_', prefix='ccle')

#
def ccle_random_control_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):
    # path to the expression data
    ccle = '/nfs/turbo/umms-csriram/daweilin/data/CCLE/CCLE_expression.csv.gz'
    # get data
    ccle_df = pd.read_csv(ccle, index_col=0)
    ccle_df.columns = pd.Series(ccle_df.columns).apply(lambda x: x.split(' ')[0])
    # single cell embryo path
    bulkRNAseq_paths = 'BulkRNAseq/'
    # find regulators
    fr = findRegulators(datasets_repo_path+bulkRNAseq_paths)
    fr.genedf = ccle_df
    fr.random_sample_regulators(rep=10)
    fr.get_top_last_genes(split_str='_', prefix=f'ccleRand')


def nci60_sigGenes_flip(datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'):

    nci60 = '/nfs/turbo/umms-csriram/daweilin/data/CellMiner/output/RNA__RNA_seq_composite_expression.xls'
    
    nci60_df = pd.read_excel(
            nci60, sheet_name='Results', header=10, index_col=0
            ).iloc[:,5:]
    nci60_df.columns = pd.Series(nci60_df.columns).apply(lambda x: re_replace(x))
    # single cell embryo path
    bulkRNAseq_paths = 'BulkRNAseq/'
    # find regulators
    fr = findRegulators(datasets_repo_path+bulkRNAseq_paths)
    fr.genedf = nci60_df.T
    fr.get_top_last_genes(split_str='_', prefix='NCI60Neg', flip=True)


def nci60_sigGenes(datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'):

    nci60 = '/nfs/turbo/umms-csriram/daweilin/data/CellMiner/output/RNA__RNA_seq_composite_expression.xls'
    
    nci60_df = pd.read_excel(
            nci60, sheet_name='Results', header=10, index_col=0
            ).iloc[:,5:]
    nci60_df.columns = pd.Series(nci60_df.columns).apply(lambda x: re_replace(x))
    # single cell embryo path
    bulkRNAseq_paths = 'BulkRNAseq/'
    # find regulators
    fr = findRegulators(datasets_repo_path+bulkRNAseq_paths)
    fr.genedf = nci60_df.T
    #fr.get_top_last_genes(split_str='_', prefix='NCI60')
    up1, dw1 = fr.get_top_last_genes(
            split_str='_', ratio=0.4, prefix_define=f'NCI60Top40', save_files=True, zscores=False, th=1
            )
    # clustermap of flux correlation
    sns.set(font_scale=2)
    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    sns.clustermap(up1.T.corr().replace([np.inf, -np.inf, np.nan], [0, 0, 0]), cmap=cmap, figsize=(40,40))
    plt.savefig('/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/nci60_z1_upgene_corrmap.png')
    
    sns.clustermap(dw1.T.corr().replace([np.inf, -np.inf, np.nan], [0, 0, 0]), cmap=cmap, figsize=(40,40))
    plt.savefig('/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/nci60_z1_dwgene_corrmap.png')
#
def nci60_random_control_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):
    # path to the expression data
    nci60 = '/nfs/turbo/umms-csriram/daweilin/data/CellMiner/output/RNA__RNA_seq_composite_expression.xls'
    # get data
    nci60_df = pd.read_excel(
            nci60, sheet_name='Results', header=10, index_col=0
            ).iloc[:,5:]
    nci60_df.columns = pd.Series(nci60_df.columns).apply(lambda x: re_replace(x))
    # single cell embryo path
    bulkRNAseq_paths = 'BulkRNAseq/'
    # find regulators
    fr = findRegulators(datasets_repo_path+bulkRNAseq_paths)
    fr.genedf = nci60_df.T
    fr.random_sample_regulators(rep=100)
    fr.get_top_last_genes(split_str='_', prefix=f'NCI60Rand')


def bulk_sigGenes(datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'):
    
    # path to the data
    hpa = '/nfs/turbo/umms-csriram/daweilin/data/HumanProteinAtlas/rna_tissue_hpa.tsv.zip'
    ccle = '/nfs/turbo/umms-csriram/daweilin/data/CCLE/CCLE_expression.csv.gz'
    nci60 = '/nfs/turbo/umms-csriram/daweilin/data/CellMiner/output/RNA__RNA_seq_composite_expression.xls'
    # humanHPA datasets
    hpa_df = pd.read_csv(hpa, sep='\t')
    hpaW_df = hpa_df.pivot(index='Gene', columns='Tissue', values='nTPM')
    hpaW_df.columns = pd.Series(hpaW_df.columns).apply(lambda x: re_replace(x))
    hpaW_df.columns = pd.Series(hpaW_df.columns).apply(lambda x: f'BulkRNAseqHumanHPA_{x}')
    hpaW_df = hpaW_df.T
    hpaW_df.columns = pd.Series(hpaW_df.columns).replace(dict(zip(hpa_df['Gene'], hpa_df['Gene name'])))
    # CCLE datasets
    ccle_df = pd.read_csv(ccle, index_col=0)
    ccle_df.columns = pd.Series(ccle_df.columns).apply(lambda x: x.split(' ')[0])
    ccle_df.index = pd.Series(ccle_df.index).apply(lambda x: f'BulkRNAseqCCLE_{x}')
    nci60_df = pd.read_excel(
            nci60, sheet_name='Results', header=10, index_col=0
            ).iloc[:,5:]
    # NCI60
    nci60_df.columns = pd.Series(nci60_df.columns).apply(lambda x: re_replace(x))
    nci60_df.columns = pd.Series(nci60_df.columns).apply(lambda x: f'BulkRNAseqNCI60_{x}')
    dfs = {'HumanHPA':hpaW_df, 'CCLE':ccle_df, 'NCI60':nci60_df.T}

    #dfs = {'HumanHPA':hpaW_df, 'NCI60':nci60_df.T}
    
    # single cell embryo path
    bulkRNAseq_paths = 'BulkRNAseq/'
    
    # find regulators
    fr = findRegulators(datasets_repo_path+bulkRNAseq_paths)
    adata = fr.get_genes_from_tables(dfs, transpose=False)


    #fr.get_regulators_by_zscores(th=1, split_str='_')
    #fr.get_transition_genes(split_str='_')

    # get top and last genes for each cell
    up1, dw1 = fr.get_top_last_genes(
            split_str='_', ratio=0.2, prefix_define=f'BulkRNAseq', save_files=True, zscores=False, th=1
            )

    #ups = pd.DataFrame({k:v.to_numpy().flatten() for k, v in ups.items()}).T

    #ups, dws = fr.get_compared_genes(
    #        alpha=0.68, split_str='_', prefix_define='', method='CI', std_num=1, save_files=False
    #        )

    ## transitions from Human to CCLE
    #ref_cells = (pd.Series(
    #        fr.genedf.index
    #        ).apply(lambda x: x.split('_')[-2]).str.contains('HumanHPA')).to_numpy()

    #exp_cells = pd.Series(
    #        fr.genedf.index
    #        ).apply(lambda x: x.split('_')[-2]).str.contains('CCLE').to_numpy()

    #up1, dw1 = fr.get_transition_genes(
    #        ref_cells, exp_cells, split_str='_', method='CI', alpha=0.68, std_num=1, save_files=True
    #        )
    #
    ## transitions from HumanHPA to NCI60
    #ref_cells = (pd.Series(
    #        fr.genedf.index
    #        ).apply(lambda x: x.split('_')[-2]).str.contains('HumanHPA')).to_numpy()

    #exp_cells = pd.Series(
    #        fr.genedf.index
    #        ).apply(lambda x: x.split('_')[-2]).str.contains('NCI60').to_numpy()

    #up2, dw2 = fr.get_transition_genes(
    #        ref_cells, exp_cells, split_str='_', method='CI', alpha=0.68, std_num=1, save_files=True
    #        )
    ## transitions from cancers (NCI60 & CCLE) to HumanHPA
    #ref_cells = (pd.Series(
    #        fr.genedf.index
    #        ).apply(lambda x: x.split('_')[-2]).str.contains('HumanHPA')==0).to_numpy()

    #exp_cells = pd.Series(
    #        fr.genedf.index
    #        ).apply(lambda x: x.split('_')[-2]).str.contains('HumanHPA').to_numpy()

    #up3, dw3 = fr.get_transition_genes(
    #        ref_cells, exp_cells, split_str='_', method='CI', alpha=0.68, std_num=1, save_files=True
    #        )
    # Get metabolic genes
    glist = pd.read_csv('/home/daweilin/StemCell/glist.csv')
    cols = dw1.columns[pd.Series(dw1.columns).isin(glist['Var1'])]
    # Merge metabolic regulator dataframes
    plot_dw = []
    for dw in [dw1]:#, dw2, dw3]:
        plot_dw.append(dw[cols].copy())
    plot_up = []
    for up in [up1]:#, up2, up3]:
        plot_up.append(up[cols].copy())

    ups = pd.concat(plot_up, axis=0)
    dws = pd.concat(plot_dw, axis=0)


    # Get labels of cell types
    labels = pd.Series(dws.index).apply(lambda x: x.split('_')[-2])
    # initiate clustering
    cf = clustering_func(
            ups.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scBulkRNAseq_upGene',
            mets_category(),
        )

    # show correlations and clustering
    cf.corr_clustermap(labels)
    cf.reduction_scatter(labels, continuous=False, func='PCA')
    #cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[50,50])
    # initiate clustering of negative regulators
    cf = clustering_func(
            dws.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scBulkRNAseq_dwGene',
            mets_category(),
        )

    # show correlations and clustering
    cf.corr_clustermap(labels)
    cf.reduction_scatter(labels, continuous=False, func='PCA')
    #cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[50,50])


def scEmbryoTransitions_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):
    # single cell embryo path
    embryoSC_paths = 'scEmbryo/GSE136714/single_cell/'
    
    # find regulators
    fr = findRegulators(datasets_repo_path+embryoSC_paths)
    adata = fr.read_scRNAseq()
    genedf = fr.get_genedf(transpose=False)

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
    fr = findRegulators(embryoSC_paths)
    adata = fr.read_scRNAseq()
    genedf = fr.get_genedf(transpose=False)
    fr.get_transition_genes(split_str='~')
    



#  
def scHuman_sigGenes(datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'):
    
    # single cell human path
    humanSC_paths = 'GSE159929/'
    # initiate the class object for converting tables into 10x format
    #fr = findRegulators(datasets_repo_path+Tcell_paths)
    #fr.table_to_10xformat(gene_cols=[0, 2], barcode_cols=[2, -1])
    
    # find regulators
    fr = findRegulators(datasets_repo_path+humanSC_paths)#+'GSM3589406_PP001swap/')
    #adata = fr.read_scRNAseq()
    fr.merge_adata(rename_cellNames=True)
    fr.get_transition_genes()

# Preparing up- and down-regulated genes from scRNAseq data of EMT transitions
def scEMT_sigGenes(datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'):

    # the path to access the expression datatables
    scEMT_path = '/scEMT/'
    # specify the cell for processing
    cellType = 'DU145-TGFB1'
    # OVCA420 # MCF7 # DU145 # A549
    # Get metadata and the day of induction of TGFB1
    metaData = pd.read_csv(
            datasets_repo_path+scEMT_path+f'/{cellType}_prep_metadata.csv'
            )
    cellIDs = dict(zip(
        metaData.iloc[:, 0], metaData.apply(lambda x: f"{x.iloc[6]}_{x.name}", axis=1)
        ))

    # Get emt scores by IDs
    emtScores = pd.DataFrame({
        'cellNames':metaData.apply(lambda x: f"{x.iloc[6]}_{x.name}", axis=1),
        'emtScores':metaData['nCount_EMT_Scores']
        })
    emtScores.to_csv('/home/daweilin/StemCell/Project_mESC_JinZhang/validation/{cellType}_emtScores.csv')
    # Get the expression level of each genes from the single cell data
    expData = pd.read_csv(
            datasets_repo_path+scEMT_path+f'/{cellType}_magic_imputed.csv', index_col=0
            )
    expData.columns = pd.Series(expData.columns).map(cellIDs).apply(
            lambda x: f'{cellType}_{x}'
            )
    # initiate the object for data processing
    fr = findRegulators(datasets_repo_path+scEMT_path)
    # assign the expression table to the objective
    fr.genedf = expData.T
    
    # get top and last genes for each cell
    #up1, dw1 = fr.get_top_last_genes(
    #        split_str='_', ratio=0.4, prefix_define=f'{cellType}', save_files=True, zscores=False, th=1
    #        )
    # Method 1:
    # Identify significant genes in each single induction day
    #fr.get_compared_genes(split_str='_')
    
    # Method 2:
    # Get gene expression dataframe

    # transitions for 0 day to 8hours
    ref_cells = (pd.Series(
            fr.genedf.index
            ).apply(lambda x: x.split('_')[-2]).str.contains('7d')).to_numpy()

    exp_cells = pd.Series(
            fr.genedf.index
            ).apply(lambda x: x.split('_')[-2]).str.contains('0d').to_numpy()

    up1, dw1 = fr.get_transition_genes(
            ref_cells, exp_cells, split_str='_', method='CI', std_num=1, alpha=0.55, save_files=False
            )

    # transitions for 8hours to 1d
    ref_cells = (pd.Series(
            fr.genedf.index
            ).apply(lambda x: x.split('_')[-2]).str.contains('0d')).to_numpy()

    exp_cells = pd.Series(
            fr.genedf.index
            ).apply(lambda x: x.split('_')[-2]).str.contains('7d').to_numpy()

    up2, dw2 = fr.get_transition_genes(
            ref_cells, exp_cells, split_str='_', method='CI', std_num=1, alpha=0.55, save_files=False
            )

    
    focus_cells = (pd.Series(
            fr.genedf.index
            ).apply(lambda x: x.split('_')[-2]).str.contains('0d|7d')).to_numpy()


    glist = pd.read_csv('/home/daweilin/StemCell/glist.csv')
    cols = dw1.columns[pd.Series(dw1.columns).isin(glist['Var1'])]
    
    plot_dw = []
    for dw in [dw1,dw2]:
        plot_dw.append(dw[cols].copy())
    plot_up = []
    for up in [up1,up2]:
        plot_up.append(up[cols].copy())
    ups = pd.concat(plot_up, axis=0)
    dws = pd.concat(plot_dw, axis=0)
    #fr.regulator_clustering(ups, dws, 'CI', True, 0.95)
    #ups = ups[focus_cells]
    #dws = dws[focus_cells]
    labels = pd.Series(dws.index).apply(lambda x: x.split('_')[-2])

    cf = clustering_func(
            ups.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scEMT_CI055_{cellType}_0d7d_upGene',
            mets_category(),
        )

    # show correlations
    cf.corr_clustermap(labels)
    cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[50,50])

    cf = clustering_func(
            dws.T,
            '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
            f'scEMT_CI055_{cellType}_0d7d_dwGene',
            mets_category(),
        )

    # show correlations
    cf.corr_clustermap(labels)
    cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[50,50])

def scCellCycleHela_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):

    # single cell cycle path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/cellCycleMetabolism/scMatrix/'
    
    # find regulators
    fr = findRegulators(SC_paths)
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



# Bulk Aging data
def yu_aging_14_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):

    # data loading function generating python generators
    def load_files_gen_edit(filenames):
        for filename in tqdm(filenames):
            yield (
                    pd.read_csv(
                        filename, sep='\t'
                        )
                    )
    # single cell cycle path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/Aging/GSE53960_yu_bulk_14/'
    # read data
    filenames = [SC_paths+f for f in os.listdir(SC_paths)]
    expData = pd.concat(load_files_gen_edit(filenames), axis=1)
    expData = expData.T.drop_duplicates().T
    expData.index = expData.iloc[:, 0]
    expData = expData.iloc[:, 1:]
    expData.columns = pd.Series(expData.columns).apply(lambda x: '-'.join(x.split('_')[:-1]))
    # get unique tissue names (separate male and female)
    tissue_types = pd.Series(expData.columns).apply(lambda x: '-'.join(x.split('-')[:-1]))
    # calculate pvalues
    for k in tissue_types.unique():
        genedf = expData.T[(tissue_types==k).to_numpy()]
        genedf = genedf[genedf.index.str.contains('002|104')]
        print(genedf.shape)
        # calculate pvalues
        pvalues = []
        foldChanges = []
        for i in range(genedf.shape[1]):

            # calculate pvalues
            _, p = ss.ttest_ind(
                    genedf.loc[f'{k}_002', genedf.columns[i]],
                    genedf.loc[f'{k}_104', genedf.columns[i]]
                    )
            # make p value not significant if nan
            p = p if np.isnan(p)==0 else 1
            # save the pvalues
            pvalues.append(p) 
            # calculate fold changes
            fc = genedf.loc[f'{k}_002', genedf.columns[i]].mean()/genedf.loc[f'{k}_104', genedf.columns[i]].mean()
            foldChanges.append(fc)
        # get up- or down-regulated genes for each cells
        dwgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)<1)
        upgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)>1)
        
        print(upgenedf.sum())
        cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/Aging/GSE53960_yu_bulk_14/sigGenes/'
        isExist = os.path.exists(cell_dir)
        if not isExist:
            # create a new dir if not existing
            os.makedirs(cell_dir)
            print(f'Create a folder for {cell_dir}')
        
        # save up-/down-regulated genes for each cells
        pd.DataFrame({
            'upgenes':genedf.columns[upgenedf].to_numpy()
            }).to_csv(cell_dir+'yu_aging_14'+f'_{k}_upgenes.csv')

        pd.DataFrame({
            'dwgenes':genedf.columns[dwgenedf].to_numpy()
            }).to_csv(cell_dir+'yu_aging_14'+f'_{k}_dwgenes.csv')

        # save regulators for reverse conditions
        pd.DataFrame({
            'dwgenes':genedf.columns[upgenedf].to_numpy()
            }).to_csv(cell_dir+'yu_aging_14_P'+f'_{k}_dwgenes.csv')

        pd.DataFrame({
            'upgenes':genedf.columns[dwgenedf].to_numpy()
            }).to_csv(cell_dir+'yu_aging_14_P'+f'_{k}_upgenes.csv')

def fujimaki_19_cells_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):

    # single cell cycle path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/fujimaki_19_quiescence_deepening.txt'
    # read data
    expData = pd.read_csv(SC_paths, index_col=0, sep='\t')
    expData.columns = pd.Series(expData.columns).apply(lambda x: x.split('_')[0])
    expData = expData.div(expData.sum(axis=0), axis=1)
    # calculate pvalues
    for k in expData.columns.unique()[1:]:
        genedf = expData[['Control', k]].T
        print(genedf.shape)
        # calculate pvalues
        pvalues = []
        foldChanges = []
        for i in range(genedf.shape[1]):

            # calculate pvalues
            _, p = ss.ttest_ind(
                    genedf.loc['Control', genedf.columns[i]],
                    genedf.loc[k, genedf.columns[i]]
                    )
            # make p value not significant if nan
            p = p if np.isnan(p)==0 else 1
            # save the pvalues
            pvalues.append(p) 
            # calculate fold changes
            fc = genedf.loc[k, genedf.columns[i]].mean()/genedf.loc['Control', genedf.columns[i]].mean()
            foldChanges.append(fc)
        # get up- or down-regulated genes for each cells
        dwgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)<1)
        upgenedf = (np.array(pvalues)<0.05) & (np.array(foldChanges)>1)
        
        print(upgenedf.sum())
        cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sigGenes/deeper_qui/'
        isExist = os.path.exists(cell_dir)
        if not isExist:
            # create a new dir if not existing
            os.makedirs(cell_dir)
            print(f'Create a folder for {cell_dir}')
        
        # save up-/down-regulated genes for each cells
        pd.DataFrame({
            'upgenes':genedf.columns[upgenedf].to_numpy()
            }).to_csv(cell_dir+'Fujimaki_19'+f'_{k}_upgenes.csv')

        pd.DataFrame({
            'dwgenes':genedf.columns[dwgenedf].to_numpy()
            }).to_csv(cell_dir+'Fujimaki_19'+f'_{k}_dwgenes.csv')

        # save regulators for reverse conditions
        pd.DataFrame({
            'dwgenes':genedf.columns[upgenedf].to_numpy()
            }).to_csv(cell_dir+'Fujimaki_19_P'+f'_{k}_dwgenes.csv')

        pd.DataFrame({
            'upgenes':genedf.columns[dwgenedf].to_numpy()
            }).to_csv(cell_dir+'Fujimaki_19_P'+f'_{k}_upgenes.csv')

def sharma_21_cells_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):

    # single cell cycle path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/sharma_21_prolif_pmid_ 34274479.csv'
    # read data
    expData = pd.read_csv(SC_paths, index_col=0)
    expData.columns = pd.Series(expData.columns).apply(lambda x: int(x.split('c')[1]))
    expData = expData[sorted(expData.columns)]
    expData = expData.div(expData.sum(axis=0), axis=1)
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

            pvalues.append(p) 
            fc = genedf.iloc[:len(df_dict[k][0]), i].mean()/genedf.iloc[len(df_dict[k][0]):, i].mean()
            foldChanges.append(fc)
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
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
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
    # p21_2E2
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

    for k in df_dict.keys():
        k='ContactIn_2E2'
        genedf = df_dict[k]
        # calculate pvalues
        pvalues = []
        foldChanges = []
        for i in range(genedf.shape[1]):
            _, p = ss.ttest_ind(genedf.iloc[len(genedf)//2:, i], genedf.iloc[:len(genedf)//2, i])
            p = p if np.isnan(p)==0 else 1

            pvalues.append(p) 
            fc = genedf.iloc[len(genedf)//2:, i].mean()/genedf.iloc[:len(genedf)//2, i].mean()
            foldChanges.append(fc)
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
    fr = findRegulators(SC_paths)
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
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):

    # single cell cycle path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/dawei_prolif_qui_datasets/johnson_18_GSE117444_prolif_qui_count.csv'


    expData = pd.read_csv(SC_paths)
    expData.index = expData.gene
    expData = expData.iloc[:, 1:7]
    

    # assign the expression table to the objective
    genedf = expData.T

    # calculate pvalues
    pvalues = []
    foldChanges = []
    for i in range(genedf.shape[1]):
        _, p = ss.ttest_ind(genedf.iloc[:3, i], genedf.iloc[3:, i])
        pvalues.append(p)
        fc = genedf.iloc[3:, i].mean()/genedf.iloc[:3, i].mean()
        foldChanges.append(fc)

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

def scCellCycleHelaTransitions_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):

    # single cell cycle path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/cellCycleMetabolism/scMatrix/'
    
    # find regulators
    fr = findRegulators(SC_paths)
    fr.merge_adata(rename_cellNames=True)
    genedf = fr.genedf

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

def NegPos2C():

    data_path = '/nfs/turbo/umms-csriram/daweilin/data/ESC_2Cpos_2Cneg/GSE121451_RAW/'

    # select files
    files = [f for f in os.listdir(data_path) if '2C' in f]

    # select files 2
    # D1 experiments
    def merge_n_identify(paths, pos_inds, neg_inds, label):
        # pos
        pos = []
        for pos_ind in pos_inds:
            df_pos_tmp = pd.read_csv(
                    data_path+paths[pos_ind],sep='\t'
                    )
            df_pos_tmp.columns = [
                    f'pos_rep{i}' for i in range(len(df_pos_tmp.columns))
                    ]
            pos.append(df_pos_tmp.iloc[3:,:])
        # merge
        df_pos = pd.concat(pos, axis=0)
        df_pos.index = df_pos.iloc[:,0]
        df_pos = df_pos.iloc[:,1:]
        # neg
        neg = []
        for neg_ind in neg_inds:
            df_neg_tmp = pd.read_csv(
                    data_path+paths[neg_ind],sep='\t'
                    )
            df_neg_tmp.columns = [
                    f'neg_rep{i}' for i in range(len(df_neg_tmp.columns))
                    ]
            neg.append(df_neg_tmp.iloc[3:,:])

    # get top and last genes for each cell
    up1, dw1 = fr.get_top_last_genes(
            split_str='_', ratio=0.4, prefix_define=f'Zhang2013_ESC', save_files=True, zscores=False, th=1
            )
        # merge
        df_neg = pd.concat(neg, axis=0)
        df_neg.index = df_neg.iloc[:,0]
        df_neg = df_neg.iloc[:,1:]
        # compare two datasets
        df_compare = pd.concat((df_neg, df_pos), axis=1).astype(float)
        # pvalues
        df_pv = df_compare.apply(
                lambda x: ss.ttest_ind(
                    x.iloc[:df_neg.shape[1]],
                    x.iloc[df_neg.shape[1]:]
                    )[1], axis=1
                )
        # fold changes
        df_fc = df_compare.apply(
                lambda x:
                    x.iloc[:df_neg.shape[1]].mean()/x.iloc[df_neg.shape[1]:].mean(), axis=1
                )
        
        # get up- or down-regulated genes
        dwgenedf = (df_pv<0.05) & (df_fc<1)
        upgenedf = (df_pv<0.05) & (df_fc>1)
        
        print(upgenedf.sum())
        print(dwgenedf.sum())

        return dwgenedf, upgenedf
        #cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/ESC_2Cpos_2Cneg/sigGenes/'
        #isExist = os.path.exists(cell_dir)
        #if not isExist:
        #    # create a new dir if not existing
        #    os.makedirs(cell_dir)
        #    print(f'Create a folder for {cell_dir}')
        #
        ## save up-/down-regulated genes for each cells
        #pd.DataFrame({
        #    'upgenes':upgenedf.index[upgenedf]
        #    }).to_csv(cell_dir+'Neg2C'+f'_{label}_upgenes.csv')

        #pd.DataFrame({
        #    'dwgenes':dwgenedf.index[dwgenedf]
        #    }).to_csv(cell_dir+'Neg2C'+f'_{label}_dwgenes.csv')

        ## save regulators for reverse conditions
        #pd.DataFrame({
        #    'upgenes':upgenedf.index[upgenedf]
        #    }).to_csv(cell_dir+'Pos2C'+f'_{label}_upgenes.csv')

        #pd.DataFrame({
        #    'dwgenes':dwgenedf.index[dwgenedf]
        #    }).to_csv(cell_dir+'Pos2C'+f'_{label}_dwgenes.csv')

    # D1Rep1
    pathsD1Rep1 = sorted([
            f for f in files if (('Rep1' in f) & ('D1' in f))
                ])
    neg_inds, pos_inds, label = [0,1], [2,3], 'D1Rep1'
    dw, up = merge_n_identify(pathsD1Rep1, pos_inds, neg_inds, label)
    
    glist = pd.read_csv('/home/daweilin/StemCell/glist.csv')
    up.index = pd.Series(up.index).apply(lambda x: x.upper())
    cols = dw[dw.isin(glist['Var1'])]
    # D1Rep2
    pathsD1Rep2 = sorted([
            f for f in files if (('Rep2' in f) & ('D1' in f))
                ])
    neg_inds, pos_inds, label = [0,1], [2,3], 'D1Rep2'
    merge_n_identify(pathsD1Rep2, pos_inds, neg_inds, label)
        

def NegPos2C():

    data_path = '/nfs/turbo/umms-csriram/daweilin/data/ESC_2Cpos_2Cneg/GSE121451_RAW/'

    # select files
    files = [f for f in os.listdir(data_path) if '2C' in f]

    # select files 2
    # D1 experiments
    def merge_n_identify(paths, pos_inds, neg_inds, label):
        # pos
        pos = []
        for pos_ind in pos_inds:
            df_pos_tmp = pd.read_csv(
                    data_path+paths[pos_ind],sep='\t'
                    )
            df_pos_tmp.columns = [
                    f'pos_rep{i}' for i in range(
                        len(df_pos_tmp.columns)
                        )
                    ]
            pos.append(df_pos_tmp.iloc[3:,:])
        # merge
        df_pos = pd.concat(pos, axis=0)
        df_pos.index = df_pos.iloc[:,0]
        df_pos = df_pos.iloc[:,1:]
        # neg
        neg = []
        for neg_ind in neg_inds:
            df_neg_tmp = pd.read_csv(
                    data_path+paths[neg_ind],sep='\t'
                    )
            df_neg_tmp.columns = [
                    f'neg_rep{i}' for i in range(
                        len(df_neg_tmp.columns)
                        )
                    ]
            neg.append(df_neg_tmp.iloc[3:,:])

        # get top and last genes for each cell
        # merge
        df_neg = pd.concat(neg, axis=0)
        df_neg.index = df_neg.iloc[:,0]
        df_neg = df_neg.iloc[:,1:]
        
        # merge two datasets
        df_compare = pd.concat(
                (df_neg, df_pos), axis=1
                ).astype(float)
        
        # find regulators
        fr = findRegulators(data_path)
        # assign the merged dataframe to the instance
        fr.genedf = df_compare
        # identify sig genes
        up1, dw1 = fr.get_top_last_genes(
                split_str='_',
                ratio=0.4,
                prefix_define=f'PosNeg2C{label}',
                save_files=False,
                zscores=False, th=1
                )
        return up1, dw1



        #cell_dir = '/nfs/turbo/umms-csriram/daweilin/data/ESC_2Cpos_2Cneg/sigGenes/'
        #isExist = os.path.exists(cell_dir)
        #if not isExist:
        #    # create a new dir if not existing
        #    os.makedirs(cell_dir)
        #    print(f'Create a folder for {cell_dir}')
        #
        ## save up-/down-regulated genes for each cells
        #pd.DataFrame({
        #    'upgenes':upgenedf.index[upgenedf]
        #    }).to_csv(cell_dir+'Neg2C'+f'_{label}_upgenes.csv')

        #pd.DataFrame({
        #    'dwgenes':dwgenedf.index[dwgenedf]
        #    }).to_csv(cell_dir+'Neg2C'+f'_{label}_dwgenes.csv')

        ## save regulators for reverse conditions
        #pd.DataFrame({
        #    'upgenes':upgenedf.index[upgenedf]
        #    }).to_csv(cell_dir+'Pos2C'+f'_{label}_upgenes.csv')

        #pd.DataFrame({
        #    'dwgenes':dwgenedf.index[dwgenedf]
        #    }).to_csv(cell_dir+'Pos2C'+f'_{label}_dwgenes.csv')

    # D1Rep1
    pathsD1Rep1 = sorted([
            f for f in files if (('Rep1' in f) & ('D1' in f))
                ])
    neg_inds, pos_inds, label = [0,1], [2,3], 'D1Rep1'
    up, dw = merge_n_identify(pathsD1Rep1, pos_inds, neg_inds, label)
    
    glist = pd.read_csv('/home/daweilin/StemCell/glist.csv')
    up.index = pd.Series(up.index).apply(lambda x: x.upper())
    dw.index = pd.Series(dw.index).apply(lambda x: x.upper())
    cols = up[up.index.isin(glist['Var1'])]
    # D1Rep2
    pathsD1Rep2 = sorted([
            f for f in files if (('Rep2' in f) & ('D1' in f))
                ])
    neg_inds, pos_inds, label = [0,1], [2,3], 'D1Rep2'
    merge_n_identify(pathsD1Rep2, pos_inds, neg_inds, label)

def Zhang2013ESC_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):

    # tcell path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/Tzelepis2016/GSE44067_esc_nsc_npc_expression.txt.gz'
    tmp = pd.read_csv(SC_paths, sep='\t', index_col=0)
    tmp = pd.DataFrame(tmp['ESC'])
    
    # find regulators
    fr = findRegulators(SC_paths)
    fr.genedf = tmp.T

    # get top and last genes for each cell
    up1, dw1 = fr.get_top_last_genes(
            split_str='_', ratio=0.4, prefix_define=f'Zhang2013_ESC', save_files=True, zscores=False, th=1
            )


def naivePluriptency_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):

    # tcell path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/GSE107060_naivePluripotency/GSE107060_processed_data_tpm.txt.gz'
    tmp = pd.read_csv(SC_paths, sep='\t', index_col=0)
    tmp = pd.DataFrame(tmp['ESC'])
    
    # find regulators
    fr = findRegulators(SC_paths)
    fr.genedf = tmp.T

    # get top and last genes for each cell
    up1, dw1 = fr.get_top_last_genes(
            split_str='_', ratio=0.4, prefix_define=f'Zhang2013_ESC', save_files=True, zscores=False, th=1
            )


def scTCell_sigGenes(
        datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/'
        ):

    # tcell path
    SC_paths = '/nfs/turbo/umms-csriram/daweilin/data/Tcell/GSE126030/lymphNode-1/'
    
    # find regulators
    fr = findRegulators(SC_paths)
    fr.merge_adata(rename_cellNames=True)
    genedf = fr.genedf
    cellType = {
            'PP006swap':'activated-1-lymphNode',
            'PP005swap':'resting-1-lymphNode',
            } 
    inds = pd.Series(genedf.index).apply(lambda x: x.split('_')[1]).replace(cellType)

    # significant genes
    #up1, dw1 = fr.get_compared_genes(
    #        alpha=0.68, split_str='_', prefix_define='', method='CI', std_num=2, save_files=False
    #        )


    # transitions from resting to activated cells
    #ref_cells = (pd.Series(
    #        inds
    #        ).str.contains('resting')).to_numpy()

    #exp_cells = pd.Series(
    #        inds
    #        ).str.contains('activated').to_numpy()

    #up1, dw1 = fr.get_transition_genes(
    #        ref_cells, exp_cells, split_str='_', method='CI', alpha=0.68, std_num=1, save_files=False
    #        )
    ## transitions from activated to resting
    #ref_cells = (pd.Series(
    #        inds
    #        ).str.contains('activated')).to_numpy()

    #exp_cells = pd.Series(
    #        inds
    #        ).str.contains('resting').to_numpy()

    #up2, dw2 = fr.get_transition_genes(
    #        ref_cells, exp_cells, split_str='_', method='CI', alpha=0.68, std_num=1, save_files=True
    #        )
    #
    #glist = pd.read_csv('/home/daweilin/StemCell/glist.csv')
    #cols = dw1.columns[pd.Series(dw1.columns).isin(glist['Var1'])]
    #
    #plot_dw = []
    #for dw in [dw1, dw2]:
    #    plot_dw.append(dw[cols].copy())
    #plot_up = []
    #for up in [up1, up2]:
    #    plot_up.append(up[cols].copy())

    #ups = pd.concat(plot_up, axis=0)
    #dws = pd.concat(plot_dw, axis=0)

    #labels = pd.Series(dws.index).apply(lambda x: x.split('_')[-2])

    #cf = clustering_func(
    #        ups.T,
    #        '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
    #        f'scTcell_upGene',
    #        mets_category(),
    #    )

    ## show correlations
    #cf.corr_clustermap(labels)
    #cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[50,50])

    #cf = clustering_func(
    #        dws.T,
    #        '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
    #        f'scTcell_dwGene',
    #        mets_category(),
    #    )

    ## show correlations
    #cf.corr_clustermap(labels)
    #cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[50,50])


print('Initiating...')
scEMT_sigGenes(datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/')
#bulk_sigGenes(datasets_repo_path='/nfs/turbo/umms-csriram/daweilin/data/')
print('Done...')
# datasets path
#datasets_repo_path = '/nfs/turbo/umms-csriram/daweilin/data/'



