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
from SCOOTI.GeneralMethods.findRegulators import findRegulators, posterior_probability, transition_expression
from SCOOTI.regressionAnalyzer import *
from statsmodels.stats.multitest import fdrcorrection
# datasets path
datasets_repo_path = '/nfs/turbo/umms-csriram/daweilin/data/'

def column_slice(chunk):

    usecols = ['gene', 'TWCM-11-103', 'TWCM-13-285']
    cols = pd.Series(
            chunk.columns
            ).apply(
                    lambda x: '_'.join(x.split('_')[:-1])
                    )
    cols[0] = 'gene'
    chunk = chunk[chunk.columns[cols.isin(usecols)]]

    return chunk

def scHeartFailure():

    # single cell embryo path
    path = '/nfs/turbo/umms-csriram/daweilin/data/GSE183852_humanHeartFailure_subset/'
    
    fr = findRegulators(path)
    fr.table_to_10xformat(
        gene_cols=[0, 1],
        barcode_cols=[],
        suffix='',
        sep=',',
        transpose=False,
        chunksize=100,
        column_slice=1,
        column_slice_func=column_slice
    )

def get_scHeartFailure_sigGenes():
    # root path to access the data
    data_path = '/nfs/turbo/umms-csriram/daweilin/data/GSE183852_humanHeartFailure_subset/GSE183852_Integrated_Counts/'
    
    #genes = pd.read_csv(data_path+'genes.tsv', '\t')
    #genes = pd.concat((genes, genes, genes), axis=1)
    #genes.to_csv(data_path+'genes.tsv', sep='\t', index=0)
    #genes = pd.read_csv(data_path+'genes.tsv', '\t')
    
    fr = findRegulators(data_path)
    adata = fr.read_scRNAseq()
    genedf = fr.get_genedf(transpose=False)
    # load data
    up1, dw1 = fr.get_top_last_genes(
            split_str='_',
            ratio=0.4,
            prefix_define=f'scHeartFailure',
            save_files=True,
            zscores=False,
            th=1
            )
    #up1 = up1.T
    #dw1 = dw1.T
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

    labels = pd.Series(dws.index).apply(lambda x: '_'.join(x.split('_')[:-1]))
    
    # +1, 0, -1 for significantly up-regulated, neutral, and down-regulated genes
    merge_df = ups.astype(int)-dws.astype(int)
    cf = clustering_func(
            merge_df.T,
            '/nfs/turbo/umms-csriram/daweilin/result_files/diseaseProject/',
            f'scHeartFailure_sigGene',
            mets_category(),
        )
    # show correlations
    cf.corr_clustermap(labels.to_numpy())
    cf.reduction_scatter(labels, continuous=False, func='UMAP', para=[50,50])

    return up1, dw1

if __name__ == "__main__":
    #scHeartFailure()
    get_scHeartFailure_sigGenes()



