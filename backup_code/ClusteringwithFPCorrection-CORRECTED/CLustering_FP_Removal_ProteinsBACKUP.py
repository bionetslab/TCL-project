import pickle
import numpy as np
import pandas as pd
from PIL import Image
import glob
import matplotlib.pyplot as plt
from skimage.morphology import convex_hull_image
from skimage import data, img_as_float
from skimage.util import invert
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from multiprocessing import Pool
import time
import math
from collections import Counter
import scanpy as sc
import networkx as nx
import pandas as pd
import numpy as np
import itertools
import random
from scipy.stats import mannwhitneyu
import os
import warnings
import pickle
import matplotlib.pyplot as plt
import pandas as pd
import gseapy as gp
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import mygene
import seaborn as sns
from gseapy import barplot, dotplot
import random
import matplotlib.pyplot as plt
import numpy as np
import decoupler as dc
from numpy.random import default_rng
from copy import deepcopy
import scanpy as sc
import squidpy as sq
import time
from neighborhood_enrichment import neighborhood_enrichment
import schist as scs
patient_id=20751

def categorize(dataframe, leiden_clusters_sorted):
    df = pd.DataFrame(columns=['group'])
    for i in leiden_clusters_sorted:
        df_=dataframe[dataframe['group']==i]
        df_=df_.sort_values(by = 'pvals_adj', ascending=True)
        df=pd.concat([df, df_], ignore_index=True)
    return df

def calculate_average_expressions_per_protein(adata_X, cells_in_cluster, proteins_in_cluster):
    n=len(cells_in_cluster)
    DF=adata_X.loc[cells_in_cluster]
    DF=DF[proteins_in_cluster]
    per_row_sum_all_proteins=DF.sum(axis=1).tolist()
    per_row_indiv_proteins=[]
    per_row_sum_all_proteins_excluding_indiv=[]
    per_row_avg_all_proteins_excluding_indiv=[]
    list_of_proteins=DF.columns.tolist()
    for per_protein in list_of_proteins:
        indiv_protein=DF[per_protein].tolist()
        per_row_indiv_proteins.append(indiv_protein)
        sub_=list(np.subtract(per_row_sum_all_proteins, indiv_protein))
        per_row_sum_all_proteins_excluding_indiv.append(sub_)
        avg_=[]
        for val in sub_:
            avg_.append(val/float(n-1))
        per_row_avg_all_proteins_excluding_indiv.append(avg_)
    return per_row_avg_all_proteins_excluding_indiv, per_row_indiv_proteins

with open('TNBC_41patients_KerenEtAl.pkl', 'rb') as f:
   pickle_= pickle.load(f)

recurrence_labels=[]
survival_days_counts=[]
survival_days_binary=[]

nhood_enrichment_zscores_recurrence_1=[]
nhood_enrichment_zscores_recurrence_0=[]

nhood_enrichment_zscores_survival_1=[]
nhood_enrichment_zscores_survival_0=[]

cnt=-1


for i in pickle_:
    recurrence_labels.append(np.unique(pickle_[i].obsm['recurrence'])[0])
    survival_days_counts.append(np.unique(pickle_[i].obsm['survival_days_capped'])[0])

survival_days_mean=np.mean(survival_days_counts)
for i in survival_days_counts:
    if i>survival_days_mean:
        survival_days_binary.append(1)
    else:
        survival_days_binary.append(0)


# protein_channels=list(pickle_[list(pickle_.keys())[0]].to_df())

# ------
pickle_items=list(pickle_.items())
first_three_items = pickle_items[7:8]
pickle_new=dict(first_three_items)
pickle_=pickle_new
# ------
   
no_of_patients=np.shape(pickle_)
error_samples=[]

for i in pickle_:
    
    print('Patient serial number' + str(cnt) + ' (patient id: ' + str(i) + ')')
    print('====================================================================')
    
    cnt+=1
    
    adata=pickle_[i]
    adata_X=adata.to_df()

    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    adata_raw=adata.copy()

    sc.tl.pca(adata_raw, svd_solver='arpack')
    
    sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=3)
    sc.tl.leiden(adata_raw)
    leiden_clusters_sorted=[str(j) for j in sorted([int(i) for i in np.unique(adata_raw.obs['leiden'])])]
    # scs.inference.nested_model(adata_raw)
    scs.inference.planted_model(adata_raw)
    ppbm_clusters_sorted=[str(j) for j in sorted([int(i) for i in np.unique(adata_raw.obs['ppbm'])])]
    
    sc.tl.umap(adata_raw)
    sc.pl.umap(adata_raw, color='leiden')
    sc.pl.umap(adata_raw, color=['ppbm'], legend_loc='on data')
    # ====================================================================================
    
    # sc.pp.scale(adata_raw)
    
    sc.tl.rank_genes_groups(adata_raw, groupby='ppbm', pts=True, method='wilcoxon')
    DEDF_ppbm=sc.get.rank_genes_groups_df(adata_raw, group=ppbm_clusters_sorted)
    DEDF_ppbm_AscPVals=categorize(DEDF_ppbm, ppbm_clusters_sorted)
    adataRaw_rankGenes_names, adataRaw_rankGenes_scores, adataRaw_rankGenes_logfoldchanges, adataRaw_rankGenes_pvals, adataRaw_rankGenes_pvals_adj, adataRaw_rankGenes_pts, adataRaw_rankGenes_pts_rest=pd.DataFrame(adata_raw.uns['rank_genes_groups']['names']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['scores']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['logfoldchanges']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['pvals']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['pvals_adj']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['pts']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['pts_rest'])
    sc.pl.rank_genes_groups(adata_raw, sharey=False, save='sc.pl.rank_genes_groups_PPBMClustering.jpg')
    
    DEDF_ppbm_AscPVals_0Point05_cutoff=DEDF_ppbm_AscPVals[DEDF_ppbm_AscPVals.pvals_adj<0.001]
    DEDF_ppbm_AscPVals_0Point05_cutoff_groups=list(np.unique(DEDF_ppbm_AscPVals_0Point05_cutoff.group))
    DEDF_ppbm_AscPVals_0Point05_cutoff_groups_sorted=[str(j) for j in sorted([int(i) for i in np.unique(DEDF_ppbm_AscPVals_0Point05_cutoff_groups)])]
    
    ppbm_no_of_clusters_dict_values=[[] for _ in range(len(DEDF_ppbm_AscPVals_0Point05_cutoff_groups_sorted))]
    ppbm_no_of_clusters_dict_keys=DEDF_ppbm_AscPVals_0Point05_cutoff_groups_sorted
    ppbm_no_of_clusters_dict=dict(zip(ppbm_no_of_clusters_dict_keys, ppbm_no_of_clusters_dict_values))
    for index, row in DEDF_ppbm_AscPVals_0Point05_cutoff.iterrows():
        ppbm_no_of_clusters_dict[row['group']].append(row['names'])
    
    df=adata_raw.obs['ppbm'].copy().to_frame()
    # ppbm_cell_lookup_keys=[str(j) for j in sorted([int(i) for i in np.unique(df.index)])]
    ppbm_cell_lookup_values=[[] for _ in range(len(ppbm_clusters_sorted))]
    # ppbm_cell_lookup_keys=[[] for _ in range(len(np.unique(ppbm_cell_lookup_keys)))]
    ppbm_cell_lookup_keys=ppbm_clusters_sorted
    ppbm_cell_lookup_dict=dict(zip(ppbm_cell_lookup_keys, ppbm_cell_lookup_values))
    for index, row in df.iterrows():
        ppbm_cell_lookup_dict[row['ppbm']].append(index)
    #////////////
    pvals_across_clusters=[]
    for key in ppbm_no_of_clusters_dict:
        print('PPBM cluster ' + str(key) + ' in patient serial number ' + str(cnt) + ' (patient id: ' + str(i) + ')')
        proteins_in_cluster=ppbm_no_of_clusters_dict[key]
        cells_in_cluster=ppbm_cell_lookup_dict[key]
        per_row_avg_all_proteins_excluding_indiv, per_row_indiv_proteins=calculate_average_expressions_per_protein(adata_X, cells_in_cluster, proteins_in_cluster)
        # all_proteins_excluding_indiv_DICT=dict(zip(ppbm_no_of_clusters_dict_keys, per_row_avg_all_proteins_excluding_indiv))
        cluster_count=-1
        pvals=[]
        for i_count in range(len(proteins_in_cluster)):
            cluster_count+=1
            U1, p = mannwhitneyu(per_row_avg_all_proteins_excluding_indiv[i_count], per_row_indiv_proteins[i_count], method="exact")
            pvals.append(p)
        proteinsAndPvalues_in_cluster_dict=dict(zip(proteins_in_cluster, pvals))
        pvals_across_clusters.append(proteinsAndPvalues_in_cluster_dict)
    
    PVALUES_dict=dict(zip(ppbm_no_of_clusters_dict_keys, pvals_across_clusters))
    
    important_proteins_ppbm=[]
    for key in PVALUES_dict:
        for key_inner in PVALUES_dict[key]:
            important_proteins_ppbm.append(key_inner)
    important_proteins_ppbm=np.unique(important_proteins_ppbm)
    
    sc.pl.dotplot(adata_raw, important_proteins_ppbm, groupby='ppbm', save='PPBM_planted_clustering_results_dotplot.jpg');
    sc.pl.stacked_violin(adata_raw, important_proteins_ppbm, groupby='ppbm', rotation=90, save='PPBM_planted_clustering_results_violinplot.jpg');
    
    
    # ========================================================================================================================
    
    
    sc.tl.rank_genes_groups(adata_raw, groupby='leiden', pts=True, method='wilcoxon')
    DEDF_leiden=sc.get.rank_genes_groups_df(adata_raw, group=leiden_clusters_sorted)
    DEDF_leiden_AscPVals=categorize(DEDF_leiden, leiden_clusters_sorted)
    adataRaw_rankGenes_names, adataRaw_rankGenes_scores, adataRaw_rankGenes_logfoldchanges, adataRaw_rankGenes_pvals, adataRaw_rankGenes_pvals_adj, adataRaw_rankGenes_pts, adataRaw_rankGenes_pts_rest=pd.DataFrame(adata_raw.uns['rank_genes_groups']['names']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['scores']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['logfoldchanges']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['pvals']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['pvals_adj']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['pts']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['pts_rest'])
    sc.pl.rank_genes_groups(adata_raw, sharey=False, save='sc.pl.rank_genes_groups_LEIDENClustering.jpg')
    
    DEDF_leiden_AscPVals_0Point05_cutoff=DEDF_leiden_AscPVals[DEDF_leiden_AscPVals.pvals_adj<0.001]
    DEDF_leiden_AscPVals_0Point05_cutoff_groups=list(np.unique(DEDF_leiden_AscPVals_0Point05_cutoff.group))
    DEDF_leiden_AscPVals_0Point05_cutoff_groups_sorted=[str(j) for j in sorted([int(i) for i in np.unique(DEDF_leiden_AscPVals_0Point05_cutoff_groups)])]
    
    leiden_no_of_clusters_dict_values=[[] for _ in range(len(DEDF_leiden_AscPVals_0Point05_cutoff_groups_sorted))]
    leiden_no_of_clusters_dict_keys=DEDF_leiden_AscPVals_0Point05_cutoff_groups_sorted
    leiden_no_of_clusters_dict=dict(zip(leiden_no_of_clusters_dict_keys, leiden_no_of_clusters_dict_values))
    for index, row in DEDF_leiden_AscPVals_0Point05_cutoff.iterrows():
        leiden_no_of_clusters_dict[row['group']].append(row['names'])
    
    df=adata_raw.obs['leiden'].copy().to_frame()
    leiden_cell_lookup_values=[[] for _ in range(len(leiden_clusters_sorted))]
    # ppbm_cell_lookup_keys=[[] for _ in range(len(np.unique(ppbm_cell_lookup_keys)))]
    leiden_cell_lookup_keys=leiden_clusters_sorted
    leiden_cell_lookup_dict=dict(zip(leiden_cell_lookup_keys, leiden_cell_lookup_values))
    for index, row in df.iterrows():
        leiden_cell_lookup_dict[row['leiden']].append(index)
    
    pvals_across_clusters=[]
    for key in leiden_no_of_clusters_dict:
        print('Leiden cluster ' + str(key) + ' in patient serial number ' + str(cnt) + ' (patient id: ' + str(i) + ')')
        proteins_in_cluster=leiden_no_of_clusters_dict[key]
        cells_in_cluster=leiden_cell_lookup_dict[key]
        per_row_avg_all_proteins_excluding_indiv, per_row_indiv_proteins=calculate_average_expressions_per_protein(adata_X, cells_in_cluster, proteins_in_cluster)
        cluster_count=-1
        pvals=[]
        for i_count in range(len(proteins_in_cluster)):
            cluster_count+=1
            U1, p = mannwhitneyu(per_row_avg_all_proteins_excluding_indiv[i_count], per_row_indiv_proteins[i_count], method="exact")
            pvals.append(p)
        proteinsAndPvalues_in_cluster_dict=dict(zip(proteins_in_cluster, pvals))
        pvals_across_clusters.append(proteinsAndPvalues_in_cluster_dict)
    
    PVALUES_dict=dict(zip(leiden_no_of_clusters_dict_keys, pvals_across_clusters))
    
    important_proteins_leiden=[]
    for key in PVALUES_dict:
        for key_inner in PVALUES_dict[key]:
            important_proteins_leiden.append(key_inner)
    important_proteins_leiden=np.unique(important_proteins_leiden)
    
    sc.pl.dotplot(adata_raw, important_proteins_leiden, groupby='leiden', save='LEIDEN_clustering_results_dotplot.jpg');
    sc.pl.stacked_violin(adata_raw, important_proteins_leiden, groupby='leiden', rotation=90, save='LEIDEN_clustering_results_violinplot.jpg');
    
    
    # ========================================================================================================================
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    