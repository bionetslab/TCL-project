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

def calculate_average_expressions_per_cell(adata_X, _count_):
    n=len(_count_)
    DF=adata_X.loc[_count_]
    per_column_sum_all_cells=DF.sum(axis=0).tolist()
    per_column_indiv_cells=[]
    per_column_sum_all_cells_excluding_indiv=[]
    per_column_avg_all_cells_excluding_indiv=[]
    for per_cell in _count_:
        indiv_cell=DF.loc[per_cell].tolist()
        per_column_indiv_cells.append(indiv_cell)
        sub_=list(np.subtract(per_column_sum_all_cells, indiv_cell))
        per_column_sum_all_cells_excluding_indiv.append(sub_)
        avg_=[]
        for val in sub_:
            avg_.append(val/float(n-1))
        per_column_avg_all_cells_excluding_indiv.append(avg_)
    return per_column_avg_all_cells_excluding_indiv, per_column_indiv_cells

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
    cnt+=1
    
    adata=pickle_[i]
    adata_X=adata.to_df()

    # --------------------------------------------------------------------------

    # =========================================================== PLOTS-2: =============================================================================

    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    adata_raw=adata.copy()

    # ======================== PCA: ========================

    sc.tl.pca(adata_raw, svd_solver='arpack')

    # # ================================================== PLOT-6: ==========================================================

    sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=3)
    sc.tl.leiden(adata_raw)
    leiden_clusters_sorted=[str(j) for j in sorted([int(i) for i in np.unique(adata_raw.obs['leiden'])])]
    # scs.inference.nested_model(adata_raw)
    scs.inference.planted_model(adata_raw)
    ppbm_clusters_sorted=[str(j) for j in sorted([int(i) for i in np.unique(adata_raw.obs['ppbm'])])]
    
    # -------------------------------------------------
    sc.tl.umap(adata_raw)
    sc.pl.umap(adata_raw, color='leiden')
    sc.pl.umap(adata_raw, color=['ppbm'], legend_loc='on data')
    # ====================================================================================

    cell_marker_df=pd.read_excel('Cell_marker_Human.xlsx')

    ### Clustering analyses:
        
    
    cell_marker_df = cell_marker_df[~cell_marker_df.duplicated(['cell_name', 'Symbol'])]
        
    dc.run_ora(
        use_raw=False,
        mat=adata_raw,
        net=cell_marker_df,
        source='cell_name',
        target='Symbol',
        verbose=True
    )

    
    sc.pp.scale(adata_raw)
    
    value_error_flag_rankGenesLeiden=1
    while (value_error_flag_rankGenesLeiden==1):
        try:
            sc.tl.rank_genes_groups(adata_raw, groupby='leiden', pts=True, method='wilcoxon')
            time.sleep(10)
            value_error_flag_rankGenesLeiden=0
        except:
            pass
    DEDF_leiden=sc.get.rank_genes_groups_df(adata_raw, group=leiden_clusters_sorted)
    DEDF_leiden_AscPVals=categorize(DEDF_leiden, leiden_clusters_sorted)
    
    value_error_flag_rankGenesPPBM=1
    while (value_error_flag_rankGenesPPBM==1):
        try:
            sc.tl.rank_genes_groups(adata_raw, groupby='ppbm', pts=True, method='wilcoxon')
            time.sleep(10)
            value_error_flag_rankGenesPPBM=0
        except:
            pass
    DEDF_ppbm=sc.get.rank_genes_groups_df(adata_raw, group=ppbm_clusters_sorted)
    DEDF_ppbm_AscPVals=categorize(DEDF_ppbm, ppbm_clusters_sorted)
    
    adataRaw_rankGenes_names, adataRaw_rankGenes_scores, adataRaw_rankGenes_logfoldchanges, adataRaw_rankGenes_pvals, adataRaw_rankGenes_pvals_adj, adataRaw_rankGenes_pts, adataRaw_rankGenes_pts_rest=pd.DataFrame(adata_raw.uns['rank_genes_groups']['names']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['scores']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['logfoldchanges']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['pvals']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['pvals_adj']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['pts']), pd.DataFrame(adata_raw.uns['rank_genes_groups']['pts_rest'])
    
    sc.pl.rank_genes_groups(adata_raw, sharey=False, save='sc.pl.rank_genes_groups.pdf')
    
    
    # ======================================================================================
    
    DEDF_ppbm_AscPVals_0Point05_cutoff=DEDF_ppbm_AscPVals[DEDF_ppbm_AscPVals.pvals_adj<10**-100]
    DEDF_ppbm_AscPVals_0Point05_cutoff_groups=list(np.unique(DEDF_ppbm_AscPVals_0Point05_cutoff.group))
    DEDF_ppbm_AscPVals_0Point05_cutoff_groups_sorted=[str(j) for j in sorted([int(i) for i in np.unique(DEDF_ppbm_AscPVals_0Point05_cutoff_groups)])]
    
    ppbm_no_of_clusters_dict_values=[[] for _ in range(len(DEDF_ppbm_AscPVals_0Point05_cutoff_groups_sorted))]
    ppbm_no_of_clusters_dict_keys=DEDF_ppbm_AscPVals_0Point05_cutoff_groups_sorted
    ppbm_no_of_clusters_dict=dict(zip(ppbm_no_of_clusters_dict_keys, ppbm_no_of_clusters_dict_values))
    for index, row in DEDF_ppbm_AscPVals_0Point05_cutoff.iterrows():
        ppbm_no_of_clusters_dict[row['group']].append(row['names'])
    
    # ======================================================================================
    
    df=adata_raw.obs['ppbm'].copy().to_frame()
    # ppbm_cell_lookup_keys=[str(j) for j in sorted([int(i) for i in np.unique(df.index)])]
    ppbm_cell_lookup_values=[[] for _ in range(len(ppbm_clusters_sorted))]
    # ppbm_cell_lookup_keys=[[] for _ in range(len(np.unique(ppbm_cell_lookup_keys)))]
    ppbm_cell_lookup_keys=ppbm_clusters_sorted
    ppbm_cell_lookup_dict=dict(zip(ppbm_cell_lookup_keys, ppbm_cell_lookup_values))
    for index, row in df.iterrows():
        ppbm_cell_lookup_dict[row['ppbm']].append(index)
        ppbm_cell_lookup_list=list(ppbm_cell_lookup_dict.values())
    cell_counts_actual_order=[]
    for cell_cnt in ppbm_cell_lookup_list:
        for indiv_cell_cnt in cell_cnt:
            cell_counts_actual_order.append(indiv_cell_cnt)
    p_val_actual=[]
    for _count_ in ppbm_cell_lookup_list:
        per_column_avg_all_cells_excluding_indiv, per_column_indiv_cells=calculate_average_expressions_per_cell(adata_X, _count_)
        pvals=[]
        cluster_count=-1
        for i_count in range(len(_count_)):
            cluster_count+=1
            U1, p = mannwhitneyu(per_column_avg_all_cells_excluding_indiv[i_count], per_column_indiv_cells[i_count], method="exact")
            p_val_actual.append({cell_counts_actual_order[cluster_count]:p})
    # p_val_actual_sorted=____
    
    total_epochs=3
    pvals=[[] for _ in range(total_epochs)]
    cell_ct_ppbm=df.index.tolist()
    cell_expression_randomizations=[]
    for epoch in range(total_epochs):
        print('Epoch: '+str(epoch))
        print('=======================================================')
        for cell_counter in cell_ct_ppbm:
            print('Cell id '+str(cell_counter)+' in epoch '+str(epoch))
            _df_=adata_X.loc[cell_counter].tolist()
            random.shuffle(_df_)
            # cell_expression_randomizations.append(cell_counter: asdf)
            _df_simulated=adata_X
            _df_simulated.loc[cell_counter]=_df_
            
            # ---------------------------------------------------
            
            adata_simulated=sc.AnnData(_df_simulated)
            sc.tl.pca(adata_simulated, svd_solver='arpack')
            sc.pp.neighbors(adata_simulated, n_neighbors=10, n_pcs=3)
            scs.inference.planted_model(adata_simulated)
            
            # ---------------------------------------------------
            
            df=adata_simulated.obs['ppbm'].copy().to_frame()
            _ppbm_clusters_sorted_=[str(j) for j in sorted([int(i) for i in np.unique(adata_simulated.obs['ppbm'])])]
            # ppbm_cell_lookup_keys=[str(j) for j in sorted([int(i) for i in np.unique(df.index)])]
            ppbm_cell_lookup_values=[[] for _ in range(len(_ppbm_clusters_sorted_))]
            # ppbm_cell_lookup_keys=[[] for _ in range(len(np.unique(ppbm_cell_lookup_keys)))]
            ppbm_cell_lookup_keys=_ppbm_clusters_sorted_
            ppbm_cell_lookup_dict=dict(zip(ppbm_cell_lookup_keys, ppbm_cell_lookup_values))
            for index, row in df.iterrows():
                ppbm_cell_lookup_dict[row['ppbm']].append(index)
                ppbm_cell_lookup_list=list(ppbm_cell_lookup_dict.values())
            cell_counts_actual_order=[]
            for cell_cnt in ppbm_cell_lookup_list:
                for indiv_cell_cnt in cell_cnt:
                    cell_counts_actual_order.append(indiv_cell_cnt)
            
            # ---------------------------------------------------
            pvals[epoch]=[]
            for _count_ in ppbm_cell_lookup_list:
                per_column_avg_all_cells_excluding_indiv, per_column_indiv_cells=calculate_average_expressions_per_cell(adata_X, _count_)
                cluster_count=-1
                for i_count in range(len(_count_)):
                    cluster_count+=1
                    U1, p = mannwhitneyu(per_column_avg_all_cells_excluding_indiv[i_count], per_column_indiv_cells[i_count], method="exact")
                    pvals.append({cell_counts_actual_order[cluster_count]:p})
    
    

    
    
    