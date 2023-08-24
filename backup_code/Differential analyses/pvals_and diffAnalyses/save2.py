
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
from statsmodels.stats.multitest import fdrcorrection

with open('TNBC_41patients_KerenEtAl_processed.pkl', 'rb') as f:
   pickle_= pickle.load(f)


leiden_clusters=[]

for i in pickle_:
    for j in np.unique(pickle_[i].obs['leiden'].tolist()):
        leiden_clusters.append(j)

leiden_clusters_sorted=list(sorted(np.unique(leiden_clusters), key=str.lower))
leiden_clusters_sorted_combinations=list(itertools.combinations(leiden_clusters_sorted, 2))

leiden_clusters_sorted_combinations_str=[]

for i in leiden_clusters_sorted_combinations:
    k_=-1
    str_=''
    for j in i:
        k_+=1
        if k_==0:
            str_+=str(j)+', '
            # leiden_clusters_sorted_combinations_str.append()
        elif k_==1:
            str_+=str(j)
    leiden_clusters_sorted_combinations_str.append(str_)



nhood_enrichment_zscores_recurrence_MWU_pvals=[]

error_samples=[]
nhood_enrichment_zscores_recurrence_1=[]
nhood_enrichment_zscores_recurrence_0=[]
_count_=-1

retrieved_leiden_clusters=[]
retrieved_celltype_clusters=[]

# recurrence_labels=[]
# for i in pickle_:
#     recurrence_labels.append(np.unique(pickle_[i].obsm['recurrence'])[0])

for i in pickle_:
    _count_+=1
    adata=pickle_[i]
    
    
    # ============================================================================
    
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    # Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells:
    sc.pp.normalize_total(adata, target_sum=1e4)
        
    # # Logarithmize the data:
    # sc.pp.log1p(adata)

    # # Identify highly-variable genes.
    # sc.pp.highly_variable_genes(adata)

    # Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.
    # You can get back an AnnData of the object in .raw by calling .raw.to_adata().
    adata.raw = adata
    adata_raw=adata.copy()

    # If you don’t proceed below with correcting the data with sc.pp.regress_out and scaling it via sc.pp.scale, you can also get away without using .raw at all.
    # The result of the previous highly-variable-genes detection is stored as an annotation in .var.highly_variable and auto-detected by PCA and hence, sc.pp.neighbors and subsequent manifold/graph tools. In that case, the step actually do the filtering below is unnecessary, too.


    # # Actually do the filtering:
    # adata = adata[:, adata.var.highly_variable]

    # Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance:
    # sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.regress_out(adata_raw, ['total_counts'])

    # Scale each gene to unit variance. Clip values exceeding standard deviation 10.
    sc.pp.scale(adata_raw, max_value=10)


    # ------------------------------ PCA: ------------------------------

    sc.tl.pca(adata_raw, svd_solver='arpack')
    
    # ============================================================================
    sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=3)
    sc.tl.umap(adata_raw)
    leiden_clustering_plot_suffix='_leidenClustering_patientID='+str(i)+'.pdf'
    sc.pl.umap(adata_raw, color='leiden', save=leiden_clustering_plot_suffix)
    leiden_clustering_withCellTypeAnnotation_plot_suffix='_leidenClusteringCellTypesAnnotated_patientID='+str(i)+'.pdf'
    sc.pl.umap(adata_raw, color='cell_type', save=leiden_clustering_withCellTypeAnnotation_plot_suffix)
    
    # ============================================================================
    retr_leiden=[]
    retr_celltype=[]
    for _j_ in adata_raw.obs['leiden']:
        retr_leiden.append(_j_)
    retr_leiden=np.unique(retr_leiden)
    retrieved_leiden_clusters.append(retr_leiden)
    retrieved_celltype_clusters.append(retr_celltype)
    # ============================================================================
    
    # sq.gr.spatial_neighbors(adata_raw)
            
    # try:
    #     leiden_nhood_enrichment, leiden_nhood_enrichment_zscore, leiden_nhood_enrichment_count, result_leiden_uns=neighborhood_enrichment(adata_raw, "leiden", True, "multiprocessing", False)
    # except EOFError:
    #     error_samples.append(i)
    
    # adata_raw.uns['leiden_nhood_enrichment']=result_leiden_uns
    
    # recurrence_label=np.unique(pickle_[i].obsm['recurrence'])[0]
    
    # if recurrence_label=='POSITIVE':
    #     nhood_enrichment_zscores_recurrence_1.append(np.array(leiden_nhood_enrichment_zscore)[np.triu_indices(np.shape(leiden_nhood_enrichment_zscore)[0], k = 1)].tolist())
    #     nhood_enrichment_zscores_recurrence_1_list=list(map(list, zip(*nhood_enrichment_zscores_recurrence_1)))
    # else:
    #     nhood_enrichment_zscores_recurrence_0.append(np.array(leiden_nhood_enrichment_zscore)[np.triu_indices(np.shape(leiden_nhood_enrichment_zscore)[0], k = 1)].tolist())
    #     nhood_enrichment_zscores_recurrence_0_list=list(map(list, zip(*nhood_enrichment_zscores_recurrence_0)))
    
common_leiden_clusters=list(set().union(*retrieved_leiden_clusters))
common_leiden_clusters_int=[int(k) for k in common_leiden_clusters]
common_leiden_clusters_sorted = sorted(common_leiden_clusters_int)
common_leiden_clusters_sorted=[str(k) for k in common_leiden_clusters_sorted]
# ----------
common_celltype_clusters=list(set().union(*retrieved_celltype_clusters))
common_celltype_clusters_sorted = sorted(common_celltype_clusters)

# for i in pickle_:
#     _df_=pickle_[i].obs['leiden'].to_df()
#     _df_=_df_[common_leiden_clusters_sorted]

for i in pickle_:
    _index_=pickle_[i].to_df().index.tolist()
    _uns_spatial_=pickle_[i].obsm['spatial'].tolist()
    index_uns_DICT=dict(zip(_index_, _uns_spatial_))
    
    leiden_df=pickle_[i].obs['leiden'].to_frame()
    leiden_df_new = leiden_df[leiden_df['leiden'].isin(common_leiden_clusters_sorted)]
    # leiden_series_new=leiden_df_new.squeeze()
    leiden_df_new_index=leiden_df_new.index.tolist()
    
    spatial_new=[index_uns_DICT[x] for x in leiden_df_new_index]
    
    adata_X=pickle_[i].to_df()
    adata_X_new=adata_X[adata_X.index.isin(leiden_df_new_index)]
    adata_temp=sc.AnnData(adata_X_new)
    # adata_temp=adata_temp.raw
    adata_temp.obsm['spatial']=np.array(spatial_new)
    # adata_temp.obs['leiden']=leiden_df_new.leiden.tolist()
    # adata_temp.obs['leiden']=pd.Series(leiden_df_new.leiden.tolist(), dtype='category')
    leids_=[]
    
    for p_ in leiden_df_new.leiden.tolist():
        leids_.append(p_)
    
    adata_temp.obs['leiden']=pd.Series(leids_, dtype='category')
    
    sq.gr.spatial_neighbors(adata_temp)
    
    try:
        leiden_nhood_enrichment, leiden_nhood_enrichment_zscore, leiden_nhood_enrichment_count, result_leiden_uns=neighborhood_enrichment(adata_temp, "leiden", True, "multiprocessing", False)
    except EOFError:
        error_samples.append(i)
    
    pickle_[i].uns['leiden_nhood_enrichment_commonClusters']=result_leiden_uns
    
    recurrence_label=np.unique(pickle_[i].obsm['recurrence'])[0]
    
    
    if recurrence_label=='POSITIVE':
        nhood_enrichment_zscores_recurrence_1.append(np.array(leiden_nhood_enrichment_zscore)[np.triu_indices(np.shape(leiden_nhood_enrichment_zscore)[0], k = 1)].tolist())
        nhood_enrichment_zscores_recurrence_1_list=list(map(list, zip(*nhood_enrichment_zscores_recurrence_1)))
    else:
        nhood_enrichment_zscores_recurrence_0.append(np.array(leiden_nhood_enrichment_zscore)[np.triu_indices(np.shape(leiden_nhood_enrichment_zscore)[0], k = 1)].tolist())
        nhood_enrichment_zscores_recurrence_0_list=list(map(list, zip(*nhood_enrichment_zscores_recurrence_0)))
    
    
    
    
    # adata_X_new=adata_X[common_leiden_clusters_sorted]
    # adata_temp=sc.AnnData(adata_X_new)
    # _df_=pickle_[i].obs['cell_type'].to_df()
    # _df_=_df_[common_leiden_clusters_sorted]

# ===============

# # Add cell type column based on annotation
# anndata_overall.obs['cell_type'] = [annotation_dict_leiden[clust] for clust in anndata_overall.obs['leiden']]

# annotation_dict_leiden_int={int(k):v for k,v in annotation_dict_leiden.items()}
# annotation_dict_leiden_sorted = dict(sorted(annotation_dict_leiden_int.items()))

# leiden_cluster_keys_sorted = list(annotation_dict_leiden_sorted.keys())

# cell_type_clusters_sorted=list(sorted(annotation_dict_leiden.values(), key=str.lower))

############# cell_type_clusters_sorted_combinations=list(itertools.combinations(cell_type_clusters_sorted, 2))

############# no_of_leiden_combinations=np.shape(cell_type_clusters_sorted_combinations)

# ===============

# for i in range(len(nhood_enrichment_zscores_recurrence_1_list)):
#     U1, p = mannwhitneyu(nhood_enrichment_zscores_recurrence_1_list[i], nhood_enrichment_zscores_recurrence_0_list[i], method="exact")
#     nhood_enrichment_zscores_recurrence_MWU_pvals.append(p)

# nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected=fdrcorrection(nhood_enrichment_zscores_recurrence_MWU_pvals, alpha=0.05, method='indep', is_sorted=False)[1].tolist()
# nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected_binary=fdrcorrection(nhood_enrichment_zscores_recurrence_MWU_pvals, alpha=0.05, method='indep', is_sorted=False)[0].tolist()

# nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected_consensus=[]
# for i in nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected_binary:
#     if i==True:
#         nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected_consensus.append('Null hypothesis upheld - no significant difference!')
#     else:
#         nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected_consensus.append('Null hypothesis rejected - significant difference found!')

# nhood_results_list = [nhood_enrichment_zscores_recurrence_MWU_pvals, nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected,
#        nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected_binary, nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected_consensus]

# df = pd.DataFrame(nhood_results_list, index =['zscores_pvals', 'zscores_pvals_multipleTestingCorrected', 'zscores_pvals_multipleTestingCorrected_binary', 'zscores_pvals_multipleTestingCorrected_consensus'],
#                                               columns =leiden_clusters_sorted_combinations_str)
