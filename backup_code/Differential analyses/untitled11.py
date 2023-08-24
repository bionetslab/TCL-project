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
first_three_items = pickle_items[0:5]
pickle_new=dict(first_three_items)
pickle_=pickle_new
# pickle_=dict(pickle_items)
# ------
   
no_of_patients=np.shape(pickle_)
error_samples=[]



_list_=[]
_patient_id_=[]
_cell_number_in_patient_=[]
_recurrence_=[]
_survival_days_binary_=[]
_obsm_spatial_=[]

no_of_observations=[]
patient_ids=[]


for i in pickle_:
    cnt+=1
    patient_ids.append(i)
    
    adata=pickle_[i]
    no_of_observations.append(adata.n_obs)
    
    if cnt==0:
        _adata_columns_=adata.to_df().columns.tolist()
    
    cnt_internal=-1
    for j in adata.X.tolist():
        cnt_internal+=1
        _list_.append(j)
        _patient_id_.append(i)
        _cell_number_in_patient_.append(adata.to_df().index[cnt_internal])
        _recurrence_.append(recurrence_labels[cnt])
        _survival_days_binary_.append(survival_days_binary[cnt])
    for k in adata.obsm['spatial']:
        _obsm_spatial_.append(k)
        # if adata.obsm['recurrence'].tolist()[cnt_internal]=='POSITIVE':
        #     _recurrence_.append(1)
        # else:
        #     _recurrence_.append(0)

_df_ = pd.DataFrame(_list_, columns=_adata_columns_) 
anndata_overall = sc.AnnData(_df_)

anndata_overall.obsm['patient_id']=np.array(_patient_id_)
anndata_overall.obsm['cell_number_in_patient']=np.array(_cell_number_in_patient_)
anndata_overall.obsm['spatial']=np.array(_obsm_spatial_)

anndata_overall.obsm['recurrence_binary']=np.array(_recurrence_)
anndata_overall.obsm['survival_days_binary']=np.array(_survival_days_binary_)

sc.pp.neighbors(anndata_overall, n_neighbors=10, n_pcs=3)
sc.tl.leiden(anndata_overall)
leiden_clusters=anndata_overall.obs['leiden']


######################################################################################################



cell_marker_df=pd.read_excel('Cell_marker_Human.xlsx')

### Clustering analyses:


cell_marker_df = cell_marker_df[~cell_marker_df.duplicated(['cell_name', 'Symbol'])]

# anndata_overall.raw=anndata_overall
    
dc.run_ora(
    use_raw=False,
    mat=anndata_overall,
    net=cell_marker_df,
    source='cell_name',
    target='Symbol',
    verbose=True
)


sc.pp.scale(anndata_overall)

sc.tl.rank_genes_groups(anndata_overall, groupby='leiden')
adataRaw_rankGenes_names=pd.DataFrame(anndata_overall.uns['rank_genes_groups']['names'])
adataRaw_rankGenes_names = adataRaw_rankGenes_names.reset_index(drop=True)
adataRaw_rankGenes_names_dict=adataRaw_rankGenes_names.to_dict()

cell_lookup_dict_leiden = {i: [adataRaw_rankGenes_names_dict[i][j] for j in adataRaw_rankGenes_names_dict[i]] for i in adataRaw_rankGenes_names_dict}

no_of_cell_types=1

annotation_dict_leiden={i: cell_lookup_dict_leiden[i][0] for i in cell_lookup_dict_leiden}

annotation_dict_leiden

# Add cell type column based on annotation
anndata_overall.obs['cell_type'] = [annotation_dict_leiden[clust] for clust in anndata_overall.obs['leiden']]

annotation_dict_leiden_int={int(k):v for k,v in annotation_dict_leiden.items()}
annotation_dict_leiden_sorted = dict(sorted(annotation_dict_leiden_int.items()))

leiden_cluster_keys_sorted = list(annotation_dict_leiden_sorted.keys())

cell_type_clusters_sorted=list(sorted(annotation_dict_leiden.values(), key=str.lower))

cell_type_clusters_sorted_combinations=list(itertools.combinations(cell_type_clusters_sorted, 2))

no_of_leiden_combinations=np.shape(cell_type_clusters_sorted_combinations)






# sns.heatmap(leiden_nhood_enrichment_zscore, cmap ='inferno', linewidths = 0.0)
# # adata_raw.uns['leiden_nhood_enrichment']=


# U-map:
sc.tl.umap(anndata_overall)
sc.pl.umap(anndata_overall, color='leiden', save='_leidenClustering_allPatients.pdf')
sc.pl.umap(anndata_overall, color='cell_type', save='_leidenClusteringCellTypesAnnotated_allPatients.pdf')

# sc.pl.umap(adata_raw, color='leiden')

adata_obs_celltype_df=anndata_overall.obs['cell_type']



######################################################################################################

patient_ids_percluster=anndata_overall.obsm['patient_id']
cellNumberInPatient_percluster=anndata_overall.obsm['cell_number_in_patient']

celltype_clusters=anndata_overall.obs['cell_type']
_sum_=-1
leiden_cluster_per_patient=[]
cell_type_per_patient=[]
for i in no_of_observations:
    leiden_cluster_=[]
    cell_type_=[]
    for j in range(1, i+1):
        _sum_+=1
        leiden_cluster_.append(leiden_clusters[_sum_])
        cell_type_.append(celltype_clusters[_sum_])
    leiden_cluster_per_patient.append(leiden_cluster_)
    cell_type_per_patient.append(cell_type_)


for i in range(len(no_of_observations)):
    pickle_[patient_ids[i]].obs['leiden']=leiden_cluster_per_patient[i]
    pickle_[patient_ids[i]].obs['cell_type']=cell_type_per_patient[i]

nhood_enrichment_zscores_recurrence_MWU_pvals=[]

error_samples=[]
nhood_enrichment_zscores_recurrence_1=[]
nhood_enrichment_zscores_recurrence_0=[]
_count_=-1
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
    sc.tl.umap(anndata_overall)
    leiden_clustering_plot_suffix='_leidenClustering_patientID='+str(i)+'.pdf'
    sc.pl.umap(anndata_overall, color='leiden', save=leiden_clustering_plot_suffix)
    leiden_clustering_withCellTypeAnnotation_plot_suffix='_leidenClusteringCellTypesAnnotated_patientID='+str(i)+'.pdf'
    sc.pl.umap(anndata_overall, color='cell_type', save=leiden_clustering_withCellTypeAnnotation_plot_suffix)
    
    # ============================================================================
    
    sq.gr.spatial_neighbors(adata_raw)
    error_flag=0
    error_=""
    
            
    try:
        leiden_nhood_enrichment, leiden_nhood_enrichment_zscore, leiden_nhood_enrichment_count, result_leiden_uns=neighborhood_enrichment(adata_raw, "leiden", True, "multiprocessing", False)
    except EOFError:
        error_samples.append(i)
    
    adata_raw.uns['leiden_nhood_enrichment']=result_leiden_uns
    
    if recurrence_labels=='POSITIVE':
        nhood_enrichment_zscores_recurrence_1.append(np.array(leiden_nhood_enrichment_zscore)[np.triu_indices(18, k = 1)].tolist())
        nhood_enrichment_zscores_recurrence_1_list=list(map(list, zip(*nhood_enrichment_zscores_recurrence_1)))
    else:
        nhood_enrichment_zscores_recurrence_0.append(np.array(leiden_nhood_enrichment_zscore)[np.triu_indices(18, k = 1)].tolist())
        nhood_enrichment_zscores_recurrence_0_list=list(map(list, zip(*nhood_enrichment_zscores_recurrence_0)))

for i in range(len(nhood_enrichment_zscores_recurrence_1)):
    U1, p = mannwhitneyu(nhood_enrichment_zscores_recurrence_1, nhood_enrichment_zscores_recurrence_0, method="exact")
    nhood_enrichment_zscores_recurrence_MWU_pvals.append(p)

nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected=fdrcorrection(nhood_enrichment_zscores_recurrence_MWU_pvals, alpha=0.05, method='indep', is_sorted=False)[1].tolist()
nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected_binary=fdrcorrection(nhood_enrichment_zscores_recurrence_MWU_pvals, alpha=0.05, method='indep', is_sorted=False)[0].tolist()

nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected_consensus=[]
for i in nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected_binary:
    if i==True:
        nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected_consensus.append('Null hypothesis upheld - no significant difference!')
    else:
        nhood_enrichment_zscores_recurrence_MWU_pvals_multipleTestingCorrected_consensus.append('Null hypothesis rejected - significant difference found!')

# # for i in range(len(leiden_clusters)):
# #     pickle_[patient_ids_percluster[i]].obs['leiden']=leiden_clusters[i]
    

# # =================================================================================================

# sc.pp.calculate_qc_metrics(anndata_overall, percent_top=None, log1p=False, inplace=True)

# # Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells:
# sc.pp.normalize_total(anndata_overall, target_sum=1e4)


# # Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.
# # You can get back an AnnData of the object in .raw by calling .raw.to_adata().
# anndata_overall.raw = anndata_overall

# # If you don’t proceed below with correcting the data with sc.pp.regress_out and scaling it via sc.pp.scale, you can also get away without using .raw at all.
# # The result of the previous highly-variable-genes detection is stored as an annotation in .var.highly_variable and auto-detected by PCA and hence, sc.pp.neighbors and subsequent manifold/graph tools. In that case, the step actually do the filtering below is unnecessary, too.


# # # Actually do the filtering:
# # adata = adata[:, adata.var.highly_variable]

# # Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance:
# # sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
# sc.pp.regress_out(anndata_overall, ['total_counts'])

# # Scale each gene to unit variance. Clip values exceeding standard deviation 10.
# sc.pp.scale(anndata_overall, max_value=10)

# # ======================== PCA: ========================

# sc.tl.pca(anndata_overall, svd_solver='arpack')

# # # ================================================== PLOT-6: ==========================================================


# sc.pp.neighbors(anndata_overall, n_neighbors=10, n_pcs=3)
# # -------------------------------------------------
# sc.tl.umap(anndata_overall)
# sc.pl.umap(anndata_overall, color='leiden')
# # ====================================================================================

# cell_marker_df=pd.read_excel('Cell_marker_Human.xlsx')

# ### Clustering analyses:
    

# cell_marker_df = cell_marker_df[~cell_marker_df.duplicated(['cell_name', 'Symbol'])]
    
# dc.run_ora(
#     mat=anndata_overall,
#     net=cell_marker_df,
#     source='cell_name',
#     target='Symbol',
#     verbose=True
# )


# sc.pp.scale(anndata_overall)

# sc.tl.rank_genes_groups(anndata_overall, groupby='leiden')
# adataRaw_rankGenes_names=pd.DataFrame(anndata_overall.uns['rank_genes_groups']['names'])
# adataRaw_rankGenes_names = adataRaw_rankGenes_names.reset_index(drop=True)
# adataRaw_rankGenes_names_dict=adataRaw_rankGenes_names.to_dict()

# cell_lookup_dict_leiden = {i: [adataRaw_rankGenes_names_dict[i][j] for j in adataRaw_rankGenes_names_dict[i]] for i in adataRaw_rankGenes_names_dict}

# no_of_cell_types=1

# annotation_dict_leiden={i: cell_lookup_dict_leiden[i][0] for i in cell_lookup_dict_leiden}

# annotation_dict_leiden

# # Add cell type column based on annotation
# anndata_overall.obs['cell_type'] = [annotation_dict_leiden[clust] for clust in anndata_overall.obs['leiden']]

# annotation_dict_leiden_int={int(k):v for k,v in annotation_dict_leiden.items()}
# annotation_dict_leiden_sorted = dict(sorted(annotation_dict_leiden_int.items()))

# leiden_cluster_keys_sorted = list(annotation_dict_leiden_sorted.keys())

# sq.gr.spatial_neighbors(anndata_overall)
# error_flag=0
# error_=""

        
# try:
#     leiden_nhood_enrichment, leiden_nhood_enrichment_zscore, leiden_nhood_enrichment_count, result_leiden_uns=neighborhood_enrichment(anndata_overall, "leiden", True, "multiprocessing", False)
# except EOFError:
#     error_samples.append(i)

# anndata_overall.uns['leiden_nhood_enrichment']=result_leiden_uns

# if recurrence_labels=='POSITIVE':
#     nhood_enrichment_zscores_recurrence_1.append(np.array(leiden_nhood_enrichment_zscore)[np.triu_indices(18, k = 1)].tolist())
# else:
#     nhood_enrichment_zscores_recurrence_0.append(np.array(leiden_nhood_enrichment_zscore)[np.triu_indices(18, k = 1)].tolist())




# # sns.heatmap(leiden_nhood_enrichment_zscore, cmap ='inferno', linewidths = 0.0)
# # # adata_raw.uns['leiden_nhood_enrichment']=


# # Visualize
# sc.pl.umap(anndata_overall, color='cell_type')
# # sc.pl.umap(adata_raw, color='leiden')

# adata_obs_celltype_df=anndata_overall.obs['cell_type']

# # # --------------------------------------------------------------------


#     # # --------------------------------------------------------------------------

#     # # =========================================================== PLOTS-2: =============================================================================

#     # sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

#     # # Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells:
#     # sc.pp.normalize_total(adata, target_sum=1e4)
        
#     # # # Logarithmize the data:
#     # # sc.pp.log1p(adata)

#     # # # Identify highly-variable genes.
#     # # sc.pp.highly_variable_genes(adata)

#     # # Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.
#     # # You can get back an AnnData of the object in .raw by calling .raw.to_adata().
#     # adata.raw = adata
#     # adata_raw=adata.copy()

#     # # If you don’t proceed below with correcting the data with sc.pp.regress_out and scaling it via sc.pp.scale, you can also get away without using .raw at all.
#     # # The result of the previous highly-variable-genes detection is stored as an annotation in .var.highly_variable and auto-detected by PCA and hence, sc.pp.neighbors and subsequent manifold/graph tools. In that case, the step actually do the filtering below is unnecessary, too.


#     # # # Actually do the filtering:
#     # # adata = adata[:, adata.var.highly_variable]

#     # # Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance:
#     # # sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
#     # sc.pp.regress_out(adata_raw, ['total_counts'])

#     # # Scale each gene to unit variance. Clip values exceeding standard deviation 10.
#     # sc.pp.scale(adata_raw, max_value=10)


#     # # ======================== PCA: ========================

#     # sc.tl.pca(adata_raw, svd_solver='arpack')

#     # # # ================================================== PLOT-6: ==========================================================

#     # sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=3)
#     # s2=sc.tl.leiden(adata_raw)
#     # sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=3)
#     # # -------------------------------------------------
#     # sc.tl.umap(adata_raw)
#     # sc.pl.umap(adata_raw, color='leiden')
#     # # ====================================================================================

#     # cell_marker_df=pd.read_excel('Cell_marker_Human.xlsx')

#     # ### Clustering analyses:
        
    
#     # cell_marker_df = cell_marker_df[~cell_marker_df.duplicated(['cell_name', 'Symbol'])]
        
#     # dc.run_ora(
#     #     mat=adata_raw,
#     #     net=cell_marker_df,
#     #     source='cell_name',
#     #     target='Symbol',
#     #     verbose=True
#     # )

    
#     # sc.pp.scale(adata_raw)
    
#     # sc.tl.rank_genes_groups(adata_raw, groupby='leiden')
#     # adataRaw_rankGenes_names=pd.DataFrame(adata_raw.uns['rank_genes_groups']['names'])
#     # adataRaw_rankGenes_names = adataRaw_rankGenes_names.reset_index(drop=True)
#     # adataRaw_rankGenes_names_dict=adataRaw_rankGenes_names.to_dict()
    
#     # cell_lookup_dict_leiden = {i: [adataRaw_rankGenes_names_dict[i][j] for j in adataRaw_rankGenes_names_dict[i]] for i in adataRaw_rankGenes_names_dict}

#     # no_of_cell_types=1
    
#     # annotation_dict_leiden={i: cell_lookup_dict_leiden[i][0] for i in cell_lookup_dict_leiden}
    
#     # annotation_dict_leiden

#     # # Add cell type column based on annotation
#     # adata_raw.obs['cell_type'] = [annotation_dict_leiden[clust] for clust in adata_raw.obs['leiden']]

#     # annotation_dict_leiden_int={int(k):v for k,v in annotation_dict_leiden.items()}
#     # annotation_dict_leiden_sorted = dict(sorted(annotation_dict_leiden_int.items()))

#     # leiden_cluster_keys_sorted = list(annotation_dict_leiden_sorted.keys())

#     # sq.gr.spatial_neighbors(adata_raw)
#     # error_flag=0
#     # error_=""
    
            
#     # try:
#     #     leiden_nhood_enrichment, leiden_nhood_enrichment_zscore, leiden_nhood_enrichment_count, result_leiden_uns=neighborhood_enrichment(adata_raw, "leiden", True, "multiprocessing", False)
#     # except EOFError:
#     #     error_samples.append(i)
    
#     # adata_raw.uns['leiden_nhood_enrichment']=result_leiden_uns
    
#     # if recurrence_labels=='POSITIVE':
#     #     nhood_enrichment_zscores_recurrence_1.append(np.array(leiden_nhood_enrichment_zscore)[np.triu_indices(18, k = 1)].tolist())
#     # else:
#     #     nhood_enrichment_zscores_recurrence_0.append(np.array(leiden_nhood_enrichment_zscore)[np.triu_indices(18, k = 1)].tolist())
    
    
    
    
#     # # # sns.heatmap(leiden_nhood_enrichment_zscore, cmap ='inferno', linewidths = 0.0)
#     # # # # adata_raw.uns['leiden_nhood_enrichment']=


#     # # Visualize
#     # sc.pl.umap(adata_raw, color='cell_type')
#     # # sc.pl.umap(adata_raw, color='leiden')

#     # adata_obs_celltype_df=adata_raw.obs['cell_type']

#     # # # --------------------------------------------------------------------
    
    

    
    
    