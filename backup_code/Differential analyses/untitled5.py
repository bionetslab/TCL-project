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
patient_id=20751

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

    # --------------------------------------------------------------------------

    # =========================================================== PLOTS-2: =============================================================================

    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    # Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells:
    sc.pp.normalize_total(adata, target_sum=1e4)
        
    # # Logarithmize the data:
    # sc.pp.log1p(adata)

    # # Identify highly-variable genes.
    # sc.pp.highly_variable_genes(adata)

    # Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.
    # You can get back an AnnData of the object in .raw by calling .raw.to_adata().
    # adata.raw = adata
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


    # ======================== PCA: ========================

    sc.tl.pca(adata_raw, svd_solver='arpack')

    # # ================================================== PLOT-6: ==========================================================

    sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=3)
    s2=sc.tl.leiden(adata_raw)
    sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=3)
    # -------------------------------------------------
    sc.tl.umap(adata_raw)
    sc.pl.umap(adata_raw, color='leiden')
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
    
    acts = dc.get_acts(adata_raw, obsm_key='ora_estimate')

    # We need to remove inf and set them to the maximum value observed
    acts_v = acts.X.ravel()
    max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
    acts.X[~np.isfinite(acts.X)] = max_e

    # We can scale the obtained activities for better visualizations
    sc.pp.scale(acts)
    acts
    
    sc.tl.rank_genes_groups(adata_raw, groupby='leiden')
    adataRaw_rankGenes_names=pd.DataFrame(adata_raw.uns['rank_genes_groups']['names'])
    adataRaw_rankGenes_names = adataRaw_rankGenes_names.reset_index(drop=True)
    adataRaw_rankGenes_names_dict=adataRaw_rankGenes_names.to_dict()
    
    cell_lookup_dict_leiden = {i: [adataRaw_rankGenes_names_dict[i][j] for j in adataRaw_rankGenes_names_dict[i]] for i in adataRaw_rankGenes_names_dict}

    no_of_cell_types=1
    
    annotation_dict_leiden={i: cell_lookup_dict_leiden[i][0] for i in cell_lookup_dict_leiden}
    
    annotation_dict_leiden

    # Add cell type column based on annotation
    adata_raw.obs['cell_type'] = [annotation_dict_leiden[clust] for clust in adata_raw.obs['leiden']]

    annotation_dict_leiden_int={int(k):v for k,v in annotation_dict_leiden.items()}
    annotation_dict_leiden_sorted = dict(sorted(annotation_dict_leiden_int.items()))

    leiden_cluster_keys_sorted = list(annotation_dict_leiden_sorted.keys())

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
    else:
        nhood_enrichment_zscores_recurrence_0.append(np.array(leiden_nhood_enrichment_zscore)[np.triu_indices(18, k = 1)].tolist())
    
    
    
    
    # # sns.heatmap(leiden_nhood_enrichment_zscore, cmap ='inferno', linewidths = 0.0)
    # # # adata_raw.uns['leiden_nhood_enrichment']=


    # Visualize
    sc.pl.umap(adata_raw, color='cell_type')
    # sc.pl.umap(adata_raw, color='leiden')

    adata_obs_celltype_df=adata_raw.obs['cell_type']

    # # --------------------------------------------------------------------
    
    

    
    
    