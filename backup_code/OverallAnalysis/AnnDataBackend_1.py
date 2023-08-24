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
patient_id=20751

# ##########################  PLOTS ####################

# =========================================================== PLOTS-1: =============================================================================
with open(f'Protein_Expression_{patient_id}.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
    Protein_Expression_, df___extracellular_expression_count, adata___extracellular_expression_count, df___relativeAreaWise_extracellular_expression_count, adata___relativeAreaWise_extracellular_expression_count, df___cytoplasm_expression_count, adata___cytoplasm_expression_count, df___relativeAreaWise_cytoplasm_expression_count, adata___relativeAreaWise_cytoplasm_expression_count, df___nuclear_expression_count, adata___nuclear_expression_count, df___relativeAreaWise_nuclear_expression_count, adata___relativeAreaWise_nuclear_expression_count, df___cell_expression_count, adata___cell_expression_count, df___relativeAreaWise_cell_expression_count, adata___relativeAreaWise_cell_expression_count = pickle.load(f)
adata=adata___relativeAreaWise_cell_expression_count.copy()
adata_df=adata.to_df()
adata_df = adata_df.iloc[1:]
adata=sc.AnnData(adata_df)

protein_channels=['IL2RA',
'SPN',
'CD163',
'KRT14',
'SELP',
'CD16a',
'FUT4',
'CD63',
'ICAM1',
'MME',
'CD4',
'THY1',
'CD5',
'CD38',
'CD276',
'CD7',
'CD11b',
'ITGB1',
'CD6',
'CD2',
'ITGAX',
'TNFRSF1A',
'MRC1',
'ADAM10',
'CD3',
'VIM',
'NOTCH1',
'CD14',
'NOTCH3',
'ITGAL',
'TNFRSF1B',
'phase',
'CD24',
'THBD',
'CD27',
'CD36',
'LY75',
'NCAM1',
'FAS',
'CD9',
'CD53',
'CD52',
'CD44',
'CD8a',
'ANPEP',
'CD45RA',
'CD45RO',
'HLA-DR',
'HLA-ABC',
'PTPRC',
'CD55',
'CD69',
'TFRC',
'BCL2',
'CD68',
'CD40',
'CD81']


import os
from skimage import io
from skimage.util import img_as_ubyte
import cv2
import imageio as iio
from skimage import img_as_ubyte

# from PIL import image
# from skimage import io
# io.use_plugin('pil')

# images=[]
# image_names=[]
# # You need to change these to valid directories on your computer
# input_dir = os.path.dirname('/Users/surya/Documents/GITHUB-REPOSITORIES/melc_databases/MELC Data - T-cell lymphoma/processed/t-cell lymphoma/20751/')
# for f in os.listdir(input_dir):
#     if f.endswith('.tif') is True:
#         image_names.append(f)
#         images.append(io.imread(os.path.join(input_dir, f), as_gray=True))
# image_names_new=[]
# for i in image_names:
#     image_names_new.append(i.replace('.tif', ''))
# counter=0
# for i in images:
#     a=image_names_new[counter]
#     # fig1, ax1 = plt.subplots()
#     # plt.imshow(i, cmap='gray')
#     plt.imsave(f'NEW_IMAGES_SAVED_NOW/{a}.png', i, cmap='gray')
#     counter+=1

images=[]
image_names=[]
# You need to change these to valid directories on your computer
input_dir = os.path.dirname('/Users/surya/Documents/GITHUB-REPOSITORIES/spatial_proteomics/2. feature_extraction/NEW_IMAGES_SAVED_NOW/')
for f in os.listdir(input_dir):
    if f.endswith('.png') is True:
        image_names.append(f)
        images.append(io.imread(os.path.join(input_dir, f), as_gray=True))

# adata.uns[spatial_key][library_id]["images"] = {}
# adata.uns[spatial_key][library_id]["images"] = {"hires": image}
# adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 1, "spot_diameter_fullres": 0.5}

# adata_melc_df=pd.read_csv('Patient_20751_relativeAreaWise_cell_expression_count.csv')
# adata_melc_df.set_index('Unnamed: 0', inplace=True)
# adata_melc_df.index.names = ['index']
# adata=sc.AnnData(adata_melc_df)



# # adata_df=adata.to_df()
# # adata_df.index.names = ['index']
# # adata=sc.AnnData(adata_df)


column_names=adata_df.columns.tolist()
column_names_new=[]

for i in column_names:
    sep = '-'
    split_=(i.split(sep))
    column_names_new.append(split_[len(split_)-2])

# # adata_df.set_axis(column_names_new, axis=1,inplace=True)
adata_df=adata_df.set_axis(column_names_new, axis=1)
adata_new=sc.AnnData(adata_df)
adata=adata_new


# --------------------------------------------------------------------------

with open('protein_and_nucleus_lookup_matrices_plus_circle.pkl', 'rb') as f:
    A=pickle.load(f)
x_min, x_max, y_min, y_max, coords_x, coords_y=A[3], A[4], A[5], A[6], A[7], A[8]

cell_names_list=adata_df.index.tolist()

x_min_dict=dict(zip(cell_names_list, x_min))
x_max_dict=dict(zip(cell_names_list, x_max))
y_min_dict=dict(zip(cell_names_list, y_min))
y_max_dict=dict(zip(cell_names_list, x_max))


del x_min[0]
del x_max[0]
del y_min[0]
del y_max[0]
del coords_x[0]
del coords_y[0]


adata.obsm['x_min']=np.array(x_min)
adata.obsm['x_max']=np.array(x_max)
adata.obsm['y_min']=np.array(y_min)
adata.obsm['y_max']=np.array(y_max)

adata.obsm['coords_x']=np.array(coords_x)
adata.obsm['coords_y']=np.array(coords_y)

coordinates=np.column_stack((coords_x, coords_y))
adata.obsm['spatial']=coordinates

# --------------------------------------------------------------------------



# adata = sc.datasets.pbmc3k()
# ----------------------------------

fig=plt.figure(layout="constrained", figsize=[12,11])
axd=fig.subplot_mosaic(
[
 
 
  ['cellular_expression']
 
 
]
)
# identify_axes(axd, fontsize=36)
# -----------------------
# fig, ax = plt.subplots(1,3, figsize=(20,6))
# sc.pl.spatial(adata, img_key="hires", color="array_row", size=1.5, ax=ax[0], show=False)
# sc.pl.spatial(bdata, img_key="hires", color="array_row", size=1.5, ax=ax[1], show=False)
# sc.pl.spatial(cdata, img_key="hires", color="array_row", size=1.5, ax=ax[2], show=False)
# plt.tight_layout(pad=3.0)
# plt.show()
g3=sc.pl.highest_expr_genes(adata, n_top=20, ax=axd['cellular_expression'], show=False)
axd['cellular_expression'].set_title("$ \mathbf{C} $", fontsize = 17, loc='left')
# axd['nuclear_expression'].set(xlabel=None, ylabel=None)
axd['cellular_expression'].set_ylabel('Protein channels', fontsize=17)
axd['cellular_expression'].set_xlabel('% of total counts', fontsize=17)
# --------------------------------
# plt.subplots_adjust(left=0.1,
#                     bottom=0.1,
#                     right=0.9,
#                     top=0.9,
#                     wspace=0.4,
#                     hspace=0.4)

# # identify_axes(axd, fontsize=36)
# # fig.suptitle("$ \mathbf{A} $" + " ROBUST-Web backend runtime", fontsize=20, x=0.069, horizontalalignment='left')
# fig.tight_layout(pad=3.0)
fig.tight_layout()
fig.savefig('fig1 - scanpy - percentage of total counts.pdf')
plt.show()
# # ============================================================================================================
# # ============================================================================================================
# # ============================================================================================================

# =========================================================== PLOTS-2: =============================================================================

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
#              jitter=0.4, multi_panel=True)

# ----------------------------------

fig=plt.figure(layout="constrained", figsize=[12,11])
axd=fig.subplot_mosaic(
[
  ['cell_expression-n_genes_by_counts','cell_expression-total_counts']
]
)
# identify_axes(axd, fontsize=36)
# -----------------------
g5=sc.pl.violin(adata, keys='n_genes_by_counts', ax=axd['cell_expression-n_genes_by_counts'], jitter=0.4, show=False)
axd['cell_expression-n_genes_by_counts'].set_title("$ \mathbf{C} $", fontsize = 17, loc='left')
# axd['nuclear_expression-n_genes_by_counts'].set(xlabel=None, ylabel=None)
axd['cell_expression-n_genes_by_counts'].set_ylabel('No. of protein\nchannels', fontsize=17)
axd['cell_expression-n_genes_by_counts'].set_xlabel('No. of genes in count matrix', fontsize=17)

g6=sc.pl.violin(adata, keys='total_counts', ax=axd['cell_expression-total_counts'], jitter=0.4, show=False)
# axd['nuclear_expression-n_genes_by_counts'].set(xlabel=None, ylabel=None)
axd['cell_expression-total_counts'].set_ylabel('Total counts', fontsize=17)
axd['cell_expression-total_counts'].set_xlabel('Total counts per cell', fontsize=17)
# # --------------------------------
fig.tight_layout()
fig.savefig('fig2 - scanpy - noOfGenes and totalCounts.pdf')
plt.show()
# # # ============================================================================================================
# # # ============================================================================================================
# # # ============================================================================================================

# =========================================================== PLOTS-3: =============================================================================
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
# sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
# ----------------------------------

fig=plt.figure(layout="constrained", figsize=[12,11])
axd=fig.subplot_mosaic(
[
 
 
  ['cellular_expression']
 
 
]
)
# identify_axes(axd, fontsize=36)
# -----------------------
# fig, ax = plt.subplots(1,3, figsize=(20,6))
# sc.pl.spatial(adata, img_key="hires", color="array_row", size=1.5, ax=ax[0], show=False)
# sc.pl.spatial(bdata, img_key="hires", color="array_row", size=1.5, ax=ax[1], show=False)
# sc.pl.spatial(cdata, img_key="hires", color="array_row", size=1.5, ax=ax[2], show=False)
# plt.tight_layout(pad=3.0)
# plt.show()
g3=sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axd['cellular_expression'], show=False)
axd['cellular_expression'].set_title("$ \mathbf{C} $", fontsize = 17, loc='left')
# axd['nuclear_expression'].set(xlabel=None, ylabel=None)
axd['cellular_expression'].set_ylabel('No. of genes expressed\nin the counts matrix', fontsize=17)
axd['cellular_expression'].set_xlabel('Total counts per cell', fontsize=17)
# --------------------------------
# plt.subplots_adjust(left=0.1,
#                     bottom=0.1,
#                     right=0.9,
#                     top=0.9,
#                     wspace=0.4,
#                     hspace=0.4)

# # identify_axes(axd, fontsize=36)
# # fig.suptitle("$ \mathbf{A} $" + " ROBUST-Web backend runtime", fontsize=20, x=0.069, horizontalalignment='left')
# fig.tight_layout(pad=3.0)
fig.tight_layout()
fig.savefig('fig3 - scanpy - scatter plot (n_genes_count vs total_counts).pdf')
plt.show()
# # ============================================================================================================
# # ============================================================================================================
# # ============================================================================================================

# Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells:
sc.pp.normalize_total(adata, target_sum=1e4)
    
# Logarithmize the data:
sc.pp.log1p(adata)

# Identify highly-variable genes.
sc.pp.highly_variable_genes(adata)


# # # ================================================== PLOT-4: ==========================================================

sc.pl.highly_variable_genes(adata)

# # =================================================================================================================



# Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.
# You can get back an AnnData of the object in .raw by calling .raw.to_adata().
# adata.raw = adata
adata_raw=adata.copy()

# If you don’t proceed below with correcting the data with sc.pp.regress_out and scaling it via sc.pp.scale, you can also get away without using .raw at all.
# The result of the previous highly-variable-genes detection is stored as an annotation in .var.highly_variable and auto-detected by PCA and hence, sc.pp.neighbors and subsequent manifold/graph tools. In that case, the step actually do the filtering below is unnecessary, too.


# Actually do the filtering:
adata = adata[:, adata.var.highly_variable]

# Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance:
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.regress_out(adata_raw, ['total_counts'])

# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata_raw, max_value=10)


# ======================== PCA: ========================
sc.tl.pca(adata_raw, svd_solver='arpack')


# # # ================================================== PLOT-5: ==========================================================

fig=plt.figure(layout="constrained", figsize=[4,4])
axd=fig.subplot_mosaic([['cellular_expression']])
g1=sc.pl.pca(adata_raw, ax=axd['cellular_expression'])
# # -----------------------
# # plt.subplots_adjust(left=0.1,
# #                     bottom=0.1,
# #                     right=0.9,
# #                     top=0.9,
# #                     wspace=0.4,
# #                     hspace=0.4)

fig.tight_layout()
fig.savefig('PCA.pdf')
plt.show()
# # # =====================================================================================================================


# # ================================================== PLOT-6: ==========================================================

# sc.pl.pca(adata, color="CD3-PE", save='PCA-CD3-PE.pdf')

# # # ============ PLOT-7: ====================
sc.pl.pca_variance_ratio(adata_raw, log=True, save='pca-percentage-variance.pdf')
# # # ====================

sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=3)
s2=sc.tl.leiden(adata_raw)
sc.pl.pca_variance_ratio(adata_raw, log=True)
# adata.write(results_file)
sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=3)

sc.tl.paga(adata_raw)
sc.pl.paga(adata_raw, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata_raw, init_pos='paga')

sc.tl.umap(adata_raw)

# sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])

adata_raw.var_names

# sc.pl.umap(adata, color=['CD45RA-PE', 'CD62P-PE', 'CD4-PE'])

# sc.pl.umap(adata, color=['CD45RA-PE', 'CD62P-PE', 'CD4-PE'], use_raw=False)

sc.tl.leiden(adata_raw)

# sc.pl.umap(adata, color=['leiden', 'CD45RA-PE', 'CD62P-PE', 'CD4-PE'])

# sc.pl.umap(adata, color=['leiden', 'CD45RA-PE', 'CD62P-PE'])

sc.tl.rank_genes_groups(adata_raw, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata_raw, n_genes=25, sharey=False)

sc.settings.verbosity = 2

sc.tl.rank_genes_groups(adata_raw, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata_raw, n_genes=25, sharey=False)

sc.tl.rank_genes_groups(adata_raw, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata_raw, n_genes=25, sharey=False)

pd.DataFrame(adata_raw.uns['rank_genes_groups']['names'])
df=pd.DataFrame(adata_raw.uns['rank_genes_groups']['names'])

result = adata_raw.uns['rank_genes_groups']

groups = result['names'].dtype.names

# pd.DataFrame(
#     {group + '_' + key[:1]: result[key][group]
#     for group in groups for key in ['names', 'pvals']}).head(5)

sc.tl.rank_genes_groups(adata_raw, 'leiden', groups=['0'], reference='1', method='wilcoxon')


sc.pl.rank_genes_groups(adata_raw, groups=['0'], n_genes=20)

sc.pl.rank_genes_groups_violin(adata_raw, groups='0', n_genes=8)

sc.pl.rank_genes_groups_violin(adata_raw, groups='0', n_genes=8)

# sc.pl.violin(adata, ['CD45RA-PE', 'CD62P-PE', 'CD4-PE'], groupby='leiden')


# ====================================================================================

cell_marker_df=pd.read_excel('Cell_marker_Human.xlsx')









# # ********************* COMPUTING THE NEIGHBORHOOD GRAPH: *********************
# # Let us compute the neighborhood graph of cells using the PCA representation of the data matrix:
# sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# # Embedding the neighborhood graph:
# # i. using umap:
# # # sc.tl.paga(adata)
# s1=sc.tl.louvain(adata)
# s2=sc.tl.leiden(adata)

# sc.tl.paga(adata, groups='louvain')
# sc.tl.paga(adata, groups='leiden')



# # # ############################################# PLOT-7: #############################################
# sc.pl.paga(adata, color=['louvain','leiden'], node_size_scale=2, edge_width_scale=2, save="louvain and leiden (default, 45 clusters).pdf")  # remove `plot=False` if you want to see the coarse-grained graph
# # # ###################################################################################################



# # sc.tl.umap(adata, init_pos='paga')
# sc.tl.umap(adata)

### ===============================================================
### ===============================================================
### ===============================================================
### ===============================================================


### Clustering analyses:

cell_marker_df = cell_marker_df[~cell_marker_df.duplicated(['cell_name', 'Symbol'])]
    
dc.run_ora(
    mat=adata_raw,
    net=cell_marker_df,
    source='cell_name',
    target='Symbol',
    verbose=True
)

acts = dc.get_acts(adata_raw, obsm_key='ora_estimate')

# We need to remove inf and set them to the maximum value observed
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

# We can scale the obtained activities for better visualizations
sc.pp.scale(acts)
acts

sc.pl.umap(acts, color=['Adipose tissue', 'leiden'], cmap='RdBu_r', vcenter=0)
sc.pl.violin(acts, keys=['Adipose tissue'], groupby='leiden')

df = dc.rank_sources_groups(acts, groupby='leiden', reference='rest', method='t-test_overestim_var')

n_ctypes = 3
ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
ctypes_dict

sc.pl.matrixplot(acts, ctypes_dict, 'leiden', dendrogram=True,
                  colorbar_title='Z-scaled scores', vmin=-2, vmax=2, cmap='RdBu_r')

sc.pl.violin(acts, keys=['Adipose tissue', 'Airway', 'Bile duct', 'Bladder'], groupby='leiden')

annotation_dict_leiden = df.groupby('group').head(1).set_index('group')['names'].to_dict()
annotation_dict_leiden

# Add cell type column based on annotation
adata_raw.obs['cell_type'] = [annotation_dict_leiden[clust] for clust in adata_raw.obs['leiden']]

annotation_dict_leiden_int={int(k):v for k,v in annotation_dict_leiden.items()}
annotation_dict_leiden_sorted = dict(sorted(annotation_dict_leiden_int.items()))

# Visualize
sc.pl.umap(adata_raw, color='cell_type')

adata_obs_celltype_df=adata_raw.obs['cell_type']

# --------------------------------------------------------------------

import squidpy as sq

tissue_image_dict = {}

for i in protein_channels:
    tissue_image_dict[i] = {}


# rng = default_rng(42)
# image = rng.uniform(0, 1, size=(10, 10, 3))  # image

spatial_key = "spatial"
# library_id = "tissue42"
adata_raw.uns[spatial_key] = tissue_image_dict

count=0
for i in protein_channels:
    adata_raw.uns[spatial_key][i]["images"] = {}
    adata_raw.uns[spatial_key][i]["images"] = {"hires": images[count]}
    adata_raw.uns[spatial_key][i]["scalefactors"] = {"tissue_hires_scalef": 1, "spot_diameter_fullres": 0.5}

sq.pl.spatial_scatter(
    adata_raw, shape=None, color="cell_type", size=0.5, library_id="spatial", figsize=(10, 10)
)

# --------------------------------------------------------------------

time.sleep(10)

sq.gr.spatial_neighbors(adata_raw)
sq.gr.nhood_enrichment(adata_raw, cluster_key="cell_type")
sq.pl.nhood_enrichment(adata_raw, cluster_key="cell_type", method="average", figsize=(5, 5))

# --------------------------------------------------------------------

sq.gr.interaction_matrix(adata_raw, cluster_key="cell_type")
sq.pl.interaction_matrix(adata_raw, cluster_key="cell_type", method="average", figsize=(5, 5))
sq.gr.co_occurrence(adata_raw, cluster_key="cell_type")
sq.pl.co_occurrence(adata_raw, cluster_key="cell_type", clusters="Adipose tissue", figsize=(8, 5))

##### ============================================================================

# import squidpy as sq

# resolution = 1.5
# print("neighbors")
# sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=20)
# print("UMAP")
# sc.tl.umap(adata_raw)
# print("Leiden")
# sc.tl.leiden(adata_raw, resolution=resolution)

# sc.set_figure_params(figsize=(10, 10))
# sc.pl.umap(adata_raw, color=["leiden"], size=5)

# sq.pl.spatial_scatter(
#     adata_raw, shape=None, color="leiden", size=0.5, library_id="spatial", figsize=(10, 10)
# )

##### ============================================================================

# # sc.pl.spatial(adata_raw, color="cell_type", spot_size=10)
sc.pl.spatial(adata_raw, color="cell_type", spot_size=10, library_id='IL2RA')
# sq.pl.spatial_scatter(adata_raw, color="cell_type", library_id='IL2RA', img_cmap="gray")

##### ============================================================================

sq.gr.centrality_scores(
    adata_raw,
    cluster_key="cell_type",
)
sq.pl.centrality_scores(adata_raw, cluster_key="cell_type", figsize=(20, 5), s=500)

#### =============================================================================

sq.gr.co_occurrence(adata_raw, cluster_key="cell_type")
sq.pl.co_occurrence(
    adata_raw,
    cluster_key="cell_type",
    clusters=["Skin", "Uterus"],
    figsize=(15, 4),
)

##### ============================================================================

sq.pl.spatial_scatter(
    adata_raw,
    color="cell_type",
    groups=[
        "Spleen",
        "Skin",
        "Soft tissue",
        "Testis",
        "Uterus",
        "Spleen",
    ],
    shape=None,
    size=2,
)

##### ============================================================================

mode = "L"
sq.gr.ripley(adata_raw, cluster_key="cell_type", mode=mode, max_dist=500)
sq.pl.ripley(adata_raw, cluster_key="cell_type", mode=mode)

sq.pl.spatial_scatter(
    adata_raw,
    color="cell_type",
    groups=["Blood vessel", "Testis"],
    size=3,
    shape=None,
)

##### ============================================================================

sq.gr.spatial_autocorr(adata_raw, mode="moran")
adata_raw.uns["moranI"].head(10)

sq.pl.spatial_scatter(
    adata_raw,
    shape=None,
    color=["CD45RA", "CD4", "TNFR2", "CD62P"],
    size=0.1,
)

##### ==============================================================================

# # sc.pp.calculate_qc_metrics(adata_raw, qc_vars=["NegPrb"], inplace=True)
# # sc.pp.calculate_qc_metrics(adata_raw, inplace=True)
# sc.pp.calculate_qc_metrics(adata_raw, percent_top=None, inplace=True)
sc.pp.calculate_qc_metrics(adata_raw, inplace=True, percent_top=None)
pd.set_option("display.max_columns", None)
# adata.obs["total_counts_NegPrb"].sum() / adata.obs["total_counts"].sum() * 100

fig, axs = plt.subplots(1, 3, figsize=(15, 4))

axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)

axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)

# ------------

sc.pp.filter_cells(adata, min_counts=100)
sc.pp.filter_genes(adata, min_cells=400)

# ------------

fig, ax = plt.subplots(1, 2, figsize=(15, 15))
sq.gr.spatial_neighbors(
    adata_raw,
    n_neighs=10,
    coord_type="generic",
)
_, idx = adata_raw.obsp["spatial_connectivities"][420, :].nonzero()
idx = np.append(idx, 420)
sq.pl.spatial_scatter(
    adata_raw[idx, :],
    library_id="CD45RA",
    color="cell_type",
    connectivity_key="spatial_connectivities",
    size=10,
    edges_width=1,
    edges_color="black",
    img=False,
    title="K-nearest neighbors",
    ax=ax[0],
)

sq.gr.spatial_neighbors(
    adata_raw,
    n_neighs=10,
    coord_type="generic",
    delaunay=True,
)
_, idx = adata_raw.obsp["spatial_connectivities"][420, :].nonzero()
idx = np.append(idx, 420)
sq.pl.spatial_scatter(
    adata_raw[idx, :],
    library_id="CD44",
    color="cell_type",
    connectivity_key="spatial_connectivities",
    size=10,
    edges_width=1,
    edges_color="black",
    img=False,
    title="Delaunay triangulation",
    ax=ax[1],
)

adata_raw.obsp["spatial_connectivities"]

adata_raw.obsp["spatial_distances"]
#---------
sq.gr.spatial_neighbors(
    adata_raw,
    radius=30,
    coord_type="generic",
)
_, idx = adata_raw.obsp["spatial_connectivities"][420, :].nonzero()
idx = np.append(idx, 420)
sq.pl.spatial_scatter(
    adata_raw[idx, :],
    library_id="CD45RA",
    color="cell_type",
    connectivity_key="spatial_connectivities",
    size=10,
    edges_width=1,
    edges_color="black",
    img=False,
)

#------------

adata_raw.obsp["spatial_connectivities"]
adata_raw.obsp["spatial_distances"]

# ========================================================

adata_spatial_neighbor = sq.gr.spatial_neighbors(
    adata_raw, coord_type="generic", delaunay=True
)

sq.gr.centrality_scores(adata_raw, cluster_key="cell_type")

# sq.pl.centrality_scores(adata_raw, cluster_key="cell_type", figsize=(10, 6))
sq.pl.centrality_scores(adata_raw, cluster_key="cell_type", figsize=(16, 5))

# ========================================================

# sc.pp.filter_cells(adata, min_counts=10)

# adata.layers["counts"] = adata.X.copy()
# sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000)
# sc.pp.normalize_total(adata, inplace=True)
# sc.pp.log1p(adata)
# sc.pp.pca(adata)
# sc.pp.neighbors(adata)
# sc.tl.umap(adata)
# sc.tl.leiden(adata)

sc.pl.umap(
    adata_raw,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "cell_type",
    ],
    wspace=0.4,
)

# -------

sq.pl.spatial_scatter(
    adata_raw,
    size=13,
    shape=None,
    color=[
        "cell_type",
    ],
    wspace=0.4,
)

# ================================

adata_subsample = sc.pp.subsample(adata_raw, fraction=0.5, copy=True)

sq.gr.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
)
sq.pl.co_occurrence(
    adata_subsample,
    cluster_key="cell_type",
    clusters=["Skin", "Soft tissue"],
    figsize=(10, 10),
)
sq.pl.spatial_scatter(
    adata_subsample,
    color="cell_type",
    shape=None,
    size=2,
)

# ================================

sq.gr.nhood_enrichment(adata_raw, cluster_key="cell_type")

fig, ax = plt.subplots(1, 2, figsize=(13, 7))
sq.pl.nhood_enrichment(
    adata_raw,
    size=10,
    cluster_key="cell_type",
    figsize=(8, 8),
    title="Neighborhood enrichment adata",
    ax=ax[0],
)
sq.pl.spatial_scatter(adata_subsample, color="cell_type", shape=None, size=2, ax=ax[1])

# ===============================

fig, ax = plt.subplots(1, 2, figsize=(15, 7))
mode = "L"

sq.gr.ripley(adata_raw, cluster_key="cell_type", mode=mode)
sq.pl.ripley(adata_raw, cluster_key="cell_type", mode=mode, ax=ax[0])

sq.pl.spatial_scatter(
    adata_subsample,
    color="cell_type",
    groups=["Skin", "Testis", "Uterus"],
    shape=None,
    size=10,
    ax=ax[1],
)

# ===============================

sq.gr.spatial_neighbors(adata_subsample, coord_type="generic", delaunay=True)
sq.gr.spatial_autocorr(
    adata_subsample,
    mode="moran",
    n_perms=100,
    n_jobs=1,
)
adata_subsample.uns["moranI"].head(10)

sq.pl.spatial_scatter(
    adata_subsample,
    color=[
        "CD45RA",
        "CD4",
    ],
    shape=None,
    size=10,
    img=False,
)

# =================================

sq.gr.spatial_neighbors(adata_raw, coord_type="generic", spatial_key="spatial")
sq.gr.nhood_enrichment(adata_raw, cluster_key="cell_type")
sq.pl.nhood_enrichment(
    adata_raw,
    cluster_key="cell_type",
    method="average",
    cmap="inferno",
    figsize=(5, 5),
)

# =================================

_counts = adata_raw.obs["cell_type"].value_counts()
_counts.name = "cell_counts"
meta_leiden = pd.DataFrame(_counts)

sq.gr.centrality_scores(adata_raw, "cell_type")
sc.set_figure_params(figsize=(20, 8))

# copy centrality data to new DataFrame
df_central = deepcopy(adata_raw.uns["cell_type_centrality_scores"])
df_central.index = meta_leiden.index.tolist()

# sort clusters based on centrality scores
################################################
# closeness centrality - measure of how close the group is to other nodes.
ser_closeness = df_central["closeness_centrality"].sort_values(ascending=False)

# degree centrality - fraction of non-group members connected to group members.
# [Networkx](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.degree_centrality.html#networkx.algorithms.centrality.degree_centrality)
# The degree centrality for a node v is the fraction of nodes it is connected to.
ser_degree = df_central["degree_centrality"].sort_values(ascending=False)

# clustering coefficient - measure of the degree to which nodes cluster together.
ser_cluster = df_central["average_clustering"].sort_values(ascending=False)


# High closeness scores:
inst_clusters = ser_closeness.index.tolist()[:5]
print(inst_clusters)
sq.pl.spatial_scatter(
    adata_raw, groups=inst_clusters, color="cell_type", size=15, img=False, figsize=(10, 10), library_id=['IL2RA']
)


# Low closeness scores:
inst_clusters = ser_closeness.index.tolist()[-5:]
print(inst_clusters)
sq.pl.spatial_scatter(
    adata_raw, groups=inst_clusters, color="cell_type", size=15, img=False, figsize=(10, 10), library_id=['IL2RA']
)


# =====================================================

# High degree centrality:

inst_clusters = ser_degree.index.tolist()[:5]
print(inst_clusters)
sq.pl.spatial_scatter(
    adata_raw, groups=inst_clusters, color="cell_type", size=15, img=False, figsize=(10, 10), library_id=['IL2RA']
)

# Low degree centrality:
inst_clusters = ser_degree.index.tolist()[-5:]
print(inst_clusters)
sq.pl.spatial_scatter(
    adata_raw, groups=inst_clusters, color="cell_type", size=15, img=False, figsize=(10, 10), library_id=['IL2RA']
)

# =====================================================

# High clustering coefficient:
inst_clusters = ser_cluster.index.tolist()[:5]
print(inst_clusters)
sq.pl.spatial_scatter(
    adata_raw, groups=inst_clusters, color="cell_type", size=15, img=False, figsize=(10, 10), library_id=['IL2RA']
)

# Low clustering coefficient:
inst_clusters = ser_cluster.index.tolist()[-5:]
print(inst_clusters)
sq.pl.spatial_scatter(
    adata_raw, groups=inst_clusters, color="cell_type", size=15, img=False, figsize=(10, 10), library_id=['IL2RA']
)

# =====================================================

# Autocorrelation: Moran’s I Score

sq.gr.spatial_autocorr(adata_raw, mode="moran")
num_view = 12
top_autocorr = (
    adata_raw.uns["moranI"]["I"].sort_values(ascending=False).head(num_view).index.tolist()
)
bot_autocorr = (
    adata_raw.uns["moranI"]["I"].sort_values(ascending=True).head(num_view).index.tolist()
)

# ---

# Genes with high spatial autocorrelation
sq.pl.spatial_scatter(
    adata_raw, color=top_autocorr, size=20, cmap="Reds", img=False, figsize=(5, 5), library_key="IL2RA"
)

# 













