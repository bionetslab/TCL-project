

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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
import os
from skimage import io
from skimage.util import img_as_ubyte
import imageio as iio
from skimage import img_as_ubyte
import re
from matplotlib.pyplot import figure
import os
from sklearn.mixture import GaussianMixture
import statistics
from bigtree import Node 

# # ------------------

# # Download the following files from HPA and unzip them:
# # - https://www.proteinatlas.org/download/rna_single_cell_type_tissue.tsv.zip
# # - https://www.proteinatlas.org/download/rna_single_cell_cluster_description.tsv.zip

# # ------------------

# # Load the data into the right format
# expr = pd.read_csv('rna_single_cell_type_tissue.tsv',sep='\t')
# meta = pd.read_csv('rna_single_cell_cluster_description.tsv', sep='\t')
# expr = expr[expr.Tissue == 'skin']
# meta = meta[meta.Tissue == 'Skin']
# expr.reset_index(inplace=True, drop=True)
# meta.reset_index(inplace=True, drop=True)
# meta.index = meta.Cluster

# # ------------------

# # Compute cell nTMPs for cell types as weighted means over corresponding cluster

# cell_types = list(set(meta['Cell type']))
# genes = list(set(expr['Gene name']))
# expr_per_cell_type = pd.DataFrame(0.0, columns=genes, index=cell_types)
# count_per_cell_type = pd.DataFrame(0, columns=genes, index=cell_types)

# for i in range(expr.shape[0]):
#     gene = expr.loc[i, 'Gene name']
#     ntpm = expr.loc[i, 'nTPM']
#     cluster = expr.loc[i, 'Cluster']
#     cell_type = meta.loc[cluster, 'Cell type']
#     cell_count = meta.loc[cluster, 'Cell count']
#     expr_per_cell_type.loc[cell_type, gene] += ntpm * cell_count
#     count_per_cell_type.loc[cell_type, gene] += cell_count

# # Cell type nTPMs
# expr_per_cell_type /= count_per_cell_type

# # Cell type nTPMs normalized by gene-wise maximum expression
# expr_per_cell_type_max_norm = expr_per_cell_type / expr_per_cell_type.max()

# # ------------------

# expr_per_cell_type.to_csv('cell_type_nTPM.csv')
# expr_per_cell_type_max_norm.to_csv('cell_type_nTPM_max_norm.csv')

# ------------------


expr_per_cell_type=pd.read_csv('cell_type_nTPM.csv')
expr_per_cell_type_max_norm=pd.read_csv('cell_type_nTPM_max_norm.csv')
expr_per_cell_type_max_norm=expr_per_cell_type_max_norm.set_index('Unnamed: 0')

def plot_heatmap(data, genes, figsize=(9,6), outfile=None, title=None):
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(data=data[genes],ax=ax)
    if not title is None:
        ax.set_title(title)
    fig.tight_layout()
    if not outfile is None:
        fig.savefig(outfile)


# ------------------


plot_heatmap(expr_per_cell_type_max_norm, ['CD4','CD8A'], outfile='heatmap.pdf', title='Gene-wise max-normalized nTPM per cell type')

# plot_heatmap(expr_per_cell_type_max_norm, ['ADAM10','LAMP1','ITGAL','ITGAX','IL3RA','SDC1','CD14','FCGR3A','CD163','CD19','CD1A','MS4A1','LY75','MRC1','CD209','CD24','IL2RA','NGFR','ITGB1','CD3', 'CD36','CD38','CD4','CD40','CD44','CD52',
# ], outfile='heatmap.pdf', title='Gene-wise max-normalized nTPM per cell type')

# plot_heatmap(expr_per_cell_type_max_norm, 
#              ['ADAM10', 'LAMP1', 'ITGAL', 'ITGAX', 'IL3RA', 'SDC1', 'CD14',
#               'FCGR3A', 'CD163', 'CD19', 'CD1A', 'MS4A1', 'LY75', 'MRC1',
#               'CD209', 'CD24', 'IL2RA', 'NGFR', 'ITGB1', 'CD3D', 'CD36', 'CD38',
#               'CD4', 'CD40', 'CD44', 'CD52', 'ICAM1', 'CD55', 'NCAM1', 'CD6',
#               'SELP', 'CD63', 'CD68', 'CD69', 'TFRC', 'CD8A', 'CD9', 'FAS', 'COL4A1',
#               'KRT14', 'CDH1', 'EGFR', 'HLA-A', 'HLA-DRA', 'CDKN1B', 'MLANA', 'NES', 
#               'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'PPARG', 'DICER1', 'ADAM17', 
#               'TP73', 'TNFRSF1A', 'TNFRSF1B', 'CTNNB1', 'TP63', 'GJA1', 
#               'BCL2', 'MME', 'ITGAM', 'ANPEP', 'THBD', 'FUT4', 'CD2', 'CD27', 'CD276', 
#               'SPN', 'CD5', 'CD53', 'CD7', 'CD81', 'CD83', 'THY1', 'MCAM', 'CD34', 'ENTPD1', 
#               'CD58', 'CEACAM1', 'NT5E', 'CD80', 'ACTB', 'CALCA', 'CSNK2A1', 
#               'CSNK2A2', 'CASP3', 'CCND1', 'CYCS', 'KRT15', 'DKK3', 'DES', 'FCER1A', 'FST', 
#               'NHEJ1', 'ITGAE', 'GFER', 'CSPG4', 'MECP2', 'MYH9', 'KLRC1', 'NRG1', 
#               'ELANE', 'ORC2', 'ORC3', 'PCNA', 'GRN', 'PRPF8', 'ETS1', 'RIMS3', 
#               'RELN', 'STAT3', 'SYT10', 'TARDBP', 'TUBB', 'TP53', 
#               'CLEC4C', 'PECAM1', 'CD46', 'SELL', 'FOXP3', 'KIT', 
#               'TUBB3', 'TSLP',
#               # -------------
#               'PTPRC' # ['CD45-PE', 'CD45RA-PE', 'CD45RO-FITC']
#               # -------------
#               ],
#              outfile='heatmap.pdf', 
#              title='Gene-wise max-normalized nTPM per cell type')


# 'ADAM10', 'LAMP1', 'ITGAL', 'ITGAX', 'IL3RA', 'SDC1', 'CD14', 
# 'FCGR3A', 'CD163', 'CD19', 'CD1A', 'MS4A1', 'LY75', 'MRC1', 
# 'CD209', 'CD24', 'IL2RA', 'NGFR', 'ITGB1', 'CD3D', 'CD36', 'CD38', 
# 'CD4', 'CD40', 'CD44', 'CD52', 'ICAM1', 'CD55', 'NCAM1', 'CD6', 
# 'SELP', 'CD63', 'CD68', 'CD69', 'TFRC', 'CD8A', 'CD9', 'FAS', 'COL4A1', 
# 'KRT14', 'CDH1', 'EGFR', 'HLA-A', 'HLA-DRA', 'CDKN1B', 'MLANA', 'NES', 
# 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'PPARG', 'DICER1', 'ADAM17', 
# 'TP73', 'TNFRSF1A', 'TNFRSF1B', 'CTNNB1', 'TP63', 'GJA1', 
# 'BCL2', 'MME', 'ITGAM', 'ANPEP', 'THBD', 'FUT4', 'CD2', 'CD27', 'CD276', 
# 'SPN', 'CD5', 'CD53', 'CD7', 'CD81', 'CD83', 'THY1', 'MCAM', 'CD34', 'ENTPD1', 
# 'CD58', 'CEACAM1', 'NT5E', 'CD80', 'ACTB', 'CALCA', 'CSNK2A1', 
# 'CSNK2A2', 'CASP3', 'CCND1', 'CYCS', 'KRT15', 'DKK3', 'DES', 'FCER1A', 'FST', 
# 'NHEJ1', 'ITGAE', 'GFER', 'CSPG4', 'MECP2', 'MYH9', 'KLRC1', 'NRG1', 
# 'ELANE', 'ORC2', 'ORC3', 'PCNA', 'GRN', 'PRPF8', 'ETS1', 'RIMS3', 
# 'RELN', 'STAT3', 'SYT10', 'TARDBP', 'TUBB', 'TP53', 
# 'CLEC4C', 'PECAM1', 'CD46', 'SELL', 'FOXP3', 'KIT', 
# 'TUBB3', 'TSLP'


plot_heatmap(expr_per_cell_type_max_norm, 
             ['ADAM10','ITGAL','ITGAX','CD14','CD163','CD24',
              'IL2RA','ITGB1','CD3D','CD36','CD38','CD4','CD40',
              'CD44','PTPRC','ICAM1','CD55','NCAM1','CD6','SELP',
              'CD63','CD68','CD69','CD8A','CD9','FAS','KRT14',
              'HLA-A','HLA-DRA','NOTCH1','NOTCH3','VIM','ACTB',
              'TSPAN7','CALCA','CEBPB','CEBPE','CSNK2A1','CSNK2A2',
              'CASP3','CCND1','CYCS','KRT15','DKK3','DROSHA','DES',
              'EBF1','KRT8','ERG','FCER1A','FLT3','HLA-DRB1','FST',
              'MKI67','MKI67','NHEJ1','ITGAE','GFER','CSPG4','MECP2',
              'MYH9','KLRC1','NRG1','ELANE','ORC2','ORC3','PCNA',
              'GRN','PRPF8','ETS1','RIMS3','RELN','STAT3','SYT10',
              'TARDBP','TUBB'],
             outfile='heatmap.pdf', 
             title='Gene-wise max-normalized nTPM per cell type')


# -------------------

# 'CD45-PE': 'PTPRC', ###
#     'CD45RA-PE': 'PTPRC (2)', ###
#     'CD45RO-FITC': 'PTPRC (3)', ###
#     'L302-FITC': 'L302', # Nothing found for L302, is it CD302?
#     'APPC-FITC': 'CD172a', # Nothing found for APPC!
#     'BOP-FITC': 'TSPAN7', # Nothing found for APPC BOP!
#     'C97-FITC': 'C97', # "Family C97 unassigned peptidases"?
#     'CEB2-FITC': 'CEBPB', # This does not look correct!
#     'CEB6-FITC': 'CEBPE', # This does not look correct!
#     'CPB3-FITC': 'CPB3', # CPB3 not found!
#     'DRO-FITC': 'DROSHA', # Not found! (rice gene: https://www.google.com/search?q=DRO+gene&sca_esv=590391945&ei=O2V5ZbLfDsaE9u8Ph9uc-A4&ved=0ahUKEwjy6KzH-YuDAxVGgv0HHYctB-8Q4dUDCBE&uact=5&oq=DRO+gene&gs_lp=Egxnd3Mtd2l6LXNlcnAiCERSTyBnZW5lMgcQABiABBgNMgcQABiABBgNMgkQABiABBgNGAoyBxAAGIAEGA0yBxAAGIAEGA0yBxAAGIAEGA0yBhAAGB4YDTIGEAAYHhgNMgYQABgeGA0yCBAAGB4YDRgPSO4QUOIDWIINcAF4AZABAJgBXaABoAOqAQE3uAEDyAEA-AEBwgIKEAAYRxjWBBiwA8ICBRAhGKABwgIHECEYoAEYCsICBhAAGBYYHsICCxAAGIAEGIoFGIYDwgIIEAAYCBgeGA3iAwQYACBBiAYBkAYI&sclient=gws-wiz-serp, https://www.google.com/search?sca_esv=590391945&q=DRO1+gene&sa=X&ved=2ahUKEwjK386v-YuDAxU2hP0HHYj_Br4Q7xYoAHoECAYQAg&cshid=1702454591935042&biw=2560&bih=1298&dpr=1)
#     'EBF-P-FITC': 'EBF1', # EBF not found!
#     'EK2-FITC': 'KRT8', # EK2 not found!
#     'ERG1-FITC': 'ERG', # ERG1 not found!
#     'FRP 13E12-FITC': '?', # This is not correct!
#     'KIS-FITC': '? (2)', # This is not correct!
#     'PBS': 'Phosphate buffered salt solution', # What to do with this?
#     'CLA-FITC': 'CLA', # Not correct!
#     'IgA-FITC': 'IGHG', # Not correct!
#     'TcR alpha': 'TRAC', # T-cell receptor (correct!)
#     'Follicular Dendritic Cells-FITC': 'CD21', # 'Mature B-cell' and 'follicular dendritic cell' marker
#     'FLT-FITC': 'FLT3', # FLT not found!
    

# expr_per_cell_type_max_norm=expr_per_cell_type_max_norm[[
#  'ADAM10','LAMP1','ITGAL','ITGAX','IL3RA','SDC1','CD14',
#  'FCGR3A','CD163','CD19','CD1A','MS4A1','LY75','MRC1',
#  'CD209','CD24','IL2RA','NGFR','ITGB1','CD3D','CD36',
#  'CD38','CD4','CD40','CD44','CD52','ICAM1','CD55','NCAM1',
#  'CD6','SELP','CD63','CD68','CD69','TFRC','CD8A','CD9',
#  'FAS','COL4A1','KRT14','CDH1','EGFR','HLA-A','HLA-DRA',
#  'CDKN1B','MLANA','NES','NOTCH1','NOTCH2','NOTCH3','NOTCH4',
#  'PPARG','DICER1','ADAM17','TP73','TNFRSF1A','TNFRSF1B',
#  'CTNNB1','TP63','GJA1','BCL2','MME','ITGAM','ANPEP','THBD',
#  'FUT4','CD2','CD27','CD276','SPN','CD5','CD53','CD7','CD81',
#  'CD83','THY1','MCAM','CD34','ENTPD1','CD58','CEACAM1','NT5E',
#  'CD80','ACTB','CALCA','CSNK2A1','CSNK2A2','CASP3','CCND1',
#  'CYCS','KRT15','DKK3','DES','FCER1A','FST','NHEJ1','ITGAE',
#  'GFER','CSPG4','MECP2','MYH9','KLRC1','NRG1','ELANE','ORC2',
#  'ORC3','PCNA','GRN','PRPF8','ETS1','RIMS3','RELN','STAT3',
#  'SYT10','TARDBP','TUBB','TP53','CLEC4C','PECAM1','CD46',
#  'SELL','FOXP3','KIT','TUBB3','TSLP','PTPRC'
#  ]]

# expr_per_cell_type_max_norm=expr_per_cell_type_max_norm[[
#   'ADAM10','ITGAL','ITGAX','CD14','CD163','CD24',
#   'IL2RA','ITGB1','CD3D','CD36','CD38','CD4','CD40',
#   'CD44','PTPRC','ICAM1','CD55','NCAM1','CD6','SELP',
#   'CD63','CD68','CD69','CD8A','CD9','FAS','KRT14',
#   'HLA-A','HLA-DRA','NOTCH1','NOTCH3','VIM','ACTB',
#   'TSPAN7','CALCA','CEBPB','CEBPE','CSNK2A1','CSNK2A2',
#   'CASP3','CCND1','CYCS','KRT15','DKK3','DROSHA','DES',
#   'EBF1','KRT8','ERG','FCER1A','FLT3','HLA-DRB1','FST',
#   'MKI67','MKI67','NHEJ1','ITGAE','GFER','CSPG4','MECP2',
#   'MYH9','KLRC1','NRG1','ELANE','ORC2','ORC3','PCNA',
#   'GRN','PRPF8','ETS1','RIMS3','RELN','STAT3','SYT10',
#   'TARDBP','TUBB'
#   ]]

expr_per_cell_type_max_norm=expr_per_cell_type_max_norm[[
  'ADAM10','LAMP1','ITGAL','ITGAX','IL3RA','SDC1','CD14',
  'FCGR3A','CD163','CD19','CD1A','MS4A1','LY75','MRC1',
  'CD209','CD24','IL2RA','NGFR','ITGB1','CD3D','CD36',
  'CD38','CD4','CD40','CD44','CD52','ICAM1','CD55','NCAM1',
  'CD6','SELP','CD63','CD68','CD69','TFRC','CD8A','CD9',
  'FAS','COL4A1','KRT14','CDH1','EGFR','HLA-A','HLA-DRA',
  'CDKN1B','MLANA','NES','NOTCH1','NOTCH2','NOTCH3','NOTCH4',
  'PPARG','DICER1','ADAM17','TP73','TNFRSF1A','TNFRSF1B',
  'CTNNB1','TP63','GJA1','BCL2','MME','ITGAM','ANPEP','THBD',
  'FUT4','CD2','CD27','CD276','SPN','CD5','CD53','CD7','CD81',
  'CD83','THY1','MCAM','CD34','ENTPD1','CD58','CEACAM1','NT5E',
  'CD80','ACTB','CALCA','CSNK2A1','CSNK2A2','CASP3','CCND1',
  'CYCS','KRT15','DKK3','DES','FCER1A','FST','NHEJ1','ITGAE',
  'GFER','CSPG4','MECP2','MYH9','KLRC1','NRG1','ELANE','ORC2',
  'ORC3','PCNA','GRN','PRPF8','ETS1','RIMS3','RELN','STAT3',
  'SYT10','TARDBP','TUBB','TP53','CLEC4C','PECAM1','CD46',
  'SELL','FOXP3','KIT','TUBB3','TSLP','PTPRC'
  ]]

# -----

clustering_tree={}
clustering_results={}

# -----
with open('/Users/surya/Documents/GITHUB-REPOSITORIES/spatial_proteomics/2. feature_extraction/anndata_allConditions.pkl', 'rb') as f:
   anndata_allConditions = pickle.load(f)
all_columns_list=list(anndata_allConditions.to_df().columns)
cell_coords_df=pd.DataFrame(list(anndata_allConditions.obsm['_cell_coords_']), columns=['x', 'y'])
# -----
anndata_allConditions_df=anndata_allConditions.to_df()
errorStatement="No_errors_all_celltypes_assigned"
# COUNT=-1
while len(list(expr_per_cell_type_max_norm.index))>0:
# while COUNT<0:
    if len(anndata_allConditions_df)==0:
        errorStatement="No_more_cells_left_to_assign_in_database"
        break
    # COUNT+=1
    anndata_allConditions_colNames=list(anndata_allConditions_df.columns)
    anndata_allConditions_colValues=[]
    for i in anndata_allConditions_colNames:
        anndata_allConditions_colValues.append(list(anndata_allConditions_df[i]))
    anndata_allConditions_dict=dict(zip(anndata_allConditions_colNames, anndata_allConditions_colValues))
    # -----
    anndata_allConditions_gmm=[]
    anndata_allConditions_gmm_labels=[]
    anndata_allConditions_gmm_upper_threshold=[]
    anndata_allConditions_gmm_lower_threshold=[]
    # -----
    good_split_genes=[]
    bad_split_genes=[]
    # -----
    for i in anndata_allConditions_dict:
        X_list=list(anndata_allConditions_dict[i])
        X=np.array(X_list)
        X=X.reshape(-1, 1)
        gm = GaussianMixture(n_components=2, random_state=0).fit(X)
        anndata_allConditions_gmm.append(gm)
        labels = gm.predict(X)
        anndata_allConditions_gmm_labels.append(labels)
        X_0=[]
        X_1=[]
        cnt=-1
        for j in labels:
            cnt+=1
            if j==0:
                X_0.append(X_list[cnt])
                X_1.append(X_list[cnt])
        X=[X_0, X_1]
        means_=gm.means_
        stddevs_0= statistics.pstdev(X_0)
        stddevs_1= statistics.pstdev(X_1)
        stddevs_=[stddevs_0, stddevs_1]
        # -----
        upper_threshold=means_[1]-1.96*stddevs_[1]
        anndata_allConditions_gmm_upper_threshold.append(upper_threshold)
        lower_threshold=means_[0]+1.96*stddevs_[0]
        anndata_allConditions_gmm_lower_threshold.append(lower_threshold)
        # -----
        if upper_threshold>lower_threshold:
            good_split_genes.append(i)
        else:
            bad_split_genes.append(i)
    
    geneName_lowerThreshold_dict=dict(zip(anndata_allConditions_colNames, anndata_allConditions_gmm_lower_threshold))
    geneName_upperThreshold_dict=dict(zip(anndata_allConditions_colNames, anndata_allConditions_gmm_upper_threshold))
    
    # ==================
    # ==================
    # ==================
    
    max_celltype_complement=[]
    gene_names=list(expr_per_cell_type_max_norm.columns)
    celltypes=list(expr_per_cell_type_max_norm.index)
    celltypes_expressions=[]
    LIST_COMP=[]
    max_impGeneValues=[]
    _count_=-1
    for i in expr_per_cell_type_max_norm.index:
        _count_+=1
        print('Count', _count_)
        _D_=expr_per_cell_type_max_norm.loc[[i]]
        celltypes_expressions.append(_D_)
        update_df = expr_per_cell_type_max_norm.drop(i)
        D_=pd.DataFrame(update_df.max()).T
        D_ = D_.rename(index={0: i})
        SUB=_D_.sub(D_)
        SUB_max=SUB.max(axis=1)
        max_impGeneValues.append(float(SUB_max))
        row=SUB.index[0]
        value=float(SUB_max)
        
        # # ---
        # list_comp = [c for c in SUB.columns if SUB[c][row] == value]
        # if len(list_comp)==1:
        #     print(list_comp[0])
        #     LIST_COMP.append(list_comp[0])
        # else:
        #     for l in list_comp:
        #         print(l)
        #         LIST_COMP.append(l)
        # # ---
        
        # # # ---
        # list_comp = [c for c in SUB.columns if SUB[c][row] == value][0]
        # print(list_comp)
        # LIST_COMP.append(list_comp)
        # # # ---
        
        
        
        # list_comp = [c for c in SUB.columns if SUB[c][row] == value]
        list_comp=list(SUB.idxmax(axis=1))
        print(list_comp)
        LIST_COMP.append(list_comp)
    
    
    # ---
    # celltypes_impGenes=dict(zip(celltypes, LIST_COMP))
    # celltypes_maxImpGeneValues=dict(zip(celltypes, max_impGeneValues))
    # ---
    
    _celltypes_impGenes_=dict(zip(celltypes, LIST_COMP))
    _celltypes_maxImpGeneValues_=dict(zip(celltypes, max_impGeneValues))
    
    _celltypes_=[]
    _impGenes_=[]
    _maxImpGeneValues_=[]
    
    for l in _celltypes_impGenes_:
        for l_ in _celltypes_impGenes_[l]:
            _celltypes_.append(l)
            _impGenes_.append(l_)
            _maxImpGeneValues_.append(_celltypes_maxImpGeneValues_[l])
    
    
    celltypes_good=[]
    impGenes_good=[]
    maxImpGeneValues_good=[]
    
    celltypes_bad=[]
    impGenes_bad=[]
    maxImpGeneValues_bad=[]
    
    celltypes_counter=-1
    for i in _celltypes_:
        celltypes_counter+=1
        if _impGenes_[celltypes_counter] in good_split_genes:
            celltypes_good.append(i)
            impGenes_good.append(_impGenes_[celltypes_counter])
            maxImpGeneValues_good.append(_maxImpGeneValues_[celltypes_counter])
        else:
            celltypes_bad.append(i)
            impGenes_bad.append(_impGenes_[celltypes_counter])
            maxImpGeneValues_bad.append(_maxImpGeneValues_[celltypes_counter])
    
    # # ---
    # for i in celltypes_impGenes:
    #     if celltypes_impGenes[i] in good_split_genes:
    #         celltypes_good.append(i)
    #         impGenes_good.append(celltypes_impGenes[i])
    #         maxImpGeneValues_good.append(celltypes_maxImpGeneValues[i])
    #     else:
    #         celltypes_bad.append(i)
    #         impGenes_bad.append(celltypes_impGenes[i])
    #         maxImpGeneValues_bad.append(celltypes_maxImpGeneValues[i])
    # # ---
    
    # ---
    if len(celltypes_good)==0:
        errorStatement="No_good_split_genes_found"
        break
    # ---
    _pos_=pd.Series(maxImpGeneValues_good).idxmax()
    print("Maximum Index position from the list: ", _pos_)  
    anndata_clustered=anndata_allConditions_df[anndata_allConditions_df[impGenes_good[_pos_]]>geneName_upperThreshold_dict[impGenes_good[_pos_]][0]]
    anndata_clustered_indices=list(anndata_clustered.index)
    # ---
    clustering_tree[impGenes_good[_pos_]]=celltypes_good[_pos_]
    clustering_results[celltypes_good[_pos_]]=anndata_clustered_indices
    # ---
    anndata_allConditions_df=anndata_allConditions_df.drop(anndata_clustered_indices)
    expr_per_cell_type_max_norm=expr_per_cell_type_max_norm.drop(celltypes_good[_pos_])
    # ---

# --------------------
# plotting_tree={}
# last_key='HPA-based celltype annotation'
# plotting_tree_list=[last_key]

# plotting_tree_count=-1
# for i in clustering_tree:
#     plotting_tree_count+=1
#     # if plotting_tree_count==0:
#     #     asdf
#     # final_plotting_tree_list=[]
#     plot_tree_dict_path='plotting_tree'
#     for j in plotting_tree_list:
#         plot_tree_dict_path=plot_tree_dict_path+'[' + '"' + j + '"' + ']'
#         # final_plotting_tree_list.append(j)
    
    
#     last_key=i
#     last_key_pos=i+'+'
#     last_key_neg=i+'-'
    
#     exec(plot_tree_dict_path+'[' + '"' + last_key_pos + '"' +']'+'='+ '"' + clustering_tree[i] + '"'  )
#     plotting_tree_list.append(last_key_neg)
# --------------------   




# --------------------
plotting_tree={}
last_key='HPA-based celltype annotation (tree):'
# plotting_tree_list=[last_key]
plotting_tree_list=[]
# ---
plotting_tree_count=-1
big_tree_count=0
exec('Node'+str(big_tree_count)+'='+'Node("' + last_key + '")')
for i in clustering_tree:
    plotting_tree_count+=1
    # if plotting_tree_count==0:
    #     asdf
    # final_plotting_tree_list=[]
    plot_tree_dict_path='plotting_tree'
    if plotting_tree_count>0:
        for j in plotting_tree_list:
            plot_tree_dict_path=plot_tree_dict_path+'[' + '"' + j + '"' + ']'
            # final_plotting_tree_list.append(j)
    else:
        pass
    # ---
    last_key=i
    last_key_pos=i+'+'
    last_key_neg=i+'-'
    # ---
    exec(plot_tree_dict_path+ '=' + '{}'  )
    exec(plot_tree_dict_path+'[' + '"' + last_key_pos + '"' +']'+'='+ '"' + clustering_tree[i] + '"'  )
    big_tree_count+=1
    exec('Node'+str(big_tree_count)+'='+'Node("' + last_key_pos + ':' + clustering_tree[i] + '"' + ',' +  'parent=' + 'Node'+str(big_tree_count-1)  + ')')
    # ---
    plotting_tree_list.append(last_key_neg)
    big_tree_count+=1
    exec('Node'+str(big_tree_count)+'='+'Node("' + last_key_neg + '"' + ',' +  'parent=' + 'Node'+str(big_tree_count-2)  + ')')

Node0.show()
# ---
print('\nHPA-based celltype annotation (dict):')
print(plotting_tree)
# --------------------




# a = Node("a", age=90)
# b = Node("b", age=65, parent=a)
# c = Node("c", age=60, parent=a)

# d = Node("x", age=90, parent=b)
# e = Node("y", age=65, parent=c)
# f = Node("z", age=60, parent=c)

# a.children
# # (Node(/a/b, age=65), Node(/a/c, age=60))

# a.depth, b.depth
# # (1, 2)

# a.max_depth
# # 2

# a.show(attr_list=["age"])
# a.show(attr_list=["age"])









