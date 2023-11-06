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
import itertools
import random
from scipy.stats import mannwhitneyu
import os
import warnings
import gseapy as gp
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import mygene
import seaborn as sns
from gseapy import barplot, dotplot
import matplotlib.pyplot as plt
import decoupler as dc
from numpy.random import default_rng
from copy import deepcopy
import squidpy as sq
import schist as scs
from pingouin import mwu
import graph_tools as gt
import squidpy
from scipy.stats import mannwhitneyu
# -------------------------------------------------------------
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
# from neighborhood_enrichment import neighborhood_enrichment
from statsmodels.stats.multitest import fdrcorrection
from itertools import chain
from sklearn.decomposition import PCA
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.inspection import permutation_importance
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split


def dasquidpy_spatial_autocorrelation(adata_pickle_path, dependent_variable_name):
    # ---
    moranI_0={}
    moranI_1={}
    # ---
    gearyC_0={}
    gearyC_1={}
    # ---
    with open(adata_pickle_path, 'rb') as f:
        pickle_= pickle.load(f)
    # ------
    Dict_NhoodEnrichment_Zscore=[]
    Dict_NhoodEnrichment_Count=[]
    # ------
    Dict_NhoodEnrichment_Zscore_1=[]
    Dict_NhoodEnrichment_Zscore_0=[]
    Dict_NhoodEnrichment_Count_1=[]
    Dict_NhoodEnrichment_Count_0=[]
    # ------
    Cluster_Pairs_List=[]
    Cluster_Pairs_List_1=[]
    Cluster_Pairs_List_0=[]
    # ------
    pickle_items=list(pickle_.items())
    first_three_items = pickle_items[:]
    pickle_new=dict(first_three_items)
    pickle_=pickle_new
    # ------
    no_of_patients=np.shape(pickle_)
    # ---
    patient_ids=[]
    recurrence_labels=[]
    # ---
    cell_types_1=[]
    cell_types_0=[]

    cell_types_1_count={}
    cell_types_0_count={}

    cnt=-1
    for i in pickle_:
        cnt+=1
        one_or_zero_flag=-1
        adata=pickle_[i]
        # ---
        patient_ids.append(i)
        rec_lab=np.unique(pickle_[i].obsm[dependent_variable_name])[0]
        recurrence_labels.append(rec_lab)
        # ---
        clusters=list(np.unique(adata.obs['celltype']))
        clusters_sorted=sorted(clusters)
        # ---
        cluster_names_matrix=[]
        for j in clusters_sorted:
            list_=[]
            for k in clusters_sorted:
                index=str(j)+'(**append_name**)'+str(k)
                list_.append(index)
            cluster_names_matrix.append(list_)
        upper_cluster_matrix = np.triu(cluster_names_matrix, 1)
        # ---
        cluster_pairs_list=[]
        for j in upper_cluster_matrix:
            for k in j:
                # print(j, k)
                if k!='' and k!=None:
                    # cluster_pairs_list.append(k)
                    cluster_pairs_list.append((k.split("(**append_name**)",1)[0], k.split("(**append_name**)",1)[1]))
        
        # ---------------------*************************---------------------
        squidpy.gr.spatial_neighbors(adata)
        
        
        
        moranI=squidpy.gr.spatial_autocorr(adata, mode='moran', copy=True, backend="multiprocessing", show_progress_bar=False).to_dict()['I']
        gearyC=squidpy.gr.spatial_autocorr(adata, mode='geary', copy=True, backend="multiprocessing", show_progress_bar=False).to_dict()['C']
        
        
        
        # ---
        if rec_lab=='POSITIVE' or rec_lab=='positive' or rec_lab=='1' or rec_lab==1:
            ct1=list(adata.obs['celltype'])
            ct1_unique=list(np.unique(ct1))
            for key in moranI:
                if key in moranI_1.keys():
                    moranI_1[key].append(moranI[key])
                else:
                    moranI_1[key]=[moranI[key]]
            # ---
            for key in gearyC:
                if key in gearyC_1.keys():
                    gearyC_1[key].append(gearyC[key])
                else:
                    gearyC_1[key]=[gearyC[key]]
            one_or_zero_flag=1
        # ---
        elif rec_lab=='NEGATIVE' or rec_lab=='negative' or rec_lab=='0' or rec_lab==0:
            ct0=list(adata.obs['celltype'])
            ct0_unique=list(np.unique(ct0))
            for key in moranI:
                if key in moranI_0.keys():
                    moranI_0[key].append(moranI[key])
                else:
                    moranI_0[key]=[moranI[key]]
            # ---
            for key in gearyC:
                if key in gearyC_0.keys():
                    gearyC_0[key].append(gearyC[key])
                else:
                    gearyC_0[key]=[gearyC[key]]
            one_or_zero_flag=0
        # ---       
        
        # ---------------------*************************---------------------
        # ---------------------*************************---------------------
        # ---------------------*************************---------------------
        
        # ---
        if one_or_zero_flag==0:
            for p in ct0_unique:
                if p in cell_types_0:
                    cell_types_0_count[p]+=ct0.count(p)
                else:
                    cell_types_0_count[p]=ct0.count(p)
            for k in ct0_unique:
                cell_types_0.append(k)
            cell_types_0=list(np.unique(cell_types_0))
        elif one_or_zero_flag==1:
            for p in ct1_unique:
                if p in cell_types_1:
                    cell_types_1_count[p]+=ct1.count(p)
                else:
                    cell_types_1_count[p]=ct1.count(p)
            for k in ct1_unique:
                cell_types_1.append(k)
            cell_types_1=list(np.unique(cell_types_1))
        # ---

    # ========================================================================================

    celltypes_1_0=list(set(cell_types_1).difference(set(cell_types_0)))
    celltypes_0_1=list(set(cell_types_0).difference(set(cell_types_1)))
    celltypes_0_and_1=list(set(cell_types_0).intersection(set(cell_types_1)))

    if len(celltypes_1_0)==0:
        print('No. of cell types in positive class not in control class = 0')
    else:
        for i in celltypes_1_0:
            cell_types_0_count[i]=0

    if len(celltypes_0_1)==0:
        print('No. of cell types in control class not in positive class = 0')
    else:
        for i in celltypes_0_1:
            cell_types_1_count[i]=0

    cell_types_count_df=pd.DataFrame({'Control class (non-recurrent TNBC)':pd.Series(cell_types_0_count),'Positive class (recurrent TNBC)':pd.Series(cell_types_1_count)})
    cell_types_count_df_transposed = cell_types_count_df.T
    cell_types_count_df_transposed['condition'] = cell_types_count_df_transposed.index
    cell_types_count_df['cell type'] = cell_types_count_df.index

    # # ========================================================================================

    moranI_list=[]
    gearyC_list=[]

    moranI_proteins=list(moranI_0.keys())
    gearyC_proteins=list(gearyC.keys())

    # ---
    for i in moranI_proteins:
        U1, p = mannwhitneyu(moranI_0[i], moranI_1[i], method="exact")
        moranI_list.append(p)
        # ---
        U1, p = mannwhitneyu(gearyC_0[i], gearyC_1[i], method="exact")
        gearyC_list.append(p)
        # ---
    # ---
    moranI_dict=dict(zip(moranI_proteins, moranI_list))
    gearyC_dict=dict(zip(gearyC_proteins, gearyC_list))

    # -------------------------

    significant_moranI_dict = dict((k, v) for k, v in moranI_dict.items() if v < 0.05)
    significant_moranI_dict_sorted = dict(sorted(significant_moranI_dict.items(), key=lambda x:x[1]))
    # ---
    significant_gearyC_dict = dict((k, v) for k, v in gearyC_dict.items() if v < 0.05)
    significant_gearyC_dict_sorted = dict(sorted(significant_gearyC_dict.items(), key=lambda x:x[1]))

    # ========================================================================================

    if len(significant_moranI_dict_sorted)==0:
        print('No significant proteins found in terms of Moran I statistic!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Moran I statistic) – proteins with significant Moran I statistic differences (p-value<0.05)')
        feature = list(significant_moranI_dict_sorted.keys())
        score = list(significant_moranI_dict_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_moranI_pvalues_significantProteins_allPatients.pdf', format='pdf')
        plt.show()

    # ---

    if len(significant_gearyC_dict_sorted)==0:
        print('No significant proteins found in terms of Geary C statistic!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Geary C statistic) – proteins with significant Geary C statistic differences (p-value<0.05)')
        feature = list(significant_gearyC_dict_sorted.keys())
        score = list(significant_gearyC_dict_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_gearyC_pvalues_significantProteins_allPatients.pdf', format='pdf')
        plt.show()

    # ---

    # ========================================================================================