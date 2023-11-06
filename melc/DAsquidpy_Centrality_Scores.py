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


def dasquidpy_centrality_scores(adata_pickle_path, dependent_variable_name):
    # ---
    degree_centrality_0={}
    degree_centrality_1={}
    # ---
    average_clustering_0={}
    average_clustering_1={}
    # ---
    closeness_centrality_0={}
    closeness_centrality_1={}
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
        nhood_enrichment=squidpy.gr.centrality_scores(adata, cluster_key='celltype', copy=True, backend="multiprocessing", show_progress_bar=False)
        # squidpy.gr.centrality_scores(adata, cluster_key='celltype', backend="multiprocessing", show_progress_bar=False)
        # sq.pl.centrality_scores(adata, cluster_key="celltype")
        nh_enrichment=nhood_enrichment.to_dict()
        # ---
        dc=nh_enrichment['degree_centrality']
        ac=nh_enrichment['average_clustering']
        cc=nh_enrichment['closeness_centrality']
        # ---
        if rec_lab=='POSITIVE' or rec_lab=='positive' or rec_lab=='1' or rec_lab==1:
            ct1=list(adata.obs['celltype'])
            ct1_unique=list(np.unique(ct1))
            for key in dc:
                if key in degree_centrality_1.keys():
                    degree_centrality_1[key].append(dc[key])
                    average_clustering_1[key].append(ac[key])
                    closeness_centrality_1[key].append(cc[key])
                else:
                    degree_centrality_1[key]=[dc[key]]
                    average_clustering_1[key]=[ac[key]]
                    closeness_centrality_1[key]=[cc[key]]
            one_or_zero_flag=1
        # ---
        elif rec_lab=='NEGATIVE' or rec_lab=='negative' or rec_lab=='0' or rec_lab==0:
            ct0=list(adata.obs['celltype'])
            ct0_unique=list(np.unique(ct0))
            for key in dc:
                if key in degree_centrality_0.keys():
                    degree_centrality_0[key].append(dc[key])
                    average_clustering_0[key].append(ac[key])
                    closeness_centrality_0[key].append(cc[key])
                else:
                    degree_centrality_0[key]=[dc[key]]
                    average_clustering_0[key]=[ac[key]]
                    closeness_centrality_0[key]=[cc[key]]
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

    degree_centrality_list=[]
    average_clustering_list=[]
    closeness_centrality_list=[]
    # ---
    for i in celltypes_0_and_1:
        U1, p = mannwhitneyu(degree_centrality_0[i], degree_centrality_1[i], method="exact")
        degree_centrality_list.append(p)
        # ---
        U1, p = mannwhitneyu(average_clustering_0[i], average_clustering_1[i], method="exact")
        average_clustering_list.append(p)
        # ---
        U1, p = mannwhitneyu(closeness_centrality_0[i], closeness_centrality_1[i], method="exact")
        closeness_centrality_list.append(p)
    # ---
    degree_centrality_dict=dict(zip(celltypes_0_and_1, degree_centrality_list))
    average_clustering_dict=dict(zip(celltypes_0_and_1, average_clustering_list))
    closeness_centrality_dict=dict(zip(celltypes_0_and_1, closeness_centrality_list))

    # -------------------------

    significant_degree_centrality_dict = dict((k, v) for k, v in degree_centrality_dict.items() if v < 0.05)
    significant_degree_centrality_dict_sorted = dict(sorted(significant_degree_centrality_dict.items(), key=lambda x:x[1]))
    # ---
    significant_average_clustering_dict = dict((k, v) for k, v in average_clustering_dict.items() if v < 0.05)
    significant_average_clustering_dict_sorted = dict(sorted(significant_average_clustering_dict.items(), key=lambda x:x[1]))
    # ---
    significant_closeness_centrality_dict = dict((k, v) for k, v in closeness_centrality_dict.items() if v < 0.05)
    significant_closeness_centrality_dict_sorted = dict(sorted(significant_closeness_centrality_dict.items(), key=lambda x:x[1]))

    # ========================================================================================

    if len(significant_degree_centrality_dict_sorted)==0:
        print('No significant celltypes found in terms of degree centrality!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Degree Centrality) – celltypes with significant degree centrality differences (p-value<0.05)')
        feature = list(significant_degree_centrality_dict_sorted.keys())
        score = list(significant_degree_centrality_dict_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_degreeCentrality_pvalues_significantCelltypes_allPatients.pdf', format='pdf')
        plt.show()

    # ---

    if len(significant_average_clustering_dict_sorted)==0:
        print('No significant celltypes found in terms of average clustering!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Average Clustering) – celltypes with significant average clustering differences (p-value<0.05)')
        feature = list(significant_average_clustering_dict_sorted.keys())
        score = list(significant_average_clustering_dict_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_averageClustering_pvalues_significantCelltypes_allPatients.pdf', format='pdf')
        plt.show()

    # ---

    if len(significant_closeness_centrality_dict_sorted)==0:
        print('No significant celltypes found in terms of closeness centrality!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Closeness Centrality) – celltypes with significant closeness centrality differences (p-value<0.05)')
        feature = list(significant_closeness_centrality_dict_sorted.keys())
        score = list(significant_closeness_centrality_dict_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_closenessCentrality_pvalues_significantCelltypes_allPatients.pdf', format='pdf')
        plt.show()    

    # ========================================================================================