


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
import schist as scs
from pingouin import mwu
import graph_tools as gt
import itertools
from itertools import permutations
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_no_of_cluster_levels_vs_no_of_clusters_per_cluster_level(adata_pickle_path, clusterings_patientLevel_dict_path):
    # method='top-down'
    # method='bottom-up' # Preferred method (this is how the algorithm works)!
    # ---
    CLUSTERING_RESULTS=[]*20
    # ---
    NO_OF_CLUSTERS_B1={}
    NO_OF_CLUSTERS_B2={}
    # ---
    SCORE=[]

    # SCORE=[None]*20

    # with open("clusterings_patientLevel_dict.pickle", "rb") as f:
    #     clusterings_patientLevel_dict = pickle.load(f)
        

    with open(adata_pickle_path, "rb") as f:
        pickle_ = pickle.load(f)

    # ------
    pickle_items=list(pickle_.items())
    first_three_items = pickle_items[0:2]
    pickle_new=dict(first_three_items)
    pickle_=pickle_new
    # ------


    with open(clusterings_patientLevel_dict_path, 'rb') as f:
        clusterings_patientLevel_dict= pickle.load(f)
    
    # ---
    clusterings_patientLevel_dict_items=list(clusterings_patientLevel_dict.items())
    first_three_items = clusterings_patientLevel_dict_items[0:2]
    clusterings_patientLevel_dict_new=dict(first_three_items)
    clusterings_patientLevel_dict=clusterings_patientLevel_dict_new
    # ---

    cnt=-1

    max_cluster_level=0
    min_cluster_level=float('inf')
    no_of_cluster_levels=[]
    cluster_size_dict=[]

    for i in pickle_:
        
        cnt+=1

        adata=pickle_[i]
        adata_X=adata.to_df()
        
        no_of_cells=len(adata_X.index)
        
        no_of_clusterings=3
        clustering_names=[]
        
        for clustering in range(no_of_clusterings):
            cl_names='Clustering-'+str(clustering)
            clustering_names.append(cl_names)
        
        # =============================================================================================================
        
        # ----------------------------------------------- Combinations: -----------------------------------------------
        for clust_count in clustering_names:
            B1_dict=clusterings_patientLevel_dict[i][0][clust_count][0]
            max_cluster_level_temp=int(max(B1_dict.keys()))
            if max_cluster_level_temp>max_cluster_level:
                max_cluster_level=max_cluster_level_temp
            if max_cluster_level_temp<min_cluster_level:
                min_cluster_level=int(max_cluster_level_temp)
            cluster_size_dict_temp={}
            for j in range(1,max_cluster_level_temp+1):
                cluster_size_dict_temp[j]=len(B1_dict[str(j)][0].keys())
            no_of_cluster_levels.append(max_cluster_level_temp)
            cluster_size_dict.append(cluster_size_dict_temp)


    clusterLevel_vs_ClusterSize={}
    for i in range(min_cluster_level, max_cluster_level+1):
        clusterLevel_vs_ClusterSize[i]={}
        for j in range(1, i+1):
            clusterLevel_vs_ClusterSize[i][j]=[]



    count=-1
    for i in no_of_cluster_levels:
        count+=1
        # cluster_dict_temp={}
        for j in range(1,i+1):
            clusterLevel_vs_ClusterSize[i][j].append(cluster_size_dict[count][j])


    # # --------

    no_of_cols=3
    no_of_rows=int(int(max_cluster_level-min_cluster_level)/no_of_cols)+1

    fig, axes = plt.subplots(no_of_rows, no_of_cols, figsize=(18, 10))
    fig.suptitle('No. of cluster levels vs. No. of clusters per cluster level')

    _cnt_=-1
    for i in clusterLevel_vs_ClusterSize:
        _cnt_+=1
        # print(_cnt_)
        row=int(_cnt_/no_of_cols)
        print(row)
        col=_cnt_%no_of_cols
        # _col_=_cnt_%no_of_cols
        # print(_col_)
        # if _col_==-1:
        #     col=no_of_cols-1
        # else:
        #     col=_col_
        print(row, col)
        df=pd.DataFrame.from_dict(clusterLevel_vs_ClusterSize[i])
        
        if no_of_rows>1:
            _ax_=axes[int(row), int(col)]
            sns.boxplot(ax=_ax_, data=df)
            _ax_.set_xlabel('No. of clustering levels')
            _ax_.set_ylabel('No. of clusters per level')
            _ax_.set_title(f'No. of clustering levels: {i}')
        elif no_of_rows==1:
            _ax_=axes[int(col)]
            sns.boxplot(ax=_ax_, data=df)
            _ax_.set_xlabel('No. of clustering levels')
            _ax_.set_ylabel('No. of clusters per level')
            _ax_.set_title(f'No. of clustering levels: {i}')
            
            
    axes[0].set_xlabel('common xlabel')
    axes[0].set_ylabel('common ylabel')


    plt.savefig("No. of cluster levels vs. No. of clusters per cluster level.jpg", format='jpg')


        
        