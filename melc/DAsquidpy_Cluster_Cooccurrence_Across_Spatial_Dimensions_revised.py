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
from sklearn.metrics import auc

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

def extract_unique_tuples_dict(input_list):
    return [tuple(x) for x in np.unique([sorted(tup) for tup in input_list], axis=0)]

def dasquidpy_cluster_cooccurrence_across_spatial_dimensions(adata_pickle_path, dependent_variable_name, interval):
    # ---
    AUC_0={}
    AUC_1={}
    # ---
    Max_0={}
    Max_1={}
    # ---
    Mean_0={}
    Mean_1={}
    # ---
    StdDev_0={}
    StdDev_1={}
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
        no_of_clusts=len(clusters_sorted)
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
        nhood_enrichment=squidpy.gr.co_occurrence(adata, cluster_key='celltype', interval=interval, copy=True, backend="multiprocessing", show_progress_bar=False)
        # ---
        clust_matrix={}
        for clust1_ in range(no_of_clusts):
            for clust2_ in range(no_of_clusts):
                clust1=clusters_sorted[clust1_]
                clust2=clusters_sorted[clust2_]
                if clust1!=clust2:
                    sorted_strs=sorted([clust1, clust2])
                    if (sorted_strs[0], sorted_strs[1]) in clust_matrix.keys():
                        clust_matrix[(sorted_strs[0], sorted_strs[1])]=(clust_matrix[(sorted_strs[0], sorted_strs[1])]+nhood_enrichment[0][clust1_,clust2_,:])/2
                    else:
                        clust_matrix[(sorted_strs[0], sorted_strs[1])]=nhood_enrichment[0][clust1_,clust2_,:]
        # ---
        no_of_bins=len(nhood_enrichment[1])
        x_axis=[]
        cnt_nbins=-1
        for n_bins in range(no_of_bins-1):
            cnt_nbins+=1
            x_axis.append((nhood_enrichment[1][n_bins]+nhood_enrichment[1][n_bins+1])/2)
        # ---
        
        # ================================== Function similarities: ==================================
        
        auc_={}
        max_={}
        mean_={}
        std_dev={} # standard_deviation = sqrt(variance)
        for key in clust_matrix:
            auc_[key]=auc(x_axis, list(clust_matrix[key]))
            max_[key]=max(list(clust_matrix[key]))
            mean_[key]=np.mean(list(clust_matrix[key]))
            std_dev[key]=np.std(list(clust_matrix[key]))
        
        
        # ---------------------*************************---------------------
        # ---------------------*************************---------------------
        # ---------------------*************************---------------------


        # ---
        if rec_lab=='POSITIVE' or rec_lab=='positive' or rec_lab=='1' or rec_lab==1:
            ct1=list(adata.obs['celltype'])
            ct1_unique=list(np.unique(ct1))
            for key in auc_:
                if key in AUC_1.keys():
                    AUC_1[key].append(auc_[key])
                else:
                    AUC_1[key]=[auc_[key]]
            # ---
            for key in max_:
                if key in Max_1.keys():
                    Max_1[key].append(max_[key])
                else:
                    Max_1[key]=[max_[key]]
            # ---
            for key in mean_:
                if key in Mean_1.keys():
                    Mean_1[key].append(mean_[key])
                else:
                    Mean_1[key]=[mean_[key]]
            # ---
            for key in std_dev:
                if key in StdDev_1.keys():
                    StdDev_1[key].append(std_dev[key])
                else:
                    StdDev_1[key]=[std_dev[key]]
            one_or_zero_flag=1
        # ---
        elif rec_lab=='NEGATIVE' or rec_lab=='negative' or rec_lab=='0' or rec_lab==0:
            ct0=list(adata.obs['celltype'])
            ct0_unique=list(np.unique(ct0))
            for key in auc_:
                if key in AUC_0.keys():
                    AUC_0[key].append(auc_[key])
                else:
                    AUC_0[key]=[auc_[key]]
            # ---
            for key in max_:
                if key in Max_0.keys():
                    Max_0[key].append(max_[key])
                else:
                    Max_0[key]=[max_[key]]
            # ---
            for key in mean_:
                if key in Mean_0.keys():
                    Mean_0[key].append(mean_[key])
                else:
                    Mean_0[key]=[mean_[key]]
            # ---
            for key in std_dev:
                if key in StdDev_0.keys():
                    StdDev_0[key].append(std_dev[key])
                else:
                    StdDev_0[key]=[std_dev[key]]
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


    AUC_list=[]
    Max_list=[]
    Mean_list=[]
    StdDev_list=[]

    AUC_celltypePairs_0=list(AUC_0.keys())
    AUC_celltypePairs_1=list(AUC_1.keys())

    AUC_celltypePairs_0and1=list(set(AUC_celltypePairs_0).intersection(set(AUC_celltypePairs_1)))

    # Max_celltypePairs=list(Max_0.keys())
    # Mean_celltypePairs=list(Mean_0.keys())
    # StdDev_celltypePairs=list(StdDev_0.keys())

    # ---
    for i in AUC_celltypePairs_0and1:
        U1, p = mannwhitneyu(AUC_0[i], AUC_1[i], method="exact")
        AUC_list.append(p)
        # ---
        U1, p = mannwhitneyu(Max_0[i], Max_1[i], method="exact")
        Max_list.append(p)
        # ---
        U1, p = mannwhitneyu(Mean_0[i], Mean_1[i], method="exact")
        Mean_list.append(p)
        # ---
        U1, p = mannwhitneyu(StdDev_0[i], StdDev_1[i], method="exact")
        StdDev_list.append(p)
        # ---
    # ---

    AUC_dict=dict(zip(celltypes_0_and_1, AUC_list))
    Max_dict=dict(zip(celltypes_0_and_1, Max_list))
    Mean_dict=dict(zip(celltypes_0_and_1, Mean_list))
    StdDev_dict=dict(zip(celltypes_0_and_1, StdDev_list))

    # -------------------------
    significant_AUC_dict = dict((k, v) for k, v in AUC_dict.items() if v < 0.05)
    significant_AUC_dict_sorted = dict(sorted(significant_AUC_dict.items(), key=lambda x:x[1]))
    # ---
    significant_Max_dict = dict((k, v) for k, v in Max_dict.items() if v < 0.05)
    significant_Max_dict_sorted = dict(sorted(significant_Max_dict.items(), key=lambda x:x[1]))
    # ---
    significant_Mean_dict = dict((k, v) for k, v in Mean_dict.items() if v < 0.05)
    significant_Mean_dict_sorted = dict(sorted(significant_Mean_dict.items(), key=lambda x:x[1]))
    # ---
    significant_StdDev_dict = dict((k, v) for k, v in StdDev_dict.items() if v < 0.05)
    significant_StdDev_dict_sorted = dict(sorted(significant_StdDev_dict.items(), key=lambda x:x[1]))

    # ========================================================================================



    # ---

    if len(significant_AUC_dict_sorted)==0:
        print('No significant celltype pairs found in terms of cluster co-occurrence across spatial dimensions (AUC)!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Cluster Co-occurrence) – celltype pairs with significant cluster co-occurrence AUC differences (p-value<0.05)')
        feature = list(significant_AUC_dict_sorted.keys())
        score = list(significant_AUC_dict_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_ClustCooccurrAUC_pvalues_significantCelltypePairs_allPatients.pdf', format='pdf')
        plt.show()

    # ---

    if len(significant_Max_dict_sorted)==0:
        print('No significant celltype pairs found in terms of cluster co-occurrence across spatial dimensions (max-value)!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Cluster Co-occurrence) – celltype pairs with significant cluster co-occurrence max-value differences (p-value<0.05)')
        feature = list(significant_Max_dict_sorted.keys())
        score = list(significant_Max_dict_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_ClustCooccurrMax_pvalues_significantCelltypePairs_allPatients.pdf', format='pdf')
        plt.show()

    # ---

    if len(significant_Mean_dict_sorted)==0:
        print('No significant celltype pairs found in terms of cluster co-occurrence across spatial dimensions (mean-value)!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Cluster Co-occurrence) – celltype pairs with significant cluster co-occurrence mean-value differences (p-value<0.05)')
        feature = list(significant_Mean_dict_sorted.keys())
        score = list(significant_Mean_dict_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_ClustCooccurrMean_pvalues_significantCelltypePairs_allPatients.pdf', format='pdf')
        plt.show()

    # ---

    if len(significant_StdDev_dict_sorted)==0:
        print('No significant celltype pairs found in terms of cluster co-occurrence across spatial dimensions (standard deviation)!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Cluster Co-occurrence) – celltype pairs with significant cluster co-occurrence standard deviation differences (p-value<0.05)')
        feature = list(significant_StdDev_dict_sorted.keys())
        score = list(significant_StdDev_dict_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_ClustCooccurrStdDev_pvalues_significantCelltypePairs_allPatients.pdf', format='pdf')
        plt.show()

    # ---