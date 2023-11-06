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
        
        nhood_enrichment=squidpy.gr.co_occurrence(adata, cluster_key='celltype', interval=interval, copy=True, backend="multiprocessing", show_progress_bar=False)
        # sq.gr.co_occurrence(adata, cluster_key="celltype", interval=interval, backend="multiprocessing", show_progress_bar=False)
        # sq.pl.co_occurrence(
        #    adata,
        #    cluster_key="celltype",
        #    # clusters=['celltype-1', 'celltype-2', 'celltype-3', ..., 'celltype-n']
        #    figsize=(10, 5),
        #    )
        nhood_enrichment_zscore=nhood_enrichment[0]
        # nhood_enrichment_zscores.append(nhood_enrichment_zscore)
        # ---
        # ---------------------*************************---------------------
        _zscore_list_=[]
        for _cnt_ in range(interval-1):
            upper_zscore_matrix = np.triu(np.array(nhood_enrichment_zscore[:,:,_cnt_]), 1)
            zscore_list=[]
            for j in upper_zscore_matrix:
                for k in j:
                    zscore_list.append(k)
            # ---
            dict_nhoodEnrichment_zscore=dict(zip(cluster_pairs_list, zscore_list))
            _zscore_list_.append(dict_nhoodEnrichment_zscore)
        # ---
        Dict_NhoodEnrichment_Zscore.append(_zscore_list_)
        # ---
        
        
        # ---------------------*************************---------------------
        # ---------------------*************************---------------------
        # ---------------------*************************---------------------


        if rec_lab=='POSITIVE' or rec_lab=='positive' or rec_lab=='1' or rec_lab==1:
            ct1=list(adata.obs['celltype'])
            ct1_unique=list(np.unique(ct1))
            one_or_zero_flag=1
        elif rec_lab=='NEGATIVE' or rec_lab=='negative' or rec_lab=='0' or rec_lab==0:
            ct0=list(adata.obs['celltype'])
            ct0_unique=list(np.unique(ct0))
            one_or_zero_flag=0
        
        # ---
        if one_or_zero_flag==0:
            Dict_NhoodEnrichment_Zscore_0.append(_zscore_list_)
        elif one_or_zero_flag==1:
            Dict_NhoodEnrichment_Zscore_1.append(_zscore_list_)
        # ---
        for j1, j2 in cluster_pairs_list:
            Cluster_Pairs_List.append((j1, j2))
            if one_or_zero_flag==0:
                Cluster_Pairs_List_0.append((j1, j2))
            elif one_or_zero_flag==1:
                Cluster_Pairs_List_1.append((j1, j2))
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




    input_list = Cluster_Pairs_List.copy()
    Cluster_Pairs_List = extract_unique_tuples_dict(input_list)
    # ---
    input_list = Cluster_Pairs_List_0.copy()
    Cluster_Pairs_List_0 = extract_unique_tuples_dict(input_list)
    # ---
    input_list = Cluster_Pairs_List_1.copy()
    Cluster_Pairs_List_1 = extract_unique_tuples_dict(input_list)
    # ---
    Cluster_Pairs_List_1_0=list(set(Cluster_Pairs_List_1).difference(set(Cluster_Pairs_List_0)))
    if len(Cluster_Pairs_List_1_0)==0:
        print('No. of cell type pairs in positive class not in control class = 0')
    else:
        pass
    # ---
    Cluster_Pairs_List_0_1=list(set(Cluster_Pairs_List_0).difference(set(Cluster_Pairs_List_1)))
    if len(Cluster_Pairs_List_0_1)==0:
        print('No. of cell type pairs in control class not in positive class = 0')
    else:
        pass
    # ---
    Cluster_Pairs_List_0_and_1=list(set(Cluster_Pairs_List_1).intersection(set(Cluster_Pairs_List_0)))


    # # ========================================================================================

    Dict_NhoodEnrichment_Zscore_aggregated_=[]
    for i in range(interval-1):
        res = {}
        for dict_ in np.array(Dict_NhoodEnrichment_Zscore)[:,i].tolist():
            for list_ in dict_:
                if list_ in res:
                    res[list_].append(dict_[list_])
                else:
                    res[list_] = [dict_[list_]]
        Dict_NhoodEnrichment_Zscore_aggregated=res
        Dict_NhoodEnrichment_Zscore_aggregated_.append(Dict_NhoodEnrichment_Zscore_aggregated)
    # ---
    Dict_NhoodEnrichment_Zscore_0_aggregated_=[]
    for i in range(interval-1):
        res = {}
        for dict_ in np.array(Dict_NhoodEnrichment_Zscore_0)[:,i].tolist():
            for list_ in dict_:
                if list_ in res:
                    res[list_].append(dict_[list_])
                else:
                    res[list_] = [dict_[list_]]
        Dict_NhoodEnrichment_Zscore_0_aggregated=res
        Dict_NhoodEnrichment_Zscore_0_aggregated_.append(Dict_NhoodEnrichment_Zscore_0_aggregated)    
    # ---
    Dict_NhoodEnrichment_Zscore_1_aggregated_=[]
    for i in range(interval-1):
        res = {}
        for dict_ in np.array(Dict_NhoodEnrichment_Zscore_1)[:,i].tolist():
            for list_ in dict_:
                if list_ in res:
                    res[list_].append(dict_[list_])
                else:
                    res[list_] = [dict_[list_]]
        Dict_NhoodEnrichment_Zscore_1_aggregated=res
        Dict_NhoodEnrichment_Zscore_1_aggregated_.append(Dict_NhoodEnrichment_Zscore_1_aggregated)    
    # ---

    # ========================================================================================

    _dict_nhood_enrichment_zscore_=[]
    _dict_nhood_enrichment_zscore_0_=[]
    _dict_nhood_enrichment_zscore_1_=[]
    for i in range(interval-1):
        # ========================================================================================
        nhood_enrichment_zscore_values=[Dict_NhoodEnrichment_Zscore_aggregated_[i][x] for x in Cluster_Pairs_List_0_and_1]
        dict_nhood_enrichment_zscore=dict(zip(Cluster_Pairs_List_0_and_1, nhood_enrichment_zscore_values))
        _dict_nhood_enrichment_zscore_.append(dict_nhood_enrichment_zscore)
        # ---
        nhood_enrichment_zscore_0_values=[Dict_NhoodEnrichment_Zscore_0_aggregated_[i][x] for x in Cluster_Pairs_List_0_and_1]
        dict_nhood_enrichment_zscore_0=dict(zip(Cluster_Pairs_List_0_and_1, nhood_enrichment_zscore_0_values))
        _dict_nhood_enrichment_zscore_0_.append(dict_nhood_enrichment_zscore_0)
        # ---
        nhood_enrichment_zscore_1_values=[Dict_NhoodEnrichment_Zscore_1_aggregated_[i][x] for x in Cluster_Pairs_List_0_and_1]
        dict_nhood_enrichment_zscore_1=dict(zip(Cluster_Pairs_List_0_and_1, nhood_enrichment_zscore_1_values))
        _dict_nhood_enrichment_zscore_1_.append(dict_nhood_enrichment_zscore_1)
        # ========================================================================================


    # ========================================================================================

    _pvals_zscore_nhoodEnrichment_=[]
    for k in range(interval-1):
        pvals_zscore_nhoodEnrichment_values=[]
        # ---
        for i in Cluster_Pairs_List_0_and_1:
            U1, p = mannwhitneyu(_dict_nhood_enrichment_zscore_0_[k][i], _dict_nhood_enrichment_zscore_1_[k][i], method="exact")
            pvals_zscore_nhoodEnrichment_values.append(p)
        # ---
        pvals_zscore_nhoodEnrichment=dict(zip(Cluster_Pairs_List_0_and_1, pvals_zscore_nhoodEnrichment_values))
        _pvals_zscore_nhoodEnrichment_.append(pvals_zscore_nhoodEnrichment)
        

    # ========================================================================================

    significant_pvals_zscore_nhoodEnrichment_=[]
    for i in range(interval-1):    
        significant_pvals_zscore_nhoodEnrichment = dict((k, v) for k, v in _pvals_zscore_nhoodEnrichment_[i].items() if v < 0.05)
        significant_pvals_zscore_nhoodEnrichment_sorted = dict(sorted(significant_pvals_zscore_nhoodEnrichment.items(), key=lambda x:x[1]))
        significant_pvals_zscore_nhoodEnrichment_.append(significant_pvals_zscore_nhoodEnrichment_sorted)


    # ========================================================================================

    significant_proteinPairs_pvals={}
    for dict_ in significant_pvals_zscore_nhoodEnrichment_:
        for list_ in dict_:
            if list_ in significant_proteinPairs_pvals:
                significant_proteinPairs_pvals[list_].append(dict_[list_])
            else:
                significant_proteinPairs_pvals[list_] = [dict_[list_]]
    Significant_ProteinPairs_Pvals=significant_proteinPairs_pvals


    # ========================================================================================

    for i in Significant_ProteinPairs_Pvals:
        Significant_ProteinPairs_Pvals[i]=min(Significant_ProteinPairs_Pvals[i])

    Significant_ProteinPairs_Pvals_sorted = dict(sorted(Significant_ProteinPairs_Pvals.items(), key=lambda x:x[1]))


    # ========================================================================================


    if len(Significant_ProteinPairs_Pvals_sorted)==0:
        print('No significant celltype pair found in terms of neighborhodd enrichment!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Celltype Cooccurrence) – celltype pairs with significant co-occurrences (p-value<0.05)')
        feature = list(Significant_ProteinPairs_Pvals_sorted.keys())
        score = list(Significant_ProteinPairs_Pvals_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_celltypeCooccurrence_pvalues_significantCelltypePairs_allPatients.pdf', format='pdf')
        plt.show()


    # ========================================================================================




























    # significant_protein_pairs={}



    # for i in range(interval-1):
    #     pvals_zscore_nhoodEnrichment_values=[]
    #     # ---
    #     for i in Cluster_Pairs_List_0_and_1:
    #         U1, p = mannwhitneyu(dict_nhood_enrichment_zscore_0[i], dict_nhood_enrichment_zscore_1[i], method="exact")
    #         pvals_zscore_nhoodEnrichment_values.append(p)
    #     # ---
    #     pvals_zscore_nhoodEnrichment=dict(zip(Cluster_Pairs_List_0_and_1, pvals_zscore_nhoodEnrichment_values))






    # ========================================================================================






    # # significant_pvals_zscore_nhoodEnrichment = dict((k, v) for k, v in pvals_zscore_nhoodEnrichment.items() if v < 0.05)
    # # significant_pvals_zscore_nhoodEnrichment_sorted = dict(sorted(significant_pvals_zscore_nhoodEnrichment.items(), key=lambda x:x[1]))

    # # significant_pvals_count_nhoodEnrichment = dict((k, v) for k, v in pvals_count_nhoodEnrichment.items() if v < 0.05)
    # # significant_pvals_count_nhoodEnrichment_sorted = dict(sorted(significant_pvals_count_nhoodEnrichment.items(), key=lambda x:x[1]))
    # # # ========================================================================================
    # # if len(significant_pvals_zscore_nhoodEnrichment)==0:
    # #     print('No significant celltype pair found in terms of neighborhodd enrichment!')
    # # else:
    # #     # ---
    # #     fig, axes = plt.subplots(figsize=(10, 4))
    # #     fig.suptitle(f'Differential Analysis (Neighborhood Enrichment) – celltype pairs with significant neighborhood enrichment differences (p-value<0.05)')
    # #     feature = list(significant_pvals_zscore_nhoodEnrichment_sorted.keys())
    # #     score = list(significant_pvals_zscore_nhoodEnrichment_sorted.values())
    # #     x_pos = np.arange(len(feature))
    # #     # ---
    # #     plt.bar(x_pos, score,align='center')
    # #     plt.xticks(x_pos, feature, rotation=90) 
    # #     plt.ylabel('p-value (MWU test)')
    # #     plt.savefig(f'DA_neighborhoodEnrichment_pvalues_significantCelltypePairs_allPatients.pdf', format='pdf')
    # #     plt.show()