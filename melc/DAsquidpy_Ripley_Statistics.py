
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

# ---------------------------------------------------------------
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

def dasquidpy_ripley_statistics(adata_pickle_path, dependent_variable_name):
    # ---
    F_AUC_0={}
    F_AUC_1={}
    # ---
    G_AUC_0={}
    G_AUC_1={}
    # ---
    L_AUC_0={}
    L_AUC_1={}
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
        
        nhood_enrichment_F=squidpy.gr.ripley(adata, mode='F', cluster_key='celltype', copy=True)['F_stat']
        nhood_enrichment_Fstats=nhood_enrichment_F.groupby('celltype')['stats'].apply(list).to_dict()
        nhood_enrichment_Fbins=nhood_enrichment_F.groupby('celltype')['bins'].apply(list).to_dict()
        F_auc={}
        for key in nhood_enrichment_Fstats:
            x=nhood_enrichment_Fbins[key]
            y=nhood_enrichment_Fstats[key]
            F_auc[key]=auc(x,y)
        sq.gr.ripley(adata, mode='F', cluster_key="celltype")
        sq.pl.ripley(
            adata,
            mode='F',
            cluster_key="celltype",
            # clusters=['celltype-1', 'celltype-2', 'celltype-3', ..., 'celltype-n']
            figsize=(10, 5),
            )
        # -----
        nhood_enrichment_G=squidpy.gr.ripley(adata, mode='G', cluster_key='celltype', copy=True)['G_stat']
        nhood_enrichment_Gstats=nhood_enrichment_G.groupby('celltype')['stats'].apply(list).to_dict()
        nhood_enrichment_Gbins=nhood_enrichment_G.groupby('celltype')['bins'].apply(list).to_dict()
        G_auc={}
        for key in nhood_enrichment_Gstats:
            x=nhood_enrichment_Gbins[key]
            y=nhood_enrichment_Gstats[key]
            G_auc[key]=auc(x,y)
        sq.gr.ripley(adata, mode='G', cluster_key="celltype")
        sq.pl.ripley(
            adata,
            mode='G',
            cluster_key="celltype",
            # clusters=['celltype-1', 'celltype-2', 'celltype-3', ..., 'celltype-n']
            figsize=(10, 5),
            )
        # -----
        nhood_enrichment_L=squidpy.gr.ripley(adata, mode='L', cluster_key='celltype', copy=True)['L_stat']
        nhood_enrichment_Lstats=nhood_enrichment_L.groupby('celltype')['stats'].apply(list).to_dict()
        nhood_enrichment_Lbins=nhood_enrichment_L.groupby('celltype')['bins'].apply(list).to_dict()
        L_auc={}
        for key in nhood_enrichment_Lstats:
            x=nhood_enrichment_Lbins[key]
            y=nhood_enrichment_Lstats[key]
            L_auc[key]=auc(x,y)
        sq.gr.ripley(adata, mode='L', cluster_key="celltype")
        sq.pl.ripley(
            adata,
            mode='L',
            cluster_key="celltype",
            # clusters=['celltype-1', 'celltype-2', 'celltype-3', ..., 'celltype-n']
            figsize=(10, 5),
            )
        
        
        # ---------------------*************************---------------------
        # ---------------------*************************---------------------
        # ---------------------*************************---------------------


        # ---
        if rec_lab=='POSITIVE' or rec_lab=='positive' or rec_lab=='1' or rec_lab==1:
            ct1=list(adata.obs['celltype'])
            ct1_unique=list(np.unique(ct1))
            for key in F_auc:
                if key in F_AUC_1.keys():
                    F_AUC_1[key].append(F_auc[key])
                else:
                    F_AUC_1[key]=[F_auc[key]]
            # ---
            for key in G_auc:
                if key in G_AUC_1.keys():
                    G_AUC_1[key].append(G_auc[key])
                else:
                    G_AUC_1[key]=[G_auc[key]]
            # ---
            for key in L_auc:
                if key in L_AUC_1.keys():
                    L_AUC_1[key].append(L_auc[key])
                else:
                    L_AUC_1[key]=[L_auc[key]]
            one_or_zero_flag=1
        # ---
        elif rec_lab=='NEGATIVE' or rec_lab=='negative' or rec_lab=='0' or rec_lab==0:
            ct0=list(adata.obs['celltype'])
            ct0_unique=list(np.unique(ct0))
            for key in F_auc:
                if key in F_AUC_0.keys():
                    F_AUC_0[key].append(F_auc[key])
                else:
                    F_AUC_0[key]=[F_auc[key]]
            # ---
            for key in G_auc:
                if key in G_AUC_0.keys():
                    G_AUC_0[key].append(G_auc[key])
                else:
                    G_AUC_0[key]=[G_auc[key]]
            # ---
            for key in L_auc:
                if key in L_AUC_0.keys():
                    L_AUC_0[key].append(L_auc[key])
                else:
                    L_AUC_0[key]=[L_auc[key]]
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

    F_AUC_list=[]
    G_AUC_list=[]
    L_AUC_list=[]

    # ---
    for i in celltypes_0_and_1:
        U1, p = mannwhitneyu(F_AUC_0[i], F_AUC_1[i], method="exact")
        F_AUC_list.append(p)
        # ---
        U1, p = mannwhitneyu(G_AUC_0[i], G_AUC_1[i], method="exact")
        G_AUC_list.append(p)
        # ---
        U1, p = mannwhitneyu(L_AUC_0[i], L_AUC_1[i], method="exact")
        L_AUC_list.append(p)
        # ---
    # ---

    F_AUC_dict=dict(zip(celltypes_0_and_1, F_AUC_list))
    G_AUC_dict=dict(zip(celltypes_0_and_1, G_AUC_list))
    L_AUC_dict=dict(zip(celltypes_0_and_1, L_AUC_list))

    # -------------------------

    significant_F_AUC_dict = dict((k, v) for k, v in F_AUC_dict.items() if v < 0.05)
    significant_F_AUC_dict_sorted = dict(sorted(significant_F_AUC_dict.items(), key=lambda x:x[1]))

    # ---

    significant_G_AUC_dict = dict((k, v) for k, v in G_AUC_dict.items() if v < 0.05)
    significant_G_AUC_dict_sorted = dict(sorted(significant_G_AUC_dict.items(), key=lambda x:x[1]))

    # ---

    significant_L_AUC_dict = dict((k, v) for k, v in L_AUC_dict.items() if v < 0.05)
    significant_L_AUC_dict_sorted = dict(sorted(significant_L_AUC_dict.items(), key=lambda x:x[1]))

    # ========================================================================================

    # ---
    if len(significant_F_AUC_dict_sorted)==0:
        print('No significant celltype pairs found in terms of Ripley F-statistic!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Ripley F-statistic) – celltype pairs with significant Ripley F-statistic differences (p-value<0.05)')
        feature = list(significant_F_AUC_dict_sorted.keys())
        score = list(significant_F_AUC_dict_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_ripleyF_pvalues_significantCelltypePairs_allPatients.pdf', format='pdf')
        plt.show()

    # ---

    if len(significant_G_AUC_dict_sorted)==0:
        print('No significant celltype pairs found in terms of Ripley G-statistic!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Ripley G-statistic) – celltype pairs with significant Ripley G-statistic differences (p-value<0.05)')
        feature = list(significant_G_AUC_dict_sorted.keys())
        score = list(significant_G_AUC_dict_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_ripleyG_pvalues_significantCelltypePairs_allPatients.pdf', format='pdf')
        plt.show()

    # ---

    if len(significant_L_AUC_dict_sorted)==0:
        print('No significant celltype pairs found in terms of Ripley L-statistic!')
    else:
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential Analysis (Ripley L-statistic) – celltype pairs with significant Ripley L-statistic differences (p-value<0.05)')
        feature = list(significant_L_AUC_dict_sorted.keys())
        score = list(significant_L_AUC_dict_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_ripleyL_pvalues_significantCelltypePairs_allPatients.pdf', format='pdf')
        plt.show()

    # ---

    # ========================================================================================