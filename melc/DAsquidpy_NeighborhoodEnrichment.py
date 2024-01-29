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

# ----------------------------------

cellTypeName_shortForm_dict={
    'Fibroblasts': 'Fib',
    'Basal keratinocytes': 'Bas',
    'Melanocytes': 'Mel',
    'Smooth muscle cells': 'SMCs',
    'Suprabasal keratinocytes': 'Sup',
    'T-cells': 'T',
    'Granulocytes': 'Gran',
    'B-cells': 'B',
    'Endothelial cells': 'Endo',
    'Langerhans cells': 'Lang',
    'Macrophages': 'Macro'
    }

cellTypeName_list=list(cellTypeName_shortForm_dict.keys())
shortForm_list=list(cellTypeName_shortForm_dict.values())

cellTypeName_shortForm_str=""
for i in cellTypeName_shortForm_dict:
    cellTypeName_shortForm_str+="'"+ cellTypeName_shortForm_dict[i] + "'" + ": " + i + "\n"
cellTypeName_shortForm_str=cellTypeName_shortForm_str[:-1]

# ----------------------------------


def extract_unique_tuples_dict(input_list):
        return [tuple(x) for x in np.unique([sorted(tup) for tup in input_list], axis=0)]


def dasquidpy_neighborhood_enrichment(adata_pickle_path, dependent_variable_name):
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

    okay_list=[]
    error_list=[]
    for i in pickle_:
        adata=pickle_[i]
        try:
            squidpy.gr.spatial_neighbors(adata)
            nhood_enrichment=squidpy.gr.nhood_enrichment(adata, cluster_key='celltype', copy=True, backend="multiprocessing", show_progress_bar=False)
            okay_list.append(i)
        except:
            error_list.append(i)    


    cnt=-1
    for i in pickle_:
        if not(i in error_list):
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
            print(i)
            squidpy.gr.spatial_neighbors(adata)
            nhood_enrichment=squidpy.gr.nhood_enrichment(adata, cluster_key='celltype', copy=True, backend="multiprocessing", show_progress_bar=False)
            # squidpy.gr.nhood_enrichment(adata, cluster_key='celltype', backend="multiprocessing", show_progress_bar=False)
            # sq.pl.nhood_enrichment(adata, cluster_key="celltype")
            nhood_enrichment_zscore=nhood_enrichment[0]
            # nhood_enrichment_zscores.append(nhood_enrichment_zscore)
            # ---
            nhood_enrichment_count=nhood_enrichment[1]
            # nhood_enrichment_counts.append(nhood_enrichment_count)
            # ---
            # if adata_raw.obsm[dependent_variable_name]=="POSITIVE":
            #     leiden_nhood_enrichment_counts_positive.append(leiden_nhood_enrichment_count)
            # elif adata_raw.obsm[dependent_variable_name]=="NEGATIVE":
            #     leiden_nhood_enrichment_counts_negative.append(leiden_nhood_enrichment_count)
            # ---
            # ---------------------*************************---------------------
            
            upper_zscore_matrix = np.triu(np.array(nhood_enrichment_zscore), 1)
            zscore_list=[]
            for j in upper_zscore_matrix:
                for k in j:
                    # print(j, k)
                    if k!=0 and k!=None:
                        # cluster_pairs_list.append(k)
                        zscore_list.append(k)
            # ---
            upper_count_matrix = np.triu(np.array(nhood_enrichment_count), 1)
            count_list=[]
            for j in upper_count_matrix:
                for k in j:
                    # print(j, k)
                    if k!=0 and k!=None:
                        # cluster_pairs_list.append(k)
                        count_list.append(k)
            # ---
            dict_nhoodEnrichment_zscore=dict(zip(cluster_pairs_list, zscore_list))
            dict_nhoodEnrichment_count=dict(zip(cluster_pairs_list, count_list))
            # ---
            Dict_NhoodEnrichment_Zscore.append(dict_nhoodEnrichment_zscore)
            Dict_NhoodEnrichment_Count.append(dict_nhoodEnrichment_count)
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
                Dict_NhoodEnrichment_Zscore_0.append(dict_nhoodEnrichment_zscore)
                Dict_NhoodEnrichment_Count_0.append(dict_nhoodEnrichment_count)
            elif one_or_zero_flag==1:
                Dict_NhoodEnrichment_Zscore_1.append(dict_nhoodEnrichment_zscore)
                Dict_NhoodEnrichment_Count_1.append(dict_nhoodEnrichment_count)
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


    # ========================================================================================


    celltypes_1_0=list(set(cell_types_1).difference(set(cell_types_0)))
    celltypes_0_1=list(set(cell_types_0).difference(set(cell_types_1)))

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
    # ---
    fig, axs = plt.subplots(figsize=(10, 4))
    # fig.suptitle(f'Cell type counts across conditions')
    # ---
    axs = cell_types_count_df.plot(x="cell type", y=["Positive class (recurrent TNBC)", "Control class (non-recurrent TNBC)"], kind="bar", rot=0, title=f'Cell type counts across conditions')                                                                                                                                                  
    # ---
    plt.ylabel('Cell count')
    plt.xlabel('Cell type')
    # ---
    plt.xticks(rotation = 45, fontsize=8)
    # ---
    plt.savefig(f'cell_counts_across_conditions_allPatients.pdf', format='pdf')
    plt.show()


    # # ================ Approach-1: =======================
    # len_=len(cell_types_0_count)
    # X = np.arange(len_)
    # # no_of_subplots='1'*len_
    # # no_of_subplots=int(no_of_subplots)
    # # ax=plt.subplot(no_of_subplots)
    # ax=plt.subplot(111)
    # # ---
    # ax.bar(X, cell_types_1_count.values(), width=0.2, color='b', align='center')
    # ax.bar(X-0.2, cell_types_0_count.values(), width=0.2, color='g', align='center')
    # # ---
    # ax.legend(('Positive class (recurrent TNBC)','Control class (non-recurrent TNBC)'))

    # plt.xticks(X, cell_types_1_count.keys())
    # plt.title("Celltype distributions across conditions", fontsize=17)
    # plt.show()
    # # =======================================

    # # ================= Approach-2: ======================
    # celltypes=list(cell_types_1_count.keys())
    # conditionWise_celltypes={
    #     'Positive class (recurrent TNBC)': [],
    #     'Control class (non-recurrent TNBC': []
    #     }

    # for i in celltypes:
    #     conditionWise_celltypes['Positive class (recurrent TNBC)'].append(cell_types_1_count[i])
    #     conditionWise_celltypes['Control class (non-recurrent TNBC'].append(cell_types_0_count[i])

    # x = np.arange(len(celltypes))  # the label locations
    # width = 0.25  # the width of the bars
    # multiplier = 0

    # fig, ax = plt.subplots(layout='constrained')

    # for attribute, measurement in conditionWise_celltypes.items():
    #     offset = width * multiplier
    #     rects = ax.bar(x + offset, measurement, width, label=attribute)
    #     ax.bar_label(rects, padding=3)
    #     multiplier += 1

    # # Add some text for labels, title and custom x-axis tick labels, etc.
    # ax.set_ylabel('Length (mm)')
    # ax.set_title('Penguin attributes by species')
    # ax.set_xticks(x + width, celltypes)
    # ax.legend(loc='upper left', ncols=3)
    # ax.set_ylim(0, 250)
    # plt.show()
    # # =======================================

    # # ========================================================================================

    input_list = Cluster_Pairs_List.copy()
    Cluster_Pairs_List = extract_unique_tuples_dict(input_list)
    # print(f"The original list : {input_list}")
    # print(f"The list after duplicated removal : {Cluster_Pairs_List}")
    # ---
    input_list = Cluster_Pairs_List_0.copy()
    Cluster_Pairs_List_0 = extract_unique_tuples_dict(input_list)
    # print(f"The original list : {input_list}")
    # print(f"The list after duplicated removal : {Cluster_Pairs_List}")
    # ---
    input_list = Cluster_Pairs_List_1.copy()
    Cluster_Pairs_List_1 = extract_unique_tuples_dict(input_list)
    # print(f"The original list : {input_list}")
    # print(f"The list after duplicated removal : {Cluster_Pairs_List}")
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


    res = {}
    for dict_ in Dict_NhoodEnrichment_Zscore:
        for list_ in dict_:
            if list_ in res:
                res[list_].append(dict_[list_])
            else:
                res[list_] = [dict_[list_]]
    Dict_NhoodEnrichment_Zscore_aggregated=res
    # ---
    res = {}
    for dict_ in Dict_NhoodEnrichment_Count:
        for list_ in dict_:
            if list_ in res:
                res[list_].append(dict_[list_])
            else:
                res[list_] = [dict_[list_]]
    Dict_NhoodEnrichment_Count_aggregated=res
    # ---
    res = {}
    for dict_ in Dict_NhoodEnrichment_Zscore_0:
        for list_ in dict_:
            if list_ in res:
                res[list_].append(dict_[list_])
            else:
                res[list_] = [dict_[list_]]
    Dict_NhoodEnrichment_Zscore_0_aggregated=res
    # ---
    res = {}
    for dict_ in Dict_NhoodEnrichment_Count_0:
        for list_ in dict_:
            if list_ in res:
                res[list_].append(dict_[list_])
            else:
                res[list_] = [dict_[list_]]
    Dict_NhoodEnrichment_Count_0_aggregated=res
    # ---
    res = {}
    for dict_ in Dict_NhoodEnrichment_Zscore_1:
        for list_ in dict_:
            if list_ in res:
                res[list_].append(dict_[list_])
            else:
                res[list_] = [dict_[list_]]
    Dict_NhoodEnrichment_Zscore_1_aggregated=res
    # ---
    res = {}
    for dict_ in Dict_NhoodEnrichment_Count_1:
        for list_ in dict_:
            if list_ in res:
                res[list_].append(dict_[list_])
            else:
                res[list_] = [dict_[list_]]
    Dict_NhoodEnrichment_Count_1_aggregated=res
    # ---


    # ========================================================================================

    # nhood_enrichment_zscore_values=[Dict_NhoodEnrichment_Zscore_aggregated[x] for x in Cluster_Pairs_List_0_and_1]
    nhood_enrichment_zscore_values=[]
    for x in Cluster_Pairs_List_0_and_1:
        try:
            nhood_enrichment_zscore_values.append(Dict_NhoodEnrichment_Zscore_aggregated[x])
        except:
            pass
    dict_nhood_enrichment_zscore=dict(zip(Cluster_Pairs_List_0_and_1, nhood_enrichment_zscore_values))
    # ---
    # nhood_enrichment_count_values=[Dict_NhoodEnrichment_Count_aggregated[x] for x in Cluster_Pairs_List_0_and_1]
    nhood_enrichment_count_values=[]
    for x in Cluster_Pairs_List_0_and_1:
        try:
            nhood_enrichment_zscore_values.append(Dict_NhoodEnrichment_Count_aggregated[x])
        except:
            pass
    dict_nhood_enrichment_count=dict(zip(Cluster_Pairs_List_0_and_1, nhood_enrichment_count_values))
    # ---------------
    # nhood_enrichment_zscore_0_values=[Dict_NhoodEnrichment_Zscore_0_aggregated[x] for x in Cluster_Pairs_List_0_and_1]
    nhood_enrichment_zscore_0_values=[]
    for x in Cluster_Pairs_List_0_and_1:
        try:
            nhood_enrichment_zscore_0_values.append(Dict_NhoodEnrichment_Zscore_0_aggregated[x])
        except:
            pass
    dict_nhood_enrichment_zscore_0=dict(zip(Cluster_Pairs_List_0_and_1, nhood_enrichment_zscore_0_values))
    # ---
    # nhood_enrichment_count_0_values=[Dict_NhoodEnrichment_Count_0_aggregated[x] for x in Cluster_Pairs_List_0_and_1]
    nhood_enrichment_count_0_values=[]
    for x in Cluster_Pairs_List_0_and_1:
        try:
            nhood_enrichment_count_0_values.append(Dict_NhoodEnrichment_Count_0_aggregated[x])
        except:
            pass
    dict_nhood_enrichment_count_0=dict(zip(Cluster_Pairs_List_0_and_1, nhood_enrichment_count_0_values))
    # ---------------
    # nhood_enrichment_zscore_1_values=[Dict_NhoodEnrichment_Zscore_1_aggregated[x] for x in Cluster_Pairs_List_0_and_1]
    nhood_enrichment_zscore_1_values=[]
    for x in Cluster_Pairs_List_0_and_1:
        try:
            nhood_enrichment_zscore_1_values.append(Dict_NhoodEnrichment_Zscore_1_aggregated[x])
        except:
            pass
    dict_nhood_enrichment_zscore_1=dict(zip(Cluster_Pairs_List_0_and_1, nhood_enrichment_zscore_1_values))
    # ---
    # nhood_enrichment_count_1_values=[Dict_NhoodEnrichment_Count_1_aggregated[x] for x in Cluster_Pairs_List_0_and_1]
    nhood_enrichment_count_1_values=[]
    for x in Cluster_Pairs_List_0_and_1:
        try:
            nhood_enrichment_count_1_values.append(Dict_NhoodEnrichment_Count_1_aggregated[x])
        except:
            pass
    dict_nhood_enrichment_count_1=dict(zip(Cluster_Pairs_List_0_and_1, nhood_enrichment_count_1_values))
    # ========================================================================================
    pvals_zscore_nhoodEnrichment_values=[]
    pvals_count_nhoodEnrichment_values=[]
    # ---
    for i in Cluster_Pairs_List_0_and_1:
        U1, p = mannwhitneyu(dict_nhood_enrichment_zscore_0[i], dict_nhood_enrichment_zscore_1[i], method="exact")
        pvals_zscore_nhoodEnrichment_values.append(p)
        # ---
        try:
            U1, p = mannwhitneyu(dict_nhood_enrichment_count_0[i], dict_nhood_enrichment_count_1[i], method="exact")
            pvals_count_nhoodEnrichment_values.append(p)
        except:
            pass
    # ---
    pvals_zscore_nhoodEnrichment=dict(zip(Cluster_Pairs_List_0_and_1, pvals_zscore_nhoodEnrichment_values))
    pvals_count_nhoodEnrichment=dict(zip(Cluster_Pairs_List_0_and_1, pvals_count_nhoodEnrichment_values))
    # ========================================================================================
    significant_pvals_zscore_nhoodEnrichment = dict((k, v) for k, v in pvals_zscore_nhoodEnrichment.items() if v < 0.05)
    significant_pvals_zscore_nhoodEnrichment_sorted = dict(sorted(significant_pvals_zscore_nhoodEnrichment.items(), key=lambda x:x[1]))

    significant_pvals_count_nhoodEnrichment = dict((k, v) for k, v in pvals_count_nhoodEnrichment.items() if v < 0.05)
    significant_pvals_count_nhoodEnrichment_sorted = dict(sorted(significant_pvals_count_nhoodEnrichment.items(), key=lambda x:x[1]))
    # ========================================================================================
    if len(significant_pvals_zscore_nhoodEnrichment)==0:
        print('No significant celltype pair found in terms of neighborhodd enrichment!')
    else:
        # ---
        from IPython.display import set_matplotlib_formats
        # set_matplotlib_formats('retina', quality=100)
        import seaborn as sns
        # sns.set_theme()
        sns.set_theme(style="whitegrid")
        # plt.rcdefaults()
        # ---
        # Set default figure size.
        plt.rcParams['figure.figsize'] = (12, 5)
        fig, axes = plt.subplots(figsize=(12, 5))
        # fig.suptitle(f'Differential Analysis (Neighborhood Enrichment) – celltype pairs with significant neighborhood enrichment differences (p-value<0.05)')
        feature = list(significant_pvals_zscore_nhoodEnrichment_sorted.keys())
        # -----
        _feature_=[]
        for i in feature:
            _feature_.append((cellTypeName_shortForm_dict[i[0]], cellTypeName_shortForm_dict[i[1]]))
        # -----
        score = list(significant_pvals_zscore_nhoodEnrichment_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        # Axis formatting.
        axes.spines['top'].set_visible(False)
        axes.spines['right'].set_visible(False)
        axes.spines['left'].set_visible(False)
        axes.spines['bottom'].set_color('#DDDDDD')
        axes.tick_params(bottom=False, left=False)
        axes.set_axisbelow(True)
        axes.yaxis.grid(True, color='#EEEEEE')
        axes.xaxis.grid(False)
        # ---
        axes.tick_params(labelsize = 14)
        # Add a footnote below and to the right side of the chart
        # axes.annotate('Footnote added below the chart with a smaller font',
        #             xy = (1.0, -0.2),
        #             xycoords='axes fraction',
        #             ha='right',
        #             va="center",
        #             fontsize=10)
        # plt.show()
        # ---
        # textstr = '\n'.join((
        # r'$\mu=%.2f$' % (mu, ),
        # r'$\mathrm{median}=%.2f$' % (median, ),
        # r'$\sigma=%.2f$' % (sigma, )))
        # ---
        
        # textstr = '\n'.join((
        # r'$\mu=%.2f$',
        # r'$\mathrm{median}=%.2f$',
        # r'$\sigma=%.2f$'))
        
        _string_ = str(cellTypeName_shortForm_dict)
        
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        
        # # axes.text(0.05, 0.95, textstr, transform=axes.transAxes, fontsize=14,
        # #     verticalalignment='top', bbox=props)
        
        # axes.text(0.05, 0.95, textstr, transform=axes.transAxes, fontsize=14,
        #     verticalalignment='top', bbox=props)
        
        my_text = ''
        my_text += cellTypeName_shortForm_str
        props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
        axes.text(1.03, 0.98, my_text, transform=axes.transAxes, fontsize=12, verticalalignment='top', bbox=props)
        
        # ---
        bars=plt.bar(x_pos, score,align='center')
        # ---
        # Grab the color of the bars so we can make the
        # text the same color.
        bar_color = bars[0].get_facecolor()

        # Add text annotations to the top of the bars.
        # Note, you'll have to adjust this slightly (the 0.3)
        # with different data.
        for bar in bars:
            axes.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0,
                round(score[int(bar.get_x() + bar.get_width() / 2)],3),
                # round(bar.get_height(), 1),
                horizontalalignment='center',
                color=bar_color,
                weight='bold'
            )
        # --- Ticks: ---
        # plt.xticks(x_pos, feature, rotation=90)
        plt.xticks(x_pos, _feature_, rotation=90)
        # --------------
        # --- Labels: ---
        # plt.ylabel('p-value (MWU test)')
        # ---------------
        # ---
        # Add labels and a title. Note the use of `labelpad` and `pad` to add some
        # extra space between the text and the tick labels.
        axes.set_xlabel('Celltype pairs', labelpad=15, color='#333333')
        axes.set_ylabel('p-value (MWU)', labelpad=15, color='#333333')
        axes.set_title('Neighborhood enrichment\n$(celltype\;pairs\;with\;p < 0.05)$', pad=15, color='#333333', weight='bold')
        # ---
        fig.tight_layout()
        # fig.constrained_layout()
        plt.savefig(f'DA_neighborhoodEnrichment_pvalues_significantCelltypePairs_allPatients.pdf', format='pdf', bbox_inches='tight')
        plt.show()