


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
import os.path
import warnings
import pandas as pd

def categorize(dataframe, leiden_clusters_sorted):
    df = pd.DataFrame(columns=['group'])
    for i in leiden_clusters_sorted:
        df_=dataframe[dataframe['group']==i]
        df_=df_.sort_values(by = 'pvals_adj', ascending=True)
        df=pd.concat([df, df_], ignore_index=True)
    return df

def calculate_average_expressions_per_protein(adata_X, cells_in_cluster, proteins_in_cluster):
    n=len(cells_in_cluster)
    DF=adata_X.loc[cells_in_cluster]
    DF=DF[proteins_in_cluster]
    per_row_sum_all_proteins=DF.sum(axis=1).tolist()
    per_row_indiv_proteins=[]
    per_row_sum_all_proteins_excluding_indiv=[]
    per_row_avg_all_proteins_excluding_indiv=[]
    list_of_proteins=DF.columns.tolist()
    for per_protein in list_of_proteins:
        indiv_protein=DF[per_protein].tolist()
        per_row_indiv_proteins.append(indiv_protein)
        sub_=list(np.subtract(per_row_sum_all_proteins, indiv_protein))
        per_row_sum_all_proteins_excluding_indiv.append(sub_)
        avg_=[]
        for val in sub_:
            avg_.append(val/float(n-1))
        per_row_avg_all_proteins_excluding_indiv.append(avg_)
    return per_row_avg_all_proteins_excluding_indiv, per_row_indiv_proteins

def calculate_average_expressions_per_protein2(adata_X, cells_in_cluster, proteins_in_cluster):
    DF=adata_X.copy()
    DF_in_cluster=DF.loc[cells_in_cluster]
    DF_rest=DF[~DF.index.isin(cells_in_cluster)]
    
    average_simulated_pvalues=[]
    actual_pvalues=[]
    for i in proteins_in_cluster:
        DF_shuffled=DF.copy()
        # -----------------------------
        DF_shuffled[i] = DF_shuffled[i].sample(frac=1).values
        DF_shuffled_in_cluster=DF_shuffled.loc[cells_in_cluster]
        DF_shuffled_rest=DF_shuffled[~DF.index.isin(cells_in_cluster)]
        
        actualSampleOne=list(DF_in_cluster[i])
        actualSampleRest=list(DF_rest[i])
        # U1_actual, p_actual=mannwhitneyu(actualSampleOne, actualSampleRest, method="exact")
        p_actual=mwu(actualSampleOne, actualSampleRest)['p-val'][0]
        actual_pvalues.append(p_actual)
        
        simulated_pvals=[]
        no_of_randomizations=3
        for j in range(no_of_randomizations):
            # print(j)
            sampleOne=list(DF_shuffled_in_cluster[i])
            sampleRest=list(DF_shuffled_rest[i])
            p=mwu(sampleOne, sampleRest)['p-val'][0]
            simulated_pvals.append(p)
        average_simulated_pval=np.mean(simulated_pvals)
        average_simulated_pvalues.append(average_simulated_pval)
    average_simulated_pvalues_dict=dict(zip(proteins_in_cluster, average_simulated_pvalues))
    actual_pvalues_dict=dict(zip(proteins_in_cluster, actual_pvalues))
    
    return actual_pvalues_dict, average_simulated_pvalues_dict
            
            
# def run(adata_pickle_path, user_selected_cluster_level):
def run(adata_pickle_path, user_selected_cluster_level):

    with open(adata_pickle_path, 'rb') as f:
        pickle_= pickle.load(f)
    res = list(pickle_.keys())[0]

    cnt=-1
    
    # ------
    pickle_items=list(pickle_.items())
    # first_three_items = pickle_items[0:3]
    # pickle_new=dict(first_three_items)
    # pickle_=pickle_new
    # ---
    all_items = pickle_items[:]
    pickle_new=dict(all_items)
    pickle_=pickle_new
    # ------
    essential_proteins_allPatients={}
    essentialProteinsPerCluster_acrossClusterLevels_forAllPatients={}
    allEssentialProteins_acrossClusterLevels_forAllPatients={}
    X_allPatients={}
    y_allPatients={}
    # -------

    for i in pickle_:
        cnt+=1
        print('Patient serial number' + str(cnt) + ' (patient id: ' + str(i) + ')')
        print('====================================================================')
        adata=pickle_[i]
        adata_X=adata.to_df()
        X_allPatients[i]=np.array(adata_X)
        adata_raw=adata.copy()
        sc.pp.neighbors(adata_raw)
        # ===========================================================================================
        scs.inference.nested_model(adata_raw)
        # ---
        uns_keys=str(adata_raw.uns_keys)
        cluster_level_max=uns_keys[uns_keys.rindex('CM_nsbm_level')+len('CM_nsbm_levelx')]
        cluster_level_max_int=int(cluster_level_max)
        actual_cluster_levels=[format(x, 'd') for x in list(range(1, cluster_level_max_int-1))]
        # ---
        cluster_level_values=[[] for _ in range(len(actual_cluster_levels))]
        cluster_level_keys=actual_cluster_levels
        cluster_level_dict=dict(zip(cluster_level_keys, cluster_level_values))
        # ===========================================================================================
        essentialProteinsPerCluster_acrossClusterLevels={}
        allEssentialProteins_acrossClusterLevels={}
    	# ===========================================================================================
        if user_selected_cluster_level=="complete":
            selected_cluster_level=cluster_level_max_int-2
        else:
            try:
                selected_cluster_level=int(user_selected_cluster_level)
            except:
                selected_cluster_level=cluster_level_max_int-2
        # ---
        cluster_level_range=list(range(1, cluster_level_max_int-1))
        if selected_cluster_level==0:
            selected_cluster_level=cluster_level_max_int-2
        elif selected_cluster_level>0:
            if not(selected_cluster_level in cluster_level_range):
                selected_cluster_level=cluster_level_max_int-2
        elif selected_cluster_level<0:
            selected_cluster_level=cluster_level_max_int-2-selected_cluster_level
            if not(selected_cluster_level in cluster_level_range):
                selected_cluster_level=cluster_level_max_int-2
        # ===========================================================================================
        
        for k in range(1, cluster_level_max_int-1):
            print('Cluster number= '+str(k)+' in '+' in patient serial number '+str(cnt)+' (patient id: ' + str(i) + ')')
            print('-------------               -----------------          ----------------------')
            # ---
            cluster_name='nsbm_level_'+str(k)
            ppbm_clusters_sorted=[str(i1) for i1 in sorted([int(i2) for i2 in np.unique(adata_raw.obs[cluster_name])])]
            # ---
            ppbm_no_of_clusters_dict_values=[[] for _ in range(len(ppbm_clusters_sorted))]
            ppbm_no_of_clusters_dict_keys=ppbm_clusters_sorted
            ppbm_no_of_clusters_dict=dict(zip(ppbm_no_of_clusters_dict_keys, ppbm_no_of_clusters_dict_values))
            # ---
            df_clusters=pd.DataFrame(adata_raw.obs[cluster_name])
            if k==selected_cluster_level:
                y_allPatients[i]=np.array(df_clusters)
            for index, row in df_clusters.iterrows():
                ppbm_no_of_clusters_dict[row[cluster_name]].append(index)
            # ---
            cluster_level_dict[str(k)].append(ppbm_no_of_clusters_dict)
            _str_='nsbm_level_'+str(k)
        	# -----
            gene_ranking_error=1
            while gene_ranking_error==1:
                try:
                    sc.tl.rank_genes_groups(adata_raw, groupby=_str_, pts=True, method='wilcoxon')
                    gene_ranking_error=0
                except:
                    pass
            # -----
            DEDF_ppbm=sc.get.rank_genes_groups_df(adata_raw, group=ppbm_clusters_sorted)
            DEDF_ppbm_AscPVals=categorize(DEDF_ppbm, ppbm_clusters_sorted)
            sc.pl.rank_genes_groups(adata_raw, sharey=False, save='sc.pl.rank_genes_groups_PPBMClustering.jpg')
            # -----
            DEDF_ppbm_AscPVals_0Point05_cutoff=DEDF_ppbm_AscPVals[DEDF_ppbm_AscPVals.pvals_adj<0.05]
            DEDF_ppbm_AscPVals_0Point05_cutoff_groups=list(np.unique(DEDF_ppbm_AscPVals_0Point05_cutoff.group))
            DEDF_ppbm_AscPVals_0Point05_cutoff_groups_sorted=[str(j) for j in sorted([int(i) for i in np.unique(DEDF_ppbm_AscPVals_0Point05_cutoff_groups)])]
            # -----
            ppbm_no_of_clusters_dict_values=[[] for _ in range(len(DEDF_ppbm_AscPVals_0Point05_cutoff_groups_sorted))]
            ppbm_no_of_clusters_dict_keys=DEDF_ppbm_AscPVals_0Point05_cutoff_groups_sorted
            ppbm_no_of_clusters_dict=dict(zip(ppbm_no_of_clusters_dict_keys, ppbm_no_of_clusters_dict_values))
            for index, row in DEDF_ppbm_AscPVals_0Point05_cutoff.iterrows():
                ppbm_no_of_clusters_dict[row['group']].append(row['names'])
            # -----
            df=adata_raw.obs[_str_].copy().to_frame()
            # ppbm_cell_lookup_keys=[str(j) for j in sorted([int(i) for i in np.unique(df.index)])]
            ppbm_cell_lookup_values=[[] for _ in range(len(ppbm_clusters_sorted))]
            # ppbm_cell_lookup_keys=[[] for _ in range(len(np.unique(ppbm_cell_lookup_keys)))]
            ppbm_cell_lookup_keys=ppbm_clusters_sorted
            ppbm_cell_lookup_dict=dict(zip(ppbm_cell_lookup_keys, ppbm_cell_lookup_values))
            for index, row in df.iterrows():
                ppbm_cell_lookup_dict[row[_str_]].append(index)
            # -----
            _all_essential_proteins_=[]
            for key in ppbm_no_of_clusters_dict:
                print('PPBM cluster ' + str(key) + ' in patient serial number ' + str(cnt) + ' (patient id: ' + str(i) + ')')
                proteins_in_cluster=ppbm_no_of_clusters_dict[key]
                cells_in_cluster=ppbm_cell_lookup_dict[key]
                # ----------------------
                actual_pvalue_in_cluster, simulated_p_value_dict_in_cluster=calculate_average_expressions_per_protein2(adata_X, cells_in_cluster, proteins_in_cluster)
                essential_proteins_in_cluster=[]
                # -----
                for pv in actual_pvalue_in_cluster:
                    if (actual_pvalue_in_cluster[pv] < simulated_p_value_dict_in_cluster[pv]):
                        essential_proteins_in_cluster.append(pv)
                _all_essential_proteins_.append(list(np.unique(essential_proteins_in_cluster)))
            
            all_essential_proteins_consolidated=[]
            for j in _all_essential_proteins_:
                for p in j:
                    all_essential_proteins_consolidated.append(p)
            all_essential_proteins_=list(np.unique(all_essential_proteins_consolidated))
            _essential_proteins_allPatients_=all_essential_proteins_
        
            # ------------------------
            
            essentialProteinsPerCluster_acrossClusterLevels[k]=_all_essential_proteins_
            allEssentialProteins_acrossClusterLevels[k]=_essential_proteins_allPatients_
        
            # ------------------------
        
            # selected_cluster_level=1
            selected_cluster_level=cluster_level_max_int-2
        
            if k==selected_cluster_level:
                all_essential_proteins=_all_essential_proteins_.copy()
                essential_proteins_allPatients[i]=_essential_proteins_allPatients_
    
        # ==============================================================================
    
        essentialProteinsPerCluster_acrossClusterLevels_forAllPatients[i]=essentialProteinsPerCluster_acrossClusterLevels
        allEssentialProteins_acrossClusterLevels_forAllPatients[i]=allEssentialProteins_acrossClusterLevels
    
        # ==============================================================================
    	
    	
        import math
    
        # ----------------------------------------------------------------------------------------------------------------------------------------------------
        # no_of_cols=3
        # no_of_rows=int(math.ceil(float((cluster_level_max_int-2)-1)/no_of_cols))
        # ---
        # # no_of_cols=1
        # # no_of_rows=cluster_level_max_int-2
        # ---
        # fig, axes = plt.subplots(no_of_rows, no_of_cols, figsize=(50, 50))
        # fig.suptitle('No. of cluster levels vs. No. of clusters per cluster level')
        # ----------------------------------------------------------------------------------------------------------------------------------------------------

        _cnt_=-1
        for k in range(1, cluster_level_max_int-1):
            print(k)
            _cnt_+=1
            print(_cnt_)
        
            # ----------------------------------------------------------------------------------------------------------------------------------------------------
            # row=int(math.floor(float(_cnt_)/no_of_cols))
            # col=int(_cnt_%no_of_cols)
            # print(row, col)
            # if k==cluster_level_max_int-2:
            #     if no_of_rows>1:
            #         _ax_=axes[int(row), int(col)]
            #         _clust_name_='nsbm_level_'+str(k)
            #         sc.pl.dotplot(adata_raw, allEssentialProteins_acrossClusterLevels_forAllPatients[i][1], groupby=_clust_name_,ax=_ax_, save="No. of cluster levels vs. No. of clusters per cluster level.jpg");
            #         _ax_.set_xlabel('No. of clustering levels')
            #         _ax_.set_ylabel('No. of clusters per level')
            #         _ax_.set_title(f'No. of clustering levels: {i}')
            #     elif no_of_rows==1:
            #         _ax_=axes[int(col)]
            #         _clust_name_='nsbm_level_'+str(k)
            #         sc.pl.dotplot(adata_raw, allEssentialProteins_acrossClusterLevels_forAllPatients[i][1], groupby=_clust_name_,ax=_ax_, save="No. of cluster levels vs. No. of clusters per cluster level.jpg");
            #         _ax_.set_xlabel('No. of clustering levels')
            #         _ax_.set_ylabel('No. of clusters per level')
            #         _ax_.set_title(f'No. of clustering levels: {i}')
            # else:
            #     if no_of_rows>1:
            #         _ax_=axes[int(row), int(col)]
            #         _clust_name_='nsbm_level_'+str(k)
            #         sc.pl.dotplot(adata_raw, allEssentialProteins_acrossClusterLevels_forAllPatients[i][1], groupby=_clust_name_,ax=_ax_, show=False);
            #         _ax_.set_xlabel('No. of clustering levels')
            #         _ax_.set_ylabel('No. of clusters per level')
            #         _ax_.set_title(f'No. of clustering levels: {i}')
            #     elif no_of_rows==1:
            #         _ax_=axes[int(col)]
            #         _clust_name_='nsbm_level_'+str(k)
            #         sc.pl.dotplot(adata_raw, allEssentialProteins_acrossClusterLevels_forAllPatients[i][1], groupby=_clust_name_,ax=_ax_, show=False);
            #         _ax_.set_xlabel('No. of clustering levels')
            #         _ax_.set_ylabel('No. of clusters per level')
            #         _ax_.set_title(f'No. of clustering levels: {i}')
            # ----------------------------------------------------------------------------------------------------------------------------------------------------
            
            _clust_name_='nsbm_level_'+str(k)
            # ---
            x_axis_len=int(len(np.unique(allEssentialProteins_acrossClusterLevels_forAllPatients[i][k])))
            x_axis_plot_len=int(x_axis_len/2)
            if x_axis_plot_len<3:
                x_axis_plot_len=3
            # ---
            y_axis_len=int(len(np.unique(adata_raw.obs[_clust_name_])))
            y_axis_plot_len=int(y_axis_len*0.75)
            if y_axis_plot_len<3:
                y_axis_plot_len=3
            # ---
            fig, axes = plt.subplots(figsize=(x_axis_plot_len, y_axis_plot_len))
            fig.suptitle(f'\n\nDotplot of potential marker genes vs. cell type clusters - cluster level {k}')
            if k==selected_cluster_level:
                sc.pl.dotplot(adata_raw, allEssentialProteins_acrossClusterLevels_forAllPatients[i][k], groupby=_clust_name_, ax=axes, save=f"____PatientID-{i}_____ClusteringLevel-{k}")
            else:
                sc.pl.dotplot(adata_raw, allEssentialProteins_acrossClusterLevels_forAllPatients[i][k], groupby=_clust_name_, ax=axes)
            
        
            # ----------------------------------------------------------------------------------------------------------------------------------------------------
            # axes.set_xlabel('No. of clustering levels')
            # axes.set_ylabel('No. of clusters per level')
            # axes.set_title(f'No. of clustering levels: {i}')
            # ----------------------------------------------------------------------------------------------------------------------------------------------------
        
        # plt.savefig("No. of cluster levels vs. No. of clusters per cluster level.jpg", format='jpg')
        # plt.show()
        
        # ==============================================================================



    with open('essential_proteins.pkl', 'wb') as f:
        pickle.dump([essential_proteins_allPatients, essentialProteinsPerCluster_acrossClusterLevels_forAllPatients, allEssentialProteins_acrossClusterLevels_forAllPatients, X_allPatients, y_allPatients], f)


    # =============================================================================================
    # =============================================================================================
    # =============================================================================================
    # =============================================================================================
    # =============================================================================================


















































