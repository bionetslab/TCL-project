


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

def plot_of_robustness_vs_intercluster_similarity(adata_pickle_path, method):
    
    with open(adata_pickle_path, 'rb') as f:
        pickle_= pickle.load(f)


    patient_ids=[]
    for i in pickle_:
        patient_ids.append(i)


    cnt=-1

    # ------
    pickle_items=list(pickle_.items())
    first_three_items = pickle_items[0:2]
    pickle_new=dict(first_three_items)
    pickle_=pickle_new
    # ------
    
    no_of_patients=np.shape(pickle_)
    error_samples=[]


    clusterings_patientLevel_values=[[] for _ in range(len(patient_ids))]
    clusterings_patientLevel_keys=patient_ids
    clusterings_patientLevel_dict=dict(zip(clusterings_patientLevel_keys, clusterings_patientLevel_values))


    for i in pickle_:
        
        print('Patient serial number' + str(cnt) + ' (patient id: ' + str(i) + ')')
        print('====================================================================')
        
        cnt+=1
        
        adata=pickle_[i]
        adata_X=adata.to_df()
        
        
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
        adata_raw=adata.copy()
        sc.tl.pca(adata_raw, svd_solver='arpack')
        sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=3)
        
        
        no_of_clusterings=3
        clustering_names=[]
        
        for clustering in range(no_of_clusterings):
            cl_names='Clustering-'+str(clustering)
            clustering_names.append(cl_names)
        
        clusterings_values=[[] for _ in range(len(clustering_names))]
        clusterings_keys=clustering_names
        clusterings_dict=dict(zip(clusterings_keys, clusterings_values))
        
        for j in range(no_of_clusterings):
            
            print('no_of_clusterings= ' + str(j) + ' in patient serial number '+str(cnt)+' (patient id: ' + str(i) + ')')
            print('-------------------------------------------------------------------------')
            
            scs.inference.nested_model(adata_raw)
            uns_keys=str(adata_raw.uns_keys)
            cluster_level_max=uns_keys[uns_keys.rindex('CM_nsbm_level')+len('CM_nsbm_levelx')]
            cluster_level_max_int=int(cluster_level_max)
            actual_cluster_levels=[format(x, 'd') for x in list(range(1, cluster_level_max_int))]
            
            
            cluster_level_values=[[] for _ in range(len(actual_cluster_levels))]
            cluster_level_keys=actual_cluster_levels     
            cluster_level_dict=dict(zip(cluster_level_keys, cluster_level_values))
            
            
            
            for k in range(1, cluster_level_max_int):
                
                print('Cluster number= '+str(k)+' in '+'no_of_clusterings= ' + str(j) + ' in patient serial number '+str(cnt)+' (patient id: ' + str(i) + ')')
                print('-------------               -----------------          ----------------------')
                
                cluster_name='nsbm_level_'+str(k)
                ppbm_clusters_sorted=[str(i1) for i1 in sorted([int(i2) for i2 in np.unique(adata_raw.obs[cluster_name])])]
                
                
                
                ppbm_no_of_clusters_dict_values=[[] for _ in range(len(ppbm_clusters_sorted))]
                ppbm_no_of_clusters_dict_keys=ppbm_clusters_sorted
                ppbm_no_of_clusters_dict=dict(zip(ppbm_no_of_clusters_dict_keys, ppbm_no_of_clusters_dict_values))
                
                df_clusters=pd.DataFrame(adata_raw.obs[cluster_name])
                
                for index, row in df_clusters.iterrows():
                    ppbm_no_of_clusters_dict[row[cluster_name]].append(index)
            
            
                cluster_level_dict[str(k)].append(ppbm_no_of_clusters_dict)
            
            
            
            clusterings_dict[clustering_names[j]].append(cluster_level_dict)
        
        clusterings_patientLevel_dict[i].append(clusterings_dict)
        
    with open('clusterings_patientLevel_dict.pkl', 'wb') as f:
        pickle.dump(clusterings_patientLevel_dict, f)    
        

    # =============================================================================================================
    # =============================================================================================================
    # =============================================================================================================
    # =============================================================================================================
    # =============================================================================================================

    # ---
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
        

    with open("clusterings_patientLevel_dict.pkl", "rb") as f:
        clusterings_patientLevel_dict = pickle.load(f)

    clusterings_patientLevel_dict_items=list(clusterings_patientLevel_dict.items())
    first_three_items = clusterings_patientLevel_dict_items[0:1]
    clusterings_patientLevel_dict_new=dict(first_three_items)
    clusterings_patientLevel_dict=clusterings_patientLevel_dict_new


    with open(adata_pickle_path, 'rb') as f:
        pickle_= pickle.load(f)

    cnt=-1

    # ------
    pickle_items=list(pickle_.items())
    first_three_items = pickle_items[0:1]
    pickle_new=dict(first_three_items)
    pickle_=pickle_new
    # ------

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
        
        
        clustering_combinations=list(itertools.combinations(clustering_names, 2))
        
        clust_comb_keys=[]
        
        
        # cluster_level_values=[[] for _ in range(len(actual_cluster_levels))]
        
        # _scores_=[]
        
        clust_comb_count=-1
        for clust_comb_iter in clustering_combinations:
            clust_comb_count+=1
            clust_comb_keys.append(clust_comb_iter)
            B1_dict=clusterings_patientLevel_dict[i][0][clust_comb_iter[0]][0]
            B1_dict_copy=B1_dict.copy()
            B2_dict=clusterings_patientLevel_dict[i][0][clust_comb_iter[1]][0]
            B2_dict_copy=B2_dict.copy()
            equal_clusters=0
            if len(set(B1_dict.keys())) == len(set(B2_dict.keys())):
                equal_clusters=1
                common_keys=list(set(B1_dict.keys()).intersection(set(B2_dict.keys())))
                common_keys=[str(i1) for i1 in sorted([int(i2) for i2 in np.unique(common_keys)])]
            
                # ------------------------------------------
            
                list_1=[]
                list_2=[]
                clustering_level=[]
                if equal_clusters==1:
                    _score_=[]
                    for j in range(len(common_keys)):
                        list_1.append(common_keys[j])
                        list_2.append(common_keys[j])
                        clustering_level.append(j)
                        # ---
                        list1=list(B1_dict[list_1[j]][0].keys())
                        list2=list(B2_dict[list_2[j]][0].keys())
                        unique_combinations=list(itertools.product(list1, list2))
                        # ---  
                        edge_weights=[]
                        B1_nodes=[]
                        B2_nodes=[]
                        B_nodes=[]
                        Edge_weights=[]
                        All_edge_attributes=[]
                        for k in common_keys:
                            b1_node='B1_'+str(k)
                            B1_nodes.append(b1_node)
                            b2_node='B2_'+str(k)
                            B2_nodes.append(b2_node)
                            B_nodes.append(b1_node)
                            B_nodes.append(b2_node)
                            # ---
                            b1_node_set=set(list(clusterings_patientLevel_dict[i][0][clust_comb_iter[0]][0][k][0].keys()))
                            b2_node_set=set(list(clusterings_patientLevel_dict[i][0][clust_comb_iter[1]][0][k][0].keys()))
                            b1_b2_intersection=b1_node_set.intersection(b2_node_set)
                            # ---
                            edge_weight=len(list(b1_b2_intersection))
                            Edge_weights.append(edge_weight)
                            all_edge_attributes=(b1_node, b2_node, {'weight': edge_weight})
                            All_edge_attributes.append(all_edge_attributes)
                            # ---
                        B1_nodes=list(np.unique(B1_nodes))
                        B2_nodes=list(np.unique(B2_nodes))
                        B_nodes=list(np.unique(B_nodes))
                        # ---
                        B=nx.Graph()
                        nx.draw(B)
                        B.add_nodes_from(B1_nodes, bipartite=0)
                        B.add_nodes_from(B2_nodes, bipartite=1)
                        B.add_edges_from(All_edge_attributes)
                        # ---
                        max_bip_matchings=sorted(nx.max_weight_matching(B))
                        score=0
                        for edges_ in max_bip_matchings:
                            score+=B[edges_[0]][edges_[1]]["weight"]
                        score/=no_of_cells
                        print(score)
                        _score_.append(score)
            else:
                if len(set(B1_dict.keys())) < len(set(B2_dict.keys())):
                    common_keys_len=len(set(B1_dict.keys()))
                    lesser_hierarchical_clusters=B1_dict.copy()
                    higher_hierarchical_clusters=B2_dict.copy()
                else:
                    common_keys_len=len(set(B2_dict.keys()))
                    lesser_hierarchical_clusters=B2_dict.copy()
                    higher_hierarchical_clusters=B1_dict.copy()
            
                # ----------------------------------------------------------------------
                # ----------------------------------------------------------------------
                # ----------------------------------------------------------------------
                
                common_keys_lesser_hierarchical=[]
                common_keys_higher_hierarchical=[]
                for p in range(len(lesser_hierarchical_clusters)):
                    if method=='top-down':
                        common_keys_lesser_hierarchical.append(list(reversed(lesser_hierarchical_clusters))[p])
                        common_keys_higher_hierarchical.append(list(reversed(higher_hierarchical_clusters))[p])
                    elif method=='bottom-up':
                        common_keys_lesser_hierarchical.append(list(lesser_hierarchical_clusters)[p])
                        common_keys_higher_hierarchical.append(list(higher_hierarchical_clusters)[p])
                
                common_keys_lesser_hierarchical=[str(i1) for i1 in sorted([int(i2) for i2 in np.unique(common_keys_lesser_hierarchical)])]
                common_keys_higher_hierarchical=[str(i1) for i1 in sorted([int(i2) for i2 in np.unique(common_keys_higher_hierarchical)])]
                
                # ------------------------------------------------------------------------    
                # ------------------------------------------------------------------------
                # ------------------------------------------------------------------------
                
                lesser_clusters=[]
                for j in common_keys_lesser_hierarchical:
                    lesser_clusters_=list(lesser_hierarchical_clusters[j][0].keys())
                    lesser_clusters.append(lesser_clusters_)
                
                
                higher_clusters=[]
                for j in common_keys_higher_hierarchical:
                    higher_clusters_=list(higher_hierarchical_clusters[j][0].keys())
                    higher_clusters.append(higher_clusters_)
                
                
                # # # lesserCluster_higherCluster_dict=dict(zip(LESSER_CLUSTERS, HIGHER_CLUSTERS))
                
                for j in range(len(common_keys_lesser_hierarchical)):
                    # print(j)
                    len_lesser_clusters=len(lesser_clusters[j])
                    len_higher_clusters=len(higher_clusters[j])
                    
                    no_of_cluster_differences=len_lesser_clusters-len_higher_clusters
                    # print(no_of_cluster_differences)
                    
                    extra_cluster_names=[]
                    lesser_clusters_new=[]
                    higher_clusters_new=[]
                    error_different_no_of_clusters=0
                        
                    if no_of_cluster_differences==0:
                        # print(no_of_cluster_differences)
                            
                        pass
                    
                    else:
                        if no_of_cluster_differences!=0:
                            if no_of_cluster_differences>0:
                                for l in range(no_of_cluster_differences):
                                    # print(l)
                                    error_different_no_of_clusters=1
                                    extra_cluster_names.append('extra_'+str(int(l)))
                                    higher_clusters_new.append('extra_'+str(int(l)))
                                    new_='extra_'+str(int(l))
                                    # higher_hierarchical_clusters[list(higher_hierarchical_clusters.keys())[j]][0][new_]=[]
                                    higher_hierarchical_clusters[common_keys_higher_hierarchical[j]][0][new_]=[]
                            else:
                                for l in range(-no_of_cluster_differences):
                                    error_different_no_of_clusters=2
                                    extra_cluster_names.append('extra_'+str(int(l)))
                                    lesser_clusters_new.append('extra_'+str(int(l)))
                                    new_='extra_'+str(int(l))
                                    # lesser_hierarchical_clusters[list(lesser_hierarchical_clusters.keys())[j]][0][new_]=[]
                                    # # common_keys_higher_hierarchical
                                    lesser_hierarchical_clusters[common_keys_lesser_hierarchical[j]][0][new_]=[]
                            
                            
                # ----------------------------  Combinations:  ----------------------------------------
            
            
            list_1=[]
            list_2=[]
            clustering_level=[]
            if equal_clusters==0:
                _score_=[]
                for j in range(len(common_keys_lesser_hierarchical)):
                    list_1.append(common_keys_lesser_hierarchical[j])
                    list_2.append(common_keys_higher_hierarchical[j])
                    clustering_level.append(j)
                    
                    list1=list(lesser_hierarchical_clusters[list_1[j]][0].keys())
                    list2=list(higher_hierarchical_clusters[list_2[j]][0].keys())
                    unique_combinations=list(itertools.product(list1, list2))
                    
                    edge_weights=[]
                    B1_nodes=[]
                    B2_nodes=[]
                    B_nodes=[]
                    Edge_weights=[]
                    All_edge_attributes=[]
                    for k in unique_combinations:
                        b1_node='B1_'+str(k[0])
                        B1_nodes.append(b1_node)
                        b2_node='B2_'+str(k[1])
                        B2_nodes.append(b2_node)
                        B_nodes.append(b1_node)
                        B_nodes.append(b2_node)
                        # ---
                        b1_node_set=set(list(lesser_hierarchical_clusters[list_1[j]][0][k[0]]))
                        b2_node_set=set(list(higher_hierarchical_clusters[list_2[j]][0][k[1]]))
                        b1_b2_intersection=b1_node_set.intersection(b2_node_set)
                        # ---
                        edge_weight=len(list(b1_b2_intersection))
                        Edge_weights.append(edge_weight)
                        all_edge_attributes=(b1_node, b2_node, {'weight': edge_weight})
                        All_edge_attributes.append(all_edge_attributes)
                        # ---
                    B1_nodes=list(np.unique(B1_nodes))
                    B2_nodes=list(np.unique(B2_nodes))
                    B_nodes=list(np.unique(B_nodes))
                    # ---
                    B=nx.Graph()
                    nx.draw(B)
                    B.add_nodes_from(B1_nodes, bipartite=0)
                    B.add_nodes_from(B2_nodes, bipartite=1)
                    B.add_edges_from(All_edge_attributes)
                    # ---
                    # score=0
                    # for u, v, d in B.edges(data=True):
                    #     # print(f"({u}, {v}) {d=}")
                    #     score+=d['weight']
                    #     print(score)
                    # score/=no_of_cells
                    # ---
                    max_bip_matchings=sorted(nx.max_weight_matching(B))
                    score=0
                    for edges_ in max_bip_matchings:
                        score+=B[edges_[0]][edges_[1]]["weight"]
                    score/=no_of_cells
                    print(score)
                    _score_.append(score)
            # ================================================
            SCORE.append(_score_)
            
    # ================================================
    # ================================================
    # ================================================

    lst = SCORE
    _, max_length=FindMaxLength(lst)

    scores_dict={}
    df_names={}
    for i in range(max_length):
        scores_dict[i]=[]
        df_names[i]='df'+str(i)

    logScores_list=[]
    labels_list=[]
    for i in SCORE:
        for j in range(len(i)):
            scores_dict[j].append(np.log10(i[j]))
            logScores_list.append(np.log10(i[j]))
            labels_list.append(int(j+1))
    df = pd.DataFrame(list(zip(labels_list, logScores_list)), columns=['Robustness (cluster level)', 'Inter-cluster similarity score'])
    clusterLevels_unique=list(np.unique(df['Robustness (cluster level)'].values.tolist()))
    # logScores_list=[]
    # labels_list=[]
    # for i in range(max_length):
    #     # eval_str='df'+str(i)+'pd.DataFrame(scores_dict[i])'
    #     # evalStr_logScores='pd.DataFrame(scores_dict[i])'
    #     # evalStr_dfNames='df'+str(i)
    #     logScores_list.append(pd.DataFrame(scores_dict[i]))

    # ---

    # fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    # fig.suptitle(f'Plot of robustness vs. inter-cluster similarity (method: {method})')

    # no_of_cols=3
    # no_of_rows=int(int(max_length)/no_of_cols)+1

    # _cnt_=-1
    # for i in range(no_of_rows):
    #     for j in range(no_of_cols):
    #         _cnt_+=1
    #         if _cnt_<max_length:
    #             sns.boxplot(ax=axes[int(i), int(j)], data=df[df['Robustness (cluster level)']==clusterLevels_unique[_cnt_]], x='Robustness (cluster level)', y='Inter-cluster similarity score')

    # plt.savefig(f"Plot of robustness vs. inter-cluster similarity (method: {method}).svg", format='svg')

    # ---

    fig, axes = plt.subplots(1, 1, figsize=(18, 10))
    fig.suptitle(f'Plot of robustness vs. inter-cluster similarity (method: {method})')

    no_of_cols=3
    no_of_rows=int(int(max_length)/no_of_cols)+1

    sns.boxplot(ax=axes, data=df, x='Robustness (cluster level)', y='Inter-cluster similarity score')
    plt.savefig(f"Plot of robustness vs. inter-cluster similarity (method: {method}).jpg", format='jpg')




