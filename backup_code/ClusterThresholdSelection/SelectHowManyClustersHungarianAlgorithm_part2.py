


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
from neighborhood_enrichment import neighborhood_enrichment
import schist as scs
from pingouin import mwu
import graph_tools as gt
import itertools
from itertools import permutations



with open("clusterings_patientLevel_dict.pickle", "rb") as f:
    clusterings_patientLevel_dict = pickle.load(f)


with open('TNBC_41patients_KerenEtAl.pkl', 'rb') as f:
   pickle_= pickle.load(f)

cnt=-1

# ------
pickle_items=list(pickle_.items())
first_three_items = pickle_items[3:4]
pickle_new=dict(first_three_items)
pickle_=pickle_new
# ------


for i in pickle_:
    
    cnt+=1
    
    adata=pickle_[i]
    adata_X=adata.to_df()
    
    
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
    
    
    for clust_comb_iter in clustering_combinations:
        clust_comb_keys.append(clust_comb_iter)
        B1_dict=clusterings_patientLevel_dict[i][0][clust_comb_iter[0]][0]
        B1_dict_copy=B1_dict.copy()
        B2_dict=clusterings_patientLevel_dict[i][0][clust_comb_iter[1]][0]
        B2_dict_copy=B2_dict.copy()
        
        common_keys=list(set(B1_dict.keys()).intersection(set(B2_dict.keys())))
        common_keys = [eval(i) for i in common_keys]
        common_keys = np.sort(common_keys)
        common_keys=[format(x, 'd') for x in common_keys]
        
        keys_count=-1
        for j in common_keys:
            keys_count+=1
            
            B1_clusters=list(B1_dict[j][0].keys())
            B1_clusters_new=[]
            # for k, v in B1_dict[j][0].items():
            #     new_key='B1_'+str(k)
            #     B1_clusters_new[new_key]=v
            # del B1_
                
            
            B2_clusters=list(B2_dict[j][0].keys())
            B2_clusters_new=[]
            # for k, v in B2_dict.items():
            #     new_key='B2_'+str(k)
            #     B2_clusters_new[new_key]=v
            # B2_dict[j][0]=B2_clusters_new
            
            no_of_cluster_differences=len(B1_clusters)-len(B2_clusters)
            extra_cluster_names=[]
            error_different_no_of_clusters=0
            if no_of_cluster_differences==0:
                pass
            else:
                if no_of_cluster_differences>0:
                    for l in range(no_of_cluster_differences):
                        error_different_no_of_clusters=1
                        extra_cluster_names.append('B1_extra_'+str(int(l)))
                        B1_clusters_new.append('B1_extra_'+str(int(l)))
                else:
                    for l in range(no_of_cluster_differences):
                        error_different_no_of_clusters=2
                        extra_cluster_names.append('B2_extra_'+str(int(l)))
                        B2_clusters_new.append('B2_extra_'+str(int(l)))
            if error_different_no_of_clusters==0:
                pass
            else:
                if error_different_no_of_clusters==1:
                    for k in extra_cluster_names:
                        B1_dict[j][0][k]=[]
                else:
                    for k in extra_cluster_names:
                        B2_dict[j][0][k]=[]
            
                
            
                # ----------------------------  Combinations:  ----------------------------------------    
        
            list_1=B1_clusters
            list_2=B2_clusters
            unique_combinations=list(itertools.product(list_1, list_2))
            
            edge_weights=[]
            
            B1_nodes=[]
            B2_nodes=[]
            B_nodes=[]
            Edge_weights=[]
            All_edge_attributes=[]
            
            for k in unique_combinations:
                b1_node='B1_'+str(k[0])
                B1_nodes.append(b1_node)
                b2_node='B1_'+str(k[1])
                B2_nodes.append(b2_node)
                B_nodes.append(b1_node)
                B_nodes.append(b2_node)
                edge_weight=set(B1_dict[j][0][k[0]]).intersection(set(B2_dict[j][0][k[1]]))
                Edge_weights.append(edge_weight)
                all_edge_attributes=(b1_node, b2_node, {'weight': edge_weight})
                All_edge_attributes.append(all_edge_attributes)
            B1_nodes=list(np.unique(B1_nodes))
            B2_nodes=list(np.unique(B2_nodes))
            B_nodes=list(np.unique(B_nodes))
            
            B=nx.Graph()
            B.add_nodes_from(B1_nodes, bipartite=0)
            B.add_nodes_from(B2_nodes, bipartite=1)
            B.add_edges_from(All_edge_attributes)
            
            
                
                
                
                
                
                
                
                # z = i[0].intersection(i[1])
                # edge_weights.append()
            
            
            
            
            
            
            
            
            
            
            
            # weighted_edge_list=[]
            
            
            
            
            
            
            