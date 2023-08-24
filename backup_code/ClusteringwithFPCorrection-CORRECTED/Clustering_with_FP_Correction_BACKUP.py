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
from community import community_louvain
patient_id=20751

def categorize(dataframe, leiden_clusters_sorted):
    df = pd.DataFrame(columns=['group'])
    for i in leiden_clusters_sorted:
        df_=dataframe[dataframe['group']==i]
        df_=df_.sort_values(by = 'pvals_adj', ascending=True)
        df=pd.concat([df, df_], ignore_index=True)
    return df

with open('TNBC_41patients_KerenEtAl.pkl', 'rb') as f:
    pickle_= pickle.load(f)

cnt=-1

# ------
pickle_items=list(pickle_.items())
first_three_items = pickle_items[7:8]
pickle_new=dict(first_three_items)
pickle_=pickle_new
# ------
   
no_of_patients=np.shape(pickle_)
error_samples=[]
weight=0.1

time_taken_per_epoch_per_patient_per_trial=[]
time_taken_per_cluster_per_epoch_per_patient_per_trial=[]
louvain_clustering_per_patient_per_trial=[]
for trial in range(10):
    print('*********** ========================== TRIAL ' + str(trial) + ' ========================== ***********' )
    time_taken_per_epoch_per_patient=[]
    time_taken_per_cluster_per_epoch_per_patient=[]
    louvain_clustering_per_patient=[]
    for i in pickle_:
        cnt+=1
        adata=pickle_[i]
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
        adata_raw=adata.copy()
        sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=3)
        # ===========================================================
        G = nx.Graph()
        weight=0.1
        G_edges_withWeights={}
        time_taken_per_epoch=[]
        time_taken_per_cluster_per_epoch=[]
        for _counter_ in range(10):
            print('Epoch: ' + str(_counter_) + ' in TRIAL '+str(trial))
            print('==========================')
            start_epoch_time = time.process_time()
            scs.inference.planted_model(adata_raw)
            ppbm_clusters_sorted=[str(j) for j in sorted([int(i) for i in np.unique(adata_raw.obs['ppbm'])])]
            # ======================================================================================
            ppbm_df=adata_raw.obs['ppbm'].to_frame()
            ppbm_df['cell_label'] = ppbm_df.index
            a = [[] for x in range(len(ppbm_clusters_sorted))]
            for index, row in ppbm_df.iterrows():
                a[int(row['ppbm'])].append(row['cell_label'])
            # ======================================================================================
            time_taken_per_cluster=[]
            cluster_count=-1
            for k in a:
                cluster_count+=1
                print('Cluster '+str(cluster_count)+' in epoch: '+str(_counter_)+'in TRIAL '+str(trial))
                print('---')
                start_cluster_time = time.process_time()
                G_=nx.complete_graph(len(k))
                # mapping = {old_label:new_label for new_label, old_label in dict(enumerate(L))}
                # G_=nx.relabel_nodes(G, mapping)
                # nx.relabel_nodes(G_,dict(enumerate(L)), copy = False)
                G_=nx.relabel_nodes(G_,dict(enumerate(k)))
                values=[weight]*len(k)
                dict_=dict(zip(k, values))
                nx.set_node_attributes(G_, dict_, "cell_number")
                
                for j in list(list(G_.edges(data=True))):
                    if (j[0], j[1]) in G_edges_withWeights.keys():
                        new_edge_weight=weight+G_edges_withWeights[(j[0], j[1])]
                        # Gdash_edges_withWeights.append((j[0], j[1], new_edge_weight))
                        G_edges_withWeights[(j[0], j[1])]=new_edge_weight
                        # print(Gdash_edges_withWeights)
                    else:
                        # Gdash_edges_withWeights.append((j[0], j[1], weight))
                        G_edges_withWeights[(j[0], j[1])]=weight
                time_taken_per_cluster.append(time.process_time()-start_cluster_time)
            time_taken_per_cluster_per_epoch.append(time_taken_per_cluster)
            time_taken_per_epoch.append(time.process_time()-start_epoch_time)
        time_taken_per_epoch_per_patient.append(time_taken_per_epoch)
        time_taken_per_cluster_per_epoch_per_patient.append(time_taken_per_cluster_per_epoch)
        for k in G_edges_withWeights:
            if k[0]==k[1]:
                print('Faulty:')
                print((k[0], k[1]))
            G.add_edge(k[0], k[1])
            G[k[0]][k[1]]['weight']=G_edges_withWeights[(k[0],k[1])]
        # ======================================================================================
        comms = community_louvain.best_partition(G)
        _dict_={}
        for i in list(sorted(np.unique(list(comms.values())))):
            _list_=[k for k, v in comms.items() if v == i]
            _dict_[i]=_list_
        
        louvain_clustering_per_patient.append(_dict_)
    louvain_clustering_per_patient_per_trial.append(louvain_clustering_per_patient)
    
    
    
    
    
    


    
    

    
    
    