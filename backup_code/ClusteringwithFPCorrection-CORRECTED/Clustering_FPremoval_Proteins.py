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

def calculate_average_expressions_per_cell(adata_X, _count_):
    n=len(_count_)
    DF=adata_X.loc[_count_]
    per_column_sum_all_cells=DF.sum(axis=0).tolist()
    per_column_indiv_cells=[]
    per_column_sum_all_cells_excluding_indiv=[]
    per_column_avg_all_cells_excluding_indiv=[]
    for per_cell in _count_:
        indiv_cell=DF.loc[per_cell].tolist()
        per_column_indiv_cells.append(indiv_cell)
        sub_=list(np.subtract(per_column_sum_all_cells, indiv_cell))
        per_column_sum_all_cells_excluding_indiv.append(sub_)
        avg_=[]
        for val in sub_:
            avg_.append(val/float(n-1))
        per_column_avg_all_cells_excluding_indiv.append(avg_)
    return per_column_avg_all_cells_excluding_indiv, per_column_indiv_cells

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
no_of_trials=3
no_of_epochs=3
for trial in range(no_of_trials):
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
        
        
        for _counter_ in range(no_of_epochs):
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


# =====================================================================================

cnt=-1

for i in pickle_:
    
    cnt+=1
    
    adata=pickle_[i]
    adata_X=adata.to_df()

    hybridClustering_clusters_sorted=list(louvain_clustering_per_patient_per_trial[0][cnt].keys())
    hybridClustering_cell_lookup_dict=louvain_clustering_per_patient_per_trial[0][cnt]
    hybridClustering_cell_lookup_list=list(hybridClustering_cell_lookup_dict.values())
    
    cell_counts_actual_order=[]
    for cell_cnt in hybridClustering_cell_lookup_list:
        for indiv_cell_cnt in cell_cnt:
            cell_counts_actual_order.append(indiv_cell_cnt)
    
    p_val_actual=[]
    for _count_ in hybridClustering_cell_lookup_list:
        per_column_avg_all_cells_excluding_indiv, per_column_indiv_cells=calculate_average_expressions_per_cell(adata_X, _count_)
        pvals=[]
        cluster_count=-1
        for i_count in range(len(_count_)):
            cluster_count+=1
            U1, p = mannwhitneyu(per_column_avg_all_cells_excluding_indiv[i_count], per_column_indiv_cells[i_count], method="exact")
            p_val_actual.append({cell_counts_actual_order[cluster_count]:p})
    
    
    total_epochs=3
    pvals=[[] for _ in range(total_epochs)]
    cell_ct_hybridClustering=adata_X.index.tolist()
    cell_expression_randomizations=[]
    for epoch in range(total_epochs):
        print('Epoch: '+str(epoch))
        print('=======================================================')
        for cell_counter in cell_ct_hybridClustering:
            print('Cell id '+str(cell_counter)+' in epoch '+str(epoch))
            _df_=adata_X.loc[cell_counter].tolist()
            random.shuffle(_df_)
            # cell_expression_randomizations.append(cell_counter: asdf)
            _df_simulated=adata_X
            _df_simulated.loc[cell_counter]=_df_
            
            
            cell_counts_actual_order=[]
            for cell_cnt in hybridClustering_cell_lookup_list:
                for indiv_cell_cnt in cell_cnt:
                    cell_counts_actual_order.append(indiv_cell_cnt)
            
            # ---------------------------------------------------
            pvals[epoch]=[]
            for _count_ in hybridClustering_cell_lookup_list:
                per_column_avg_all_cells_excluding_indiv, per_column_indiv_cells=calculate_average_expressions_per_cell(adata_X, _count_)
                cluster_count=-1
                for i_count in range(len(_count_)):
                    cluster_count+=1
                    U1, p = mannwhitneyu(per_column_avg_all_cells_excluding_indiv[i_count], per_column_indiv_cells[i_count], method="exact")
                    pvals.append({cell_counts_actual_order[cluster_count]:p})
    
    
    
    
    
    


    
    

    
    
    