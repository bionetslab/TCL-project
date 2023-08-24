


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


with open('TNBC_41patients_KerenEtAl.pkl', 'rb') as f:
   pickle_= pickle.load(f)


patient_ids=[]
for i in pickle_:
    patient_ids.append(i)


cnt=-1

# ------
pickle_items=list(pickle_.items())
first_three_items = pickle_items[3:4]
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
    
    
    
    # =============================================================================================================
    
    # ----------------------------------------------- Combinations: -----------------------------------------------
    
    
    
    
    
    
    
    
    
    
    
    
    
