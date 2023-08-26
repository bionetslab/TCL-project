
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

def da_multiple_protein_profiles(adata_pickle_path, dependent_variable_name, N, threshold):

    with open(adata_pickle_path, 'rb') as f:
        pickle_= pickle.load(f)
    res = list(pickle_.keys())[0]
    proteins_list=list(pickle_[res].to_df().columns)
    no_of_proteins=len(proteins_list)

    n_element_combinations_proteinNames = list(itertools.combinations(proteins_list, N))

    # ===================================

    labels=[]

    for i in pickle_:
        labels.append(np.unique(pickle_[i].obsm[dependent_variable_name])[0])


    # ====================================

    n_element_combination_counts_1=[]
    n_element_combination_counts_0=[]


    cnt=-1
    for i in pickle_:
        
        cnt+=1
        
        print('Patient serial number' + str(cnt) + ' (patient id: ' + str(i) + ')')
        print('====================================================================')
        
        adata=pickle_[i]
        adata_X=adata.to_df()
        
        adata_X_filtered=adata_X.copy()
        
        protein_comb_counts_temp=[]
        for j in n_element_combinations_proteinNames: 
            for k in j:
                adata_X_filtered=adata_X_filtered[adata_X_filtered[k]>threshold]
            protein_comb_counts_temp.append(len(adata_X_filtered))
        
        label=np.unique(pickle_[i].obsm[dependent_variable_name])[0]
        
        if label.lower()=='positive' or label==1 or label=='1':
            n_element_combination_counts_1.append(protein_comb_counts_temp)
        else:
            n_element_combination_counts_0.append(protein_comb_counts_temp)
        

    protein_profile_wise_cell_count_df_1=pd.DataFrame(n_element_combination_counts_1, columns=n_element_combinations_proteinNames)
    protein_profile_wise_cell_count_df_0=pd.DataFrame(n_element_combination_counts_0, columns=n_element_combinations_proteinNames)


    imp_proteins_names=[]
    imp_proteins_pValues=[]
    

    cnt=-1
    for i in n_element_combinations_proteinNames:
        cnt+=1
        U1, p = mannwhitneyu(list(protein_profile_wise_cell_count_df_1[i]), list(protein_profile_wise_cell_count_df_0[i]), method="exact")
        if p<0.05:
            imp_proteins_names.append(i)
            imp_proteins_pValues.append(p)
        

    important_proteins=dict(zip(imp_proteins_names, imp_proteins_pValues))



    # =====================================================================================================================



    important_proteins_sorted={}
    if len(list(important_proteins.keys()))!=0:
        sorted_list =sorted(important_proteins.items(), key=lambda x:x[1])
        for i in range(len(sorted_list)):
            important_proteins_sorted[sorted_list[i][0]]=sorted_list[i][1]

        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential analysis ({dependent_variable_name}) – significant protein coexpression profiles (p-value<0.05)')
        feature = list(important_proteins_sorted.keys())
        score = list(important_proteins_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'multi_protein_coexpression_pvalues_impProteins_allPatients.pdf', format='pdf')
        plt.show()



































