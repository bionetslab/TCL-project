
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

def da_protein_correlation(adata_pickle_path, dependent_variable_name):
    with open(adata_pickle_path, 'rb') as f:
        pickle_= pickle.load(f)
    res = list(pickle_.keys())[0]
    proteins_list=list(pickle_[res].to_df().columns)
    no_of_proteins=len(proteins_list)

    proteins_list_matrix=[]
    for i in proteins_list:
        proteins_list_matrix_temp=[]
        for j in proteins_list:
            proteins_list_matrix_temp.append((i,j))
        proteins_list_matrix.append(proteins_list_matrix_temp)

    proteins_list_matrix_array=np.array(proteins_list_matrix)
    protein_combination_names=proteins_list_matrix_array[np.triu_indices(no_of_proteins, k = 1)]
    res = tuple(tuple(sub) for sub in protein_combination_names)
    protein_combination_names=res


    protein_correlation_1=[]
    protein_correlation_0=[]


    # ===================================

    labels=[]

    for i in pickle_:
        labels.append(np.unique(pickle_[i].obsm[dependent_variable_name])[0])


    cnt=-1
    for i in pickle_:
        
        
        print('Patient serial number' + str(cnt) + ' (patient id: ' + str(i) + ')')
        print('====================================================================')
        
        cnt+=1
        
        adata=pickle_[i]
        adata_X=adata.to_df()
        
        adata_protein_correlation=np.array(adata_X.corr())
        protein_combination_values=list(adata_protein_correlation[np.triu_indices(no_of_proteins, k = 1)])
        
        label=np.unique(pickle_[i].obsm[dependent_variable_name])[0]
        
        
        if label.lower()=='positive' or label==1 or label=='1':
            protein_correlation_1.append(protein_combination_values)
        else:
            protein_correlation_0.append(protein_combination_values)

        
    protein_correlation_df_1=pd.DataFrame(protein_correlation_1, columns=protein_combination_names)
    protein_correlation_df_0=pd.DataFrame(protein_correlation_0, columns=protein_combination_names)

    imp_proteins_names=[]
    imp_proteins_pValues=[]

    cnt=-1
    for i in protein_combination_names:
        cnt+=1
        U1, p = mannwhitneyu(list(protein_correlation_df_1[i]), list(protein_correlation_df_0[i]), method="exact")
        if p<0.05:
            imp_proteins_names.append(i)
            imp_proteins_pValues.append(p)   

    important_proteins=dict(zip(imp_proteins_names, imp_proteins_pValues))

    important_proteins_sorted={}
    if len(list(important_proteins.keys()))!=0:
        sorted_list =sorted(important_proteins.items(), key=lambda x:x[1])
        for i in range(len(sorted_list)):
            important_proteins_sorted[sorted_list[i][0]]=sorted_list[i][1]

        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential analysis ({dependent_variable_name}) – significant protein correlations (p-value<0.05)')
        feature = list(important_proteins_sorted.keys())
        score = list(important_proteins_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'protein_correlation_pvalues_impProteins_allPatients.pdf', format='pdf')
        plt.show()
        
    










