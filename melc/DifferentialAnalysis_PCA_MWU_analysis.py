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






def da_pca_mwu_analysis(adata_pickle_path, dependent_variable_name):

    with open(adata_pickle_path, 'rb') as f:
       pickle_= pickle.load(f)
    res = list(pickle_.keys())[0]
    proteins_list=list(pickle_[res].to_df().columns)
    
    # ===================================
    
    labels=[]
    
    for i in pickle_:
        labels.append(np.unique(pickle_[i].obsm[dependent_variable_name])[0])
    
    # ====================================
    pc1_1=[]
    pc1_0=[]
    # ---
    information_content_OR_variance_ratio_1=[]
    information_content_OR_variance_ratio_0=[]
    # ---
    singular_values_1=[]
    singular_values_0=[]
    # ====================================
    overall_patient_dataframe=pd.DataFrame(columns=proteins_list)
    cnt=-1
    for i in pickle_:
        cnt+=1
        adata=pickle_[i]
        adata_X=adata.to_df()
        overall_patient_dataframe=pd.concat([overall_patient_dataframe, adata_X], axis=0)
    
    # ====================================
    
    # randomized_overall_patient_dataframe = dict(zip(proteins_list, [None]*len(proteins_list)))
    randomized_overall_patient_dataframe = {}
    
    # ====================================
    randomized_pc1_1={}
    randomized_pc1_0={}
    # ---
    randomized_information_content_OR_variance_ratio_1={}
    randomized_information_content_OR_variance_ratio_0={}
    # ---
    randomized_singular_values_1={}
    randomized_singular_values_0={}
    # ====================================
    
    for i in proteins_list:
        shuffled_column=list(overall_patient_dataframe[i].sample(frac=1).values)
        randomized_df=overall_patient_dataframe.copy()
        randomized_df[i]=shuffled_column
        randomized_overall_patient_dataframe[i]=randomized_df
        # --------------------------------------------------------------------
        randomized_pc1_1[i], randomized_pc1_0[i]=[], []
        randomized_information_content_OR_variance_ratio_1[i], randomized_information_content_OR_variance_ratio_0[i]=[], []
        randomized_singular_values_1[i], randomized_singular_values_0[i]=[], []
        
    # ====================================
    
    count_rows_begin=-1
    count_rows_end=-1
    # ---
    cnt=-1
    for i in pickle_:
        
        cnt+=1
        
        if cnt==0:
            count_rows_begin+=1
        else:
            count_rows_begin=count_rows_end+1
        
        print('Patient serial number' + str(cnt) + ' (patient id: ' + str(i) + ')')
        print('====================================================================')
        
        adata=pickle_[i]
        adata_X=adata.to_df()
        
        print(len(adata_X.index))
        
        X=np.array(adata_X)
        pca = PCA(n_components=1)
        principal_components=pca.fit_transform(X)
        principal_components_list=list(chain.from_iterable([l.tolist() for l in principal_components]))
        principal_components=principal_components_list
        # print(len(principal_components))
        # -----
        principal_components_df = pd.DataFrame(data = principal_components, columns = ['principal component 1'])
        pca_explained_variance_ratio=pca.explained_variance_ratio_[0]
        pca_singular_values=pca.singular_values_[0]
        # ---
        label=np.unique(pickle_[i].obsm[dependent_variable_name])[0]
        # ---
        if label.lower()=='positive' or label==1 or label=='1':
            pc1_1.append(list(principal_components))
            information_content_OR_variance_ratio_1.append(pca_explained_variance_ratio)
            singular_values_1.append(pca_singular_values)
        else:
            pc1_0.append(list(principal_components))
            information_content_OR_variance_ratio_0.append(pca_explained_variance_ratio)
            singular_values_0.append(pca_singular_values)

        # ====================================================================================================
        count_rows_end=count_rows_begin+len(list(adata_X.index))-1
        # print(count_rows_begin)
        # print(count_rows_end)
        for protein_count in proteins_list:
            random_X=randomized_overall_patient_dataframe[protein_count].iloc[count_rows_begin:count_rows_end+1]
            # print(len(random_X.index))
            X=np.array(random_X)
            pca = PCA(n_components=1)
            principal_components=pca.fit_transform(X)
            principal_components_list=list(chain.from_iterable([l.tolist() for l in principal_components]))
            principal_components=principal_components_list
            # print(len(principal_components))
            # -----
            principal_components_df = pd.DataFrame(data = principal_components, columns = ['principal component 1'])
            pca_explained_variance_ratio=pca.explained_variance_ratio_[0]
            pca_singular_values=pca.singular_values_[0]
            # ------------------------------------------------------------------------------------------------
            if label.lower()=='positive' or label==1 or label=='1':
                randomized_pc1_1[protein_count].append(list(principal_components))
                randomized_information_content_OR_variance_ratio_1[protein_count].append(pca_explained_variance_ratio)
                randomized_singular_values_1[protein_count].append(pca_singular_values)
            else:
                randomized_pc1_0[protein_count].append(list(principal_components))
                randomized_information_content_OR_variance_ratio_0[protein_count].append(pca_explained_variance_ratio)
                randomized_singular_values_0[protein_count].append(pca_singular_values)
    
    pc1_1=list(chain.from_iterable(pc1_1))
    pc1_0=list(chain.from_iterable(pc1_0))
    pValues_pca={}
    
    print('\n--------------------- X ---------------------\n')
    
    accuracy_allproteins={}
    cnt=-1
    for protein_count in proteins_list:
        print(f'Running MWU-test for protein channel: {protein_count}')
        cnt+=1
        randomized_pc1_1[protein_count], randomized_pc1_0[protein_count]=list(chain.from_iterable(randomized_pc1_1[protein_count])), list(chain.from_iterable(randomized_pc1_0[protein_count]))                                                                                 
        U1, p = mannwhitneyu(randomized_pc1_1[protein_count], randomized_pc1_0[protein_count], method="exact")
        pValues_pca[protein_count]=p
        # # ---
        # if cnt==0:
        #     y_1, y_0=[1]*len(randomized_pc1_1[protein_count]), [0]*len(randomized_pc1_0[protein_count])
        #     y=y_1.copy()
        #     for j in y_0:
        #         y.append(j)
        #     y=np.array(y)
        #     # y.append(y_0)
        #     # y=list(chain.from_iterable(y))
        # X_1=randomized_pc1_1[protein_count]
        # X_0=randomized_pc1_0[protein_count]
        # X=X_1.copy()
        # for j in X_0:
        #     X.append(j)
        # X=np.array(X).reshape(-1, 1)
        # # X=list(chain.from_iterable(X))
        # # ---
        # X_train,X_test,y_train,y_test=train_test_split(X,y,test_size=0.25)
        # # ---
        # # clf = OneVsRestClassifier(SVC()).fit(X, y)
        # # clf.fit(X_train, y_train)
        # # y_pred = clf.predict(X_test)
        # # ---
        randomized_information_content_OR_variance_ratio_1[protein_count], randomized_information_content_OR_variance_ratio_0[protein_count]=randomized_information_content_OR_variance_ratio_1[protein_count], randomized_information_content_OR_variance_ratio_0[protein_count]
        # ---
        randomized_singular_values_1[protein_count], randomized_singular_values_0[protein_count]=randomized_singular_values_1[protein_count], randomized_singular_values_0[protein_count]                                                                               


    # # =============================================================
    # X_allproteins=randomized_pc1_1
    # # =============================================================


    U1, p = mannwhitneyu(pc1_1, pc1_0, method="exact")
    pValues_pca_important=pValues_pca.copy()
    for i in pValues_pca:
        if pValues_pca[i]<p:
            pValues_pca_important.pop(i)
    # ====================
    important_proteins_sorted={}
    if len(list(pValues_pca_important.keys()))!=0:
        sorted_list =sorted(pValues_pca_important.items(), key=lambda x:x[1])
        for i in range(len(sorted_list)):
            important_proteins_sorted[sorted_list[i][0]]=sorted_list[i][1]
    
        # ---
        fig, axes = plt.subplots(figsize=(10, 4))
        fig.suptitle(f'Differential analysis ({dependent_variable_name}) – significant protein expressions (p-value<0.05)')
        feature = list(important_proteins_sorted.keys())
        score = list(important_proteins_sorted.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('p-value (MWU test)')
        plt.savefig(f'DA_significant_protein_expressions_pvalues.pdf', format='pdf')
        plt.show()
    
    
    
    
    
    