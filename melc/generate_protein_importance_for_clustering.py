


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
import itertools
from itertools import permutations
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.inspection import permutation_importance
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from matplotlib import pyplot as plt


# =================================================== Different importance metrics on different classifiers: ============================================================
    
def generate_protein_importance_for_clustering(adata_pickle_path, essential_proteins_per_patient_pickle_path):
    with open(essential_proteins_per_patient_pickle_path, 'rb') as f:
        essential_proteins_allPatients, essentialProteinsPerCluster_acrossClusterLevels_forAllPatients, allEssentialProteins_acrossClusterLevels_forAllPatients, X_allPatients, y_allPatients=pickle.load(f)

    # -----------

    essential_proteins_allPatients_items=list(essential_proteins_allPatients.items())
    first_three_items = essential_proteins_allPatients_items[0:1]
    essential_proteins_allPatient_new=dict(first_three_items)
    essential_proteins_allPatients=essential_proteins_allPatient_new

    # -----------

    essentialProteinsPerCluster_acrossClusterLevels_forAllPatients_items=list(essentialProteinsPerCluster_acrossClusterLevels_forAllPatients.items())
    first_three_items = essentialProteinsPerCluster_acrossClusterLevels_forAllPatients_items[0:1]
    essentialProteinsPerCluster_acrossClusterLevels_forAllPatients_new=dict(first_three_items)
    essentialProteinsPerCluster_acrossClusterLevels_forAllPatients=essentialProteinsPerCluster_acrossClusterLevels_forAllPatients_new

    # -----------

    allEssentialProteins_acrossClusterLevels_forAllPatients_items=list(allEssentialProteins_acrossClusterLevels_forAllPatients.items())
    first_three_items = allEssentialProteins_acrossClusterLevels_forAllPatients_items[0:1]
    allEssentialProteins_acrossClusterLevels_forAllPatients_new=dict(first_three_items)
    allEssentialProteins_acrossClusterLevels_forAllPatients=allEssentialProteins_acrossClusterLevels_forAllPatients_new

    # -----------

    X_allPatients_items=list(X_allPatients.items())
    first_three_items = X_allPatients_items[0:1]
    X_allPatients_new=dict(first_three_items)
    X_allPatients=X_allPatients_new

    # -----------

    y_allPatients_items=list(y_allPatients.items())
    first_three_items = y_allPatients_items[0:1]
    y_allPatients_new=dict(first_three_items)
    y_allPatients=y_allPatients_new

    # -----------


    with open(adata_pickle_path, 'rb') as f:
        pickle_= pickle.load(f)
    res = list(pickle_.keys())[0]
    all_columns=pickle_[res].to_df().columns



    cnt=-1

    # ------
    pickle_items=list(pickle_.items())
    first_three_items = pickle_items[0:1]
    pickle_new=dict(first_three_items)
    pickle_=pickle_new
    # ------

    for i in essential_proteins_allPatients:
        protein_names=essential_proteins_allPatients[i]
        
        cnt+=1
        essential_proteins=essential_proteins_allPatients[i]
        selected_proteins_index=[]
        column_count=-1
        
        for j in all_columns:
            column_count+=1
            if j in essential_proteins:
                selected_proteins_index.append(column_count)
        
        X=X_allPatients[i][:,selected_proteins_index]
        y=y_allPatients[i]

        X=X[0:10]
        y=y[0:10]
        
        # --------------------- Permutation importance on one-vs-all classifier: -----------------------
        clf = OneVsRestClassifier(SVC()).fit(X, y)
        # to get permutation: 
        results = permutation_importance(clf, X, y, scoring='accuracy')
        # get important features:
        important_features = results.importances_mean
        # protein_names=essential_proteins_allPatients[i]
        proteinNames_featureImportances_dict=dict(zip(protein_names, list(important_features)))
        sorted_list = sorted(proteinNames_featureImportances_dict.items(), key = lambda x:x[1], reverse = True)
        sorted_list = dict(sorted_list)
        # # list all features:
        # for k,v in enumerate(important_features):
        #  print('Feature: %0d, Score: %.5f' % (k,v))
        for k in  sorted_list:
            print(f'Feature: {k}, Score: {sorted_list[k]}')
        # ---
        # save the names and their respective scores separately
        # reverse the tuples to go from most frequent to least frequent 
        feature = list(sorted_list.keys())
        score = list(sorted_list.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('Score')
        plt.savefig(f'feature_imp_scores_patient{i}_PermImportance_OneVsAllClassifier.pdf', format='pdf')
        plt.show()
        # ------------------------------------------------------------------------------------------------
        
        # --------------------- Gini index on RF classifier: -----------------------
        rf = RandomForestRegressor(n_estimators=100)
        rf.fit(X, y)
        important_features=rf.feature_importances_
        proteinNames_featureImportances_dict=dict(zip(protein_names, list(important_features)))
        sorted_list = sorted(proteinNames_featureImportances_dict.items(), key = lambda x:x[1], reverse = True)
        sorted_list = dict(sorted_list)
        # # list all features:
        # for k,v in enumerate(important_features):
        #  print('Feature: %0d, Score: %.5f' % (k,v))
        for k in  sorted_list:
            print(f'Feature: {k}, Score: {sorted_list[k]}')
        
        
        feature = list(sorted_list.keys())
        score = list(sorted_list.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('Score')
        plt.savefig(f'feature_imp_scores_patient{i}_Gini_RFClassifier.pdf', format='pdf')
        plt.show()
        
        # ------------------------------------------------------------------------------------------------
        
        # --------------------- Permutation importance index on RF classifier: -----------------------
        rf = RandomForestRegressor(n_estimators=100)
        rf.fit(X, y)
        # to get permutation: 
        results = permutation_importance(rf, X, y)
        # get important features:
        important_features = results.importances_mean
        # protein_names=essential_proteins_allPatients[i]
        proteinNames_featureImportances_dict=dict(zip(protein_names, list(important_features)))
        sorted_list = sorted(proteinNames_featureImportances_dict.items(), key = lambda x:x[1], reverse = True)
        sorted_list = dict(sorted_list)
        # # list all features:
        # for k,v in enumerate(important_features):
        #  print('Feature: %0d, Score: %.5f' % (k,v))
        for k in  sorted_list:
            print(f'Feature: {k}, Score: {sorted_list[k]}')
        # ---
        # save the names and their respective scores separately
        # reverse the tuples to go from most frequent to least frequent 
        feature = list(sorted_list.keys())
        score = list(sorted_list.values())
        x_pos = np.arange(len(feature))
        # ---
        plt.bar(x_pos, score,align='center')
        plt.xticks(x_pos, feature, rotation=90) 
        plt.ylabel('Score')
        plt.savefig(f'feature_imp_scores_patient{i}_PermImportance_RFClassifier.pdf', format='pdf')
        plt.show()
        
        # ------------------------------------------------------------------------------------------------    
    
# How to run:    
# generate_protein_importance_for_clustering('<path to anndata object as a pickle file>', '<path to essential_proteins_for_clustering dict from a previous step, as a pickle file>')    
    
    





































