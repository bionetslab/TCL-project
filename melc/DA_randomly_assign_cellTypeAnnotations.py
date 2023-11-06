


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
import itertools
import random
from scipy.stats import mannwhitneyu
import os
import warnings
import gseapy as gp
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import mygene
import seaborn as sns
from gseapy import barplot, dotplot
import matplotlib.pyplot as plt
import decoupler as dc
from numpy.random import default_rng
from copy import deepcopy
import squidpy as sq
# from neighborhood_enrichment import neighborhood_enrichment
import schist as scs
from pingouin import mwu
import graph_tools as gt


def da_randomly_assign_cell_type_annotations(adata_pickle_path):
    with open(adata_pickle_path, 'rb') as f:
        pickle_= pickle.load(f)


    patient_ids=[]
    for i in pickle_:
        patient_ids.append(i)


    cnt=-1

    # ------
    pickle_items=list(pickle_.items())
    first_three_items = pickle_items[:]
    pickle_new=dict(first_three_items)
    pickle_=pickle_new
    # ------
    
    no_of_patients=np.shape(pickle_)
    no_of_celltypes_limits=[5, 20]

    # def calculate_


    for i in pickle_:
        
        print('Patient serial number' + str(cnt) + ' (patient id: ' + str(i) + ')')
        print('====================================================================')
        
        cnt+=1
        
        adata=pickle_[i]
        adata_X=adata.to_df()
        
        no_of_samples_X=np.shape(adata_X)[0]
        no_of_columns_X=np.shape(adata_X)[1]
        
        no_of_celltypes=random.randint(no_of_celltypes_limits[0], no_of_celltypes_limits[1])
        print(no_of_celltypes)
        
        celltypes=[]
        for samples_cnt in range(no_of_samples_X):
            celltypes.append(random.randint(1, no_of_celltypes))
        
        cell_types_cnt=-1
        for _celltypes_ in celltypes:
            cell_types_cnt+=1
            celltypes[cell_types_cnt]='CT-'+str(_celltypes_)
        
        celltypes=pd.Categorical(celltypes)
        adata.obs['celltype']=celltypes
        
        # celltypes=pd.Series(celltypes, dtype="category")
        # celltypes=pd.DataFrame(celltypes, columns=['abc'])
        # adata.obs['celltype']=celltypes['abc']
        

    with open('TNBC_41patients_KerenEtAl_with_celltype_annotations.pkl', 'wb') as f:
        pickle.dump(pickle_, f)


