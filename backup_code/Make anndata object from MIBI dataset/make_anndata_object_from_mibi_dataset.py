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

import os
from skimage import io
from skimage.util import img_as_ubyte
import cv2
import imageio as iio
from skimage import img_as_ubyte
import re

images=[]
image_names=[]
# You need to change these to valid directories on your computer
input_dir = os.path.dirname('MIBI datasets/4. TNBC_shareCellData-1/')
for f in os.listdir(input_dir):
    if f.endswith('.tiff') is True:
        image_names.append(f)
        images.append(io.imread(os.path.join(input_dir, f), as_gray=True))
image_names_new=[]
for i in image_names:
    image_names_new.append(i.replace('.tiff', ''))

count=0
for i in image_names_new:
    image_names_new[count]=re.sub("[^0-9]", "", i)
    count+=1

image_names_new_=[eval(i) for i in image_names_new]
image_names_new_dict=dict(zip(image_names_new_, images))

# -------------------------------------------------------------------------------

cellData1=pd.read_csv('MIBI datasets/4. TNBC_shareCellData-1/cellData.csv', nrows=150091)
cellData_column_names=cellData1.columns.tolist()
cellData2=pd.read_csv('MIBI datasets/4. TNBC_shareCellData-1/cellData.csv', skiprows=150092, nrows=1000000, names=cellData_column_names, header=None)
cellData=pd.concat([cellData1, cellData2], axis=0)
cellData = cellData.reset_index(drop=True)
cellData.dropna()

patient_class=pd.read_csv('MIBI datasets/4. TNBC_shareCellData-1/patient_class.csv')

mmc2_supp=pd.read_excel('MIBI datasets/4. TNBC_shareCellData-1/mmc2.xlsx', usecols = ['InternalId', 'YEAR', 'AGE_AT_DX', 'STAGE', 'LATERAL', 'GRADE', 'CS_TUM_SIZE', 'RECURRENCE_LABEL', 'Survival_days_capped*'])

mmc2_supp_imp_columns=['YEAR', 'AGE_AT_DX', 'STAGE', 'GRADE', 'RECURRENCE_LABEL', 'Survival_days_capped*']

for i in mmc2_supp_imp_columns:
    mmc2_supp = mmc2_supp[mmc2_supp[i].notna()]

patient_labels_mmc2_supp=np.unique(mmc2_supp.InternalId).tolist()
cellData=cellData[cellData['SampleID'].isin(patient_labels_mmc2_supp)]

# -----------------------------------------------------------------------------------------------------

patient_labels_cellData=np.unique(cellData.SampleID).tolist()
mmc2_supp=mmc2_supp[mmc2_supp['InternalId'].isin(patient_labels_cellData)]

patient_labels=np.unique(cellData.SampleID).tolist()

count=0
anndatas_list=[]
for i in patient_labels:
    index_of_i=patient_labels.index(i)
    grouped = cellData.groupby(['SampleID'])
    df_=grouped.get_group(i)
    no_of_rows=df_.shape[0]
    df_X = df_.iloc[:,4:52]
    df_X['cellLabelInImage']=df_['cellLabelInImage'].tolist()
    df_X.set_index('cellLabelInImage', inplace=True)
    # # # ----------------------------------------------------------------------------------------------
    anndata_name='anndata'+str(i)
    exec(anndata_name + "=sc.AnnData(df_X)")
    exec(anndata_name + ".obsm['cellLabelInImage']=np.array(df_['cellLabelInImage'].tolist())")
    exec(anndata_name + ".obsm['cellSize']=np.array(df_['cellSize'].tolist())")
    exec(anndata_name + ".obsm['tumorYN']=np.array(df_['tumorYN'].tolist())")
    exec(anndata_name + ".obsm['tumorCluster']=np.array(df_['tumorCluster'].tolist())")
    exec(anndata_name + ".obsm['Group']=np.array(df_['Group'].tolist())")
    exec(anndata_name + ".obsm['immuneCluster']=np.array(df_['immuneCluster'].tolist())")
    exec(anndata_name + ".obsm['immuneGroup']=np.array(df_['immuneGroup'].tolist())")
    patient_label, year_list, age_list, stage_list, grade_list, recurrence_label_list, survival_days_capped_list=np.array([i]*no_of_rows), np.array([mmc2_supp['YEAR'].iloc[index_of_i]]*no_of_rows), np.array([mmc2_supp['AGE_AT_DX'].iloc[index_of_i]]*no_of_rows), np.array([mmc2_supp['STAGE'].iloc[index_of_i]]*no_of_rows), np.array([mmc2_supp['GRADE'].iloc[index_of_i]]*no_of_rows), np.array([mmc2_supp['RECURRENCE_LABEL'].iloc[index_of_i]]*no_of_rows), np.array([mmc2_supp['Survival_days_capped*'].iloc[index_of_i]]*no_of_rows)
    exec(anndata_name + ".obsm['patient_label']=patient_label")
    exec(anndata_name + ".obsm['year']=year_list")
    exec(anndata_name + ".obsm['stage']=stage_list")
    exec(anndata_name + ".obsm['grade']=grade_list")
    exec(anndata_name + ".obsm['recurrence']=recurrence_label_list")
    exec(anndata_name + ".obsm['survival_days_capped']=survival_days_capped_list")
    # exec(anndata_name + ".uns['spatial']['images']['cell_segments']=image_names_new_dict[i]")
    cell_label_indices=[]
    cell_coordinates=[]
    coords_x=[]
    coords_y=[]
    for j in np.unique(df_.cellLabelInImage).tolist():
        cell_label_indices.append(j)
        x_=list(np.where(df_==j))[0]
        x_avg=sum(x_)/len(x_)
        y_=list(np.where(df_==j))[1]
        y_avg=sum(y_)/len(y_)
        cell_coordinates.append(np.column_stack((x_, y_)))
        coords_x.append(x_avg)
        coords_y.append(y_avg)
    centroids=np.column_stack((coords_x, coords_y))
    exec(anndata_name + ".obsm['spatial']=centroids")
    # exec(anndata_name + ".obsm['cell_coordinates']=np.array(cell_coordinates)")
    cell_coordinates_dict=dict(zip(cell_label_indices, cell_coordinates))
    exec(anndata_name + ".uns['cell_coordinates']=cell_coordinates_dict")
    # ------
    exec(anndata_name + ".uns['patient_id']=i")
    # ------
    cell_segment_dict = {}
    cell_segment_dict['segmentation_mask'] = {}
    exec(anndata_name + ".uns['spatial']=cell_segment_dict")
    exec(anndata_name + ".uns['spatial']['segmentation_mask']['images']={}")
    exec(anndata_name + ".uns['spatial']['segmentation_mask']['images']={'hires': images[count]}")
    exec(anndata_name + ".uns['spatial']['segmentation_mask']['scalefactors']={'tissue_hires_scalef': 1, 'spot_diameter_fullres': 0.5}")
    # -----
    exec('anndatas_list.append('+anndata_name+')')
    # -----
    count+=1

anndatas_dict=dict(zip(patient_labels, anndatas_list))

import pickle
with open('TNBC_41patients_KerenEtAl.pkl', 'wb') as f:
    pickle.dump(anndatas_dict, f)
    
    
        
        
        



# from PIL import image
# from skimage import io
# io.use_plugin('pil')

# images=[]
# image_names=[]
# # You need to change these to valid directories on your computer
# input_dir = os.path.dirname('/Users/surya/Documents/GITHUB-REPOSITORIES/melc_databases/MELC Data - T-cell lymphoma/processed/t-cell lymphoma/20751/')
# for f in os.listdir(input_dir):
#     if f.endswith('.tif') is True:
#         image_names.append(f)
#         images.append(io.imread(os.path.join(input_dir, f), as_gray=True))
# image_names_new=[]
# for i in image_names:
#     image_names_new.append(i.replace('.tif', ''))
    




