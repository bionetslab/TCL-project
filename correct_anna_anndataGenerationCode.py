


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

import imageio as iio
from skimage import img_as_ubyte
import re

# ----------------------------------------

antibodies={'ADAM10-FITC': 'ADAM10',
    'ADAM10-PE': 'ADAM10',
    'CD107a-FITC': 'LAMP1',
    'CD11a-PE': 'ITGAL',
    'CD11c-PE': 'ITGAX',
    'CD123-FITC': 'IL3RA',
    'CD138-FITC': 'SDC1',
    'CD14-PE': 'LBP',
    'CD16-PE': 'FCGR3A',
    'CD163-PE': 'CD163',
    'CD19-PE': 'CD19',
    'CD1a-PE': 'CD1A',
    'CD20-FITC': 'MS4A1',
    'CD205-PE': 'LY75',
    'CD206-FITC': 'MRC1',
    'CD209-FITC': 'CD209',
    'CD24-FITC': 'CD24',
    'CD25-PE': 'IL2RA',
    'CD271-FITC': 'NGFR',
    'CD29-FITC': 'ITGB1',
    'CD3-PE': 'CD3E',
    'CD36-FITC': 'CD36',
    'CD38-PE': 'CD38',
    'CD4-PE': 'CD4',
    'CD40-PE': 'CD40',
    'CD44-PE': 'CD44',
    'CD45-PE': 'PTPRC',
    'CD45RA-PE': 'PTPRC',
    'CD45RO-FITC': 'PTPRC',
    'CD52-FITC': 'CD52',
    'CD54-FITC': 'ICAM1',
    'CD55-PE': 'CD55',
    'CD56-PE': 'NCAM1',
    'CD6-PE': 'CD6',
    'CD62P-PE': 'SELP',
    'CD63-FITC': 'CD63',
    'CD68-FITC': 'CD68',
    'CD69-PE': 'CD69',
    'CD71-FITC': 'TFRC',
    'CD8-PE': 'CD8A',
    'CD9-PE': 'CD9',
    'CD95-PE': 'TNFRSF6',
    'Collagen IV-FITC': 'COL4A1',
    'Cytokeratin-14-FITC': 'KRT14',
    'E-Cadherin-FITC': 'CDH1',
    'EGFR-AF488': 'EGFR',
    'HLA-ABC-PE': 'HLA-A, HLA-B, HLA-C',
    'HLA-DR-PE': 'HLA-DRB1',
    'KIP1-FITC': 'CDKN1B',
    'Melan-A-FITC': 'MLANA',
    'Nestin-AF488': 'NES',
    'Notch-1-FITC': 'NOTCH1',
    'Notch-2-PE': 'NOTCH2',
    'Notch-3-PE': 'NOTCH3',
    'Notch-4-PE': 'NOTCH4',
    'PPARgamma-FITC': 'PPARG',
    'PPB-FITC': 'VIM',
    'TACE-FITC': 'ADAM17',
    'TAP73-FITC': 'TP73',
    'TNFR1-PE': 'TNFRSF1A',
    'TNFR2-PE': 'TNFRSF1B',
    'Vimentin-FITC': 'VIM',
    'beta-Catenin-FITC': 'CTNNB1',
    'p63-FITC': 'TP63',
    'phospho-Connexin-FITC': 'GJA1',
    'Melan-A-AF488': 'MLANA',
    'Bcl-2-FITC': 'BCL2',
    'CD10-FITC': 'CD10',
    'CD11b-FITC': 'ITGAM',
    'CD13-FITC': 'ANPEP',
    'CD141-PE': 'THBD',
    'CD15-FITC': 'FUT4',
    'CD2-FITC': 'CD2',
    'CD27-PE': 'CD27',
    'CD276-FITC': 'B7-H3',
    'CD43-PE': 'SPN',
    'CD5-PE': 'CD5',
    'CD53-FITC': 'F11R',
    'CD7-PE': 'CD7',
    'CD81-FITC': 'CD81',
    'CD83-FITC': 'CD83',
    'CD90-FITC': 'THY1',
    'CD146-FITC': 'MCAM',
    'CD34-FITC': 'CD34',
    'CD39-FITC': 'ENTPD1',
    'CD56-FITC': 'NCAM1',
    'CD58-FITC': 'CD58',
    'CD66abce-FITC': 'CEACAM',
    'CD73-FITC': 'NT5E',
    'CD80-PE': 'CD80',
    'L302-FITC': 'CD302',
    'CD141-FITC': 'THBD',
    'APPC-FITC': 'CD172a',
    'Actin-FITC': 'ACTB',
    'BOP-FITC': 'TSPAN7',
    'C97-FITC': 'CD97',
    'CALCA-FITC': 'CALCA',
    'CEB2-FITC': 'CEBPB',
    'CEB6-FITC': 'CEBPE',
    'CK2-A1-FITC': 'CSNK2A1',
    'CK2-A2-FITC': 'CSNK2A2',
    'CPB3-FITC': 'CPB3',
    'Caspase 3 active-FITC': 'CASP3',
    'Cyclin D1-FITC': 'CCND1',
    'Cytochrome C-FITC': 'CYCS',
    'Cytokeratin-15-DyLight488': 'KRT15',
    'DKK3-FITC': 'DKK3',
    'DRO-FITC': 'DROSHA',
    'Desmin-FITC': 'DES',
    'EBF-P-FITC': 'EBF1',
    'EK2-FITC': 'KRT8',
    'ERG1-FITC': 'ERG',
    'FCepsilonRIa-FITC': 'FCER1A',
    'FLT-FITC': 'FLT3',
    'FRP 13E12-FITC': 'HLA-DRB1',
    'FST-FITC': 'FST',
    'KIS-FITC': 'MKI67',
    'Ki67-FITC': 'MKI67',
    'LOC-HE-FITC': 'NHEJ1',
    'LPAM-1-FITC': 'ITGAE',
    'LSDP-FITC': 'GFER',
    'MCSP-FITC': 'CSPG4',
    'MECP2-FITC': 'MECP2',
    'MPV-FITC': 'MYH9',
    'NKp80-FITC': 'KLRC1',
    'NRG9-FITC': 'NRG1',
    'Neutrophil Elastase-FITC': 'ELANE',
    'ORC2-FITC': 'ORC2',
    'ORC3-FITC': 'ORC3',
    'PBS': 'Phosphate buffered salt solution',
    'PCNA-FITC': 'PCNA',
    'PGRN-FITC': 'GRN',
    'PRPB-FITC': 'PRPF8',
    'RIK-2-FITC': 'ETS1',
    'RIM3-FITC': 'RIMS3',
    'Reelin-FITC': 'RELN',
    'STAT3-FITC': 'STAT3',
    'SYT10-FITC': 'SYT10',
    'TDP-FITC': 'TARDBP',
    'beta-Tubulin-FITC': 'TUBB',
    'p53-FITC': 'TP53',
    'CD163-FITC': 'CD163',
    'Ki67-AF488': 'MKI67',
    'CD20-PE': 'MS4A1',
    'CD3-FITC': 'CD3E',
    'CD303-PE': 'BDCA-2 (CD303)',
    'CD31-FITC': 'PECAM1 (CD31)',
    'CD46-FITC': 'CD46',
    'CD62L-FITC': 'SELL (CD62L)',
    'CD71-PE': 'TFRC',
    'CLA-FITC': 'CLA',
    'FoxP3-FITC': 'FOXP3',
    'IgA-FITC': 'IGHG',
    'TcR alpha': 'TRAC',
    'CD117-PE': 'KIT',
    'CD2-PE': 'CD2',
    'CD31-PE': 'PECAM1 (CD31)',
    'CD34-PE': 'CD34',
    'CD68-PE': 'CD68',
    'CD73-PE': 'NT5E',
    'Follicular Dendritic Cells-FITC': 'CD21',
    'FoxP3-PE': 'FOXP3',
    'III beta-Tubulin-FITC': 'TUBB3',
    'TSLP-PE': 'TSLP'}

antibodies_list=list(antibodies.keys())

anndatas_list=[]

patient_labels=[]

with open('/Users/surya/Documents/preprocessed_anndata_files_ctcl/adata_cell.pickle', 'rb') as f:
   pickle_= pickle.load(f)

# ----------------------------------------

for i in pickle_:
    
    df_X=pickle_[i].to_df()
    df_X.columns=antibodies_list
    
    # ---
    
    _anndata_=sc.AnnData(df_X)

    # ===============================

    _anndata_.obsm['cellLabelInImage']=pickle_[i].obsm['cellLabelInImage']
    _anndata_.obsm['cellSize']=pickle_[i].obsm['cellSize']
    _anndata_.obsm['Group']=pickle_[i].obsm['Group']
    _anndata_.obsm['patient_label']=pickle_[i].obsm['patient_label']
    _anndata_.obsm['field_of_view']=pickle_[i].obsm['field_of_view']
    _anndata_.obsm['control_mean_expression']=pickle_[i].obsm['control_mean_expression']
    _anndata_.obsm['control_std_expression']=pickle_[i].obsm['control_std_expression']

    # ===============================

    # _anndata_.uns['patient_id']=pickle_[i].uns['patient_id']
    _anndata_.uns['patient_id']=i
    _anndata_.uns['control_mean_expression']=pickle_[i].uns['control_mean_expression']
    _anndata_.uns['control_std_expression']=pickle_[i].uns['control_std_expression']
    _anndata_.uns['cell_coordinates']=pickle_[i].uns['cell_coordinates']

    # ---

    cell_segment_dict = {}
    cell_segment_dict['segmentation_mask'] = {}
    _anndata_.uns['spatial']=cell_segment_dict
    _anndata_.uns['spatial']['segmentation_mask']['images']={}
    _anndata_.uns['spatial']['segmentation_mask']['images']={'hires': pickle_[i].uns['spatial']['segmentation']}

    # ---

    _anndata_.uns['spatial']['segmentation_mask']['protein_channels']={}
    _anndata_.uns['spatial']['segmentation_mask']['protein_channels']=pickle_[i].uns['spatial']['images']

    # ===============================
    
    anndatas_list.append(_anndata_)
    patient_labels.append(int(i))


anndatas_dict=dict(zip(patient_labels, anndatas_list))


with open('CTCL_preprocessed_FAUuniclinic.pkl', 'wb') as f:
    pickle.dump(anndatas_dict, f)
















