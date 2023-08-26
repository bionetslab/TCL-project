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

def neighborhood_enrichment(anndata_object, cluster_key_option, copy_option, backend_option, show_progress_bar_option):
    leiden_nhood_enrichment=sq.gr.nhood_enrichment(anndata_object, cluster_key=cluster_key_option, copy=copy_option, backend=backend_option, show_progress_bar=show_progress_bar_option)
    leiden_nhood_enrichment_zscore=pd.DataFrame(leiden_nhood_enrichment[0])
    leiden_nhood_enrichment_count=pd.DataFrame(leiden_nhood_enrichment[1])
    result_leiden_uns = {'zscore': leiden_nhood_enrichment[0], 'count': leiden_nhood_enrichment[1]}
    return leiden_nhood_enrichment, leiden_nhood_enrichment_zscore, leiden_nhood_enrichment_count, result_leiden_uns