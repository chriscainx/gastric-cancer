import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
import bbknn
from matplotlib import rcParams

adata = sc.read_h5ad("02.ICA_corrected.h5ad")

# Dimensionality reduction
# Incorporate bbknn to contruct a chemistry-balanced KNN

sc.tl.pca(adata, svd_solver='arpack', n_comps=100)
bbknn.bbknn(adata, batch_key='Platform', save_knn=True)
sc.tl.umap(adata, min_dist=0.4, spread=1.7)
sc.tl.leiden(adata)

# First level of clustering (major cell types) complete.
adata.write("03.major_cell_type.h5ad")