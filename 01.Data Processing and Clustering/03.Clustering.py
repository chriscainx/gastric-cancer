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
# Cell clusters are then merged based on the expression of canonical markers into major cell groups.
adata.write("03.major_cell_type.h5ad")

pears_res = np.load("residual.corrected.npy")

# For each of the major cell types, cells are isolated and clusters to acquire sub celltypes.
# 600 highly variable genes are selected and then used for downstream clustering.
adata_cd8 = adata[adata.obs['Annotation'] == "CD8 T", :]
pears_res_cd8 = pears_res[:, adata.obs['Annotation'] == "CD8 T"]
adata_cd8.var['variance'] = np.var(pears_res_cd8, 1)
adata_cd8 = adata_cd8[:, adata_cd8.var.nlargest(600, 'variance').index]

# Perform dimensionality reduction, then construct a chemistry-balanced kNN using top 25 pcs.
# Clusters are then calculated from the reduced dimensions, and annotated based on their DEGs.
sc.tl.pca(adata_cd8, n_comps=100, svd_solver='arpack')
bbknn.bbknn(adata_cd8, batch_key='Platform', n_pcs=25)

# Since the constructed kNN is chemistry-balanced, we perform Leiden clustering on cells to identify stable clusters.
sc.tl.leiden(adata_cd8)
sc.tl.rank_genes_groups(adata_cd8, 'leiden', method='wilcoxon')

# The following parameters are used for generating the UMAP visualization.
sc.tl.umap(adata_cd8, min_dist=0.4, spread=1.7)

# Second level of clustering (cell subtypes) complete.
adata_cd8.write_h5ad("processed_data_CD8.h5ad")

# The same process and parameters are applied to other major cell types, including 
# CD4 T cells, Fibroblasts, Myeloid cells, Endothelial cells, Epithelial cells, Plasma cells, etc.
# Number of variable genes: 600.
# Number of Principal Components for construction of BBKNN: 25.
# Leiden clustering: default parameters on the chemistry-balanced kNN.
# UMAP visualization: min_dist = 0.4, spread = 1.7.
