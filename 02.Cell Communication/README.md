# Code for Cell Communication Analysis

This folder contains the source code for constructing the cell communication network, and the ligand/receptor pairs.

- `bbann` is Python/C++ module for constructing cell-wise affinity network. 
- `LR_pairs` is the ligand/receptor pairs used for the analysis.

The C++ module is bridged to Python using [Pybind11](https://github.com/pybind/pybind11). To compile it for calling in python, one can run:

```bash
c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) bbann.cpp -o bbann$(python3-config --extension-suffix)
```

In the python environment, a simple usage is as follows:

```python
import pandas as pd
import scanpy as sc
import ann.bbann as ann

# Read AnnData of the dataset
adata = sc.read_h5ad('adata.h5ad')

# Read ligand/receptor pairs
lrp = pd.read_csv("LR_pairs.txt", sep='\t')

# Compute cell-wise affinity network
# matrix is the cell-gene expression matrix, with cells as rows
# genes are the gene names of the matrix
# lrp.Receptor.values is a list of receptors
# lrp.Ligand.values is a list of ligands
# lrp.Weight.values is the weight of receptor-ligand pairing, by default it's all 1
# in_memory controls whether the cell-cell affinity matrix will be stored in memory
# if not, a matrix file will be required for the on-disk computation
affy = ann.compute_affy(matrix = adata.X, genes, lrp.Receptor.values, lrp.Ligand.values, lrp.Weight.values, in_memory=False, matrix_file="01.affy.uncentered.mat.memmap")

# Next, we construct a k-most affinity network, similar to kNN, from the cell-wise network
# batch_list is a list of batch identifiers, for dataset containing multiple chemistry or batches that needs blocking
# neighbors_within_batch is the number of nearest cells in each batch, adopted from BBKNN
# set_op_mix_ratio, local_connectivity are parameters for constructing the fuzzy simplicial set
connectivities, distances, affinities = ann.bbann_matrix(affy, batch_list=adata.obs['Platform'].values, neighbors_within_batch=6, 1, 1)

# We store the cell-wise network into default slots of AnnData, and perform the downstream analyses 
adata.uns['neighbors'] = {}
adata.uns['neighbors']['params'] = {'n_neighbors': 12, 'method': 'bbann'}
adata.uns['neighbors']['distances'] = distances
adata.uns['neighbors']['connectivities'] = connectivities[0]
adata.uns['neighbors']['affinities'] = affinities

# Finally, PAGA is run to summarize the cell-wise comminication intensity to cluster level
sc.tl.paga(adata, groups='Subtype', model="v1.0")
sc.pl.paga(adata, layout= 'fa', frameon=False)

# Exporting PAGA network
cluster_names = adata.obs['Subtype'].cat.categories.tolist()
paga_con = adata.uns['paga']['connectivities'].tocoo()
paga_network = pd.DataFrame({'class_a': [cluster_names[i] for i in paga_con.row], 'class_b': [cluster_names[i] for i in paga_con.col], 'connectivity': paga_con.data})
paga_network.to_excel("paga_network.xlsx")
```

The final PAGA network can be ingested to R/Python/Cytoscape for downstream network analysis.
