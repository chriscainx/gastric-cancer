import numpy as np
import pandas as pd
import dask.array as da
from scipy.sparse import coo_matrix
from umap.umap_ import fuzzy_simplicial_set

import sys
sys.path.insert(0, '/data2/gastric/running_tests/final_interaction/ann')
import bbann_internals as bi

def compute_affy(exprs, genes, receptors, ligands, weights, in_memory, matrix_file):
    '''
    generate batch-balanced affinity network
    '''
    # Tidy the exprs matrix
    n_obs = exprs.shape[0]
    # exprs should be cell * gene
    lr_pairs = pd.DataFrame({"r": receptors, "l": ligands, "w": weights})
    rec_keep = [Rec in genes for Rec in receptors]
    lig_keep = [Lig in genes for Lig in ligands]
    lr_pairs = lr_pairs[np.logical_and(rec_keep, lig_keep)]
    # generate receiver * emitter matrix
    rec_ind = [genes.index(Rec) for Rec in lr_pairs.r]
    lig_ind = [genes.index(Lig) for Lig in lr_pairs.l]
    if in_memory:
        affy = np.matmul(exprs[:, rec_ind] * np.array(lr_pairs.w), exprs[:, lig_ind].T)
        affy = affy + affy.T
    else:
        a = da.from_array(exprs[:, rec_ind] * np.array(lr_pairs.w), chunks=1000)
        b = da.from_array(exprs[:, lig_ind].T, chunks=1000)
        da_affy = da.matmul(a, b)
        da_affy = da_affy + da_affy.T
        affy = np.memmap(matrix_file, dtype=np.single, mode='w+', shape=(n_obs, n_obs), order='C')
        da.store(da_affy, affy)
    np.fill_diagonal(affy, 0)
    return affy

def bbann_matrix(affy, batch_list, neighbors_within_batch, set_op_mix_ratio, local_connectivity, flavor="umap"):
    # Encode batch var
    batch_names, batch_ind = np.unique(batch_list, return_inverse=True)
    n_batches = len(batch_names)
    n_obs = affy.shape[0]
    # and generate ann matrix. ann_dist = 1 / ann_affinity
    ann_indices, ann_affy = bi.get_graph(affy, batch_ind, n_batches, neighbors_within_batch)
    ann_dists = 1 / ann_affy
    n_neighbors = n_batches * neighbors_within_batch
    distances = coo_matrix((ann_dists, (np.repeat(np.arange(n_obs), n_neighbors), ann_indices)), shape=(n_obs, n_obs)).tocsr()
    affinities = coo_matrix((ann_affy, (np.repeat(np.arange(n_obs), n_neighbors), ann_indices)), shape=(n_obs, n_obs)).tocsr()
    ann_indices.shape = (n_obs, n_neighbors)
    ann_dists.shape = (n_obs, n_neighbors)
    if flavor == "umap":
    # compute connectivities 
        X = coo_matrix(([], ([], [])), shape=(n_obs, 1))
        connectivities = fuzzy_simplicial_set(X, n_neighbors, None, None, knn_indices=ann_indices, knn_dists=ann_dists, 
                                              set_op_mix_ratio=set_op_mix_ratio, local_connectivity=local_connectivity)
    else:
        connectivities = None
    return connectivities, distances, affinities

def bbann(adata, batch_key, genes=None, receptors=None, ligands=None, weights=None, in_memory=True, matrix_file="mat.memmap", 
          compute=True, neighbors_within_batch=6, set_op_mix_ratio=1, local_connectivity=1):
    batch_list = adata.obs[batch_key].values
    if compute:
        affy = compute_affy(adata.raw.X, adata.raw.var.index.tolist(), receptors, ligands, weights, in_memory, matrix_file)
    else:
        affy = np.memmap(matrix_file, dtype=np.single, mode='r', shape=(adata.shape[0], adata.shape[0]), order='C')
    connectivities, distances, affinities = bbann_matrix(affy, batch_list, neighbors_within_batch, set_op_mix_ratio, local_connectivity)
    adata.uns['neighbors'] = {}
    adata.uns['neighbors']['params'] = {'n_neighbors': distances.shape[1], 'method': 'bbann'}
    adata.uns['neighbors']['distances'] = distances
    adata.uns['neighbors']['connectivities'] = connectivities
    adata.uns['neighbors']['affinities'] = affinities
    return None

def contrib(exprs, genes, receptors, ligands, weights, mode, label, a, b, affinities):
    '''
    if mode == "ann" only use ann connected cells
    if mode == "all" use all cells
    '''
    # Tidy the exprs matrix
    n_obs = exprs.shape[0]
    # exprs should be cell * gene
    lr_pairs = pd.DataFrame({"r": receptors, "l": ligands, "w": weights})
    rec_keep = [Rec in genes for Rec in receptors]
    lig_keep = [Lig in genes for Lig in ligands]
    lr_pairs = lr_pairs[np.logical_and(rec_keep, lig_keep)]
    # generate receiver * emitter matrix
    rec_ind = [genes.index(Rec) for Rec in lr_pairs.r]
    lig_ind = [genes.index(Lig) for Lig in lr_pairs.l]
    gene1_ind = rec_ind + lig_ind
    gene2_ind = lig_ind + rec_ind
    mask_a = [x == a for x in label]
    mask_b = [x == b for x in label]
    if mode == "all":
        contributions, affinity = bi.calc_contrib_all(exprs[mask_a, :][:, gene1_ind], exprs[mask_b, :][:, gene2_ind])
    if mode == "ann":
        ab_affinities = affinities[np.array(mask_a), :][:, np.array(mask_b)].tocoo()
        contributions, affinity = bi.calc_contrib_ann(ab_affinities.row, ab_affinities.col, exprs[mask_a, :][:, gene1_ind], exprs[mask_b, :][:, gene2_ind])
    df_contrib = pd.DataFrame({"PartA":a, "PartB":b, "GeneA":[genes[i] for i in gene1_ind], "GeneB":[genes[i] for i in gene2_ind], "Contrib":contributions, "Affinity":affinity})
    return df_contrib