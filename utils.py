# Application-specific routines for working with (smallish matrices of) data.

import base64, io, os, time, json
import numpy as np, scipy as sp, pandas as pd


def calc_nbrs_exact(raw_data, k=10):
    """
    Calculate list of `k` exact Euclidean nearest neighbors for each point.
    
    Parameters
    ----------
    raw_data: array of shape (n_samples, n_features)
        Input dataset.
    Returns
    -------
    nbr_list_sorted: array of shape (n_samples, n_neighbors)
        Indices of the `n_neighbors` nearest neighbors in the dataset, for each data point.
    """
    a = sklearn.metrics.pairwise_distances(raw_data)
    nbr_list_sorted = np.argsort(a, axis=1)[:, 1:]
    return nbr_list_sorted[:, :k]


def make_nn_graph(raw_data, k=10, symmetrize_type='fuzzy'):
    nbrs_list = calc_nbrs_exact(raw_data, k=k)
    samt = sp.sparse.lil_matrix((nbrs_list.shape[0], nbrs_list.shape[0]))
    for i in range(nbrs_list.shape[0]):
        samt[i, nbrs_list[i,:]] = 1
    samt = samt.tocsr()
    if (symmetrize_type == 'mutual'):
        sparse_adj = samt.minimum(samt.transpose())
    elif (symmetrize_type == 'inclusive'):
        sparse_adj = samt.maximum(samt.transpose())
    elif (symmetrize_type == 'fuzzy'):
        sparse_adj = samt + samt.transpose() - samt.multiply(samt.transpose())
    return sparse_adj


"""
Functions for manipulating neighborhoods on the graph.
"""

def get_neighbors(query_expt_index, embeddings, experiments):
    neighbors = experiments[numpy.where(embeddings[query_expt_index,:] != 0)]
    print("Query: {} in {}".format(experiments[query_expt_index][1], cell_type_map[experiments[query_expt_index][0]]))
    print("Neighbors: (Cell-type ==> Assay)")
    for n in neighbors:
        #print('{}\t\t\t{}'.format(cell_type_map[n[0]], n[1]))
        print('{:>12}\t\t{:>12}'.format(cell_type_map[n[0]], n[1]))
    #Summarize
    print("===========SUMMARY=============")
    print("Neighbor Assay\t\t\tNeighbor Count")
    uniq_assays = numpy.unique(neighbors[:,1], return_counts=True)
    for i in range(len(uniq_assays[0])):
        print('{}\t\t{}'.format(uniq_assays[0][i], uniq_assays[1][i]))
    return neighbors


def get_neighbor_mask(experiment, neighbor):
    return experiment[0] == neighbor[0] and experiment[1] == neighbor[1]


def plot_query_neighbors(query_expt_index, neighbors, experiments):
    neighbor_mask = numpy.squeeze([numpy.where(numpy.apply_along_axis(get_neighbor_mask, 1, experiments, neighbors[i])) for i in range(neighbors.shape[0])])
    plt.figure(figsize=(20, 20), facecolor='w')
    plt.scatter(*X_umap.T, s=0.1, color='0.7')
    plt.scatter(*X_umap[neighbor_mask].T, s=40, label="Neighbors", color='r', alpha=0.5) #Neighbors
    plt.scatter(*X_umap[query_expt_index].T, s=40, label="Query", color='b', alpha=0.5) #Query
    plt.title("{} in {}".format(experiments[query_expt_index][1], cell_type_map[experiments[query_expt_index][0]]), fontsize=14)
    plt.axis('off')
    plt.legend(fontsize=14, loc=(1.01, 0.3), markerscale=5)
    plt.show()