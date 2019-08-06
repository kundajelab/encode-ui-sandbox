# Application-specific routines for working with (smallish matrices of) data.
# Author: Akshay Balsubramani

import base64, io, os, time, json, numpy as np, scipy as sp, pandas as pd, diffmap as dm
import app_config, building_block_divs
# import goterm_caller
# from goatools.associations import read_ncbi_gene2go
# from goatools.test_data.genes_NCBI_9606_ProteinCoding import GENEID2NT
import re
import umap

# =======================================================
# ================== Utility functions ==================
# =======================================================


# Input: a sparse matrix
def deviance_residuals(rawdata, sqrt=True, binomial=False):
    rawdata = sp.sparse.csr_matrix(rawdata)
    signal_per_cell = np.ravel(rawdata.sum(axis=1))
    null_gene_abund = np.ravel(rawdata.sum(axis=0))/rawdata.sum()
    logmat = rawdata.T.multiply(1.0/signal_per_cell).T.multiply(1.0/null_gene_abund)
    logmat.data = np.log(logmat.data)
    dev_contribs = 2*logmat.multiply(rawdata)
    if binomial:
        rho = rawdata.T.multiply(1.0/signal_per_cell).T.toarray()
        comp_xent = np.multiply(1 - rho, np.log(1 - rho) - np.log(1 - null_gene_abund))
        comp_dev_contribs = np.multiply(comp_xent.T, 2*signal_per_cell).T
        contribs = dev_contribs + comp_dev_contribs
        if sqrt:
            sgn_logmat = np.sign(contribs)
            return np.multiply(np.sqrt(contribs), sgn_logmat)
        else:
            return contribs
    else:
        if sqrt:
            sgn_logmat = logmat._with_data(np.sign(logmat.data), copy=True)
            return dev_contribs.multiply(sgn_logmat).sqrt().multiply(sgn_logmat)
        else:
            return dev_contribs


# Normalize data matrix with cells as rows.
def normalize_data(datamat, scale_per_point=None, mode='log'):
    if scale_per_point is None:
        if mode == 'deviance':
            scale_per_point = 1
        else:
            scale_per_point = 10000    # np.percentile(signal_per_cell, q=90)
    signal_per_cell = np.ravel(datamat.sum(axis=1))
    scale_factors = scale_per_point/signal_per_cell
    if mode == 'log':
        return datamat.T.multiply(scale_factors).T.log1p()
    elif mode == 'deviance':
        return deviance_residuals(datamat)


# Fit_data is a sparse matrix.
def interesting_feat_ndces(fit_data, num_feats_todisplay=500, mode='var'):
    num_feats_todisplay = min(fit_data.shape[1], num_feats_todisplay)
    if ((fit_data is None) or (np.prod(fit_data.shape) == 0)):
        return np.arange(num_feats_todisplay)
    if mode == 'var':
        interestingnesses = dm.sparse_variance(fit_data, axis=0)  # np.std(fit_data, axis=0)
    elif mode == 'deviance':     # Each entry is the sqrt of its contribution to deviance.
        interestingnesses = np.ravel(fit_data.power(2).sum(axis=0))
    feat_ndces = np.argsort(interestingnesses)[::-1][:num_feats_todisplay]
    return feat_ndces


def quantile_norm(dtd):
    qtiles = np.zeros(len(dtd))
    nnz_ndces = np.nonzero(dtd)[0]
    qtiles[nnz_ndces] = sp.stats.rankdata(dtd[nnz_ndces])/len(nnz_ndces)
    return qtiles


def reviz_embed_data(fit_data, alg='UMAP'):
#     feat_ndces = interesting_feat_ndces(fit_data)
#     fit_data = fit_data[:, feat_ndces]
    if alg == 'UMAP':
        reducer = umap.UMAP(n_neighbors=15, n_components=2, random_state=42)
        embedding = reducer.fit_transform(fit_data)
    return embedding



# ========================================================
# =================== Raw data heatmap ===================
# ========================================================


def hm_row_scatter(fit_data, scatter_fig, hm_point_names, view_cocluster, row_clustIDs=None):
    row_scat_traces = []
    all_hm_point_names = []
    hmscat_mode = 'markers'
    # Decide if few enough points are around to display row labels
    if len(hm_point_names) <= 35:
        hmscat_mode = 'markers+text'
    if (scatter_fig is not None) and ('data' in scatter_fig):
        pts_so_far = 0
        # Re-sort rows or not? should be >=1 trace, so can check if first is continuous
        resort_rows = (not view_cocluster)
        is_continuous_color = None
        hm_row_ndces = []
        for trace in scatter_fig['data']:
            trace_markers = trace['marker']
            if 'line' in trace_markers:
                trace_markers['line']['width'] = 0.0
            if is_continuous_color is None:
                is_continuous_color = isinstance(trace_markers['color'], (list, tuple, np.ndarray))
            # Of the point names, choose the ones in this trace and get their indices...
            hm_point_names_this_trace = np.intersect1d(hm_point_names, trace['text'])           # point IDs in this trace
            num_in_trace = len(hm_point_names_this_trace)
            hm_point_ndces_this_trace = np.where(np.isin(hm_point_names, hm_point_names_this_trace))[0]        # this trace's row indices in heatmap
            y_coords_this_trace = np.arange(len(hm_point_names))[hm_point_ndces_this_trace]
            
            # At this point, rows are sorted in order of co-clustering. 
            # Now sort points by color within each trace. 
            # This does nothing if the colors are discrete (many traces), and is just for continuous plotting.
            if resort_rows:   # Row order determined by sorting, and continuous variable being plotted
                y_coords_this_trace = np.arange(pts_so_far, pts_so_far+num_in_trace)
                if is_continuous_color:
                    # Extract subset of rows in heatmap
                    resorted_ndces_this_trace = np.argsort(
                        np.array(trace_markers['color'])[hm_point_ndces_this_trace])
                    trace_markers['color'] = np.array(trace_markers['color'])[resorted_ndces_this_trace]
                    # TODO y_coords_this_trace = y_coords_this_trace[resorted_ndces_this_trace]
                else:
                    pts_so_far += num_in_trace
                    resorted_ndces_this_trace = np.arange(len(hm_point_names_this_trace))
                hm_point_names_this_trace = hm_point_names_this_trace[resorted_ndces_this_trace]
                hm_point_ndces_this_trace = hm_point_ndces_this_trace[resorted_ndces_this_trace]
            
            hm_row_ndces.extend(hm_point_ndces_this_trace)
            all_hm_point_names.extend(hm_point_names_this_trace)
            new_trace = {
                'name': trace['name'], 
                'x': np.zeros(num_in_trace), 
                'y': y_coords_this_trace, # all_hm_point_names, 
                'xaxis': 'x2', 
                'hoverinfo': 'text+name', 
                'text': hm_point_names_this_trace, 
                'mode': hmscat_mode, 
                'textposition': 'center left', 
                'textfont': building_block_divs.hm_font_macro, 
                'marker': trace_markers, 
                'selected': building_block_divs.style_selected, 
                'type': 'scatter'
            }
            row_scat_traces.append(new_trace) 
        # reorganize matrix if things were re-sorted.
        if resort_rows:
            fit_data = np.array([]) if len(hm_row_ndces) == 0 else fit_data[np.array(hm_row_ndces), :]
    return row_scat_traces, fit_data, all_hm_point_names


def hm_col_plot(
    fit_data, 
    reordered_groups=None, reordered_featnames=None, col_clustIDs=None, clustersep_coords=[]
):
    if reordered_groups is None:
        reordered_groups = reordered_featnames
    transform_col_order = np.arange(fit_data.shape[1])
    # Sort genes by color within each cluster, editing ordered_cols appropriately.
    for cid in np.unique(col_clustIDs):
        ndcesc = np.where(col_clustIDs == cid)[0]
        clust_colors = reordered_groups[ndcesc]
        new_order_perclust = np.argsort(clust_colors)
        transform_col_order[ndcesc] = transform_col_order[ndcesc][new_order_perclust]
    fit_data = fit_data[:, transform_col_order]
    reordered_groups = reordered_groups[transform_col_order]
    reordered_featnames = reordered_featnames[transform_col_order]
    
    col_scat_traces = []
    pt_text = []
    num_clusters = len(np.unique(col_clustIDs))
    widths = np.diff(clustersep_coords)    # len(widths) == num_clusters
    
#     new_trace = {
#         'name': 'sdfjak', 
#         'x': reordered_featnames, 
#         'y': np.zeros(len(reordered_featnames)), 
#         'yaxis': 'y2', 
#         'hoverinfo': 'text+name', 
#         'text': reordered_featnames, 
#         'mode': 'markers', 
#         'textposition': 'top center', 
#         'textfont': building_block_divs.hm_font_macro, 
#         'marker': { 'color': col_clustIDs }, 
#         'selected': building_block_divs.style_selected, 
#         'type': 'scatter'
#     }
#     col_scat_traces.append(new_trace)
    col_scat_traces.append({
        'y': np.ones(len(col_clustIDs)), 'yaxis': 'y2', 
        'x': reordered_featnames, 
        'width': 1, 'base': -0.5, 'offset': -0.5, 
        'hoverinfo': 'text', 'hovertext': reordered_featnames, 
        'text': reordered_featnames, 
        'hoverdistance': 40, 
        'colorbar': {
            'len': 0.3, 'thickness': 20, 
            'xanchor': 'left', 'yanchor': 'top', 
            'title': 'Expression', 'titleside': 'top', 'ticks': 'outside', 
            'titlefont': building_block_divs.colorbar_font_macro, 
            'tickfont': building_block_divs.colorbar_font_macro
        }, 
        'marker': { 'color': col_clustIDs }, 
        'type': 'bar'
    })
    return col_scat_traces, fit_data


def hm_hovertext(data, rownames, colnames):
    pt_text = []
    # First the rows, then the cols
    for r in range(data.shape[0]):
        trow = np.ravel(data[r].toarray())
        pt_text.append(["Cell {}".format(str(rownames[r])) for k in trow])
        for c in range(data.shape[1]):
            pt_text[r][c] += "<br>Gene: {}<br>Log expression: {}".format(str(colnames[c]), str(round(float(trow[c]), 2)))
    return pt_text

    
def display_heatmap_cb(
    hm_raw_data,    # 2D numpy array of selected data
    feat_names,     # col labels of hm_raw_data
    hm_point_names,    # (unique!) row labels of hm_raw_data
    scatter_fig,    # Scatterplot panel which this is mirroring.
    view_cocluster, 
    feat_group_names=None, 
    scatter_frac_domain=0.10, 
    scatter_frac_range=0.08, 
    show_legend=False, 
    normalize_mode='log'
):
    if np.prod(hm_raw_data.shape) == 0:
        return { 
            'data': [], 
            'layout': building_block_divs.create_hm_layout(
                scatter_frac_domain=scatter_frac_domain, scatter_frac_range=scatter_frac_range, 
                show_legend=show_legend
            )
        }
    fit_data = hm_raw_data
    if (normalize_mode == 'log') and (not app_config.params['hm_diverging']):
        fit_data = normalize_data(fit_data)
        intrst_mode = 'var'
        colorscale = app_config.cmap_custom_blackbody
        hm_colorvar_name = 'Log expr.'
    elif (normalize_mode == 'mult_dev'):
        fit_data = deviance_residuals(fit_data, binomial=False, sqrt=True)
        intrst_mode = 'deviance'
        colorscale = app_config.cmap_custom_ylbu_diverging
        hm_colorvar_name = 'Dev. residual'
    elif (normalize_mode == 'bin_dev'):
        fit_data = sp.sparse.csr_matrix(deviance_residuals(fit_data, binomial=True, sqrt=True))
        intrst_mode = 'deviance'
        colorscale = app_config.cmap_cubehelix
        hm_colorvar_name = 'Dev. residual'
    # Identify (interesting) genes to plot. Currently: high-variance genes
    feat_ndces = interesting_feat_ndces(fit_data, mode=intrst_mode)
    absc_labels = feat_names[feat_ndces]
    absc_group_labels = np.array(['Genes']*len(absc_labels)) if feat_group_names is None else feat_group_names[feat_ndces]
    fit_data = fit_data.tocsc()[:, feat_ndces]
    # Quantile normalize the data if necessary to better detect patterns.
    if app_config.params['hm_qnorm_plot']:
        qtiles = np.zeros_like(fit_data)
        nnz_ndces = np.nonzero(fit_data)
        qtiles[nnz_ndces] = sp.stats.rankdata(fit_data[nnz_ndces]) / len(fit_data[nnz_ndces])
        fit_data = qtiles
    # Spectral coclustering to cluster the heatmap. We always order rows (points) by spectral projection, but cols (features) can have different orderings for different viewing options.
    row_clustIDs = np.zeros(fit_data.shape[0])
    col_clustIDs = np.zeros(fit_data.shape[1])
    
    if (fit_data.shape[0] > 1):
        ordered_rows, ordered_cols, row_clustIDs, col_clustIDs = dm.compute_coclustering(fit_data.toarray())
        fit_data = fit_data[ordered_rows, :]
        hm_point_names = hm_point_names[ordered_rows]
    else:
        ordered_cols = np.arange(fit_data.shape[1])
    # assemble coordinates of lines adumbrating clusters.
    clustersep_line_coords = []
    for cid in np.unique(col_clustIDs):
        ndcesc = np.where(col_clustIDs == cid)[0]
        clustersep_line_coords.append(np.min(ndcesc) - 0.5)
    if len(fit_data.shape) > 1:
        clustersep_line_coords.append(fit_data.shape[1] - 0.5)
    
    fit_data = fit_data[:, ordered_cols]
    absc_labels = absc_labels[ordered_cols]
    if absc_group_labels is not None:
        absc_group_labels = absc_group_labels[ordered_cols]
    # Copy trace metadata from scatter_fig, in order of hm_point_names, to preserve colors etc.
    row_scat_traces, fit_data, hm_point_names = hm_row_scatter(
        fit_data, scatter_fig, hm_point_names, view_cocluster, row_clustIDs=row_clustIDs
    )
    col_scat_traces, fit_data = hm_col_plot(
        fit_data, 
        reordered_groups=absc_group_labels, reordered_featnames=absc_labels, 
        col_clustIDs=col_clustIDs, clustersep_coords=clustersep_line_coords
    )
    pt_text = hm_hovertext(fit_data, hm_point_names, absc_labels)
    hm_trace = {
        'z': fit_data.toarray(), 
        'x': absc_labels, 
        # 'y': hm_point_names, 
        'hoverinfo': 'text',
        'text': pt_text, 
        'colorscale': colorscale, 
        'zmin': 0, 
        'colorbar': {
            'len': 0.3, 
            'thickness': 20, 
            'xanchor': 'left', 
            'yanchor': 'top', 
            'title': hm_colorvar_name,
            'titleside': 'top',
            'ticks': 'outside', 
            'titlefont': building_block_divs.colorbar_font_macro, 
            'tickfont': building_block_divs.colorbar_font_macro
        }, 
        'type': 'heatmap'
    }
    if normalize_mode == 'mult_dev':
        max_magnitude = np.percentile(np.abs(fit_data.toarray()), 99) if fit_data.shape[0] > 0 else 2
        hm_trace['zmin'] = -max_magnitude
        hm_trace['zmax'] = max_magnitude
    elif normalize_mode == 'bin_dev':
        max_magnitude = np.percentile(np.abs(fit_data.toarray()), 99) if fit_data.shape[0] > 0 else 2
        hm_trace['zmax'] = max_magnitude
    
    return {
        'data': [ hm_trace ] + row_scat_traces + col_scat_traces, 
        'layout': building_block_divs.create_hm_layout(
            scatter_frac_domain=scatter_frac_domain, scatter_frac_range=scatter_frac_range, 
            show_legend=show_legend, clustersep_coords=clustersep_line_coords
        )
    }
