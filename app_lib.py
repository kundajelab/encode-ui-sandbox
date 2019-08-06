# Application-specific routines for working with (smallish matrices of) data.
# Author: Akshay Balsubramani

import base64, io, os, time, json, numpy as np, scipy as sp, pandas as pd, diffmap as dm
import app_config, building_block_divs
import goterm_caller
# from goatools.associations import read_ncbi_gene2go
# from goatools.test_data.genes_NCBI_9606_ProteinCoding import GENEID2NT
import re
import umap



# =======================================================
# ================== Utility functions ==================
# =======================================================

# Loads anndata from 3 files: {data, rowdata, coldata}
def read_as_anndata(path_pfx, transp=True):
    m = scipy.io.mmread(path_pfx + '_data.mtx')
    if transp:
        m = m.transpose()
    c_annot = pd.read_csv(path_pfx + '_mdata_points.csv', index_col=None, header=0)
    f_annot = pd.read_csv(path_pfx + '_mdata_feats.csv', index_col=None, header=0)
    return anndata.AnnData(X=m, obs=c_annot, var=f_annot)


# Writes anndata in the form of 3 dataframes: {data, rowdata, coldata}
def write_anndata(path_pfx, adata, transp=True):
    if transp:
        adata = adata.transpose()
    adata.obs.to_csv(path_pfx + '_mdata_points.csv', index=False)
    adata.var.to_csv(path_pfx + '_mdata_feats.csv', index=False)
    scipy.io.mmwrite(path_pfx + '_data.mtx', adata.X)


# Analysis mini-workflow for RNA-seq data
def fastfilt_rna(fdata, min_points_per_gene=25, max_expr_qtile_per_gene=99.5, scale_mrna_per_point=10000.0):
    mrna_per_point = np.ravel(fdata.sum(axis=1))
    mrna_per_gene = np.ravel(fdata.sum(axis=0))
    min_mrna_per_point = np.percentile(mrna_per_point, q=2)
    max_mrna_per_gene = np.percentile(mrna_per_gene, q=max_expr_qtile_per_gene)
    npoints, ngenes = fdata.shape
    point_mask = np.ones(npoints, dtype=bool)
    gene_mask = np.ones(ngenes, dtype=bool)
    point_mask = point_mask & (np.ravel(fdata.sum(axis=1)) >= min_mrna_per_point)
    gene_mask = gene_mask & (np.ravel((fdata > 0).sum(axis=0)) >= min_points_per_gene) 
    gene_mask = gene_mask & (np.ravel(fdata.sum(axis=0)) <= max_mrna_per_gene)
    print("Fraction of points kept:\t" + str(np.mean(point_mask)))
    print("Fraction of genes kept:\t" + str(np.mean(gene_mask)))
    return point_mask, gene_mask


def quantile_norm(dtd):
    qtiles = np.zeros(len(dtd))
    nnz_ndces = np.nonzero(dtd)[0]
    qtiles[nnz_ndces] = sp.stats.rankdata(dtd[nnz_ndces])/len(nnz_ndces)
    return qtiles


def reviz_embed_data(fit_data, alg='UMAP'):
    if alg == 'UMAP':
        reducer = umap.UMAP(n_neighbors=15, n_components=2, random_state=42)
        embedding = reducer.fit_transform(fit_data)
    return embedding



# =========================================================
# =================== Main scatter plot ===================
# =========================================================

"""
(Data, layout) for the main graph panel.
Color_var is either a field of the plotting df, or a numpy array.
"""


# TODO: Add legend groups as applicable, to bunch colors within a group

# Here selected_point_ids is a list of unique string IDs of points. 
def traces_scatter(
    data_df, 
    color_var, 
    colorscale, 
    selected_point_ids, 
    marker_size=app_config.params['marker_size'], 
    style_selected=building_block_divs.style_selected
):
    traces_list = []
    display_ndces = app_config.params['display_coordinates']
    cumu_color_dict = {}
    # Check to see if color_var is continuous or discrete and plot points accordingly
    if isinstance(color_var, (list, tuple, np.ndarray)):     # Color_var is an array, not a col index.
        # print('continuous: {}'.format(str(color_var)))
        continuous_color_var = color_var
        point_names = list(data_df[app_config.params['display_ID_var']])
        spoints = np.where(np.isin(point_names, selected_point_ids))[0]
        colorbar_title = app_config.params['hm_colorvar_name']
        if app_config.params['qnorm_plot']:
            continuous_color_var = quantile_norm(continuous_color_var)
            colorbar_title = 'Percentile'
        pt_text = ["{}<br>Quantile: {}".format(point_names[i], round(continuous_color_var[i], 3)) for i in range(len(point_names))]
        traces_list.append({ 
            'name': 'Data', 
            'x': data_df[display_ndces['x']], 
            'y': data_df[display_ndces['y']], 
            'selectedpoints': spoints, 
            'hoverinfo': 'text', 
            'text': pt_text, 
            'mode': 'markers', 
            'marker': {
                'size': marker_size, 
                'opacity': app_config.params['marker_opacity'], 
                'symbol': 'circle', 
                'showscale': True, 
                'colorbar': {
                    'len': 0.3, 
                    'thickness': 20, 
                    'xanchor': 'right', 
                    'yanchor': 'top', 
                    'title': colorbar_title,
                    'titleside': 'top',
                    'ticks': 'outside', 
                    'titlefont': building_block_divs.colorbar_font_macro, 
                    'tickfont': building_block_divs.colorbar_font_macro
                }, 
                'color': continuous_color_var, 
                'colorscale': colorscale
            }, 
            'selected': style_selected, 
            'type': 'scattergl'
        })
    else:    # Categorical color scheme, one trace per color
        cnt = 0
        for idx, val in data_df.groupby(color_var):
            point_ids_this_trace = list(val[app_config.params['display_ID_var']])
            spoint_ndces_this_trace = np.where(np.isin(point_ids_this_trace, selected_point_ids))[0]
            if app_config.params['legendgroup']:
                legendgroup = idx.split('_')[0]
            else:
                legendgroup = idx
            if legendgroup not in cumu_color_dict:
                trace_color = colorscale[cnt]
                cnt += 1
                cumu_color_dict[legendgroup] = trace_color
            trace_opacity = 1.0
            trace_info = {
                'name': str(idx), 
                'x': val[display_ndces['x']], 
                'y': val[display_ndces['y']], 
                'selectedpoints': spoint_ndces_this_trace, 
                'hoverinfo': 'text+name', 
                'text': point_ids_this_trace, 
                'mode': 'markers', 
                'opacity': trace_opacity, 
                'marker': {
                    'size': marker_size, 
                    'opacity': app_config.params['marker_opacity'] if str(idx) != 'Other' else app_config.params['bg_marker_opacity'], 
                    'symbol': 'circle', 
                    'color': trace_color
                }, 
                'legendgroup': legendgroup, 
                'selected': style_selected
            }
            if not app_config.params['three_dims']:
                trace_info.update({'type': 'scattergl'})
            else:
                trace_info.update({ 'type': 'scatter3d', 'z': val[display_ndces['z']] })
            traces_list.append(trace_info)
    return traces_list


def layout_scatter(annots):
    display_ndces = app_config.params['display_coordinates']
    new_layout = building_block_divs.create_scatter_layout(annots)
    return new_layout


def build_main_scatter(data_df, color_var, colorscale, highlight=False, 
                       marker_size=app_config.params['marker_size'], 
                       annots=[], selected_point_ids=[], 
                       style_selected = building_block_divs.style_selected
                      ):
    if highlight:
        style_selected['marker']['color'] = 'white'
    else:
        style_selected['marker'].pop('color', None)    # Remove color if exists
    trace_list = traces_scatter(
        data_df, 
        color_var, 
        colorscale, 
        selected_point_ids, 
        marker_size=marker_size, 
        style_selected=style_selected
    )
    return { 
        'data': trace_list, 
        'layout': layout_scatter(annots)
    }


"""
# TODO: Finish this function, which plots mean feature values over each heatmap cluster.

def generate_percluster_viz(raw_data, cell_cluster_list, cell_color_list, featID='Gene'):
    cluster_IDs = np.unique(cell_cluster_list)
    plot_values = np.zeros((len(cluster_IDs), raw_data.shape[1]))
    for i in range(len(cluster_IDs)):
        plot_values[i, :] = np.mean(raw_data[cell_cluster_list == cluster_IDs[i], :], axis=0)
    panel_layout = {
        'margin': { 'l': 0, 'r': 0, 'b': 0, 't': 0}, 
        'hovermode': 'closest', 
        'autosize': True,
        'xaxis': {
            'title': {'text': featID, 'font': building_block_divs.legend_font_macro }, 
            'tickfont': building_block_divs.legend_font_macro 
        }, 
        'yaxis': {
            'showticklabels': False,
            'automargin': True, 
            'ticks': 'outside', 
            'tickcolor': app_config.params['legend_font_color']
        },
        'plot_bgcolor': app_config.params['bg_color'],
        'paper_bgcolor': app_config.params['bg_color'], 
        'showlegend': True, 
        'legend': {
            'font': building_block_divs.legend_font_macro
        }
    }
    go_results = np.array(gp.gprofile(selected_genes))
    top_go_logpvals = np.array([])
    top_go_termnames = np.array([])
    top_go_dbIDs = np.array([])
    if go_results.shape[0] > 0:
        go_results = go_results[:topk, :]
        top_go_logpvals = -np.log10(go_results[:,2].astype(float))
        top_go_dbIDs = go_results[:,9]
        top_go_termnames = go_results[:,11]
    database_colors = { 'MF': '#CB3C19'}
    database_IDs = {'MF': 'Molecular function'}
    bar_colors = np.array([database_colors[x] for x in top_go_dbIDs])
    panel_data = []
    ordi = np.arange(len(top_go_dbIDs))[::-1]
    return {'data': panel_data, 'layout': panel_layout }
"""


#==========================================================================================================
#= Differential scatterplot analysis ======================================================================
#==========================================================================================================

def create_diff_scatter_panel(
    pointset1, 
    pointset2, 
    raw_data, 
    feat_names
):
    graph_layout = create_scatter_layout([])
    if (scatter_fig is None) or ('data' not in scatter_fig):
        return
    for trace in scatter_fig['data']:
        trace_markers = trace['marker']
        hm_point_names_this_trace = np.intersect1d(feat_names, trace['text'])
    graph_data = (pointset1, pointset2)
    return { 'data': graph_data, 'layout': graph_layout }



"""
Update GO enrichment panel.
https://biit.cs.ut.ee/gprofiler/page/apis. or 
g:GOSt API (in class header of gprofiler.py).

* ``all_results`` - (*Boolean*) All results, including those deemed not significant.
* ``ordered`` - (*Boolean*) Ordered query.
* ``exclude_iea`` - (*Boolean*) Exclude electronic GO annotations.
* ``underrep`` - (*Boolean*) Measure underrepresentation.
* ``evcodes`` - (*Boolean*) Request evidence codes in output as the
  final column.
* ``hier_sorting`` - (*Boolean*) Sort output into subgraphs.
* ``hier_filtering`` - (*Boolean*) Hierarchical filtering.
* ``max_p_value`` - (*Float*) Custom p-value threshold.
* ``min_set_size`` - (*Int*) Minimum size of functional category.
* ``max_set_size`` - (*Int*) Maximum size of functional category.
* ``min_isect_size`` - (*Int*) Minimum size of query / functional
  category intersection.
* ``max_isect_size`` - (*Int*) Maximum size of query / functional
  category intersection.
* ``correction_method`` - Algorithm used for multiple testing correction, one of:
  - ``GProfiler.THR_GSCS`` **Default** g:SCS.
  - ``GProfiler.THR_FDR`` Benjamini-Hochberg FDR.
  - ``GProfiler.THR_BONFERRONI`` Bonferroni.
* ``domain_size`` - Statistical domain size, one of:
  - ``GProfiler.DOMAIN_ANNOTATED`` - **Default** Only annotated genes.
  - ``GProfiler.DOMAIN_KNOWN`` - All known genes.
* ``custom_bg`` - (*String* | *List*) Custom statistical background
* ``src_filter`` - (*List*) A list of data source ID strings, e.g.
  ``["GO:BP", "KEGG"]``. These currently include GO (GO:BP, GO:MF,
  GO:CC to select a particular GO branch), KEGG, REAC, TF, MI, CORUM,
  HP, HPA, OMIM.
"""
def display_goenrich_panel_func(selected_genes, topk=20):
    panel_layout = {
        'margin': { 'l': 0, 'r': 0, 'b': 0, 't': 25}, 
        'hovermode': 'closest',
        'orientation': 90,
        'autosize': True,
        'xaxis': {
            'title': {'text': '-log(p)', 'font': building_block_divs.legend_font_macro }, 
            'tickfont': building_block_divs.legend_font_macro 
        }, 
        'yaxis': {
            'showticklabels': False,
            'automargin': True, 
            'ticks': 'outside', 
            'tickcolor': app_config.params['legend_font_color']
        },
        'plot_bgcolor': app_config.params['bg_color'],
        'paper_bgcolor': app_config.params['bg_color'], 
        'showlegend': True, 
        'legend': {
            'font': building_block_divs.legend_font_macro
        }
    }
    if len(selected_genes) == 0:
        return {'data': [], 'layout': panel_layout }
    go_results = goterm_caller.gprofiler(selected_genes)
    top_go_logpvals = np.array([])
    top_go_termnames = np.array([])
    top_go_dbIDs = np.array([])
    if (go_results is not None) and (go_results.shape[0] > 0):
        x = np.array(np.argsort(go_results['p.value']))
        go_results = go_results.iloc[x[:topk], :]
        top_go_logpvals = -np.log10(go_results['p.value'].astype(float))
        top_go_dbIDs = go_results['domain']
        top_go_termnames = go_results['term.name']
        top_go_termIDs = go_results['term.id']
        top_go_queryhits = go_results['intersection']
    database_colors = { 'MF': '#CB3C19', 'BP': '#FA9D18', 'CC': '#198520', 'keg': '#D9869B', 'rea': '#335ACC', 'tf': '#4F6E9B', 'mir': '#53B8AD', 'hpa': '#542DB1', 'cor': '#6AAB19', 'hp': '#8C0784'}
    database_IDs = {'MF': 'Molecular function', 'BP': 'Biological process', 'CC': 'Cellular component', 'keg': 'KEGG', 'rea': 'REAC', 'tf': 'TF', 'mir': 'MIRNA', 'hpa': 'HPA', 'cor': 'CORUM', 'hp': 'HP'}
    bar_colors = np.array([database_colors[x] for x in top_go_dbIDs])
    panel_data = []
    ordi = np.arange(len(top_go_dbIDs))[::-1]
    for c in np.unique(bar_colors):
        trace_locs = np.where(bar_colors == c)[0]
        trace_locs = trace_locs[::-1]      # Reverse because it's better.
        panel_data.append({
            'name': database_IDs[top_go_dbIDs[trace_locs[0]]], 
            'x': top_go_logpvals[trace_locs],
            # 'y': ordi[trace_locs], 
            'y': top_go_termIDs[trace_locs],
            'hovertext': [
                "-log10(p): {}<br>{}<br>{}".format(
                    str(round(top_go_logpvals[t], 2)), 
                    top_go_termnames[t], 
                    top_go_termIDs[t]
                ) for t in trace_locs], 
            'text': top_go_termnames[trace_locs], 
            'hoverinfo': 'text', 
            'insidetextfont': { 'family': 'sans-serif', 'color': 'white' }, 
            'outsidetextfont': { 'family': 'sans-serif', 'color': 'white' }, 
            'marker': { 'color': c },      
            'textposition': 'auto', 
            'orientation': 'h', 
            'type': 'bar'
        })
    return {'data': panel_data, 'layout': panel_layout }
