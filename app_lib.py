# Application-specific routines for working with (smallish matrices of) data.
# Author: Akshay Balsubramani

import base64, io, os, time, json, numpy as np, scipy as sp, pandas as pd, diffmap as dm
import app_config, building_block_divs
import re
import umap
import dash_table


# =======================================================
# ================== Utility functions ==================
# =======================================================

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
Here selected_point_ids is a list of unique string IDs of points.
"""
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
            if idx == 0:
                marker_size = 3.0
            else:
                marker_size = 9.0
            point_ids_this_trace = list(val[app_config.params['display_ID_var']])
            spoint_ndces_this_trace = np.where(np.isin(point_ids_this_trace, selected_point_ids))[0]
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
                    'opacity': 1.0,
                    'symbol': 'circle'#, 'color': trace_color
                },
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


def build_main_scatter(
    data_df, color_var, colorscale, highlight=False,
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
