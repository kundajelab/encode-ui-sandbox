#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runner for an interactive browser / analyzer.
Dash application (Dash [Python] <- Plotly <- React.js <- D3.js)
"""

import base64, io, os, time, json
import numpy as np, scipy as sp, pandas as pd, dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import app_config, app_lib, building_block_divs, interact_heatmap
import dash_html_components as html
# import matplotlib, matplotlib.pyplot as plt, matplotlib.colors as colors
"""
For more on jobs that take a while: set up workers https://github.com/WileyIntelligentSolutions/wiley-boilerplate-dash-app
"""



# =========================================================
# ================== Initialize Dash app ==================
# =========================================================

# Load gene embedded coordinates.
plot_data_df = pd.read_csv(app_config.params['plot_data_df_path'], sep="\t", index_col=False)
graph_adj = sp.sparse.load_npz(app_config.params['adj_mat_path'])
# data_ess = pd.read_csv(app_config.params['raw_data_path'], index_col=0, sep='\t')

point_names = np.array(plot_data_df[app_config.params['display_ID_var']])
additional_colorvars = []

assay_names = [x.split('_signal')[0] for x in plot_data_df['Assay'] if '_signal'  in x]
cell_types = np.array(plot_data_df['Biosample'])

app = dash.Dash(__name__)
if not app_config._DEPLOY_LOCALLY:
    app.config.update({'routes_pathname_prefix':'/ENCODE_databrowser/', 'requests_pathname_prefix':'/ENCODE_databrowser/'})

server=app.server
app.layout = building_block_divs.create_div_mainapp(
    point_names,
    assay_names=np.unique(assay_names),
    celltype_names=np.unique(cell_types),
    more_colorvars=additional_colorvars,
    align_options_list=['Unaligned', 'Aligned']
)



# ============================================================
# ===================== Callback helpers =====================
# ============================================================

"""
Utility function to take union of a list of selected subsets that have been stored.
Returns: dictionary. Keys: point IDs in the union. Values: {'pointIndex': p, 'curveNumber': v }.
"""
def union_of_selections(selected_subsets, subset_store):
    points_so_far = np.array([])
    dict_to_return = {}
    if selected_subsets is not None and isinstance(selected_subsets, (list, tuple)):
        for v in selected_subsets:
            new_point_ids = np.setdiff1d(list(subset_store[v].keys()), points_so_far)
            points_so_far = np.union1d(points_so_far, new_point_ids)
            for c in new_point_ids:
                dict_to_return[c] = subset_store[v][c]
    return dict_to_return


def make_store_points(selectedData_points):
    if (selectedData_points is not None) and ('points' in selectedData_points):
        toret = {}
        for p in selectedData_points['points']:
            toret[p['text']] = {
                # 'pointIndex': p['pointIndex'], 'curveNumber': p['curveNumber']
            }
        return toret
    else:
        return {}


# Inverse function of the above.
def make_selected(stored_dict):
    toret = { 'range': None }
    toret['points'] = [
        {
#             'pointIndex': stored_dict[k]['pointIndex'], 'curveNumber': stored_dict[k]['curveNumber'], 'pointNumber': stored_dict[k]['pointIndex'],
            'text': k
        } for k in stored_dict
    ]
    return toret


def get_cell_subsets(cell_IDs, cell_colors):
    return {ct: list(cell_IDs[cell_colors == ct]) for ct in np.unique(cell_colors)}


"""
Update the main graph panel with selected points annotated, using the given dataset.
"""
def highlight_landscape_func(
    annotated_points,
    data_df,
    point_names_to_use,
    marker_size=app_config.params['marker_size'],
    style_selected=building_block_divs.style_selected,
    color_var=app_config.params['default_color_var'], # Could be an array of continuous colors!
    colorscale=app_config.params['colorscale'],
    selectedpoint_ids=[],
    highlight_selected=False,
    absc_arr=None,
    ordi_arr=None
):
    annots = []
    looked_up_ndces = np.where(np.in1d(point_names_to_use, annotated_points))[0]
    for point_ndx in looked_up_ndces:
        absc = absc_arr[point_ndx]
        ordi = ordi_arr[point_ndx]
        cname = point_names_to_use[point_ndx]
        annots.append({
            'x': absc,
            'y': ordi,
            'xref': 'x', 'yref': 'y',
            # 'text': '<b>Cell {}</b>'.format(cname),
            'font': {
                'color': 'white',
                'size': 15
            },
            'arrowcolor': 'white',
            'showarrow': True,
            'arrowhead': 2, 'arrowwidth': 2, 'arrowsize': 2,
            'ax': 0,
            'ay': -50
        })
    toret = app_lib.build_main_scatter(
        data_df,
        color_var,
        colorscale,
        highlight=highlight_selected,
        annots=annots,
        selected_point_ids=selectedpoint_ids,
        marker_size=marker_size,
        style_selected=style_selected
    )
    return toret


def run_update_landscape(
    annotated_points,      # Selected points annotated
    subset_store,       # Store of selected point subsets.
    data_df,
    point_names,
    marker_size,
    style_selected,
    highlighted_points=[],
    nbr_points=[],
    diffmap_embed=False,
    color_scheme='Dummy'          # Feature(s) selected to plot as color.
):
    pointIDs_to_select = highlighted_points if (len(highlighted_points) > 0) else list(subset_store['_current_selected_data'].keys())
    if annotated_points is None:
        annotated_points = []
    if ~diffmap_embed:
        absc_arr = data_df[app_config.params['display_coordinates']['x']]
        ordi_arr = data_df[app_config.params['display_coordinates']['y']]
    else:
        absc_arr = data_df[app_config.params['display_coordinates_diffmap']['x']]
        ordi_arr = data_df[app_config.params['display_coordinates_diffmap']['y']]
    data_df[color_scheme][np.isin(point_names, nbr_points)] = 1

    print('Color scheme: {}'.format(color_scheme))
    # Check if a continuous feature is chosen to be plotted.
    if ((color_scheme != app_config.params['default_color_var']) and
        (color_scheme not in additional_colorvars) and
        (len(color_scheme) > 0)
       ):
        if not isinstance(color_scheme, (list, tuple)):
            color_scheme = [color_scheme]
        new_colors = 0
        return highlight_landscape_func(
            annotated_points,
            data_df,
            point_names,
            color_var=new_colors,
            colorscale=app_config.params['colorscale_continuous'],
            selectedpoint_ids=pointIDs_to_select,
            highlight_selected=True,
            absc_arr=absc_arr,
            ordi_arr=ordi_arr,
            marker_size=marker_size,
            style_selected=style_selected
        )
    else:    # color_scheme is a col ID indexing a discrete column.
        return highlight_landscape_func(
            annotated_points,
            data_df,
            point_names,
            color_var=color_scheme,
            colorscale=app_config.params['colorscale_discrete'],
            selectedpoint_ids=pointIDs_to_select,
            highlight_selected=False,
            absc_arr=absc_arr,
            ordi_arr=ordi_arr,
            marker_size=marker_size,
            style_selected=style_selected
        )


def parse_upload_contents(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    return decoded.decode('utf-8').splitlines()
    if 'json' in filename:
        return json.loads(decoded)
    elif 'csv' in filename:
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    elif 'xls' in filename:
        df = pd.read_excel(io.BytesIO(decoded))


#def get_contents_to_df(celltypes, assays):
#    data = {'celltypes': celltypes, 'assays':assays}
#    dataframe = pd.DataFrame(data)
#    return html.Table(
#        # Header
#        [html.Tr([html.Th(col) for col in dataframe.columns]) ] +
#        # Body
#        [html.Tr([
#            html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
#        ]) for i in range(min(len(dataframe), max_rows))]
#    )

# =================================================================
# =========================== Callbacks ===========================
# =================================================================


#@app.callback(
#    Output('datatable', 'data'),
#    [Input('stored-pointsets', 'data'),
#     Input('landscape-plot', 'clickData')]
#)
#def update_graph(
#    data_store, clickData
#):
#    df = get_contents_to_df(data_store)
#    return df.to_json()
#    toret = ""
#    for setname in data_store:
#        toret = toret + "{}\t{}\n".format(data_store[setname], setname)
#    return ""
#     "***SELECTED DATA***\n{}\nCLICK DATA: \n{}".format(
#         toret,
#         json.dumps(clickData, indent=2)
#     )


# https://community.plot.ly/t/download-raw-data/4700/7
@app.callback(
    Output('download-set-link', 'href'),
    [Input('stored-pointsets', 'data')]
)
def save_selection(subset_store):
    save_contents = '\n'.join(list(subset_store['_current_selected_data'].keys()))
    return "data:text/csv;charset=utf-8," + save_contents
#     save_contents = json.dumps(subset_store['_current_selected_data'], indent=4)
#     return "data:text/json;charset=utf-8," + save_contents
    """
    with open('tmp.txt', 'w+') as f:
        f.write('\n'.join(tosavecells))
    """


import urllib
# Downloads the currently selected dataframe before extracting and exporting the desired columns.
@app.callback(
    Output('download-layout-link', 'href'),
    [Input('sourcedata-select', 'value')]
)
def update_download_layout_link(sourcedata_select):
    dataset_names = app_config.params['dataset_options']
    ndx_selected = dataset_names.index(sourcedata_select) if sourcedata_select in dataset_names else 0
    data_df = pd.read_csv(app_config.params['plot_data_df_path'], sep="\t", index_col=False)
    coords_to_load = list(app_config.params['display_coordinates'].values()) + ['gene_names']
    for c in coords_to_load:
        if c not in data_df:
            return ""
    csvString = data_df[coords_to_load].to_csv(sep="\t", index=False, encoding='utf-8')
    return "data:text/csv;charset=utf-8," + csvString


@app.callback(
    Output('display-selected-marker-size-factor', 'children'),
    [Input('slider-selected-marker-size-factor', 'value')]
)
def update_bgmarker_size(sel_size):
    return 'Selected marker size: {}'.format(sel_size)


@app.callback(
    Output('display-marker-size-factor', 'children'),
    [Input('slider-marker-size-factor', 'value')]
)
def update_marker_size(marker_size):
    return 'Marker size: {}'.format(marker_size)


### Updating datatable
@app.callback(
    Output('datatable', 'children'),
    [Input('landscape-plot', 'clickData')]
)
def display_results(plot_click_data):
    if plot_click_data is None:
        return ""
    expname = plot_click_data['points'][0]['text']
    exp_ndx = np.where(point_names == expname)[0][0]
    nbrs = point_names[np.where(np.ravel(graph_adj[exp_ndx,:].toarray()) != 0)[0]]
    assays = [s.rsplit('_', 1)[0] for s in nbrs]
    celltypes = [s.rsplit('_', 1)[1] for s in nbrs]
    data = {'celltype': celltypes, 'assays':assays}
    df = pd.DataFrame(data)
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in df.columns]) ] +
        # Body
        [html.Tr([
            html.Td(df.iloc[i][col], style={'color':'white'}) for col in df.columns
        ]) for i in range(len(df))]
    )
    #return get_contents_to_df(celltypes, assays)
    #return df
    #return "\n".join(nbrs)


"""
Update the main graph panel.
"""
@app.callback(
    Output('landscape-plot', 'figure'),
    [Input('points_annot', 'value'),
     Input('stored-pointsets', 'data'),
     Input('sourcedata-select', 'value'),
     Input('slider-marker-size-factor', 'value'),
     Input('slider-selected-marker-size-factor', 'value'),
     Input('assay_select', 'value'),
     Input('celltype_select', 'value'),
     Input('landscape-plot', 'clickData')]
)
def update_landscape(
    annotated_points,      # Selected points annotated
    subset_store,          # Store of selected point subsets.
    sourcedata_select,
    marker_size,
    sel_marker_size,
    assay_selected,
    celltype_selected,
    plot_click_data
):
    dataset_names = app_config.params['dataset_options']
    ndx_selected = dataset_names.index(sourcedata_select) if sourcedata_select in dataset_names else 0
    data_df = pd.read_csv(app_config.params['plot_data_df_path'], sep="\t", index_col=False)
    style_selected = building_block_divs.style_selected
    style_selected['marker']['size'] = sel_marker_size

    if plot_click_data is not None:
        expname = plot_click_data['points'][0]['text']
        exp_ndx = np.where(point_names == expname)[0][0]
        # Now get top few neighbors
        nbrs = point_names[np.where(np.ravel(graph_adj[exp_ndx,:].toarray()) != 0)[0]]
        data_df['Dummy'][np.isin(point_names, point_names[exp_ndx])] = 1
        data_df['Dummy'][np.isin(point_names, nbrs)] = 2
    else:
        nbrs = []

    curr_highlighted = []
    #curr_highlighted = np.where(np.logical_or(assay_names == assay_selected, cell_types == celltype_selected))[0]
    data_df['Dummy'][np.array(assay_names) == assay_selected] = 3
    data_df['Dummy'][cell_types == celltype_selected] = 4

    return run_update_landscape(
        annotated_points,
        subset_store,
        data_df,
        point_names,
        marker_size,
        style_selected,
        nbr_points=nbrs,
        highlighted_points=curr_highlighted    # List of most recently highlighted cells.
    )



# =======================================================
# ===================== Run the app =====================
# =======================================================

if __name__ == '__main__':
    app.run_server(port=8052, debug=True)
