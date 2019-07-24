# -*- coding: utf-8 -*-
"""
Runner for an interactive browser / analyzer.
Dash application (Dash [Python] <- Plotly <- React.js <- D3.js)
Author: Akshay Balsubramani
"""

import base64, io, os, time, json
import numpy as np, scipy as sp, pandas as pd, dash
from dash.dependencies import Input, Output, State
import app_config, app_lib, building_block_divs, learning_algs, interact_heatmap, goterm_caller



# =========================================================
# ================== Initialize Dash app ==================
# =========================================================

# Load gene embedded coordinates.
plot_data_df = pd.read_csv(app_config.params['plot_data_df_path'][0], sep="\t", index_col=False)
# graph_adj = sp.sparse.load_npz(app_config.params['adj_mat_path'])

point_names = np.array(plot_data_df[app_config.params['display_ID_var']])
additional_colorvars = []

feat_names = np.load(app_config.params['feat_names_path'], allow_pickle=True)

# raw_data = np.load(app_config.params['raw_datamat_path'])
raw_data = sp.sparse.load_npz(app_config.params['raw_datamat_path'])


app = dash.Dash(__name__)
if not app_config._DEPLOY_LOCALLY:
    app.config.update({'routes_pathname_prefix':'/encode_ui/', 'requests_pathname_prefix':'/encode_ui/'})

server=app.server
app.layout = building_block_divs.create_div_mainapp(
    point_names, 
    feat_names, 
    more_colorvars=additional_colorvars
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
    toret['points'] = [{ 'text': k } for k in stored_dict]
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
    three_dim_plot=False, 
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
        point_names, 
        colorscale, 
        highlight=highlight_selected, 
        annots=annots, 
        three_dim_plot=three_dim_plot, 
        selected_point_ids=selectedpoint_ids, 
        marker_size=marker_size, 
        style_selected=style_selected
    )
    return toret


def run_update_landscape(
    color_scheme,          # Feature(s) selected to plot as color.
    annotated_points_list,      # Selected points annotated
    subset_store,       # Store of selected point subsets.
    data_df, 
    point_names, 
    raw_data_to_use, 
    marker_size, 
    style_selected, 
    highlighted_points=[], 
    three_dim_plot=False
):
    # Here is where what is rendered as selected in the landscape differs from _current_selected_data
    pointIDs_to_select = highlighted_points if (len(highlighted_points) > 0) else list(subset_store['_current_selected_data'].keys())
    absc_arr = data_df[app_config.params['display_coordinates']['x']]
    ordi_arr = data_df[app_config.params['display_coordinates']['y']]
    
    print('Color scheme: {}'.format(color_scheme))
    # Check if a continuous feature is chosen to be plotted.
    if ((color_scheme != app_config.params['default_color_var']) and 
        (color_scheme not in additional_colorvars) and 
        (len(color_scheme) > 0)
       ):
        if not isinstance(color_scheme, (list, tuple)):
            color_scheme = [color_scheme]
        feat_ndces = np.isin(feat_names, color_scheme)
        # If there are multiple continuous features given, plot their mean; otherwise the mean's a no-op anyway so it's safe to use.
        if sp.sparse.issparse(raw_data_to_use):
            new_colors = np.squeeze(np.array(raw_data_to_use[:, feat_ndces].mean(axis=1)))
        else:
            new_colors = np.mean(raw_data_to_use[:, feat_ndces], axis=1)
        return highlight_landscape_func(
            annotated_points_list, 
            data_df, 
            point_names, 
            color_var=new_colors, 
            colorscale=app_config.params['colorscale_continuous'], 
            selectedpoint_ids=pointIDs_to_select, 
            highlight_selected=True, 
            absc_arr=absc_arr, 
            ordi_arr=ordi_arr, 
            three_dim_plot=three_dim_plot, 
            marker_size=marker_size, 
            style_selected=style_selected
        )
    else:    # color_scheme is a col ID indexing a discrete column.
        return highlight_landscape_func(
            annotated_points_list, 
            data_df, 
            point_names, 
            color_var=color_scheme, 
            colorscale=app_config.params['colorscale_discrete'], 
            selectedpoint_ids=pointIDs_to_select, 
            highlight_selected=True, 
            absc_arr=absc_arr, 
            ordi_arr=ordi_arr, 
            three_dim_plot=three_dim_plot, 
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

        

# =================================================================
# =========================== Callbacks ===========================
# =================================================================


@app.callback(
    Output('test-select-data', 'children'),
    [Input('stored-panel-settings', 'data'), 
     Input('stored-landscape-selected', 'data'), 
     Input('stored-pointsets', 'data'), 
     Input('main-heatmap', 'selectedData'), 
     Input('stored-most-recently-highlighted', 'data')]
)
def display_test(
    panel_data, 
    sel_data, 
    data_store, 
    hmsel_data, 
    hlight_store
):
    see_hm = "0" if hmsel_data is None else str(len(hmsel_data['points']))
    see_sel = "0" if sel_data is None else str(len(sel_data))
    see_hlight = "{}\t{}".format(hlight_store['_last_panel_highlighted'], len(hlight_store.keys()) - 1)
    toret = ""
    for setname in data_store:
        toret = toret + "{}\t{}\n".format(data_store[setname], setname)
    if panel_data['debug_panel']:
        return "***STORED SELECTED DATA***\n{}\n***Landscape SELECTED DATA***\n{}\n***Heatmap SELECTED DATA***\n{}\n***Most recently used panel:\t{}".format(
            toret, 
            see_sel, 
            see_hm, 
            see_hlight
        )
    else:
        return ""


@app.callback(
    Output('goenrich-panel', 'figure'),
    [Input('stored-featsets', 'data'), 
     Input('select-topk-goterms', 'n_submit'),
     Input('select-topk-goterms', 'n_blur')],
    [State('select-topk-goterms', 'value')]
)
def display_goenrich_panel(featset_store, dummy1, dummy2, topk):
    return goterm_caller.display_goenrich_panel_func(
        list(featset_store['_current_selected_feats'].keys()), topk=int(topk), 
        legend_font_style=building_block_divs.legend_font_macro
    )


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


# Render selectable point subsets.
@app.callback(
    Output('list-pointsets', 'options'), 
    [Input('stored-pointsets', 'data')]
)
def update_subset_options(stored_setlist):
    toret = []    # This could be done less readably with a nested list comprehension.
    blacklist = ['_current_subselected_data', '_current_selected_data']
    if stored_setlist is not None:
        for s in stored_setlist.keys():
            if s not in blacklist:
                toret.append({'label': s, 'value': s})
    return toret


# Render selectable point subsets.
@app.callback(
    Output('list-featsets', 'options'), 
    [Input('stored-featsets', 'data')]
)
def update_featset_options(stored_setlist):
    toret = []    # This could be done less readably with a nested list comprehension.
    blacklist = ['_current_subselected_feats', '_current_selected_feats']
    if stored_setlist is not None:
        for s in []:
            if s not in blacklist:
                toret.append({'label': s, 'value': s})
    return toret


import urllib
# Downloads the currently selected dataframe before extracting and exporting the desired columns.
@app.callback(
    Output('download-layout-link', 'href'), 
    [Input('sourcedata-select', 'value')]
)
def update_download_layout_link(sourcedata_select):
    dataset_names = app_config.params['dataset_options']
    ndx_selected = dataset_names.index(sourcedata_select) if sourcedata_select in dataset_names else 0
    data_df = pd.read_csv(app_config.params['plot_data_df_path'][ndx_selected], sep="\t", index_col=False)
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


@app.callback(
    Output('stored-landscape-selected', 'data'), 
    [Input('landscape-plot', 'selectedData')]
)
def update_stored_landscape_data(landscape_data):
    return make_store_points(landscape_data)


@app.callback(
    Output('stored-heatmap-selected', 'data'), 
    [Input('main-heatmap', 'selectedData')], 
    [State('stored-heatmap-selected', 'data')]
)
def update_stored_heatmap_data(hm_selection, old_hm_selcells):
    hmd = make_store_points(hm_selection)
    return hmd if ((len(hmd) > 0) and (list(hmd.keys())[0] in point_names)) else old_hm_selcells


@app.callback(
    Output('landscape-plot', 'selectedData'), 
    [Input('stored-pointsets', 'data')]
)
def update_selected_landscape_data(subset_store):
    return make_selected(subset_store['_current_selected_data'])


@app.callback(
    Output('last-loaded', 'data'), 
    [Input('list-pointsets', 'value')]
)
def update_last_loaded_subsets(selected_subsets):
    return selected_subsets


@app.callback(
    Output('last-loaded-feats', 'data'), 
    [Input('list-featsets', 'value')]
)
def update_last_loaded_subsets(selected_subsets):
    return selected_subsets


@app.callback(
    Output('last-stored-subsetname', 'data'), 
    [Input('pointset-name', 'value')]
)
def update_last_stored_subset_name(stored_subset_name):
    return stored_subset_name


@app.callback(
    Output('last-stored-featsname', 'data'), 
    [Input('featset-name', 'value')]
)
def update_last_stored_featset_name(stored_subset_name):
    return stored_subset_name


# The following determines what was highlighted most recently in an auxiliary panel.
@app.callback(
    Output('stored-most-recently-highlighted', 'data'), 
    [Input('stored-landscape-selected', 'modified_timestamp'), 
     Input('stored-heatmap-selected', 'modified_timestamp')], 
    [State('stored-landscape-selected', 'data'), 
     State('stored-heatmap-selected', 'data')]
)
def update_most_recently_highlighted_cells(
    landscape_time, 
    hm_time, 
    landscape_selection, 
    hm_selection
):
    recent_times = np.array([ int(landscape_time), int(hm_time), 0 ])
    most_recent_time = max(recent_times)
    if most_recent_time == int(landscape_time):
        rec_panel = 'landscape'
        data_selection = landscape_selection
    elif most_recent_time == int(hm_time):
        rec_panel = 'heatmap'
        data_selection = hm_selection
    # TODO add more aux panels here. keys should be _last_panel_highlighted as well as cell IDs in highlighted subset.
    data_selection['_last_panel_highlighted'] = rec_panel
    return data_selection


@app.callback(
    Output('diff-heatmap-panel', 'figure'), 
    [Input('diff-foreset-select', 'value'), 
     Input('diff-backset-select', 'value'), 
     Input('main-heatmap', 'figure')],  
    [State('stored-featsets', 'data'), 
     State('stored-pointsets', 'data')]
)
def update_diff_heatmap_panel(
    subsets_fg, 
    subsets_bg, 
    fig_hm_main, 
    featset_store, 
    subset_store
):
    fg_ids = list(union_of_selections(subsets_fg, subset_store).keys())
    bg_ids = list(union_of_selections(subsets_bg, subset_store).keys())
    panel_data = []
    hmfig_data = fig_hm_main['data']
    if (len(fg_ids) > 0) and (len(bg_ids) > 0) and (len(hmfig_data) > 0):
        raw_data_to_use = raw_data
        fg_ndces = np.isin(point_names, fg_ids)
        fg_data = raw_data_to_use[fg_ndces, :]
        fg_point_names = point_names[fg_ndces]
        bg_ndces = np.isin(point_names, bg_ids)
        bg_data = raw_data_to_use[bg_ndces, :]
        bg_point_names = point_names[bg_ndces]
        if sp.sparse.issparse(raw_data_to_use):
            fg_data = fg_data.toarray()
            bg_data = bg_data.toarray()
        # Now just isolate selected genes
        selgenes = list(featset_store['_current_selected_feats'].keys())
        selected_feats = np.isin(feat_names, selgenes)
        fg_data = fg_data[:, selected_feats]
        bg_data = bg_data[:, selected_feats]
        fg_trace = {
            'z': fg_data, 
            'x': selgenes, 'yaxis': 'y2', 
            'hoverinfo': 'text', 
            'hovertext': interact_heatmap.hm_hovertext(fg_data, fg_point_names, selgenes),
            'colorscale': hmfig_data[0]['colorscale'], 
            'zmin': hmfig_data[0]['zmin'], 'zmax': hmfig_data[0]['zmax'], 
            'colorbar': {
                'len': 0.3, 'thickness': 20, 
                'xanchor': 'left', 'yanchor': 'top', 
                'titleside': 'top',
                'ticks': 'outside', 
                'titlefont': building_block_divs.colorbar_font_macro, 
                'tickfont': building_block_divs.colorbar_font_macro
            }, 
            'type': 'heatmap'
        }
        panel_data.append(fg_trace)
        bg_trace = {
            'z': bg_data, 
            'x': selgenes, 
            'hoverinfo': 'text', 
            'hovertext': interact_heatmap.hm_hovertext(bg_data, bg_point_names, selgenes),
            'colorscale': hmfig_data[0]['colorscale'], 
            'zmin': hmfig_data[0]['zmin'], 'zmax': hmfig_data[0]['zmax'], 
            'colorbar': {
                'len': 0.3, 'thickness': 20, 
                'xanchor': 'left', 'yanchor': 'top', 
                'titleside': 'top',
                'ticks': 'outside', 
                'titlefont': building_block_divs.colorbar_font_macro, 
                'tickfont': building_block_divs.colorbar_font_macro
            }, 
            'type': 'heatmap'
        }
        panel_data.append(bg_trace)
    scatter_frac_domain = 0.10
    panel_layout = building_block_divs.create_diff_hm_layout(scatter_frac_domain=scatter_frac_domain)
    return {'data': panel_data, 'layout': panel_layout }


@app.callback(
    Output('num-selected-counter', 'children'), 
    [Input('stored-pointsets', 'data')]
)
def update_numselected_counter(subset_store):
    num_selected = len(subset_store['_current_selected_data'])
    return 'Cells selected: {}'.format(num_selected)


@app.callback(
    Output('numfeats-selected-counter', 'children'), 
    [Input('stored-featsets', 'data')]
)
def update_numfeatselected_counter(subset_store):
    num_selected = len(subset_store['_current_selected_feats'])
    return 'Genes selected: {}'.format(num_selected)


@app.callback(
    Output('display-genelist', 'value'), 
    [Input('stored-featsets', 'data')]
)
def update_gogenes_display(subset_store):
    toret_list = []
    for x in subset_store['_current_selected_feats']:
        toret_list.append(x)
    return ', '.join(toret_list)


"""
Updates the stored dictionary of saved subsets. 
Contains control logic for subset selection and storage.
"""
@app.callback(
    Output('stored-pointsets', 'data'), 
    [Input('last-loaded', 'modified_timestamp'), 
     Input('stored-most-recently-highlighted', 'modified_timestamp'), 
     Input('last-stored-subsetname', 'modified_timestamp'), 
     Input('upload-pointsets', 'contents')], 
    [State('last-loaded', 'data'), 
     State('stored-most-recently-highlighted', 'data'), 
     State('upload-pointsets', 'filename'), 
     State('last-stored-subsetname', 'data'), 
     State('stored-pointsets', 'data')]
)
def update_subset_storage(
    load_time, 
    highlighted_time, 
    memstore_time, 
    file_contents, 
    load_value_subsetIDs, 
    highlighted_data, 
    file_paths, 
    newset_name, 
    subset_store
):
    new_sets_dict = {} if subset_store is None else subset_store
    recent_times = np.array([ int(load_time), int(highlighted_time), int(memstore_time), 0 ])
    most_recent_time = max(recent_times)
    if most_recent_time == int(load_time):     # Update _current_selected_data by loading subsets.
        new_sets_dict['_current_selected_data'] = union_of_selections(load_value_subsetIDs, subset_store)
    elif most_recent_time == int(highlighted_time):   # Update _current_selected_data from the main plot / heatmap.
        # Logic to display points as selected from an auxplot (most recently used for selection).
        last_hlight_panel = highlighted_data.pop('_last_panel_highlighted', None)
        if last_hlight_panel == 'landscape':
            new_sets_dict['_current_selected_data'] = highlighted_data
        elif last_hlight_panel == 'heatmap':
            new_sets_dict['_current_selected_data'] = highlighted_data
    elif most_recent_time == int(memstore_time):        # Store current selected data as a new set
        if ((newset_name is not None) and 
            (newset_name not in new_sets_dict) and 
            (newset_name != '')):
            new_sets_dict[newset_name] = new_sets_dict['_current_selected_data']
    # By default load all color annotations as well.
    for color in [app_config.params['default_color_var']] + additional_colorvars:
        new_subsets = get_cell_subsets(point_names, plot_data_df[color])
        for name in new_subsets:
            new_sets_dict[name] = { x: {} for x in new_subsets[name]}
    # Load a bunch of cell sets with names derived from their filenames.
    if file_contents is not None and len(file_contents) > 0:
        for contents, fname in zip(file_contents, file_paths):   # fname here is a relative (NOT an absolute) file path
            fname_root = fname.split('/')[-1].split('.')[0]
            new_sets_dict[fname_root] = { x: {} for x in parse_upload_contents(contents, fname) }
    return new_sets_dict


@app.callback(
    Output('heatmap-arrays', 'data'), 
    [Input('stored-pointsets', 'data')], 
    [State('heatmap-arrays', 'data')]
)
def update_stored_hm_cells(
    subset_store, 
    old_heatmap_stored, 
    num_points_to_sample=900
):
    toret = old_heatmap_stored
    point_names_to_use = point_names
    pointIDs_to_display = list(subset_store['_current_selected_data'].keys())
    # Subsample down to a max #points, for smoothly interactive heatmap display.
    if len(pointIDs_to_display) > num_points_to_sample:
        pointIDs_to_display = np.random.choice(pointIDs_to_display, num_points_to_sample, replace=False)
    point_ndces_to_display = np.where(np.isin(point_names_to_use, pointIDs_to_display))[0]
    toret['ndces_to_select'] = list(point_ndces_to_display)
    return toret


"""
Update the main graph panel.
"""
@app.callback(
    Output('landscape-plot', 'figure'), 
    [Input('landscape-color', 'value'), 
     Input('points_annot', 'n_submit'), 
     Input('stored-pointsets', 'data'), 
     Input('sourcedata-select', 'value'), 
     Input('stored-most-recently-highlighted', 'data'), 
     Input('slider-marker-size-factor', 'value'), 
     Input('slider-selected-marker-size-factor', 'value'), 
     Input('main-landscape-dimension', 'value')], 
    [State('points_annot', 'value')]
)
def update_landscape(
    color_scheme,          # Feature(s) selected to plot as color.
    annotated_points_submit, 
    subset_store,          # Store of selected point subsets.
    sourcedata_select, 
    aux_highlighted, 
    marker_size, 
    sel_marker_size, 
    dim_landscape_plot, 
    annotated_points,      # Selected points annotated
):
    dataset_names = app_config.params['dataset_options']
    ndx_selected = dataset_names.index(sourcedata_select) if sourcedata_select in dataset_names else 0
    data_df = pd.read_csv(app_config.params['plot_data_df_path'][ndx_selected], sep="\t", index_col=False)
    style_selected = building_block_divs.style_selected
    style_selected['marker']['size'] = sel_marker_size
    annotated_points_list = annotated_points.split(', ')

    recently_highlighted = [x for x in aux_highlighted.keys()]
    recently_highlighted.remove('_last_panel_highlighted')    
    return run_update_landscape(
        color_scheme, 
        annotated_points_list, 
        subset_store, 
        data_df, 
        point_names, 
        raw_data, 
        marker_size, 
        style_selected, 
        highlighted_points=recently_highlighted,     # List of most recently highlighted cells.
        three_dim_plot=(dim_landscape_plot == '3D')
    )



# =======================================================
# ===================== Run the app =====================
# =======================================================

if __name__ == '__main__':
    app.run_server(port=8052, debug=True)
