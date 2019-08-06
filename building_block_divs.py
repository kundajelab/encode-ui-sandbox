# -*- coding: utf-8 -*-
# Generic front-end code for making biological dataset browsing interfaces.
# Author: Akshay Balsubramani
"""
This module contains the following parts:
- Component styles: for panels
- Component divs: save/load
- Component divs: heatmap
- Component divs
- Component divs: differential signal
- Component divs: cosmetic panel (debugging)
- Aggregation of component divs
"""


import dash_core_components as dcc, dash_html_components as html
import app_config, numpy as np


# ================================================================
# =================== Component styles/layouts ===================
# ================================================================


legend_font_macro = {
    'family': 'sans-serif', 
    'size': app_config.params['legend_font_size'], 
    'color': app_config.params['legend_font_color'] 
}

colorbar_font_macro = {
    'family': 'sans-serif', 
    'size': 8, 
    'color': app_config.params['legend_font_color'] 
}

hm_font_macro = {
    'family': 'sans-serif', 
    'size': 8, 
    'color': app_config.params['legend_font_color'] 
}

style_unselected = {
    'marker': {
        'size': 3.0, 
        'opacity': 0.7
    }
}

style_selected = {
    'marker': {
        'size': 6.0, 
        'opacity': 122.2
    }
}

style_outer_dialog_box = {
    # 'user-select': 'none', '-moz-user-select': 'none', '-webkit-user-select': 'none', '-ms-user-select': 'none', 
    'padding': 5, 
    # 'margin': 5, 
    # 'borderRadius': 5, 
    'border': 'thin lightgrey solid'
}

style_invis_dialog_box = {
    # 'user-select': 'none', '-moz-user-select': 'none', '-webkit-user-select': 'none', '-ms-user-select': 'none', 
    'padding': 5, 
    'margin': 5
}

style_hm_colorbar = {
    'len': 0.3, 
    'thickness': 20, 
    'xanchor': 'left', 
    'yanchor': 'top', 
    'title': app_config.params['hm_colorvar_name'], 
    'titleside': 'top', 
    'ticks': 'outside', 
    'titlefont': legend_font_macro, 
    'tickfont': legend_font_macro
}

style_text_box = {
    'textAlign': 'center', 
    'width': '100%', 
    'color': app_config.params['font_color']
}

style_upload = {
    'width': '100%', 
    'border': 'thin lightgrey solid',
    'textAlign': 'center', 
    'color': app_config.params['font_color'], 
    'padding-top': '5px', 
    'padding-bottom': '5px'
}

style_legend = {
    'font': legend_font_macro, 
    # bgcolor=app_config.params['legend_bgcolor'], 
    # 'borderwidth': app_config.params['legend_borderwidth'], 
    'padding': 0, 
    'margin': 0, 
    'border': 'thin lightgrey solid', 
    'traceorder': 'normal', 
    'orientation': 'h'
}


def create_hm_layout(
    scatter_frac_domain=0.10, scatter_frac_range=0.08, 
    show_legend=False, clustersep_coords=[]
):
    shape_list = []
    for x in clustersep_coords:
        shape_list.append({
            'type': 'line',
            'x0': x, 'x1': x, 'y0': -1.0, 'y1': 1.0, 'yref': 'y2', 
            'line': { 'color': 'white', 'width': 3 }
        })
    hm_layout = {
        'annotations': [{
                'x': 0.5, 'y': 1.05, 'showarrow': False, 
                'font': { 'family': 'sans-serif', 'size': 15, 'color': app_config.params['legend_font_color'] }, 
                'text': 'Genes',
                'xref': 'paper', 'yref': 'paper'
            }, 
            {
                'x': 0.0, 'y': 0.5, 'showarrow': False, 
                'font': { 'family': 'sans-serif', 'size': 15, 'color': app_config.params['legend_font_color'] }, 
                'text': 'Cells', 
                'textangle': -90, 
                'xref': 'paper', 'yref': 'paper'
            }
        ], 
        'margin': { 'l': 0, 'r': 0, 'b': 0, 't': 30 }, 
        'clickmode': 'event+select',  # https://github.com/plotly/plotly.js/pull/2944/
        'hovermode': 'closest', 
        'uirevision': 'Default dataset', 
        'xaxis': {
            'showticklabels': True, 'side': 'top', 
            'tickcolor': app_config.params['legend_bgcolor'], 
            'tickfont': { 'family': 'sans-serif', 'size': app_config.params['hm_font_size'], 'color': app_config.params['legend_font_color'] }, # 'dtick': 1, 
            'showgrid': False, 'showline': False, 'zeroline': False, 'visible': False, 
            'domain': [scatter_frac_domain, 1]
        }, 
        'xaxis2': {
            'showgrid': False, 'showline': False, 'zeroline': False, 'visible': False, 
            'domain': [0, scatter_frac_domain], 
            'range': [-1, 0.2]
        }, 
        'yaxis': {
            'automargin': True, 
            'showticklabels': False, 
            'showgrid': False, 'showline': False, 'zeroline': False, 'visible': False, 
            'domain': [0, 1-scatter_frac_range] 
        }, 
        'yaxis2': {
            'showgrid': False, 'showline': False, 'zeroline': False, 'visible': False, 
            'domain': [1-scatter_frac_range, 1]
        }, 
        'legend': style_legend, 
        'showlegend': show_legend, 
        'plot_bgcolor': app_config.params['bg_color'], 
        'paper_bgcolor': app_config.params['bg_color'], 
        'shapes': shape_list
    }
    return hm_layout


def create_scatter_layout(annotations):
    return {
        'margin': { 'l': 0, 'r': 0, 'b': 0, 't': 20}, 
        'clickmode': 'event+select',  # https://github.com/plotly/plotly.js/pull/2944/
        'hovermode': 'closest', 
        'uirevision': 'Default dataset',     # https://github.com/plotly/plotly.js/pull/3236
        'xaxis': {
            'automargin': True, 
            'showticklabels': False, 
            'showgrid': False, 'showline': False, 'zeroline': False, #'visible': False, 
            'style': {'display': 'none'}
        }, 
        'yaxis': {
            'automargin': True, 
            'showticklabels': False, 
            'showgrid': False, 'showline': False, 'zeroline': False, #'visible': False, 
            'style': {'display': 'none'}
        }, 
        'legend': style_legend, 
        'annotations': annotations, 
        'plot_bgcolor': app_config.params['bg_color'], 
        'paper_bgcolor': app_config.params['bg_color']
    }



# =================================================================
# =================== Component divs: save/load ===================
# =================================================================


div_pointset_select = html.Div(
    className='row', 
    children=[
        html.Div(
            className='row', 
            children=[
                html.Div(
                    id='num-selected-counter', 
                    className='four columns', 
                    children='Cells selected: 0', 
                    style={
                        'textAlign': 'center', 
                        'color': app_config.params['font_color'], 
                        'padding-top': '10px'
                    }
                ), 
                html.Div(
                    className='eight columns', 
                    children=[
                        dcc.Input(
                            id='pointset-name', 
                            type='text', 
                            placeholder="Store cell set with name...", 
                            value='', 
                            debounce=True, 
                            n_submit_timestamp='0', 
                            style={'width': '100%'}
                        )], 
                    style={'padding-top': '5px'}
                )]
        ), 
        html.Div(
            className='row', 
            children=[
                dcc.Dropdown(
                    id='list-pointsets', 
                    placeholder="Select stored cell sets to load", 
                    multi=True, 
                    value=[]
                )], 
            style={'padding-top': '5px'}
        )]
)


div_featset_select = html.Div(
    className='row', 
    children=[
        html.Div(
            className='row', 
            children=[
                html.Div(
                    id='numfeats-selected-counter', 
                    className='four columns', 
                    children='Genes selected: 0', 
                    style={
                        'textAlign': 'center', 
                        'color': app_config.params['font_color'], 
                        'padding-top': '10px'
                    }
                ), 
                html.Div(
                    className='eight columns', 
                    children=[
                        dcc.Input(
                            id='featset-name', 
                            type='text', 
                            placeholder="Store gene set with name...", 
                            value='', 
                            debounce=True, 
                            n_submit_timestamp='0', 
                            style={'width': '100%'}
                        )], 
                    style={'padding-top': '5px'}
                )]
        )]
)



# ======================================================
# =================== Component divs ===================
# ======================================================


def create_div_align_selection(options_list):
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='four columns', 
                children=[ 
                    html.Button(
                        id='align-button', 
                        children='Display alignment', 
                        n_clicks=0, 
                        n_clicks_timestamp=0,
                        style=style_text_box
                    )], 
                style={'padding-top': '5px'}
            ), 
            html.Div(
                className='eight columns', 
                children=[
                    dcc.RadioItems(
                        id='align-method-selection', 
                        options=[ {'label': v, 'value': v} for v in options_list ], 
                        style=legend_font_macro, 
                        labelStyle={
                            'display': 'inline-block', 
                            'margin-right': '5px'
                        }, 
                        value='Unaligned'
                    )], 
                style={'padding-top': '10px'}
            )]
    )


# Default dataset first in the given list of dataset options.
def create_div_select_dataset(dataset_options):
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='four columns', 
                children=[
                    html.P(
                        "Browse dataset: ", 
                        style=style_text_box
                    )], 
                style={'padding-top': '10px'}
            ), 
            html.Div(
                className='eight columns', 
                children=[
                    dcc.Dropdown(
                        id='sourcedata-select', 
                        options = [ {'value': dn, 'label': dn} for dn in dataset_options ], # style={'height': '30px'}, 
                        value=dataset_options[0], 
                        clearable=False
                    )]
            )], 
        style=style_outer_dialog_box
    )


def create_div_go_ctrl(gene_names):
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='row', 
                children=[
                    dcc.Dropdown(
                        id='list-featsets', 
                        placeholder="Select stored gene(s) to load", 
                        multi=True, 
                        value=[]
                    )], 
                style={'padding-top': '5px'}
            ), 
            html.Div(
                className='row', 
                children=[
                    html.Div(
                        className='six columns', 
                        children=[
                            html.P(
                                "Gene set enrichment",
                                style={ 'textAlign': 'center', 'width': '100%', 'color': app_config.params['font_color'] }
                            )], 
                        style={'padding-top': '10px'}
                    ), 
                    html.Div(
                        className='three columns',
                        children=[
                            html.P(
                                "# terms to display: ",
                                style={ 'textAlign': 'right', 'width': '100%', 'color': app_config.params['font_color'] }
                            )], 
                        style={'padding-top': '7px'}
                    ),
                    html.Div(
                        className='three columns', 
                        children=[
                            dcc.Input(
                                id='select-topk-goterms', 
                                type='text', 
                                value='20', 
                                style={'textAlign': 'center', 'width': '100%'}
                            )]
                    )]
            )], 
        style=style_outer_dialog_box
    )


div_reviz_scatter = html.Div(
    className='row', 
    children=[
        html.Div(
            className='four columns', 
            children=[ 
                dcc.Checklist(
                    id='reviz-status', 
                    options=[
                        {'label': 'Visualize selection', 'value': 'viz'}
                    ],
                    values=[], 
                    style={
                        'textAlign': 'center', 
                        # 'width': '80%', 
                        'color': app_config.params['font_color']
                    }
                )], 
            style={'padding-top': '10px'}
        ), 
        html.Div(
            className='eight columns', 
            children=[
                dcc.RadioItems(
                    id='reviz-method-selection', 
                    options=[ {'label': v, 'value': v} for v in ['dummy', 'UMAP'] ], 
                    style=legend_font_macro, 
                    labelStyle={
                        'display': 'inline-block', 
                        'margin-right': '5px'
                    }, 
                    value='dummy'
                )], 
            style={'padding-top': '10px'}
        )]
)


def create_div_hm_panel():
    return html.Div(
        className="row",
        children = [
            html.Div(
                id='hm-feat-control', 
                children=[]
            ), 
            dcc.Graph(
                id='main-heatmap', 
                style={ 'height': '55vh' }, 
                config={'displaylogo': False, 'displayModeBar': True}
            )]
    )



# ================================================================
# =================== Component divs: cosmetic ===================
# ================================================================


# Default dataset first in the given list of dataset options.
def create_div_cosmetic_panel():
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='two columns', 
                children=[
                    html.P(
                        "Interface options: ", 
                        style=style_text_box
                    )], 
                style={'padding-top': '10px'}
            ), 
            html.Div(
                className='three columns', 
                children=[
                    dcc.Slider(
                        id='slider-marker-size-factor',
                        min=0, max=8, step=0.2, 
                        value=app_config.params['marker_size']
                    ), 
                    html.Div(
                        id='display-marker-size-factor', 
                        style={
                            'textAlign': 'center', 
                            'color': app_config.params['font_color'], 
                            'padding-top': '10px'
                        }
                    )]
            ), 
            html.Div(
                className='three columns', 
                children=[
                    dcc.Slider(
                        id='slider-selected-marker-size-factor', # marks={ 4: {'label': '4'} }
                        min=0, max=15, step=0.2, 
                        value=app_config.params['sel_marker_size']
                    ), 
                    html.Div(
                        id='display-selected-marker-size-factor', 
                        style={
                            'textAlign': 'center', 
                            'color': app_config.params['font_color'], 
                            'padding-top': '10px'
                        }
                    )]
            ), 
            html.Div(
                className='four columns', 
                children=[
                    dcc.Checklist(
                        id='toggle-hm-feat-panels', 
                        options=[
                            {'label': 'Cell line selection', 'value': 'bars'}, 
                            {'label': 'Per-cluster', 'value': 'percluster'}, 
                            {'label': 'Dendrogram', 'value': 'dendrogram'}, 
                            {'label': 'Heatmap', 'value': 'heatmap'}
                        ],
                        values=[], 
                        style={
                            'textAlign': 'left', 
                            'width': '80%', 
                            'color': app_config.params['font_color']
                        }, 
                        labelStyle={
                            'display': 'inline-block', 
                            'margin-right': '5px'
                        }
                    )], 
                style={'padding-top': '0px'}
            ), 
            html.Div(
                className='nine columns', 
                children=[
                    dcc.Checklist(
                        id='toggle-future-panels', 
                        options=[
                            {'label': 'Diff. feature sets', 'value': 'diff_features'}
                        ],
                        values=['diff_features'], 
                        style={
                            'textAlign': 'left', 
                            'width': '80%', 
                            'color': app_config.params['font_color']
                        }, 
                        labelStyle={
                            'display': 'inline-block', 
                            'margin-right': '5px'
                        }
                    ), 
                    dcc.Checklist(
                        id='toggle-debug-panels', 
                        options=[
                            {'label': 'Debug panel', 'value': 'debug-panel'}
                        ],
                        values=[], 
                        style={
                            'textAlign': 'left', 
                            'width': '80%', 
                            'color': app_config.params['font_color']
                        }, 
                        labelStyle={
                            'display': 'inline-block', 
                            'margin-right': '5px'
                        }
                    )], 
                style={'padding-top': '0px'}
            ), 
            html.Div(
                className='two columns', 
                children=[
                    html.A(
                        html.Button(
                            id='download-layout-button', 
                            children='Get CSV', 
                            style=style_text_box, 
                            n_clicks='0', 
                            n_clicks_timestamp='0'
                        ), 
                        id='download-layout-link',
                        download="selected_layout.csv", 
                        href="",
                        target="_blank", 
                        style={
                            'width': '100%', 
                            'textAlign': 'center', 
                            'color': app_config.params['font_color']
                        }
                    )], 
                style={ 'padding-top': '10px' }
            )
        ], 
        style=style_outer_dialog_box
    )



# ==================================================================
# =================== Aggregating component divs ===================
# ==================================================================


def create_div_mainctrl(point_names, feat_names, more_colorvars):
#     download_image = app_config.params['download_img_path'],
#     encoded_image = base64.b64encode(open(download_image, 'rb').read())
    return html.Div(
        className='row', 
        children=[
            html.Div(
                className='three columns', 
                children=[
                    dcc.Dropdown(
                        id='points_annot', 
                        options = [ {'value': gn, 'label': gn} for gn in point_names ], 
                        placeholder="Cell(s)...", multi=True
                    )], 
                style={ 'padding': 2, 'margin': 2}
            ), 
            html.Div(
                className='four columns', 
                children=[
                    dcc.Dropdown(
                        id='landscape-color', 
                        options = [{
                            'value': app_config.params['default_color_var'], 
                            'label': app_config.params['default_color_var']
                        }] + [
                            {'value': n, 'label': n} for n in more_colorvars
                        ] + [
                            {'value': gn, 'label': gn} for gn in feat_names 
                        ], 
                        value=app_config.params['default_color_var'], 
                        placeholder="Select colors to plot", 
                        clearable=False
                    )], 
                style={ 'padding': 2, 'margin': 2}
            ), 
            html.Div(
                className='five columns', 
                children=[
                    html.Div(
                        className='six columns', 
                        children=[
                            dcc.RadioItems(
                                id='main-heatmap-roworder', 
                                options=[ 
                                    {'label': 'Cocluster', 'value': 'Cocluster'}, 
                                    {'label': 'Sort by color', 'value': 'Sort by color'}], 
                                style=legend_font_macro, 
                                labelStyle={
                                    'margin-right': '5px'
                                }, 
                                value='Sort by color'
                            )]
                    ), 
                    html.Div(
                        className='six columns', 
                        children=[
                            dcc.RadioItems(
                                id='main-heatmap-normalize', 
                                options=[ 
                                    {'label': 'Log-transform', 'value': 'log'}, 
                                    {'label': 'Multinomial residual', 'value': 'mult_dev'}, 
                                    {'label': 'Binomial residual', 'value': 'bin_dev'}
                                ], 
                                style=legend_font_macro, 
                                labelStyle={ 'margin-right': '5px' }, 
                                value='log'
                            )]
                    )]
            )]
    )


def create_div_landscapes(feat_names):
    return html.Div(
        className="seven columns",
        children=[
            dcc.Graph(
                id='landscape-plot',
                config={'displaylogo': False, 'displayModeBar': True}, 
                style={ 'height': '100vh'}
            ), 
            html.Div(
                className="row", 
                children=[
                    html.Div(
                        className="four columns", 
                        children=[
                            div_pointset_select, 
                            html.Div(
                                className='row', 
                                children=[
                                    html.Div(
                                        className='row', 
                                        children='Cell set download/upload:', 
                                        style={
                                            'textAlign': 'center', 
                                            'color': app_config.params['font_color'], 
                                            'padding-top': '0px'
                                        }
                                    ), 
                                    html.Div(
                                        className='row', 
                                        children=[
                                            html.Div(
                                                className='six columns', 
                                                children=[
                                                    html.A(
                                                        html.Button(
                                                            id='download-button', 
                                                            children='Save', 
                                                            style=style_text_box, 
                                                            n_clicks='0', 
                                                            n_clicks_timestamp='0'
                                                        ), 
                                                        id='download-set-link',
                                                        download="selected_set.csv", 
                                                        href="",
                                                        target="_blank", 
                                                        style={
                                                            'width': '100%', 
                                                            'textAlign': 'center', 
                                                            'color': app_config.params['font_color']
                                                        }
                                                    )], 
                                                style={'padding-top': '0px'}
                                            ), 
                                            html.Div(
                                                className='six columns', 
                                                children=[
                                                    dcc.Upload(
                                                        id='upload-pointsets',
                                                        children=html.Div([
                                                            html.Button(
                                                                id='upload-button', 
                                                                children='Load', 
                                                                style=style_text_box, 
                                                                n_clicks='0', 
                                                                n_clicks_timestamp='0'
                                                            )]
                                                        ),
                                                        style={
                                                            'width': '100%', 
                                                            'textAlign': 'center', 
                                                            'color': app_config.params['font_color']
                                                        }, 
                                                        multiple=True
                                                    )]
                                            )]
                                    )], style={ 'border': 'thin lightgrey solid', 'padding': 1, 'margin': 1 }
                            ), 
                            create_div_select_dataset(app_config.params['dataset_options']), 
                            div_featset_select, 
                            html.Div(
                                className='row', 
                                children=[
                                    html.Div(
                                        className='row', 
                                        children='Gene set download/upload:', 
                                        style={
                                            'textAlign': 'center', 
                                            'color': app_config.params['font_color'], 
                                            'padding-top': '0px'
                                        }
                                    ), 
                                    html.Div(
                                        className='row', 
                                        children=[
                                            html.Div(
                                                className='six columns', 
                                                children=[
                                                    html.A(
                                                        html.Button(
                                                            id='feat-download-button', 
                                                            children='Save', 
                                                            style=style_text_box, 
                                                            n_clicks='0', 
                                                            n_clicks_timestamp='0'
                                                        ), 
                                                        id='feat-download-set-link',
                                                        download="selected_feats.csv", 
                                                        href="",
                                                        target="_blank", 
                                                        style={
                                                            'width': '100%', 
                                                            'textAlign': 'center', 
                                                            'color': app_config.params['font_color']
                                                        }
                                                    )], 
                                                style={'padding-top': '0px'}
                                            ), 
                                            html.Div(
                                                className='six columns', 
                                                children=[
                                                    dcc.Upload(
                                                        id='upload-featsets',
                                                        children=html.Div([
                                                            html.Button(
                                                                id='feat-upload-button', 
                                                                children='Load', 
                                                                style=style_text_box, 
                                                                n_clicks='0', 
                                                                n_clicks_timestamp='0'
                                                            )]
                                                        ),
                                                        style={
                                                            'width': '100%', 
                                                            'textAlign': 'center', 
                                                            'color': app_config.params['font_color']
                                                        }, 
                                                        multiple=True
                                                    )]
                                            )]
                                    )], style={ 'border': 'thin lightgrey solid', 'padding': 1, 'margin': 1 }
                            )
                        ], 
                        style=style_invis_dialog_box
                    ), 
                    html.Div(
                        className="seven columns", 
                        id='hm-future-panels', 
                        children=[
                            html.Div(
                                className='row', 
                                children=[
                                    create_div_go_ctrl(feat_names), 
                                    html.Div(
                                        className="row",
                                        children = [
                                            dcc.Graph(
                                                id='goenrich-panel', 
                                                style={ 'height': '30vh' }, 
                                                config={'displaylogo': False, 'displayModeBar': True}
                                            )]
                                    )]
                            )]#, style=style_outer_dialog_box
                    )]
            )]
    )


def create_div_sidepanels(point_names, more_colorvars, align_options_list):
    return html.Div(
        className='five columns', 
        children=[
            create_div_hm_panel(), 
            html.Div(
                className="row", 
                children=[
                    html.Div(
                        className='row', 
                        children=[
                            html.Div(
                                className='four columns', 
                                children=[
                                    dcc.Dropdown(
                                        id='diff-foreset-select', 
                                        value=[], multi=True, 
                                        placeholder="Foreground..."
                                    )], 
                                style={'padding-top': '0px'}
                            ), 
                            html.Div(
                                className='four columns', 
                                children=[
                                    dcc.Dropdown(
                                        id='diff-backset-select', 
                                        value=[], multi=True, 
                                        placeholder="Background..."
                                    )], 
                                style={'padding-top': '0px'}
                            ), 
                            html.Div(
                                className='two columns',
                                children=[
                                    html.P(
                                        "# displayed: ",
                                        style={ 'textAlign': 'right', 'width': '100%', 'color': app_config.params['font_color'] }
                                    )], 
                                style={'padding-top': '7px'}
                            ),
                            html.Div(
                                className='two columns', 
                                children=[
                                    dcc.Input(
                                        id='select-numgenes-disc', 
                                        type='text', 
                                        value='20', 
                                        style={'textAlign': 'center', 'width': '100%'}
                                    )]
                            )]
                    ), 
                    html.Div(
                        className='row', 
                        children=[
                            html.Div(
                                className='three columns', 
                                children=[
                                    html.P(
                                        "Select marker genes (logistic regression)",
                                        style={ 'textAlign': 'center', 'padding-top': '5px', 'padding-left': '5px', 'width': '100%', 'color': app_config.params['font_color'] }
                                    )]
                            ), 
                            html.Div(
                                className='nine columns', 
                                children=[
                                    dcc.Textarea(
                                        id='display-genelist', 
                                        wrap='True', value = '', 
                                        rows=1, placeholder="Selected genes", 
                                        style={'width': '100%'}
                                    )]
                            )], 
                        style={'padding-top': '5px'}
                    ), 
                    html.Div(
                        className='row', 
                        children=[
                            html.Div(
                                className='six columns', 
                                children=[
                                    dcc.Graph(
                                        id='diff-signal-panel', 
                                        config={ 'displaylogo': False, 'displayModeBar': True }
                                    )]
                            ), 
                            html.Div(
                                className='six columns', 
                                children=[
                                    dcc.Graph(
                                        id='diff-heatmap-panel', 
                                        config={ 'displaylogo': False, 'displayModeBar': False }
                                    )]
                            )]
                    ), 
                    html.Div(
                        className='row', 
                        children=[
                            html.Div(
                                className='four columns', 
                                children=[
                                    html.P(
                                        "Importance: ",
                                        style={ 'textAlign': 'center', 'width': '100%', 'color': app_config.params['font_color'] }
                                    )]
                            ), 
                            html.Div(
                                className='eight columns', 
                                children=[
                                    dcc.RadioItems(
                                        id='select-impts-display', 
                                        options=[ 
                                            {'label': 'Mean weight (\u03BC)', 'value': 'Mean'}, 
                                            {'label': 'Dispersion (\u03BC / \u03C3)', 'value': 'Dispersion'}, 
                                            {'label': 'Combinatorial potential', 'value': 'Combinatorial'}
                                        ], 
                                        style=legend_font_macro, 
                                        labelStyle={
                                            'display': 'inline-block', 
                                            'margin-right': '5px'
                                        }, 
                                        value='Mean'
                                    )], 
                                style={'padding-top': '0px'}
                            )], 
                        style={'padding-top': '10px'}
                    ), 
                    html.Div(
                        className='row', 
                        children=[
                            html.Div(
                                className='four columns', 
                                children=[
                                    html.P(
                                        "Color by: ",
                                        style={ 'textAlign': 'center', 'width': '100%', 'color': app_config.params['font_color'] }
                                    )]
                            ), 
                            html.Div(
                                className='four columns', 
                                children=[
                                    dcc.RadioItems(
                                        id='select-impts-color', 
                                        options=[ 
                                            {'label': 'Log fold change', 'value': 'Logfc'}, 
                                            {'label': 'Importance', 'value': 'Importance'}
                                        ], 
                                        style=legend_font_macro, 
                                        labelStyle={
                                            'display': 'inline-block', 
                                            'margin-right': '5px'
                                        }, 
                                        value='Logfc'
                                    )], 
                                style={'padding-top': '0px'}
                            ), 
                            html.Div(
                                className='four columns', 
                                children=[
                                    html.P(
                                        id='disc-err', 
                                        children="AUC: ", 
                                        style={ 'textAlign': 'center', 'width': '100%', 'color': app_config.params['font_color'] }
                                    )], 
                                style=style_outer_dialog_box
                            )], 
                        style={'padding-bottom': '0px'}
                    )], 
                style=style_outer_dialog_box
            )
            # create_div_align_selection(align_options_list), 
            # div_reviz_scatter, 
        ],
        style=style_invis_dialog_box
    )


"""
Main layout.
"""

def create_div_mainapp(point_names, feat_names, more_colorvars=[], align_options_list=['Unaligned', 'Aligned']):
    return html.Div(
        className="container", 
        children=[
            html.Div(
                className='row', 
                children=[
                    html.H1(
                        id='title', 
                        children=app_config.params['title'], 
                        style=style_text_box
                    )]
            ), 
            create_div_mainctrl(point_names, feat_names, more_colorvars), 
            html.Div(
                className="row", 
                children=[
                    create_div_landscapes(feat_names), 
                    create_div_sidepanels(point_names, more_colorvars, align_options_list)
                ]
            ), 
            create_div_cosmetic_panel(), 
            html.Div([ html.Pre(id='test-select-data', style={ 'color': app_config.params['font_color'], 'overflowX': 'scroll' } ) ]),     # For testing purposes only!
            html.Div(
                className='row', 
                children=[ 
                    dcc.Markdown(
                        """Queries? Requests? Contact [Akshay Balsubramani](abalsubr@stanford.edu). """ 
                        + """Source [repository](https://github.com/kundajelab/hodos)."""
                        )], 
                style={
                    'textAlign': 'center', 
                    'color': app_config.params['font_color'], 
                    'padding-bottom': '10px'
                }
            ), 
            dcc.Store(
                id='stored-pointsets', 
                data={ '_current_selected_data': {} }    # Maintained as the short-term state of a point subset.
            ), 
            dcc.Store(
                id='stored-featsets', 
                data={ '_current_selected_feats': {} }    # Maintained as the short-term state of a point subset.
            ), 
            dcc.Store(
                id='stored-landscape-selected', 
                data={ }, 
                modified_timestamp='0'
            ), 
            dcc.Store(
                id='stored-heatmap-selected', 
                data={ }, 
                modified_timestamp='0'
            ), 
            dcc.Store(
                id='stored-most-recently-highlighted', 
                data={ '_last_panel_highlighted': 'landscape' }, 
                modified_timestamp='0'
            ), 
            dcc.Store(
                id='stored-panel-settings', 
                data={
                    'debug_panel': False
                }, 
                modified_timestamp='0'
            ), 
            dcc.Store(
                id='stored-selected-feats', 
                data={ }, 
                modified_timestamp='0'
            ), 
            dcc.Store(
                id='last-stored-subsetname', 
                data='', 
                modified_timestamp='0'
            ), 
            dcc.Store(
                id='last-stored-featsname', 
                data='', 
                modified_timestamp='0'
            ), 
            dcc.Store(
                id='last-loaded', 
                data=[], 
                modified_timestamp='0'
            ), 
            dcc.Store(
                id='last-loaded-feats', 
                data=[], 
                modified_timestamp='0'
            ), 
            dcc.Store(
                id='stored-classifier-info', 
                data={ 'losses': [], 'importances': [], 'log_fc': [], 'feat_names': [] }
            )
        ],
        style={ 
            'width': '100vw', 
            'max-width': 'none'
        }
    )
