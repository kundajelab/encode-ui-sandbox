import dash
import dash_core_components as dcc
import dash_html_components as html

from dash.dependencies import State, Input, Output


external_scripts = [
    {'src': "https://igv.org/web/release/2.0.1/dist/igv.min.js"},
    {'src': "https://igv.org/web/jb/release/1.0.0/dist/juicebox.min.js"},
    {'src': "https://aframe.io/releases/0.8.0/aframe.min.js"},
    {'src': "https://code.jquery.com/jquery-3.2.1.slim.min.js",
     'integrity': "sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN",
     'crossorigin': "anonymous"},
    {'src': "https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js",
     'integrity': "sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q",
     'crossorigin': "anonymous"},
    {'src': "https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js",
     'integrity': "sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl",
     'crossorigin': "anonymous"},
    {'src': "https://unpkg.com/epgg@latest/umd/epgg.js"}
]


def custom_index():
    return '''
    <!DOCTYPE html>
    <html>

        <head>
            <meta charset="UTF-8">
            <meta charset="utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
            <meta name="theme-color" content="#000000">

            <title>ENCODE Comparison Tool</title>
            <link ref="stylesheet" href="/static/my-custom-css-normalize.css">
            {dash_renderer_css_bundle}
            {dash_component_css_bundles}
            <link ref="stylesheet" href="/static/my-component-css-override.css">

        </head>

        <body>
            <div id="react-entry-point">
                <div class="_dash-loading">
                    Loading...
                </div>
            </div>
        </body>

        <footer>
            <script type="javascript" src="/static/my-custom-javascript-bundle.js">
            {dash_renderer_javascript_bundle}
            {dash_component_javascript_bundles}
            <script type="javascript" src="/static/my-custom-javascript-override.js">
        </footer>

    </html>
    '''

app = dash.Dash(__name__, external_scripts = external_scripts)
app.index = custom_index()

app.layout = html.Div([
    html.Meta(charSet = 'utf-8'),
    html.Meta(name = 'viewport', content = 'width=device-width, initial-scale=1, shrink-to-fit=no'),
    html.Meta(name = 'theme-color', content = '#000000'),
    html.Div([
        html.H2(
            children = 'ENCODE Comparison Tool',
            style = {'margin':'0',
                     'padding-top':'10px',
                     'font-family':'"Open Sans", "HelveticaNeue", "Helvetica Neue", Helvetica, Arial, sans-serif'
                    }),
        html.H5(
            children = 'Display the difference between Google and ENCODE searches',
            style = {'margin':'0',
                     'padding':'0',
                     'font-family':'"Open Sans", "HelveticaNeue", "Helvetica Neue", Helvetica, Arial, sans-serif'
                    }),
        dcc.Input(
            id = 'query',
            placeholder = 'What would you like to search?',
            type = 'text',
            value = '',
            style = {'width':'400px',
                     'margin':'15px'
                    }
        ),
        html.Button(id='submit-button', n_clicks=0, children='Submit', style={'padding': '0 30px'})
    ],
    style={'color':'white',
           'text-align':'center',
           'background-color':'#1e2c45',
           'margin':'0',
    }),
    html.Div([
        html.Iframe(
            id='google-frame',
            src='',
            style={'width':'48%',
                   'float':'left',
                   'height':'625px',
                   'margin-left':'25px',
                   'margin-top':'15px',
                   'frameBorder':'0'
                  }),
        html.Iframe(
            id='encode-frame',
            src='',
            style={'width':'48%',
                   'float':'left',
                   'height':'625px',
                   'margin-right':'23px',
                   'margin-top':'15px',
                   'frameBorder':'0'
                  })
    ],
    style={'display':'flex'}),
    html.Div([
        html.Div(
                children = 'test',
                id = 'embed',
                style = {'width':'1000px',
                         'margin':'25px',
                         'background-color':'white',
                         'padding':'10px'}
        ),
        html.Div([
            dcc.Input(
                id = 'address',
                placeholder = 'Enter a genomic address',
                type = 'text',
                value = '',
                style = {
                    'margin-top':'25px',
                    'margin-bottom':'25px',
                    'width':'100%'
                }
            ),
            html.Button(id = 'address-submit', n_clicks=0, children='Submit', style={'padding': '0 30px', 'text-align':'center'})
        ])
    ],
    style = {
        'display':'flex',
        'text-align':'center',
        'align-items':'center'
    }),
    html.Div([
        html.Button(id='feedback-button',
                    n_clicks=0,
                    children='Give us Feedback',
                    style={'text-align':'center', 'margin-bottom':'100px', 'padding': '0 30px'}),
        html.Iframe(
            id='feedback-form',
            src='https://docs.google.com/forms/d/e/1FAIpQLScE5B3bJM8-L3Gj8kvCEMcqkV60VN6jvWQ8FLQ_3HKu0u5x7w/viewform?embedded=true',
            width='600' ,
            height='471',
            style={'margin-top':'5px', 'frameBorder':'0'}
        )
    ],
    style={
        'margin-top':'5px',
        'display':'flex',
        'flex-direction':'column',
        'flex':'1 1 20em',
        'align-items':'center',
        'justify-content':'center',
    }),
])


@app.callback(Output('google-frame', 'src'),
             [Input('submit-button', 'n_clicks')],
             [State('query', 'value')])
def updateGoogle(n_clicks, query):
    queryList = query.split()
    if n_clicks > 0:
        queryList.append('encode')
    googleURL = 'https://www.google.com/search?igu=1&ei=&q='
    googleURL = googleURL + '+'.join(queryList)
    return googleURL

@app.callback(Output('encode-frame', 'src'),
             [Input('submit-button', 'n_clicks')],
             [State('query', 'value')])
def updateEncode(n_clicks, query):
    queryList = query.split()
    encodeURL = 'https://www.encodeproject.org/search/?searchTerm='
    encodeURL = encodeURL + '+'.join(queryList)
    return encodeURL

@app.callback(Output('feedback-form', 'style'),
             [Input('feedback-button', 'n_clicks')])
def showFeedback(n_clicks):
    if(n_clicks > 0):
        return {'display':'block'}
    else:
        return {'display':'none'}

@app.callback(Output('feedback-button', 'style'),
             [Input('feedback-button', 'n_clicks')])
def hideButton(n_clicks):
    if(n_clicks > 0):
        return {'display':'none'}


if __name__ == '__main__':
    app.run_server(debug=True)
