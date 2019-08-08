import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import State, Input, Output


app = dash.Dash(__name__)

app.layout = html.Div([
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
        html.Button(id='submit-button', n_clicks=0, children='Submit')
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
                   'margin-top':'15px'
                  })
    ])
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

if __name__ == '__main__':
    app.run_server(debug=True)
