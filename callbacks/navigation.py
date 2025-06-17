import time
from dash import Input, Output, html
import dash
from layouts import main, converter, viewer 
from app import app 

@app.callback(
    [Output('page-content', 'children'),
     Output('main-page-loaded-signal', 'data')],
    Input('url', 'pathname')
)
def display_page(pathname):
    """Renders the layout based on the URL pathname and signals main page load."""
    page_layout = None
    signal_data = dash.no_update # Default to no_update for the signal

    if pathname == '/convert':
        page_layout = converter.create_layout()
        # signal_data remains dash.no_update for non-main pages
    elif pathname == '/file-viewer':
        page_layout = viewer.create_layout()
        # signal_data remains dash.no_update
    elif pathname == '/':
        page_layout = main.create_layout()
        signal_data = time.time() # Set signal only when main page is loaded
    else:
        # Handle 404
        page_layout = html.Div([
            html.H1("404: Page Not Found"),
            html.P(f"The pathname '{pathname}' was not recognised.")
        ], className="text-center text-danger mt-5")
        # signal_data remains dash.no_update

    return page_layout, signal_data