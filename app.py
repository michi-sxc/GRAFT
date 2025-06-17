import dash
import dash_bootstrap_components as dbc
from flask import Flask
from dash import DiskcacheManager 
import diskcache

from config import config_data

# --- Flask Server Initialization ---
server = Flask(__name__)
# Increased session timeout (for background tasks)
server.config['PERMANENT_SESSION_LIFETIME'] = config_data.get('session_lifetime', 3600)

cache = diskcache.Cache("./graft_cache") # Ensure this directory is writable
background_callback_manager = DiskcacheManager(cache)

# --- Dash App Initialization ---
app = dash.Dash(
    __name__,
    server=server,
    external_stylesheets=[
        dbc.themes.SOLAR, # base theme, modified in assets/custom.css
        "/assets/custom.css" 
    ],
    suppress_callback_exceptions=True, # Necessary because callbacks are defined in other files
    title="GRAFT",
    update_title='Updating...',
    background_callback_manager=background_callback_manager # Text shown in browser tab during callbacks

)

# --- Application Layout ---
# Layout is defined in index.py using components from layouts/


# --- Custom Index String (Minimal) ---
# Most styling is in assets/custom.css
app.index_string = '''
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
'''