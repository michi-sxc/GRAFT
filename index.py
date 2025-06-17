from dash import dcc, html

# --- 1. Initialize App ---
# This import sequence is important: config -> app -> layouts -> callbacks -> processing/plotting/utils
import config # Load config and constants first
from app import app, server # Import app and server instance after config is set

# layout modules
from layouts import main as main_layout # Alias to avoid conflict with callbacks.main
from layouts import converter as converter_layout
from layouts import viewer as viewer_layout
from layouts import shared as shared_layout
from layouts.main import alignment_stats_offcanvas # Import specific components if needed

# --- 2. Register Callbacks ---
from callbacks import navigation, main, viewer, converter, export, merging, settings


# --- 3. App Layout ---
# The main layout structure, referencing shared components and the page content container
app.layout = html.Div([
    dcc.Location(id='url', refresh=False), # URL routing component
    shared_layout.create_navbar(),             # Navbar from shared layout
    shared_layout.create_merge_modal(),        # Merge modal (initially hidden)
    shared_layout.create_settings_offcanvas(), # Settings offcanvas (initially hidden)
    alignment_stats_offcanvas,          # Stats offcanvas from main layout module

    # Container where page content is rendered by the navigation callback
    html.Div(id='page-content'),

    # --- Global Stores ---
    # Store for uploaded file data
    dcc.Store(id='main-page-loaded-signal', storage_type='memory'),
    dcc.Store(id='file-store', storage_type='memory', data={'files': []}), # memory storage should be sufficient
    # Store for processed data from main filters
    dcc.Store(id='processed-data-store', storage_type='memory'),
    # Store for alignment stats
    dcc.Store(id='alignment-stats-store', storage_type='memory'),
    # Store for mismatch frequencies
    dcc.Store(id='mismatch-frequencies-store', storage_type='memory'),
    # Store for selected lengths in the read length histogram
    dcc.Store(id='selected-lengths', storage_type='memory', data=[]),
    # Store for processed data in the viewer tab
    dcc.Store(id='viewer-processed-data-store', storage_type='memory'),
     # Store for merged file data before download
    dcc.Store(id='merged-file-data-store', storage_type='memory'),

    dcc.Store(id='persistent-settings-store', storage_type='local'), # filter values, selected file name
    dcc.Store(id='persistent-file-metadata-store', storage_type='local'), # file selector options

    
    # --- Download Components ---
    dcc.Download(id="download-selected-reads"),
    dcc.Download(id="download-plots-zip"),
    dcc.Download(id="convert-download-file"),
    dcc.Download(id="download-merged-file"),

], style={'backgroundColor': config.colors['background']}) # Apply background globally (avoids some styling issues)

# --- 4. Validation Layout ---
app.validation_layout = html.Div([
    app.layout, 
    main_layout.create_layout(),        
    converter_layout.create_layout(),  
    viewer_layout.create_layout(),     
    
])


# --- 5. Run Server ---
if __name__ == '__main__':
    app.run_server(debug=True, host='127.0.0.1', port=8050)