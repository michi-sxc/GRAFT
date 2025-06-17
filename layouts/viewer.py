import dash_bootstrap_components as dbc
from dash import dcc, html
from dash_iconify import DashIconify

from config import colors

def create_layout():
    return dbc.Container([
        dbc.Row([
            dbc.Col([
                dbc.Card([
                     dbc.CardHeader(html.H4("File Viewer", className="card-title mb-0")), # Use CardHeader
                     dbc.CardBody([
                        # --- File Selection ---
                        dbc.Label("Select File to View:", html_for="viewer-file-selector", className="fw-bold"),
                        dcc.Dropdown(
                            id='viewer-file-selector', options=[], value=None,
                            clearable=False, placeholder="Select a file...",
                            className="mb-3",
                        ),

                        # --- Filtering Accordion ---
                        dbc.Accordion([
                            dbc.AccordionItem([
                                dbc.Row([
                                    dbc.Col([
                                        dbc.Label("Mapping Quality (MAPQ):", className="mb-1 fw-bold"),
                                        dcc.RangeSlider(
                                            id='viewer-mapq-range', min=0, max=255, step=1, value=[0, 255],
                                            marks={i: str(i) for i in range(0, 256, 50)},
                                            tooltip={"placement": "bottom", "always_visible": False},
                                            updatemode='drag', className="mb-3 dbc",
                                        ),
                                    ], md=6),
                                     dbc.Col([
                                        dbc.Label("Min Base Quality:", className="mb-1 fw-bold"),
                                        dbc.Input(id='viewer-min-base-quality', type='number', value=20, min=0, max=93, step=1, className="mb-3", size="sm"),
                                    ], md=6),
                                ]),
                                dbc.Row([
                                     dbc.Col([
                                        dbc.Label("Mismatch (NM) Filter:", className="mb-1 fw-bold"),
                                        dcc.Dropdown(
                                            id='viewer-nm-dropdown',
                                            options=[{'label': f'NM <= {i}', 'value': i} for i in range(6)] + [{'label': 'All Reads', 'value': 'all'}],
                                            value='all', clearable=False, className="mb-3",
                                        ),
                                    ], md=6),
                                    dbc.Col([
                                        dbc.Label("C>T Change Count:", className="mb-1 fw-bold"),
                                        dcc.Dropdown(
                                            id='viewer-ct-count-dropdown',
                                            options=[{'label': 'Any', 'value': 'any'}] + [{'label': str(i), 'value': i} for i in range(1, 6)],
                                            value='any', clearable=False, className="mb-3",
                                        ),
                                    ], md=6),
                                ]),
                                dbc.Label("C>T Filter Logic:", className="mb-1 fw-bold"),
                                dcc.Checklist(
                                    id='viewer-ct-checklist',
                                    options=[ # Same options as main page
                                        {'label': ' Require C>T/G>A', 'value': 'only_ct'},
                                        {'label': ' Require ONLY C>T/G>A', 'value': 'exclusively_ct'},
                                        {'label': ' Exclude C>T/G>A', 'value': 'exclude_ct'},
                                        {'label': ' Subtract C>T/G>A from NM', 'value': 'subtract_ct'}
                                    ], value=[], className="mb-3", inline=True, labelStyle={'margin-right':'10px'}
                                ),
                                # --- Display Options ---
                                html.Hr(),
                                html.H6("Display Options", className="mt-3 mb-2"),
                                dbc.Label("Soft-clip Display:", className="mb-1 fw-bold"),
                                dbc.RadioItems(
                                    id='viewer-softclip-display',
                                    options=[
                                        {'label': ' Show all', 'value': 'show_all', 'title': 'Display all bases, including soft-clipped'},
                                        {'label': ' Highlight', 'value': 'highlight', 'title': 'Visually highlight soft-clipped bases'},
                                        {'label': ' Hide Clips', 'value': 'hide', 'title': 'Do not render soft-clipped bases'},
                                        {'label': ' Exclude Reads', 'value': 'exclude', 'title': 'Filter out reads containing any soft-clipping'}
                                    ], value='highlight', inline=True, className="mb-3"
                                ),
                                html.Label("Sequence Display Legend:", className="mb-1 fw-bold"),
                                html.Div([
                                    html.Span("Match", className="me-2 badge bg-secondary"),
                                    html.Span("C>T/G>A", className="me-2 badge bg-warning text-dark"),
                                    html.Span("Mismatch", className="me-2 badge bg-danger"),
                                    html.Span("Soft Clip", className="me-2 badge bg-light text-dark"),
                                ], className="small mb-3")

                            ], title="Filtering & Display Options", item_id="viewer-filter-item")
                        ], start_collapsed=True, flush=True, className="mb-4"),

                        # --- Navigation & Search ---
                        dbc.Card([
                             dbc.CardBody([
                                dbc.Row([
                                    # Search Inputs
                                    dbc.Col(dbc.InputGroup([
                                        dbc.Input(id='viewer-search-input', placeholder='Search read name...', size="sm"),
                                        dbc.Button(DashIconify(icon="mdi:magnify"), id='viewer-search-button', n_clicks=0, size="sm", color="info", outline=True)
                                    ]), lg=4, md=6, className="mb-2 mb-lg-0"),
                                    dbc.Col(dbc.InputGroup([
                                        dbc.Input(id='viewer-sequence-search-input', placeholder='Search sequence content...', size="sm"),
                                        dbc.Button(DashIconify(icon="mdi:magnify"), id='viewer-sequence-search-button', n_clicks=0, size="sm", color="info", outline=True)
                                    ]), lg=4, md=6, className="mb-2 mb-lg-0"),

                                    # Sorting
                                    dbc.Col(dbc.InputGroup([
                                        dbc.InputGroupText("Sort:", style={"fontSize": "0.8em"}),
                                        dbc.Select(id='viewer-sort-field', options=[
                                            {'label': 'Name', 'value': 'name'}, {'label': 'Length', 'value': 'length'}, {'label': 'Mismatches (NM)', 'value': 'mismatches'}, {'label': 'MAPQ', 'value': 'mapq'} # Added MAPQ
                                        ], value='name', size="sm"),
                                        dbc.Select(id='viewer-sort-order', options=[
                                            {'label': 'Asc', 'value': 'asc'}, {'label': 'Desc', 'value': 'desc'}
                                        ], value='asc', size="sm"),
                                    ]), lg=4, md=12), # Full width on medium
                                ], className="g-2 align-items-center mb-2"),

                                # Pagination
                                dbc.Row([
                                    dbc.Col(html.Div(id='viewer-page-info', className="text-muted small text-start my-auto"), width="auto"),
                                    dbc.Col(dbc.Pagination(id='viewer-pagination', max_value=1, fully_expanded=False, first_last=True, previous_next=True, size="sm"), width="auto", className="ms-auto")
                                ], className="g-2 align-items-center")
                             ], className="p-2") # Reduced padding
                        ], className="mb-3", color="light", inverse=True), # Light background card


                        # --- File Viewer Content Area ---
                        dcc.Loading(
                            id="loading-viewer", type="circle", color=colors['primary'],
                            children=[
                                html.Div(
                                    id='file-viewer-content',
                                    style={
                                        "maxHeight": "70vh", # Relative height
                                        "overflowY": "auto",
                                        "padding": "15px",
                                        "backgroundColor": colors['surface'],
                                        "border": f"1px solid {colors['grid']}",
                                        "borderRadius": "4px",
                                        "fontFamily": "'Courier New', Courier, monospace", # Monospace font
                                        "fontSize": "0.9em",
                                        "color": colors['on_surface']
                                    }
                                ),
                                # Store to hold the total number of pages
                                dcc.Store(id='viewer-total-pages', data=1)
                             ]
                        ),
                    ])
                ], className="mb-4")
            ], width=12), # Full width
        ], justify='center'),
    ], fluid=True, className="px-4 py-3", style={"backgroundColor": colors['background']})