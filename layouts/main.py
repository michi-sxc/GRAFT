import dash_bootstrap_components as dbc
from dash import dcc, html
from dash_iconify import DashIconify
from config import colors

def create_layout():
    return dbc.Container([
        dbc.Row([
            # --- Left Sidebar (Filters & Actions) ---
            dbc.Col([
                # File Selection Card
                dbc.Card([
                    dbc.CardHeader([
                        DashIconify(icon="fluent:document-data-24-filled", width=20, className="me-2"),
                        html.Span("File Selection", className="fw-semibold")
                    ], className="d-flex align-items-center py-2"),
                    dbc.CardBody([
                        dcc.Dropdown(
                            id='file-selector',
                            options=[],
                            value=None,
                            clearable=False,
                            placeholder="Select or Upload a File",
                            className="mb-3",
                            style={"fontSize": "14px"}
                        ),
                        dbc.Button(
                            [DashIconify(icon="fluent:folder-link-24-regular", width=18, className="me-2"), 
                             "Merge Files"],
                            id="open-merge-modal-button",
                            color="primary",
                            outline=True,
                            size="sm",
                            className="w-100"
                        ),
                    ], className="py-3")
                ], className="mb-3 shadow-sm border-0"),

                # Filters Card
                dbc.Card([
                    dbc.CardHeader([
                        DashIconify(icon="fluent:filter-24-filled", width=20, className="me-2"),
                        html.Span("Analysis Filters", className="fw-semibold")
                    ], className="d-flex align-items-center py-2"),
                    dbc.CardBody([
                        # Sequencing Type Toggle
                        html.Div([
                            html.Label("Sequencing Type", className="small text-muted mb-2"),
                            dbc.ButtonGroup([
                                dbc.Button("Single", id="seq-single", size="sm", active=True),
                                dbc.Button("Double", id="seq-double", size="sm")
                            ], className="w-100 mb-3 sequencing-button-group"),
                            # Hidden radio items for compatibility
                            dcc.RadioItems(id='sequencing-type', options=[
                                {'label': 'Single-stranded', 'value': 'single'},
                                {'label': 'Double-stranded', 'value': 'double'}
                            ], value='single', style={'display': 'none'}),
                        ]),

                        html.Hr(className="my-3"),

                        # Quality Filters Section
                        html.Div([
                            html.H6("Quality Filters", className="text-muted small mb-3"),
                            
                            # MAPQ Range
                            html.Label("MAPQ Range", className="small fw-semibold mb-1"),
                            dcc.RangeSlider(
                                id='mapq-range',
                                min=0, max=255, step=1,
                                value=[0, 255],
                                marks={i: {'label': str(i), 'style': {'fontSize': '10px'}} for i in [0, 50, 100, 150, 200, 255]},
                                tooltip={"placement": "bottom", "always_visible": False},
                                className="mb-4"
                            ),

                            # Base Quality
                            html.Label("Min Base Quality", className="small fw-semibold mb-1"),
                            dbc.InputGroup([
                                dbc.Input(id='min-base-quality', type='number', value=20, min=0, max=93, step=1, size="sm"),
                                dbc.InputGroupText("Phred", className="small")
                            ], size="sm", className="mb-3"),

                            # # Deduplication Switch
                            # dbc.Row([
                            #     dbc.Col([
                            #         dbc.Switch(
                            #             id='apply-sequence-deduplication-switch',
                            #             label="Remove Duplicates",
                            #             value=False,
                            #             className="mt-1"
                            #         ),
                            #     ], width=12)
                            # ], className="mb-3"),
                            dbc.Switch(
                                id='apply-sequence-deduplication-switch',
                                label="Remove Duplicates",
                                value=False,
                                className="mt-1 mb-3",
                                style={'marginLeft': '2.6rem'} # Adjust '1rem' (approx 16px) as needed
                            ),                          
                        ]),

                        html.Hr(className="my-3"),

                        # Alignment Filters Section
                        html.Div([
                            html.H6("Alignment Filters", className="text-muted small mb-3"),
                            
                            # Soft Clipping
                            html.Label("Soft Clipping", className="small fw-semibold mb-1"),
                            dbc.Select(
                                id='soft-clip-handling',
                                options=[
                                    {'label': 'Show all bases', 'value': 'show'},
                                    {'label': 'Exclude clipped regions', 'value': 'exclude_regions'},
                                    {'label': 'Exclude reads with clips', 'value': 'exclude_all'},
                                ],
                                value='show',
                                size="sm",
                                className="mb-3"
                            ),

                            # NM Tag Filter
                            html.Label("Mismatch Filter (NM)", className="small fw-semibold mb-1"),
                            dbc.Select(
                                id='nm-dropdown',
                                options=[{'label': f'NM ≤ {i}', 'value': i} for i in range(6)] + 
                                       [{'label': 'All Reads', 'value': 'all'}],
                                value='all',
                                size="sm",
                                className="mb-3"
                            ),
                        ]),

                        html.Hr(className="my-3"),

                        # Damage Pattern Filters
                        html.Div([
                            html.H6("Damage Pattern Filters", className="text-muted small mb-3"),
                            
                            # C>T Options with better visual grouping
                            dbc.Card([
                                dbc.CardBody([
                                    dcc.Checklist(
                                        id='ct-checklist',
                                        options=[
                                            {
                                                'label': html.Span([
                                                    DashIconify(icon="mdi:filter-plus", width=16, className="me-2 text-primary"),
                                                    "Require C>T/G>A"
                                                ], className="d-flex align-items-center"),
                                                'value': 'only_ct'
                                            },
                                            {
                                                'label': html.Span([
                                                    DashIconify(icon="mdi:filter-check", width=16, className="me-2 text-info"),
                                                    "Only C>T/G>A"
                                                ], className="d-flex align-items-center"),
                                                'value': 'exclusively_ct'
                                            },
                                            {
                                                'label': html.Span([
                                                    DashIconify(icon="mdi:filter-minus", width=16, className="me-2 text-danger"),
                                                    "Exclude C>T/G>A"
                                                ], className="d-flex align-items-center"),
                                                'value': 'exclude_ct'
                                            },
                                            {
                                                'label': html.Span([
                                                    DashIconify(icon="mdi:minus-circle", width=16, className="me-2 text-warning"),
                                                    "Subtract from NM"
                                                ], className="d-flex align-items-center"),
                                                'value': 'subtract_ct'
                                            }
                                        ],
                                        value=[],
                                        className="damage-pattern-checklist",
                                        labelStyle={'display': 'block', 'marginBottom': '6px'}
                                    ),
                                ], className="py-2")
                            ], className="border-0 mb-3", style={"backgroundColor": "#3a3f44"}),

                            # C>T Count
                            html.Label("C>T/G>A Count", className="small fw-semibold mb-1"),
                            dbc.Select(
                                id='ct-count-dropdown',
                                options=[{'label': 'Any', 'value': 'any'}] + 
                                       [{'label': str(i), 'value': i} for i in range(1, 6)],
                                value='any',
                                size="sm"
                            ),
                        ]),
                    ], className="py-3")
                ], className="mb-3 shadow-sm border-0"),

                # Action Buttons Card
                dbc.Card([
                    dbc.CardBody([
                        dbc.Button(
                            [DashIconify(icon="fluent:arrow-export-24-filled", width=18, className="me-2"),
                             "Export Selected"],
                            id="export-button",
                            color="success",
                            className="w-100 mb-2"
                        ),
                        dbc.Button(
                            [DashIconify(icon="fluent:data-histogram-24-filled", width=18, className="me-2"),
                             "View Statistics"],
                            id="stats-button",
                            color="info",
                            className="w-100 mb-2"
                        ),
                        dbc.Button(
                            [DashIconify(icon="fluent:dismiss-square-24-regular", width=18, className="me-2"),
                             "Clear Selection"],
                            id="clear-selection-button",
                            color="secondary",
                            outline=True,
                            className="w-100"
                        ),
                    ], className="py-2")
                ], className="shadow-sm border-0"),

            ], lg=3, md=4, className="mb-4 sidebar-column"),

            # --- Main Content Area ---
            dbc.Col([
                # Upload Card
                dbc.Card([
                    dbc.CardBody([
                        dbc.Row([
                            dbc.Col([
                                html.H5("Upload Sequence Files", className="mb-3 fw-bold"),
                                html.P("Drag and drop your files or click to browse", 
                                      className="text-muted small mb-0")
                            ], width=8),
                            dbc.Col([
                                dbc.Badge("Supports: BAM, SAM, FASTA, FASTQ", 
                                         color="primary", 
                                         pill=True, 
                                         className="float-end mt-2")
                            ], width=4)
                        ]),
                        html.Hr(className="my-3"),
                        dcc.Upload(
                            id='upload-data',
                            children=html.Div([
                                html.Div([
                                    DashIconify(icon="fluent:cloud-arrow-up-48-filled", 
                                               width=48, 
                                               className="mb-3",
                                               style={"color": colors['primary']}),
                                    html.H6("Drop files here or click to browse", 
                                           className="mb-2"),
                                    html.P("Multiple files can be uploaded at once", 
                                          className="text-muted small mb-0")
                                ], className="text-center py-4")
                            ]),
                            multiple=True,
                            className="upload-area"
                        ),
                        html.Div(id="upload-status-messages", className="mt-3")
                    ])
                ], className="mb-4 shadow-sm border-0"),

                # Analysis Results Card
                dbc.Card([
                    dbc.CardHeader([
                        html.Div([
                            html.H5("Analysis Results", className="mb-0 fw-bold"),
                            html.Div(id="main-loader-status", className="small text-muted")
                        ]),
                        dbc.Progress(
                            id="main-progress-bar",
                            value=0,
                            striped=True,
                            animated=True,
                            style={'display': 'none', 'height': '4px'},
                            className="mt-2"
                        )
                    ], className="py-3"),
                    dbc.CardBody([
                        dcc.Loading(
                            id="loading-graphs",
                            type="circle",
                            color=colors['primary'],
                            children=dbc.Tabs([
                                dbc.Tab(
                                    dcc.Graph(id='read-length-histogram', config={'displaylogo': False}),
                                    label="Length Distribution",
                                    tab_id="tab-length",
                                    tab_class_name="custom-tab",
                                    active_tab_class_name="custom-tab-active"
                                ),
                                dbc.Tab(
                                    dcc.Graph(id='overall-cg-histogram', config={'displaylogo': False}),
                                    label="Overall CG",
                                    tab_id="tab-cg",
                                    tab_class_name="custom-tab",
                                    active_tab_class_name="custom-tab-active"
                                ),
                                dbc.Tab(
                                    dcc.Graph(id='cg-content-histogram', config={'displaylogo': False}),
                                    label="CG vs Length",
                                    tab_id="tab-cg-length",
                                    tab_class_name="custom-tab",
                                    active_tab_class_name="custom-tab-active"
                                ),
                                dbc.Tab(
                                    dcc.Graph(id='damage-patterns', config={'displaylogo': False}),
                                    label="End Patterns",
                                    tab_id="tab-damage",
                                    tab_class_name="custom-tab",
                                    active_tab_class_name="custom-tab-active"
                                ),
                                dbc.Tab(
                                    dcc.Graph(id='mismatch-frequency-plot', config={'displaylogo': False}),
                                    label="Damage Freq",
                                    tab_id="tab-freq",
                                    tab_class_name="custom-tab",
                                    active_tab_class_name="custom-tab-active"
                                ),
                                dbc.Tab(
                                    dcc.Graph(id='mismatch-type-bar-chart', config={'displaylogo': False}),
                                    label="Mismatch Types",
                                    tab_id="tab-types",
                                    tab_class_name="custom-tab",
                                    active_tab_class_name="custom-tab-active"
                                ),
                                dbc.Tab(
                                    dcc.Graph(id='damage-pattern-plot', config={'displaylogo': False}),
                                    label="Damage Heatmap",
                                    tab_id="tab-heatmap",
                                    tab_class_name="custom-tab",
                                    active_tab_class_name="custom-tab-active"
                                ),
                                dbc.Tab(
                                    dcc.Graph(id='mapq-histogram', config={'displaylogo': False}),
                                    label="MAPQ",
                                    tab_id="tab-mapq",
                                    tab_class_name="custom-tab",
                                    active_tab_class_name="custom-tab-active"
                                ),
                                dbc.Tab(
                                    dcc.Loading(
                                        html.Pre(id='data-summary', 
                                               className="p-3 bg-dark text-light",
                                               style={'fontSize': '12px', 'maxHeight': '600px', 'overflowY': 'auto'})
                                    ),
                                    label="Summary",
                                    tab_id="tab-summary",
                                    tab_class_name="custom-tab",
                                    active_tab_class_name="custom-tab-active"
                                ),
                            ], id="analysis-tabs", active_tab="tab-length", className="custom-tabs mt-0")
                        ),
                    ], className="pt-0")
                ], className="shadow-sm border-0"),
            ], lg=9, md=8),
        ]),
    ], fluid=True, className="px-3 py-4")

# --- Redesigned Alignment Stats Offcanvas ---
alignment_stats_offcanvas = dbc.Offcanvas([
    # Add the alignment-stats-content div that the callback expects
    html.Div(id="alignment-stats-content", style={'display': 'none'}),
    
    dbc.Container([
        # Summary Cards
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            DashIconify(icon="fluent:data-area-24-filled", width=32, className="text-primary mb-2 stat-icon"),
                            html.H3(id="total-reads-stat", className="mb-0 fw-bold"),
                            html.P("Total Reads", className="text-muted small mb-0")
                        ], className="text-center")
                    ])
                ], className="border-0 shadow-sm h-100 stats-card")
            ], lg=6, className="mb-3"),
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            DashIconify(icon="fluent:location-24-filled", width=32, className="text-success mb-2 stat-icon"),
                            html.H3(id="mapped-reads-stat", className="mb-0 fw-bold"),
                            html.P("Mapped", className="text-muted small mb-0")
                        ], className="text-center")
                    ])
                ], className="border-0 shadow-sm h-100 stats-card")
            ], lg=6, className="mb-3"),
        ]),
        
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            DashIconify(icon="fluent:copy-24-filled", width=32, className="text-warning mb-2 stat-icon"),
                            html.H3(id="duplicate-reads-stat", className="mb-0 fw-bold"),
                            html.P("Duplicates", className="text-muted small mb-0")
                        ], className="text-center")
                    ])
                ], className="border-0 shadow-sm h-100 stats-card")
            ], lg=6, className="mb-3"),
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.Div([
                            DashIconify(icon="fluent:arrow-swap-24-filled", width=32, className="text-info mb-2 stat-icon"),
                            html.H3(id="avg-mismatches-stat", className="mb-0 fw-bold"),
                            html.P("Avg Mismatches", className="text-muted small mb-0")
                        ], className="text-center")
                    ])
                ], className="border-0 shadow-sm h-100 stats-card")
            ], lg=6, className="mb-3"),
        ]),

        html.Hr(className="my-4"),

        # Detailed Stats
        html.H6("Mismatch Analysis", className="text-muted mb-3"),
        html.Div(id="mismatch-details-content"),
        
        html.Hr(className="my-4"),
        
        # Pie Chart
        html.H6("Mismatch Distribution", className="text-muted mb-3"),
        dcc.Loading(
            dcc.Graph(id="stats-pie-chart", config={'displaylogo': False, 'responsive': True}),
            type="circle",
            color=colors['primary']
        )
    ], fluid=True, className="py-3")
], id="stats-offcanvas", title="Alignment Statistics", placement="end", is_open=False, 
   style={"width": "600px"}, className="custom-offcanvas")