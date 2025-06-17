import dash_bootstrap_components as dbc
from dash import dcc, html
from dash_iconify import DashIconify
from config import colors

# --- Navbar ---
def create_navbar():
    return dbc.Navbar(
        dbc.Container([
            dcc.Link(
                dbc.Row([
                    dbc.Col(
                        html.Div([
                            DashIconify(icon="fluent:dna-24-filled", width=36, 
                                      style={"color": colors['secondary']}),
                            html.Span("GRAFT", 
                                    className="ms-2 fw-bold", 
                                    style={"fontSize": "20px", "letterSpacing": "1px"})
                        ], className="d-flex align-items-center")
                    ),
                ], align="center", className="g-0"),
                href="/",
                style={"textDecoration": "none", "color": "white"},
            ),
            dbc.NavbarToggler(id="navbar-toggler", className="border-0"),
            dbc.Collapse(
                dbc.Nav([
                    dbc.NavItem(
                        dcc.Link(
                            html.Div([
                                DashIconify(icon="fluent:home-24-regular", width=18, className="me-1"),
                                "Home"
                            ], className="d-flex align-items-center"),
                            href="/", 
                            className="nav-link nav-item-custom"
                        )
                    ),
                    dbc.NavItem(
                        dcc.Link(
                            html.Div([
                                DashIconify(icon="fluent:arrow-sync-24-regular", width=18, className="me-1"),
                                "Convert"
                            ], className="d-flex align-items-center"),
                            href="/convert", 
                            className="nav-link nav-item-custom"
                        )
                    ),
                    dbc.NavItem(
                        dcc.Link(
                            html.Div([
                                DashIconify(icon="fluent:document-search-24-regular", width=18, className="me-1"),
                                "Viewer"
                            ], className="d-flex align-items-center"),
                            href="/file-viewer", 
                            className="nav-link nav-item-custom"
                        )
                    ),
                    dbc.NavItem(
                        dbc.Button(
                            [DashIconify(icon="fluent:settings-24-regular", width=18, className="me-1"), 
                             "Settings"],
                            id="settings-button", 
                            color="light", 
                            outline=True, 
                            size="sm",
                            className="ms-3"
                        )
                    ),
                ], className="ms-auto", navbar=True),
                id="navbar-collapse",
                navbar=True,
            ),
        ], fluid=True),
        color="dark",
        dark=True,
        className="navbar-custom shadow-sm",
        sticky="top",
        style={"backgroundColor": "#212529", "borderBottom": "2px solid #343a40"}
    )

# --- Merge Modal ---
def create_merge_modal():
    return dbc.Modal([
        dbc.ModalHeader(
            dbc.ModalTitle([
                DashIconify(icon="fluent:folder-link-24-filled", width=24, className="me-2"),
                "Merge Sequence Files"
            ], className="d-flex align-items-center"),
            close_button=True
        ),
        dbc.ModalBody([
            dbc.Alert(
                [DashIconify(icon="fluent:info-24-regular", width=20, className="me-2"),
                 "Select multiple files of the same type to merge them into a single file."],
                color="info",
                className="d-flex align-items-center",
                dismissable=True,
                id="merge-validation-alert",
                is_open=False
            ),
            html.Div([
                html.H6("Available Files", className="mb-3"),
                dbc.Card([
                    dbc.CardBody([
                        dcc.Checklist(
                            id='merge-file-checklist',
                            options=[],
                            value=[],
                            className="merge-checklist",
                            labelStyle={'display': 'block', 'padding': '8px', 
                                       'marginBottom': '4px', 'cursor': 'pointer',
                                       'borderRadius': '4px', 'transition': 'background 0.2s'}
                        ),
                    ], className="p-3", style={"maxHeight": "300px", "overflowY": "auto"})
                ], className="border-0 bg-light")
            ]),
            dbc.Progress(
                id="merge-progress",
                value=0,
                striped=True,
                animated=True,
                className="mt-4 mb-2",
                style={"height": "6px", "display": "none"}
            ),
            html.Div(id="merge-status", className="text-center small mt-2"),
        ]),
        dbc.ModalFooter([
            dbc.Button(
                [DashIconify(icon="fluent:checkmark-24-filled", width=18, className="me-1"), 
                 "Merge Selected"],
                id="execute-merge-button", 
                color="primary",
                disabled=True
            ),
            dbc.Button(
                "Cancel", 
                id="close-merge-modal", 
                color="secondary", 
                outline=True,
                className="ms-2"
            )
        ])
    ], id="merge-modal", size="lg", centered=True, backdrop="static", 
      className="modal-custom")

# --- Settings Offcanvas ---
def create_settings_offcanvas():
    return dbc.Offcanvas([
        dbc.Container([
            # Header with icon
            html.Div([
                DashIconify(icon="fluent:settings-48-filled", width=48, 
                          className="text-primary mb-3"),
                html.H4("Application Settings", className="fw-bold mb-1"),
                html.P("Configure export and performance options", 
                      className="text-muted small")
            ], className="text-center mb-4"),
            
            html.Hr(),

            # Performance Settings
            dbc.Card([
                dbc.CardHeader([
                    DashIconify(icon="fluent:rocket-24-filled", width=20, className="me-2"),
                    "Performance"
                ], className="d-flex align-items-center py-2"),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            dbc.Switch(
                                id="show-performance-stats",
                                label=" Show Performance Metrics",
                                value=False,
                                className="mb-2",
                                disabled=True,
                                style={'marginLeft': '2.6rem'}
                            ),
                            html.Small("(Coming soon)", className="text-muted")
                        ], width=12),
                    ]),
                ])
            ], className="mb-3 border-0"),

            # Export Settings
            dbc.Card([
                dbc.CardHeader([
                    DashIconify(icon="fluent:arrow-export-24-filled", width=20, className="me-2"),
                    "Export Settings"
                ], className="d-flex align-items-center py-2"),
                dbc.CardBody([
                    # Format Selection
                    html.Label("Export Format", className="small fw-semibold mb-2"),
                    dbc.RadioItems(
                        id="plot-format-select",
                        options=[
                            {"label": html.Span([
                                DashIconify(icon="fluent:image-24-regular", width=16, className="me-1"),
                                "PNG"
                            ], className="d-flex align-items-center"), "value": "png"},
                            {"label": html.Span([
                                DashIconify(icon="fluent:shapes-24-regular", width=16, className="me-1"),
                                "SVG"
                            ], className="d-flex align-items-center"), "value": "svg"},
                            {"label": html.Span([
                                DashIconify(icon="fluent:document-24-regular", width=16, className="me-1"),
                                "PDF"
                            ], className="d-flex align-items-center"), "value": "pdf"},
                            {"label": html.Span([
                                DashIconify(icon="fluent:code-24-regular", width=16, className="me-1"),
                                "HTML"
                            ], className="d-flex align-items-center"), "value": "html"}
                        ],
                        value="png",
                        inline=True,
                        className="mb-4 custom-radio"
                    ),
                    
                    # Dimensions
                    dbc.Row([
                        dbc.Col([
                            html.Label("Width (px)", className="small fw-semibold mb-1"),
                            dbc.Input(
                                id="plot-width-input", 
                                type="number", 
                                value=1200,
                                min=600, 
                                max=3000, 
                                step=100,
                                size="sm"
                            ),
                        ], width=6),
                        dbc.Col([
                            html.Label("Height (px)", className="small fw-semibold mb-1"),
                            dbc.Input(
                                id="plot-height-input", 
                                type="number", 
                                value=800,
                                min=400, 
                                max=2000, 
                                step=100,
                                size="sm"
                            ),
                        ], width=6),
                    ], className="mb-3"),
                    
                    # Scale Factor
                    html.Label("Scale Factor", className="small fw-semibold mb-1"),
                    dcc.Slider(
                        id="plot-scale-input",
                        min=1,
                        max=5,
                        step=0.5,
                        value=2,
                        marks={i: str(i) for i in range(1, 6)},
                        className="mb-4"
                    ),
                    html.Small("Higher scale = better quality, larger file size", 
                             className="text-muted d-block mb-3"),
                    
                    # Export Button
                    dbc.Button(
                        [DashIconify(icon="fluent:arrow-download-24-filled", 
                                   width=18, className="me-2"),
                         "Export All Plots"],
                        id="export-plots-button",
                        color="primary",
                        className="w-100"
                    ),
                ])
            ], className="mb-3 border-0"),

            # Action Buttons
            html.Div([
                dbc.Button(
                    [DashIconify(icon="fluent:checkmark-24-filled", width=18, className="me-1"),
                     "Apply"],
                    id="apply-settings-button",
                    color="success",
                    className="me-2",
                    disabled=True
                ),
                dbc.Button(
                    [DashIconify(icon="fluent:arrow-reset-24-regular", width=18, className="me-1"),
                     "Reset"],
                    id="reset-settings-button",
                    color="secondary",
                    outline=True
                ),
            ], className="d-flex justify-content-center mt-4")
        ], fluid=True, className="py-3")
    ],
    id="settings-offcanvas",
    placement="end",
    is_open=False,
    style={"width": "400px"},
    className="custom-offcanvas"
)