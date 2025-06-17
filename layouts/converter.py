import dash_bootstrap_components as dbc
from dash import dcc, html
from dash_iconify import DashIconify
from config import colors 

def create_layout():
    return dbc.Container([
        dbc.Row(dbc.Col(html.H2("File Conversion", className="mb-4 text-center"))),
        
        dbc.Row([
            # --- Left Column: Upload and Main Conversion Settings ---
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        dbc.Row([
                            dbc.Col(DashIconify(icon="mdi:file-document-multiple-outline", width=24), width="auto", className="me-2"),
                            dbc.Col(html.H5("1. Upload & Select Output", className="mb-0")),
                        ], align="center")
                    ),
                    dbc.CardBody([
                        dcc.Upload(
                            id='convert-upload-data',
                            children=html.Div([
                                DashIconify(icon="ph:upload-simple-duotone", width=50, color=colors.get('secondary')),
                                html.Div("Drag & Drop or Click to Select Files", className="mt-2"),
                                html.Small("Supports: BAM, SAM, FASTA, FASTQ", className="text-muted"),
                            ]),
                            multiple=True,
                            className="upload-box mb-3 p-4",
                        ),
                        html.Div(id='convert-upload-filenames', className='mt-2 mb-3 small text-muted text-center'),

                        dbc.Label("Select Output Format:", html_for='convert-output-format', className="fw-bold"),
                        dcc.Dropdown(
                            id='convert-output-format',
                            options=[
                                {'label': 'SAM (Sequence Alignment/Map)', 'value': 'sam'},
                                {'label': 'BAM (Binary Alignment/Map)', 'value': 'bam'},
                                {'label': 'FASTA (Sequence Only)', 'value': 'fasta'},
                                {'label': 'FASTQ (Sequence + Quality)', 'value': 'fastq'},
                            ],
                            value='sam',
                            clearable=False,
                            className="mb-3",
                        ),
                    ])
                ], style={"overflow": "visible"}),
            ], md=6, className="mb-3 mb-md-0"),


            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        dbc.Row([
                            dbc.Col(DashIconify(icon="mdi:tune-variant", width=24), width="auto", className="me-2"),
                            dbc.Col(html.H5("2. Format-Specific Options", className="mb-0")),
                        ], align="center")
                    ),
                    dbc.CardBody([
                        html.Div(id='converter-advanced-options-container', children=[
                            # Section for BAM/SAM -> FASTA/FASTQ (initially hidden)
                            html.Div(id='adv-opts-bam-to-seq-div', style={'display': 'none'}, children=[
                                html.H6("FASTA/FASTQ Output Options (from BAM/SAM):", className="mb-3 text-muted"),

                                dbc.Checklist(
                                    options=[
                                        {"label": "Exclude unmapped reads", "value": "exclude_unmapped"},
                                        {"label": "Exclude secondary alignments", "value": "exclude_secondary"},
                                        {"label": "Exclude supplementary alignments", "value": "exclude_supplementary"},
                                    ],
                                    value=["exclude_unmapped", "exclude_secondary", "exclude_supplementary"],
                                    id="converter-bam-filter-checklist",
                                    switch=True,
                                    className="mb-3",
                                    inputStyle={"marginLeft": "0"} 
                                ),
                                # -------------------------

                                dbc.Label("Minimum Mapping Quality (MAPQ) to include:", html_for="converter-min-mapq", className="fw-bold"),
                                dbc.Input(type="number", id="converter-min-mapq", value=0, min=0, max=255, step=1, className="mb-3", size="sm"),
                                dbc.Label("Phred offset for dummy qualities (if needed):", html_for="converter-dummy-quality-phred", className="fw-bold"),
                                dcc.Dropdown(
                                    id="converter-dummy-quality-phred",
                                    options=[{'label': f"Phred+{val} (ASCII: {chr(val+33)})", 'value': val} for val in [33, 64]],
                                    value=33, clearable=False, className="mb-3"
                                ),
                                html.Small("Note: For BAM->FASTQ, original qualities used if present. Dummy qualities if missing or for FASTA->FASTQ.", className="text-muted d-block mb-3")
                            ]),

                            html.Div(id='adv-opts-seq-to-bam-div', style={'display': 'none'}, children=[
                                html.H6("BAM/SAM Output Options (from FASTA/FASTQ):", className="mb-3 text-muted"),
                                dbc.Label("Header Source:", className="fw-bold"),
                                dcc.RadioItems(
                                    id='converter-header-source',
                                    options=[{'label': 'Minimal Unaligned Header', 'value': 'minimal_unaligned'}],
                                    value='minimal_unaligned', className="mb-3"
                                ),
                            ]),
                            html.Div(id='adv-opts-bam-sam-div', style={'display': 'none'}, children=[
                                html.H6("BAM <-> SAM Options:", className="mb-3 text-muted"),
                                dbc.Checklist(
                                    options=[{"label": "Include header in output", "value": "include_header"}],
                                    value=["include_header"],
                                    id="converter-header-include-checklist",
                                    switch=True, className="mb-3",
                                    inputStyle={"marginLeft": "0"}
                                ),
                                dbc.Label("BAM Compression Level (0-9):", className="fw-bold"),
                                dcc.Slider(id="converter-bam-compression", min=0, max=9, step=1, value=6, marks={i: str(i) for i in range(10)}, className="mb-3 dbc"),
                            ]),
                            html.Div(id='adv-opts-fasta-to-fastq-div', style={'display': 'none'}, children=[
                                html.H6("FASTA -> FASTQ Options:", className="mb-3 text-muted"),
                                dbc.Label("Assign Dummy Quality Score (Phred value):", className="fw-bold"),
                                dbc.Input(type="number", id="converter-dummy-quality-char-code", value=40, min=0, max=93, step=1, className="mb-3", size="sm"),
                                html.Small("This Phred score (e.g., 40) will be assigned to each base. Phred+33 offset is assumed for ASCII conversion.", className="text-muted d-block mb-3")
                            ]),
                            
                            dbc.Alert("Select an output format to see relevant advanced options.", color="info", id="converter-adv-options-initial-message")
                        ]),
                    ])

                ], style={"overflow": "visible"}), 

                dbc.Card([
                    dbc.CardHeader(
                        dbc.Row([
                            dbc.Col(DashIconify(icon="mdi:play-box-outline", width=24), width="auto", className="me-2"),
                            dbc.Col(html.H5("3. Execute Conversion", className="mb-0")),
                        ], align="center")
                    ),
                    dbc.CardBody([
                         dbc.Button(
                            [DashIconify(icon="mdi:cog-transfer-outline", className="me-1"), "Convert & Download Files"],
                            id="convert-button", color="primary", className="w-100 mb-3", disabled=True
                        ),
                        dcc.Download(id="convert-download-file"),
                        html.Div(id='convert-status-progress-wrapper', children=[
                            html.Div(id='convert-status', className='mt-2 small'),
                            # Future: Could add a dbc.Progress bar here if conversion is long
                        ], className="mt-2")
                    ])
                ], className="mt-3") # top margin
            ], md=6)
        ], className="mb-4"),

    ], fluid=True, className="px-4 py-5", style={"backgroundColor": colors.get('background')})
