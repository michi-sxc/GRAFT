import dash
from dash import html, dcc, callback, Input, Output, State
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import plotly.express as px
import pysam
import numpy as np
from pathlib import Path
import datetime
from Bio.Seq import Seq
from Bio import SeqIO
import random
from scipy import stats
import io
import base64
from scipy.stats import gamma, beta
from app import app
import warnings
warnings.filterwarnings('ignore')

class DamagePattern:
    def __init__(self):
        self.patterns = {
            'terminal_ct': {'rate': 0.45, 'length': 10, 'decay': 'exponential'},
            'internal_ct': {'rate': 0.02},
            'terminal_ga': {'rate': 0.35, 'length': 8, 'decay': 'exponential'},
            'internal_ga': {'rate': 0.02},
            'transitions': {'rate': 0.01},  # C->T, G->A, A->G, T->C
            'transversions': {'rate': 0.005},  # All other substitutions
            'insertions': {'rate': 0.002, 'max_length': 3},
            'deletions': {'rate': 0.003, 'max_length': 3},
            'homopolymer_errors': {'rate': 0.01, 'max_extension': 2}
        }
        self.damage_decay_patterns = ['exponential', 'linear', 'sigmoid']
        
    def apply_damage(self, sequence, strand='positive', temperature_factor=1.0):
        """Apply damage patterns to sequence with safe index handling"""
        damaged_seq = list(sequence)
        original_length = len(sequence)
        
        # Create a copy of patterns with temperature factor applied
        temp_patterns = {}
        for key, pattern in self.patterns.items():
            temp_patterns[key] = pattern.copy()
            if 'rate' in temp_patterns[key]:
                temp_patterns[key]['rate'] *= temperature_factor
        
        # Phase 1: Apply terminal damage (no length changes)
        pattern = 'terminal_ct' if strand == 'positive' else 'terminal_ga'
        base_rate = temp_patterns[pattern]['rate']
        
        positions = range(original_length)
        if strand == 'negative':
            positions = reversed(positions)
        
        for i in positions:
            current_rate = self._calculate_position_specific_rate(
                base_rate, 
                i if strand == 'positive' else original_length-1-i,
                original_length,
                temp_patterns[pattern]
            )
            
            base = damaged_seq[i]
            if ((strand == 'positive' and base == 'C') or 
                (strand == 'negative' and base == 'G')):
                if random.random() < current_rate:
                    damaged_seq[i] = 'T' if strand == 'positive' else 'A'
        
        # Phase 2: Apply substitutions (no length changes)
        for i in range(len(damaged_seq)):
            base = damaged_seq[i]
            
            # Transitions
            if random.random() < temp_patterns['transitions']['rate']:
                transitions = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
                if base in transitions:
                    damaged_seq[i] = transitions[base]
            
            # Transversions
            elif random.random() < temp_patterns['transversions']['rate']:
                transversions = {
                    'A': ['C', 'T'], 'T': ['G', 'A'],
                    'C': ['G', 'A'], 'G': ['C', 'T']
                }
                if base in transversions:
                    damaged_seq[i] = random.choice(transversions[base])
        
        # Phase 3: Apply length-changing mutations
        # First collect all modifications, then apply them safely
        modifications = []
        
        # Collect homopolymer modifications
        for i in range(1, len(damaged_seq)):
            if damaged_seq[i] == damaged_seq[i-1]:
                if random.random() < temp_patterns['homopolymer_errors']['rate']:
                    if random.random() < 0.5:  # Extension
                        extension = random.randint(1, temp_patterns['homopolymer_errors']['max_extension'])
                        modifications.append(('insert', i, [damaged_seq[i]] * extension))
                    else:  # Contraction
                        modifications.append(('delete', i, 1))
        
        # Collect insertions and deletions
        i = 0
        while i < len(damaged_seq):
            if random.random() < temp_patterns['insertions']['rate']:
                length = random.randint(1, temp_patterns['insertions']['max_length'])
                insertion = [random.choice('ACGT') for _ in range(length)]
                modifications.append(('insert', i, insertion))
            
            if random.random() < temp_patterns['deletions']['rate']:
                length = random.randint(1, temp_patterns['deletions']['max_length'])
                modifications.append(('delete', i, length))
            i += 1
        
        # Apply modifications in reverse order (from end to start)
        for mod_type, pos, mod_data in sorted(modifications, key=lambda x: x[1], reverse=True):
            if mod_type == 'insert':
                damaged_seq[pos:pos] = mod_data
            elif mod_type == 'delete':
                del damaged_seq[pos:pos + mod_data]
        
        return ''.join(damaged_seq)

    def _calculate_position_specific_rate(self, base_rate, position, total_length, pattern_params):
        """Calculate position-specific damage rate based on decay pattern"""
        if position >= pattern_params.get('length', 0):
            return self.patterns['internal_ct']['rate']
            
        decay_type = pattern_params.get('decay', 'exponential')
        length = pattern_params.get('length', 10)
        
        if decay_type == 'exponential':
            return base_rate * np.exp(-position / (length / 3))
        elif decay_type == 'linear':
            return base_rate * (1 - position / length)
        elif decay_type == 'sigmoid':
            midpoint = length / 2
            steepness = 1
            return base_rate / (1 + np.exp(steepness * (position - midpoint)))
        return base_rate

class QualityProfile:
    def __init__(self, mean_quality=30, terminal_decay=5, quality_shape='normal'):
        self.mean_quality = mean_quality
        self.terminal_decay = terminal_decay
        self.quality_shape = quality_shape
        self.quality_shapes = ['normal', 'gamma', 'beta']
        
    def generate_qualities(self, length):
        if self.quality_shape == 'normal':
            qualities = np.random.normal(self.mean_quality, 2, length)
        elif self.quality_shape == 'gamma':
            shape = self.mean_quality / 4
            scale = 4
            qualities = gamma.rvs(shape, scale=scale, size=length)
        elif self.quality_shape == 'beta':
            a = self.mean_quality / 2
            b = (40 - self.mean_quality) / 2
            qualities = beta.rvs(a, b, size=length) * 40
        
        # Apply terminal quality decay with smoothing
        window = self.terminal_decay
        for i in range(window):
            decay_factor = (window - i) / window
            smooth_factor = 0.5 * (1 - np.cos(np.pi * i / window))
            decay = decay_factor * smooth_factor
            
            qualities[i] *= (1 - decay * 0.3)
            qualities[-(i+1)] *= (1 - decay * 0.3)
        
        return np.clip(qualities, 0, 40).astype(int)


def generate_fragment_lengths(n_reads, params):
    """Generate fragment lengths with enhanced distribution control"""
    # Set default values if not provided
    default_params = {
        'mean_length': 150,
        'length_std': 0.4,
        'min_length': 35,
        'max_length': 500,
        'distribution': 'lognormal',
        'length_skew': 0,
        'length_kurtosis': 0,
        'bimodal_proportion': 0.7
    }
    
    # Update defaults with provided params, excluding None values
    params = {k: v for k, v in params.items() if v is not None}
    params = {**default_params, **params}
    
    distribution = params['distribution']
    
    if distribution == 'lognormal':
        mu = np.log(float(params['mean_length']))
        sigma = float(params['length_std'])
        lengths = np.random.lognormal(mu, sigma, n_reads)
    elif distribution == 'gamma':
        shape = (float(params['mean_length']) / float(params['length_std'])) ** 2
        scale = float(params['length_std']) ** 2 / float(params['mean_length'])
        lengths = gamma.rvs(shape, scale=scale, size=n_reads)
    elif distribution == 'bimodal':
        # Mix of two normal distributions
        prop = float(params['bimodal_proportion'])
        mean1 = float(params['mean_length'])
        mean2 = mean1 * 1.5
        std1 = float(params['length_std']) * mean1
        std2 = float(params['length_std']) * mean2
        
        mask = np.random.random(n_reads) < prop
        lengths = np.zeros(n_reads)
        lengths[mask] = np.random.normal(mean1, std1, mask.sum())
        lengths[~mask] = np.random.normal(mean2, std2, (~mask).sum())
    
    # Apply length modifiers
    if params['length_skew'] != 0:
        lengths = np.exp(np.log(np.maximum(lengths, 1e-10)) + float(params['length_skew']))
    
    if params['length_kurtosis'] != 0:
        lengths = np.power(lengths, 1 + float(params['length_kurtosis']))
    
    return np.clip(lengths, params['min_length'], params['max_length']).astype(int)

# Add at the top with other imports
colors = {
    'primary': '#6200EA',    # Deep Purple
    'secondary': '#03DAC6',  # Teal
    'background': '#343a40', # Dark background
    'surface': '#1E1E1E',    # Slightly lighter surface
    'stats_tab_color': '#002B36',
    'on_surface': '#FFFFFF', # White text on surface
    'muted': '#c4c4c4',      # Muted text color
    'highlight': '#FF4081',  # Highlight color
    'highlight2': '#FFC107', # Amber color for second highlight
    'plot_bg': '#1E1E1E',    # Plot background
    'grid': 'rgba(255, 255, 255, 0.2)', # Grid color
    'line': '#1ca8ff',       # Line color
    'marker': '#ff7f0e',     # Marker color
}

# Update app initialization
app = dash.Dash(__name__, 
                external_stylesheets=[dbc.themes.SOLAR, "/assets/custom.css"],
                suppress_callback_exceptions=True)



def create_bam_header(reference_genome):
    """Create BAM header with reference information"""
    ref_lengths = {
        'hg38': {'chr1': 248956422, 'chrM': 16569},
        'mm10': {'chr1': 195471971, 'chrM': 16299},
        'custom': {'chr1': 1000000}  # Default for custom reference
    }
    
    lengths = ref_lengths.get(reference_genome, ref_lengths['custom'])
    
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [{'LN': length, 'SN': chrom} for chrom, length in lengths.items()],
        'RG': [
            {
                'ID': 'synthetic1',
                'PL': 'ILLUMINA',
                'SM': f'synthetic_{reference_genome}',
                'LB': 'synthetic_lib',
                'DT': datetime.datetime.now().isoformat()
            }
        ],
        'PG': [
            {
                'ID': 'SyntheticADNA',
                'PN': 'SyntheticADNA',
                'VN': '1.0',
                'CL': 'Synthetic ancient DNA generator'
            }
        ]
    }
    return header

def load_reference_sequence(reference_genome, custom_fasta=None):
    """Load reference sequence from built-in or custom FASTA"""
    if reference_genome == 'custom' and custom_fasta:
        content_type, content_string = custom_fasta.split(',')
        decoded = base64.b64decode(content_string)
        fasta = io.StringIO(decoded.decode('utf-8'))
        return str(next(SeqIO.parse(fasta, 'fasta')).seq)
    
    # Generate a random DNA sequence
    else:
        random.seed(42)  # For reproducibility, if desired
        return ''.join(random.choices('ACGT', k=1000000))
    


def create_synthetic_read(reference_genome, length, damage_pattern, quality_profile, strand, ref_seq, temperature_factor=1.0):
    """
    Create a synthetic read with specified parameters
    
    Parameters:
    -----------
    reference_genome : str
        Reference genome identifier
    length : int
        Desired read length
    damage_pattern : DamagePattern
        Damage pattern object with mutation parameters
    quality_profile : QualityProfile
        Quality profile object for generating quality scores
    strand : str
        'positive' or 'negative' strand
    ref_seq : str
        Reference sequence to generate read from
    temperature_factor : float, optional
        Factor to modify damage rates (default: 1.0)
    """
    # Random start position
    start = random.randint(0, len(ref_seq) - length)

    # Extract sequence
    read_seq = ref_seq[start:start + length]

    # Apply damage pattern with temperature factor
    damaged_seq = damage_pattern.apply_damage(
        read_seq, 
        strand=strand, 
        temperature_factor=temperature_factor
    )

    # Generate quality scores based on the new length
    qualities = quality_profile.generate_qualities(len(damaged_seq))

    # Convert quality scores to ASCII
    qual = ''.join(chr(q + 33) for q in qualities)

    # Handle reverse strand
    if strand == 'negative':
        damaged_seq = str(Seq(damaged_seq).reverse_complement())
        qual = qual[::-1]

    return {
        'sequence': damaged_seq,
        'quality': qual,
        'position': start,
        'reference': 'chr1',
        'strand': strand
    }



def write_bam_file(reads, output_path, header, reference_sequence):
    """Write reads to BAM file with proper headers, flags, MD, and NM tags"""
    with pysam.AlignmentFile(output_path, 'wb', header=header) as outf:
        for i, read in enumerate(reads):
            a = pysam.AlignedSegment()
            a.query_name = f"synthetic_{i}"
            a.query_sequence = read['sequence']
            a.reference_id = 0  # chr1
            a.reference_start = read['position']
            a.query_qualities = pysam.qualitystring_to_array(read['quality'])
            
            # Set flags based on strand
            if read['strand'] == 'negative':
                a.flag = 16  # reverse strand
            else:
                a.flag = 0   # forward strand
            
            a.mapping_quality = 60
            a.cigar = [(0, len(read['sequence']))]  # Match all bases
            
            # Calculate and set the MD and NM tags
            ref_seq = reference_sequence[a.reference_start:a.reference_start + len(a.query_sequence)]
            
            # Handle reverse strand
            if a.is_reverse:
                ref_seq = str(Seq(ref_seq).reverse_complement())
            
            md_tag, nm_tag = calculate_md_and_nm_tags(a.query_sequence, ref_seq)
            a.set_tag("MD", md_tag)
            a.set_tag("NM", nm_tag)
            
            outf.write(a)


def calculate_md_and_nm_tags(query_seq, ref_seq):
    """
    Calculate the MD tag string and NM tag value by comparing the query sequence to the reference sequence.
    """
    md_string = ""
    match_count = 0
    nm_count = 0  # Number of mismatches
    for q_base, r_base in zip(query_seq, ref_seq):
        if q_base == r_base:
            match_count += 1
        else:
            nm_count += 1
            if match_count > 0:
                md_string += str(match_count)
                match_count = 0
            md_string += r_base
    if match_count > 0:
        md_string += str(match_count)
    return md_string, nm_count



layout = html.Div([
    dbc.Container([
        html.H1("Synthetic Ancient DNA Generator", className="text-center my-4"),
        
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H3("Basic Parameters", className="card-title"),
                        html.Hr(),
                        
                        dbc.Form([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Number of Reads"),
                                    dbc.Input(id='num-reads', type='number', value=1000, min=1, max=1000000),
                                ], width=12),
                            ], className="mb-3"),
                            
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Reference Genome"),
                                    dcc.Dropdown(
                                        id='reference-genome',
                                        options=[
                                            {'label': 'Human (hg38)', 'value': 'hg38'},
                                            {'label': 'Mouse (mm10)', 'value': 'mm10'},
                                            {'label': 'Custom FASTA', 'value': 'custom'}
                                        ],
                                        value='hg38',
                                        clearable=False,
                                        style={
                                            "background-color": colors['surface'],
                                            "color": colors['on_surface'],
                                            "border": f"1px solid {colors['muted']}",
                                            "width": "100%",
                                            "margin-bottom": "10px"
                                        }
                                    ),
                                ], width=12),
                            ], className="mb-3"),
                            
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Upload Custom FASTA"),
                                    dcc.Upload(
                                        id='custom-fasta',
                                        children=dbc.Button("Upload File", color="secondary", outline=True),
                                        style={'width': '100%'},
                                        multiple=False
                                    ),
                                    html.Div(id='custom-fasta-status', className="mt-2"),
                                ], width=12),
                            ], className="mb-3"),
                        ]),
                    ])
                ], className="mb-4"),
                
                dbc.Card([
                    dbc.CardBody([
                    html.H3("Advanced Parameter Controls", className="card-title"),
                    html.Hr(),
                        dbc.Row([
                            dbc.Col([
                                html.Label("Length Distribution Type"),
                                dcc.Dropdown(
                                    id='length-distribution-type',
                                    options=[
                                        {'label': 'Log-normal', 'value': 'lognormal'},
                                        {'label': 'Gamma', 'value': 'gamma'},
                                        {'label': 'Bimodal', 'value': 'bimodal'}
                                    ],
                                    value='lognormal',
                                    clearable=False,
                                    style={
                                        "background-color": colors['surface'],
                                        "color": colors['on_surface'],
                                        "border": f"1px solid {colors['muted']}",
                                        "width": "100%",
                                        "margin-bottom": "10px"
                                    }
                                )
                            ]),
                            dbc.Col([
                                html.Label("Length Distribution Skew"),
                                dcc.Slider(
                                    id='length-skew',
                                    min=-2, max=2, step=0.01,
                                    value=0,
                                    marks={i: str(i) for i in range(-2, 3)}
                                )
                            ]),
                            dbc.Col([
                                html.Label("Length Distribution Kurtosis"),
                                dcc.Slider(
                                    id='length-kurtosis',
                                    min=-0.5, max=0.5, step=0.01,
                                    value=0,
                                    marks={i*0.5: str(i) for i in range(-1, 2)}
                                )
                            ])
                        ]),
                        dbc.Row([
                            dbc.Col([
                                html.Label("Damage Decay Pattern"),
                                dcc.Dropdown(
                                    id='damage-decay-pattern',
                                    options=[
                                        {'label': 'Exponential', 'value': 'exponential'},
                                        {'label': 'Linear', 'value': 'linear'},
                                        {'label': 'Sigmoid', 'value': 'sigmoid'}
                                    ],
                                    value='exponential',
                                    clearable=False,
                                    style={
                                        "background-color": colors['surface'],
                                        "color": colors['on_surface'],
                                        "border": f"1px solid {colors['muted']}",
                                        "width": "100%",
                                        "margin-bottom": "10px"
                                    }                                    
                                )
                            ]),
                            dbc.Col([
                                html.Label("Quality Score Distribution"),
                                dcc.Dropdown(
                                    id='quality-distribution',
                                    options=[
                                        {'label': 'Normal', 'value': 'normal'},
                                        {'label': 'Gamma', 'value': 'gamma'},
                                        {'label': 'Beta', 'value': 'beta'}
                                    ],
                                    value='normal',
                                    clearable=False,
                                    style={
                                        "background-color": colors['surface'],
                                        "color": colors['on_surface'],
                                        "border": f"1px solid {colors['muted']}",
                                        "width": "100%",
                                        "margin-bottom": "10px"
                                    }                                      
                                )
                            ])
                        ]),
                        dbc.Row([
                            dbc.Col([
                                html.Label("Temperature Factor"),
                                dcc.Slider(
                                    id='temperature-factor',
                                    min=0.5, max=2, step=0.1,
                                    value=1,
                                    marks={i/2: str(i/2) for i in range(1, 5)}
                                )
                            ]),
                            dbc.Col([
                                html.Label("Homopolymer Error Rate"),
                                dcc.Slider(
                                    id='homopolymer-error-rate',
                                    min=0, max=5, step=0.1,
                                    value=1,
                                    marks={i: f"{i}%" for i in range(0, 6)}
                                )
                            ])
                        ])
                    ])
                ], className="mb-4"),

                dbc.Card([
                    dbc.CardBody([
                        html.H3("Fragment Length Distribution", className="card-title"),
                        html.Hr(),
                        
                        dbc.Form([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Mean Length"),
                                    dbc.Input(id='mean-length', type='number', value=150, min=20, max=1000),
                                ], width=6),
                                dbc.Col([
                                    dbc.Label("Standard Deviation"),
                                    dbc.Input(id='length-std', type='number', value=0.4, min=0.1, max=2, step=0.1),
                                ], width=6),
                            ], className="mb-3"),
                            
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Minimum Length"),
                                    dbc.Input(id='min-length', type='number', value=35, min=20, max=500),
                                ], width=6),
                                dbc.Col([
                                    dbc.Label("Maximum Length"),
                                    dbc.Input(id='max-length', type='number', value=500, min=50, max=1000),
                                ], width=6),
                            ], className="mb-3"),
                            
                            dcc.Graph(id='length-distribution-preview')
                        ]),
                    ])
                ]),
            ], md=4),
            
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H3("Damage Patterns", className="card-title"),
                        html.Hr(),
                        
                        dbc.Form([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Terminal C>T Damage"),
                                    dcc.Slider(id='terminal-ct-rate', min=0, max=100, value=45,
                                               marks={i: f'{i}%' for i in range(0, 101, 20)}, 
                                               tooltip={'placement': 'bottom', 'always_visible': True}),
                                ], width=12),
                            ], className="mb-3"),
                            
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Terminal G>A Damage"),
                                    dcc.Slider(id='terminal-ga-rate', min=0, max=100, value=35,
                                               marks={i: f'{i}%' for i in range(0, 101, 20)},
                                               tooltip={'placement': 'bottom', 'always_visible': True}),
                                ], width=12),
                            ], className="mb-3"),
                            
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Internal Damage Rate"),
                                    dcc.Slider(id='internal-damage', min=0, max=10, value=2,
                                               marks={i: f'{i}%' for i in range(0, 11, 2)},
                                               tooltip={'placement': 'bottom', 'always_visible': True}),
                                ], width=12),
                            ], className="mb-3"),
                            
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Other Mutations Rate"),
                                    dcc.Slider(id='other-mutations', min=0, max=10, value=1,
                                               marks={i: f'{i}%' for i in range(0, 11, 2)},
                                               tooltip={'placement': 'bottom', 'always_visible': True}),
                                ], width=12),
                            ], className="mb-3"),
                            
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Strand Bias"),
                                    dcc.Slider(id='strand-bias', min=0, max=100, value=50,
                                               marks={i: f'{i}%' for i in range(0, 101, 20)},
                                               tooltip={'placement': 'bottom', 'always_visible': True}),
                                ], width=12),
                            ], className="mb-3"),
                            
                            dcc.Graph(id='damage-pattern-preview')
                        ]),
                    ])
                ], className="mb-4"),
                
                dbc.Card([
                    dbc.CardBody([
                        html.H3("Quality Score Profile", className="card-title"),
                        html.Hr(),
                        
                        dbc.Form([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Mean Quality Score"),
                                    dbc.Input(id='mean-quality', type='number', value=30, min=0, max=40),
                                ], width=6),
                                dbc.Col([
                                    dbc.Label("Terminal Quality Decay"),
                                    dbc.Input(id='terminal-decay', type='number', value=5, min=0, max=20),
                                ], width=6),
                            ], className="mb-3"),
                            
                            dcc.Graph(id='quality-profile-preview')
                        ]),
                    ])
                ]),
            ], md=8)
        ]),
        
        dbc.Row([
            dbc.Col([
                dbc.Button("Generate BAM", id='generate-button', color="primary", size="lg", className="mt-4"),
                dcc.Download(id='download-bam')
            ], className="text-center")
        ])
    ], fluid=True)
])


@callback(
    Output('custom-fasta-status', 'children'),
    Input('custom-fasta', 'contents'),
    State('custom-fasta', 'filename')
)
def update_custom_fasta(contents, filename):
    if contents is None:
        return "No file uploaded"
    return f"Uploaded: {filename}"

@callback(
    [Output('length-distribution-preview', 'figure'),
     Output('damage-pattern-preview', 'figure'),
     Output('quality-profile-preview', 'figure')],
    [Input('num-reads', 'value'),
     Input('mean-length', 'value'),
     Input('length-std', 'value'),
     Input('min-length', 'value'),
     Input('max-length', 'value'),
     Input('terminal-ct-rate', 'value'),
     Input('terminal-ga-rate', 'value'),
     Input('internal-damage', 'value'),
     Input('other-mutations', 'value'),
     Input('strand-bias', 'value'),
     Input('mean-quality', 'value'),
     Input('terminal-decay', 'value'),
     Input('length-distribution-type', 'value'),
     Input('length-skew', 'value'),
     Input('length-kurtosis', 'value'),
     Input('damage-decay-pattern', 'value'),
     Input('quality-distribution', 'value'),
     Input('temperature-factor', 'value'),
     Input('homopolymer-error-rate', 'value')]
)
def update_previews(num_reads, mean_length, length_std, min_length, max_length,
                   terminal_ct, terminal_ga, internal_damage, other_mutations,
                   strand_bias, mean_quality, terminal_decay,
                   length_dist_type, length_skew, length_kurtosis,
                   damage_decay_pattern, quality_dist, temperature_factor, 
                   homopolymer_rate):
    # Set default values for all parameters
    defaults = {
        'num_reads': 1000,
        'mean_length': 150,
        'length_std': 0.4,
        'min_length': 35,
        'max_length': 500,
        'terminal_ct': 45,
        'terminal_ga': 35,
        'internal_damage': 2,
        'other_mutations': 1,
        'strand_bias': 50,
        'mean_quality': 30,
        'terminal_decay': 5,
        'length_distribution_type': 'lognormal',
        'length_skew': 0,
        'length_kurtosis': 0,
        'damage_decay_pattern': 'exponential',
        'quality_distribution': 'normal',
        'temperature_factor': 1.0,
        'homopolymer_rate': 1.0
    }
    
    # Create a dictionary of provided values, excluding None
    params = {
        'num_reads': num_reads,
        'mean_length': mean_length,
        'length_std': length_std,
        'min_length': min_length,
        'max_length': max_length,
        'distribution': length_dist_type,
        'length_skew': length_skew,
        'length_kurtosis': length_kurtosis
    }
    
    # Remove None values
    params = {k: v for k, v in params.items() if v is not None}
    
    # Update with defaults for missing values
    for k, v in defaults.items():
        if k not in params or params[k] is None:
            params[k] = v
    
    # Limit preview reads
    preview_reads = min(params['num_reads'], 10000)
    
    # Generate length distribution
    try:
        lengths = generate_fragment_lengths(preview_reads, params)
        
        # Create length distribution preview
        length_fig = go.Figure()
        length_fig.add_trace(go.Histogram(
            x=lengths,
            nbinsx=50,
            name='Fragment Lengths',
            marker_color=colors['line']
        ))
        length_fig.update_layout(
            title='Fragment Length Distribution',
            template='plotly_dark',
            xaxis_title='Length (bp)',
            yaxis_title='Count',
            plot_bgcolor=colors['plot_bg'],
            paper_bgcolor=colors['surface']
        )
    except Exception as e:
        print(f"Error generating length distribution: {e}")
        # Create empty figure if there's an error
        length_fig = go.Figure()
        length_fig.update_layout(title='Error generating length distribution')
    
    # Enhanced damage pattern preview
    x = list(range(100))
    damage_calc = DamagePattern()
    damage_calc.patterns['terminal_ct']['decay'] = damage_decay_pattern
    damage_calc.patterns['terminal_ga']['decay'] = damage_decay_pattern
    
    y_ct = []
    y_ga = []
    y_transitions = []
    y_transversions = []
    
    for i in x:
        # Calculate terminal damage rates with decay pattern
        ct_rate = damage_calc._calculate_position_specific_rate(
            terminal_ct/100 * temperature_factor,
            i, 100,
            {'length': 10, 'decay': damage_decay_pattern}
        )
        ga_rate = damage_calc._calculate_position_specific_rate(
            terminal_ga/100 * temperature_factor,
            99-i, 100,
            {'length': 10, 'decay': damage_decay_pattern}
        )
        
        y_ct.append(ct_rate)
        y_ga.append(ga_rate)
        y_transitions.append(other_mutations/100 * temperature_factor)
        y_transversions.append(other_mutations/200 * temperature_factor)
    
    damage_fig = go.Figure()
    damage_fig.add_trace(go.Scatter(x=x, y=y_ct, name='C→T frequency', line={'color': colors['line']}))
    damage_fig.add_trace(go.Scatter(x=x, y=y_ga, name='G→A frequency', line={'color': colors['highlight']}))
    damage_fig.add_trace(go.Scatter(x=x, y=y_transitions, name='Other transitions', line={'dash': 'dash', 'color': colors['highlight2']}))
    damage_fig.add_trace(go.Scatter(x=x, y=y_transversions, name='Transversions', line={'dash': 'dot', 'color': colors['muted']}))
    
    damage_fig.update_layout(
        title='Expected Damage Pattern',
        template='plotly_dark',
        xaxis_title='Position',
        yaxis_title='Mutation Frequency',
        plot_bgcolor=colors['plot_bg'],
        paper_bgcolor=colors['surface']
    )
    
    # Enhanced quality profile preview
    qual_profile = QualityProfile(
        mean_quality=mean_quality,
        terminal_decay=terminal_decay,
        quality_shape=quality_dist
    )
    
    # Generate multiple quality profiles to show variation
    n_profiles = 10
    qual_fig = go.Figure()
    
    for i in range(n_profiles):
        qualities = qual_profile.generate_qualities(100)
        qual_fig.add_trace(go.Scatter(
            x=list(range(100)),
            y=qualities,
            name=f'Profile {i+1}' if i == 0 else None,
            line={'color': colors['line'], 'width': 1 if i > 0 else 2},
            opacity=0.3 if i > 0 else 1
        ))
    
    qual_fig.update_layout(
        title='Quality Score Profile',
        template='plotly_dark',
        xaxis_title='Position',
        yaxis_title='Quality Score',
        plot_bgcolor=colors['plot_bg'],
        paper_bgcolor=colors['surface'],
        showlegend=False
    )
    
    return length_fig, damage_fig, qual_fig

@callback(
    Output('download-bam', 'data'),
    [Input('generate-button', 'n_clicks')],
    [State('num-reads', 'value'),
     State('reference-genome', 'value'),
     State('mean-length', 'value'),
     State('length-std', 'value'),
     State('min-length', 'value'),
     State('max-length', 'value'),
     State('terminal-ct-rate', 'value'),
     State('terminal-ga-rate', 'value'),
     State('internal-damage', 'value'),
     State('other-mutations', 'value'),
     State('strand-bias', 'value'),
     State('mean-quality', 'value'),
     State('terminal-decay', 'value'),
     State('length-distribution-type', 'value'),
     State('length-skew', 'value'),
     State('length-kurtosis', 'value'),
     State('damage-decay-pattern', 'value'),
     State('quality-distribution', 'value'),
     State('temperature-factor', 'value'),
     State('homopolymer-error-rate', 'value')],
    prevent_initial_call=True
)
def generate_bam(n_clicks, num_reads, reference_genome, mean_length, length_std,
                min_length, max_length, terminal_ct, terminal_ga, internal_damage,
                other_mutations, strand_bias, mean_quality, terminal_decay,
                length_dist_type, length_skew, length_kurtosis,
                damage_decay_pattern, quality_dist, temperature_factor,
                homopolymer_rate):
    if not n_clicks:
        return
    
    # Set default values
    if temperature_factor is None:
        temperature_factor = 1.0
    if homopolymer_rate is None:
        homopolymer_rate = 1.0
    if damage_decay_pattern is None:
        damage_decay_pattern = 'exponential'
    if quality_dist is None:
        quality_dist = 'normal'
    
    # Initialize enhanced damage pattern and quality profile
    damage_pattern = DamagePattern()
    quality_profile = QualityProfile(
        mean_quality=mean_quality if mean_quality is not None else 30,
        terminal_decay=terminal_decay if terminal_decay is not None else 5,
        quality_shape=quality_dist
    )
    
    # Update damage pattern values with all new parameters
    damage_pattern.patterns.update({
        'terminal_ct': {
            'rate': (terminal_ct if terminal_ct is not None else 45)/100,
            'length': 10,
            'decay': damage_decay_pattern
        },
        'terminal_ga': {
            'rate': (terminal_ga if terminal_ga is not None else 35)/100,
            'length': 10,
            'decay': damage_decay_pattern
        },
        'internal_ct': {'rate': (internal_damage if internal_damage is not None else 2)/100},
        'internal_ga': {'rate': (internal_damage if internal_damage is not None else 2)/100},
        'transitions': {'rate': (other_mutations if other_mutations is not None else 1)/100},
        'transversions': {'rate': (other_mutations if other_mutations is not None else 1)/200},
        'homopolymer_errors': {
            'rate': homopolymer_rate/100,
            'max_extension': 2
        }
    })
    
    # Generate fragment lengths with enhanced parameters
    lengths = generate_fragment_lengths(
        num_reads if num_reads is not None else 1000,
        {
            'mean_length': mean_length,
            'length_std': length_std,
            'min_length': min_length,
            'max_length': max_length,
            'distribution': length_dist_type,
            'length_skew': length_skew,
            'length_kurtosis': length_kurtosis,
            'bimodal_proportion': 0.7
        }
    )
    
    # Load reference sequence
    ref_seq = load_reference_sequence(reference_genome)
    
    # Generate reads with enhanced parameters
    reads = []
    for length in lengths:
        strand = 'positive' if random.random() < (strand_bias if strand_bias is not None else 50)/100 else 'negative'
        read = create_synthetic_read(
            reference_genome,
            length,
            damage_pattern,
            quality_profile,
            strand,
            ref_seq,
            temperature_factor=temperature_factor
        )
        reads.append(read)
    
    # Create BAM file
    output_path = Path('temp_synthetic.bam')
    write_bam_file(reads, output_path, create_bam_header(reference_genome), ref_seq)
    
    return dcc.send_file(
        output_path,
        filename=f'synthetic_adna_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.bam'
    )

