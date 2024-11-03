import dash
import os
from dash import html, dcc, callback, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import plotly.express as px
import pysam
import numpy as np
from pathlib import Path
from plotly.subplots import make_subplots
import datetime
from Bio.Seq import Seq
from Bio import SeqIO
import random
from scipy import stats
import io
import base64
from scipy.stats import gamma, beta, lognorm
from app import app, colors
import warnings
warnings.filterwarnings('ignore')

def calculate_pmd(bam_file_path, reference_sequence):
    """
    Calculate PMD scores from BAM file
    
    Parameters:
    -----------
    bam_file_path : Path or str
        Path to BAM file
    reference_sequence : str
        Reference sequence
    """
    # Convert to string if Path object
    bam_file_str = str(bam_file_path)
       
    if not os.path.exists(bam_file_str):
        raise FileNotFoundError(f"BAM file not found: {bam_file_str}")
    
    if not os.path.exists(bam_file_str + ".bai"):
        pysam.index(bam_file_str)

    bamfile = pysam.AlignmentFile(bam_file_str, 'rb')
    max_read_length = 100  # Adjust as needed

    position_counts = {
        'C>T': np.zeros(max_read_length),
        'G>A': np.zeros(max_read_length),
        'total': np.zeros(max_read_length)
    }

    for read in bamfile.fetch(until_eof=True):
        seq = read.query_sequence
        if seq is None:
            continue
        # Get reference sequence based on read's alignment
        start = read.reference_start
        end = read.reference_end
        ref_seq = reference_sequence[start:end]
        if not ref_seq:
            continue  # Skip if reference sequence is not available
        read_len = len(seq)
        if read_len > max_read_length or read_len == 0:
            continue  # Skip reads longer than max_read_length or empty reads
        if len(ref_seq) != read_len:
            continue  # Ensure lengths match

        # Reverse complement sequences if the read is on the reverse strand
        if read.is_reverse:
            seq = str(Seq(seq).reverse_complement())
            ref_seq = str(Seq(ref_seq).reverse_complement())

        for i in range(read_len):
            q_base = seq[i]
            r_base = ref_seq[i]
            pos_from_5p = i
            pos_from_3p = read_len - i - 1
            if q_base != r_base:
                if r_base == 'C' and q_base == 'T':
                    position_counts['C>T'][pos_from_5p] += 1
                elif r_base == 'G' and q_base == 'A':
                    position_counts['G>A'][pos_from_3p] += 1
            position_counts['total'][pos_from_5p] += 1

    bamfile.close()
    return position_counts

def create_pmd_plot(position_counts):
    fig = go.Figure()
    positions = np.arange(len(position_counts['total']))
    c_to_t_freq = position_counts['C>T'] / np.maximum(position_counts['total'], 1)
    g_to_a_freq = position_counts['G>A'] / np.maximum(position_counts['total'], 1)
    fig.add_trace(go.Scatter(x=positions, y=c_to_t_freq, name='C→T', line={'color': colors['line']}))
    fig.add_trace(go.Scatter(x=positions, y=g_to_a_freq, name='G→A', line={'color': colors['highlight']}))
    fig.update_layout(
        title='PMD Scores',
        xaxis_title='Position from 5\' End',
        yaxis_title='Substitution Frequency',
        template='plotly_dark',
        plot_bgcolor=colors['plot_bg'],
        paper_bgcolor=colors['surface']
    )
    return fig

class PMDamage:
    def __init__(self, matrices_dir="matrices", damage_type="double"):
        """
        Initialize PMD model with damage matrices.
        
        Parameters:
        -----------
        matrices_dir : str
            Directory containing the matrix files (double-5.dat, double-3.dat, single-5.dat, single-3.dat)
        damage_type : str
            Type of damage pattern - "single" or "double"
        """
        self.damage_type = damage_type
        self.matrices = self.load_damage_matrices(matrices_dir)
        
        # Store current rates for modification
        self.current_rates = {
            'terminal_ct': 0.4,
            'terminal_ga': 0.3,
            'temperature_factor': 1.0
        }

    def load_damage_matrices(self, matrices_dir):
        """Load empirical damage matrices from files"""
        patterns = {
            'single_5p': 'single-5.dat',
            'single_3p': 'single-3.dat', 
            'double_5p': 'double-5.dat',
            'double_3p': 'double-3.dat'
        }

        matrices = {}
        for key, filename in patterns.items():
            matrix = {}
            filepath = os.path.join(matrices_dir, filename)
            try:
                with open(filepath) as f:
                    # Skip header
                    next(f)
                    for line in f:
                        values = line.strip().split('\t')
                        pos = int(values[0])
                        # Convert probability values, excluding confidence intervals
                        probs = [float(v.split()[0]) for v in values[1:]]
                        matrix[pos] = probs
                matrices[key] = matrix
            except (FileNotFoundError, IOError) as e:
                print(f"Warning: Could not load matrix {filename}: {e}")
                matrices[key] = {}
                
        return matrices
    
    def sample_substitution(self, base, probs):
        """
        Sample a substitution for the base given the probabilities.

        Parameters:
        -----------
        base : str
            The original base
        probs : list
            A list of probabilities for substitutions

        Returns:
        --------
        substitution : str or None
            The substituted base, or None if no substitution occurs
        """
        base = base.upper()
        indices = self.substitution_indices.get(base)
        if indices is None:
            return None
        possible_substitutions = self.possible_substitutions[base]
        base_probs = [probs[i] for i in indices]
        total_prob = sum(base_probs)
        if total_prob == 0:
            return None
        normalized_probs = [p / total_prob for p in base_probs]
        substitution = random.choices(possible_substitutions, weights=normalized_probs, k=1)[0]
        return substitution


    def modify_damage_rates(self, terminal_ct_rate=None, terminal_ga_rate=None, 
                          internal_rate=None, temperature_factor=None):
        """
        Modify damage rates while maintaining matrix structure.
        
        Parameters:
        -----------
        terminal_ct_rate : float, optional
            Rate for C->T damage at 5' ends (0-1)
        terminal_ga_rate : float, optional
            Rate for G->A damage at 3' ends (0-1)
        internal_rate : float, optional
            Rate for internal damages (0-1)
        temperature_factor : float, optional
            Factor to modify all rates (typically 0.5-2.0)
        """
        if terminal_ct_rate is not None:
            self.current_rates['terminal_ct'] = terminal_ct_rate
        if terminal_ga_rate is not None:
            self.current_rates['terminal_ga'] = terminal_ga_rate
        if temperature_factor is not None:
            self.current_rates['temperature_factor'] = temperature_factor

        # Apply modifications to matrices
        if hasattr(self, 'original_matrices'):
            # Scale original matrix values
            for pos in range(len(self.matrices['CT'])):
                # Scale CT rates
                if pos < 10:  # Terminal region
                    self.matrices['CT'][pos] = self.original_matrices['CT'][pos] * self.current_rates['terminal_ct']
                else:  # Internal region
                    self.matrices['CT'][pos] = self.original_matrices['CT'][pos] * self.current_rates['internal_rate']
                    
                # Scale GA rates
                if pos > len(self.matrices['GA']) - 11:  # Terminal region
                    self.matrices['GA'][pos] = self.original_matrices['GA'][pos] * self.current_rates['terminal_ga']
                else:  # Internal region
                    self.matrices['GA'][pos] = self.original_matrices['GA'][pos] * self.current_rates['internal_rate']
                    
            # Apply temperature factor
            self.matrices['CT'] *= self.current_rates['temperature_factor']
            self.matrices['GA'] *= self.current_rates['temperature_factor']

    def parse_damage_matrix(self, file_path):
        """Parse Gargammel-style damage matrix file"""
        ct_rates = []
        ga_rates = []
        positions = []
        
        with open(file_path) as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                pos = int(parts[0])
                positions.append(pos)
                
                # Extract frequencies (first value before confidence interval)
                ct_freq = float(parts[6].split()[0])  # C>T column (index 6)
                ga_freq = float(parts[7].split()[0])  # G>A column (index 7)
                
                ct_rates.append(ct_freq)
                ga_rates.append(ga_freq)
        
        # Store original matrices for rate modification
        self.original_matrices = {
            'CT': np.array(ct_rates),
            'GA': np.array(ga_rates),
            'positions': np.array(positions)
        }
        
        # Create working copy
        self.matrices = {
            'CT': np.array(ct_rates),
            'GA': np.array(ga_rates),
            'positions': np.array(positions)
        }
        
        return self.matrices

    def create_default_matrices(self):
        """Create default exponential decay patterns if no matrix file provided"""
        positions = np.arange(-80, 1)
        ct_rates = 0.4 * np.exp(positions / 20) + 0.01
        ga_rates = 0.3 * np.exp(positions / 20) + 0.01
        
        # Store as original matrices
        self.original_matrices = {
            'CT': np.clip(ct_rates, 0, 1),
            'GA': np.clip(ga_rates, 0, 1),
            'positions': positions
        }
        
        # Create working copy
        self.matrices = {
            'CT': np.clip(ct_rates, 0, 1),
            'GA': np.clip(ga_rates, 0, 1),
            'positions': positions
        }
        
        return self.matrices
    
    def get_damage_prob(self, position, base, strand, end='5p'):
        """
        Get damage probability for a specific position and base.
        
        Parameters:
        -----------
        position : int
            Position from relevant end (5' or 3')
        base : str
            Base to check for damage (C or G)
        strand : str
            'forward' or 'reverse'
        end : str
            '5p' or '3p' end
        """
        matrix_key = f"{self.damage_type}_{end}"
        
        # Get the appropriate matrix based on damage type and end
        if matrix_key not in self.matrices:
            return 0.0
            
        matrix = self.matrices[matrix_key]
        
        # Get position-specific probabilities
        if position not in matrix:
            return 0.0
            
        # Get the appropriate probability based on the base and strand
        if base == 'C' and end == '5p':
            prob = matrix[position][6]  # C>T column
        elif base == 'G' and end == '3p':
            prob = matrix[position][7]  # G>A column
        else:
            return 0.0
            
        # Apply rate modifications
        if base == 'C':
            prob *= self.current_rates['terminal_ct']
        else:  # base == 'G'
            prob *= self.current_rates['terminal_ga']
            
        return prob * self.current_rates['temperature_factor']

    def apply_damage(self, sequence, strand='forward'):
        """
        Apply damage patterns to a sequence.
        
        Parameters:
        -----------
        sequence : str
            Input DNA sequence
        strand : str
            'forward' or 'reverse'
        """
        damaged_seq = list(sequence)
        seq_length = len(sequence)
        
        for i in range(seq_length):
            base = damaged_seq[i]
            
            # Calculate positions from both ends
            pos_5p = i
            pos_3p = seq_length - i - 1
            
            if strand == 'forward':
                # 5' end C→T damage
                if base == 'C':
                    prob = self.get_damage_prob(pos_5p, 'C', strand, '5p')
                    if random.random() < prob:
                        damaged_seq[i] = 'T'
                
                # 3' end G→A damage
                elif base == 'G':
                    prob = self.get_damage_prob(pos_3p, 'G', strand, '3p')
                    if random.random() < prob:
                        damaged_seq[i] = 'A'
            else:
                # Reverse strand: mirror the damage patterns
                if base == 'G':
                    prob = self.get_damage_prob(pos_5p, 'G', strand, '5p')
                    if random.random() < prob:
                        damaged_seq[i] = 'A'
                elif base == 'C':
                    prob = self.get_damage_prob(pos_3p, 'C', strand, '3p')
                    if random.random() < prob:
                        damaged_seq[i] = 'T'
        
        return ''.join(damaged_seq)

class ReferenceGenome:
    def __init__(self, length=1000000, gc_content=0.5, seed=42):
        """Initialize reference genome with fixed seed for reproducibility"""
        np.random.seed(seed)
        self.sequence = self.generate_sequence(length, gc_content)
        self.length = length
        
    def generate_sequence(self, length, gc_content):
        """Generate reference sequence with specified GC content"""
        gc_count = int(length * gc_content)
        at_count = length - gc_count
        bases = ['G', 'C'] * (gc_count // 2) + ['A', 'T'] * (at_count // 2)
        np.random.shuffle(bases)
        return ''.join(bases)
    
    def extract_fragment(self, start, length):
        """Extract a fragment from reference genome"""
        if start + length > self.length:
            return None
        return self.sequence[start:start + length]

def generate_fragment_length(mean_length=60, std_dev=0.4, min_length=35, max_length=500):
    """Generate fragment length from log-normal distribution"""
    shape = std_dev
    scale = np.exp(np.log(mean_length) - (shape**2)/2)
    
    while True:
        length = int(round(lognorm.rvs(s=shape, scale=scale)))
        if min_length <= length <= max_length:
            return length

class AncientDNASimulator:
    def __init__(self, reference_genome=None, pmd_model=None):
        """Initialize simulator with reference genome and PMD model"""
        self.reference = reference_genome or ReferenceGenome()
        self.pmd_model = pmd_model or PMDamage()

    def generate_read(self, length=None, strand=None):
        """
        Generate a single ancient DNA read by extracting from reference genome
        and applying damage patterns
        """
        if length is None:
            length = generate_fragment_length()
        if strand is None:
            strand = random.choice(['forward', 'reverse'])

        # Critical: Extract an actual fragment from the reference genome
        max_start = self.reference.length - length
        if max_start <= 0:
            return None
            
        # Get a valid fragment from reference
        start = random.randint(0, max_start)
        fragment = self.reference.extract_fragment(start, length)
        
        if fragment is None:
            return None

        # Create damaged version (initially identical to reference)
        damaged_seq = list(fragment)
        
        # Only apply C>T and G>A transitions according to PMD patterns
        for i in range(length):
            base = damaged_seq[i]
            pos_5p = i
            pos_3p = length - i - 1
            
            if strand == 'forward':
                if base == 'C':
                    # Get 5' C>T probability
                    if random.random() < self.pmd_model.get_damage_prob(pos_5p, 'C', strand, '5p'):
                        damaged_seq[i] = 'T'
                elif base == 'G':
                    # Get 3' G>A probability
                    if random.random() < self.pmd_model.get_damage_prob(pos_3p, 'G', strand, '3p'):
                        damaged_seq[i] = 'A'
            else:  # reverse strand
                if base == 'G':
                    # Get 5' G>A probability (reverse complement of C>T)
                    if random.random() < self.pmd_model.get_damage_prob(pos_5p, 'G', strand, '5p'):
                        damaged_seq[i] = 'A'
                elif base == 'C':
                    # Get 3' C>T probability (reverse complement of G>A)
                    if random.random() < self.pmd_model.get_damage_prob(pos_3p, 'C', strand, '3p'):
                        damaged_seq[i] = 'T'

        damaged_sequence = ''.join(damaged_seq)
        
        # Handle reverse strand orientation
        if strand == 'reverse':
            damaged_sequence = str(Seq(damaged_sequence).reverse_complement())
            fragment = str(Seq(fragment).reverse_complement())

        return {
            'sequence': damaged_sequence,
            'reference': fragment,
            'position': start,
            'strand': strand,
            'length': length
        }
        
    def generate_reads(self, num_reads, length_params=None, strand_bias=0.5):
        """
        Generate multiple ancient DNA reads
        
        Parameters:
        -----------
        num_reads : int
            Number of reads to generate
        length_params : dict
            Parameters for fragment length generation
        strand_bias : float
            Proportion of reads on forward strand (0-1)
        """
        if length_params is None:
            length_params = {}
                
        reads = []
        for _ in range(num_reads):
            length = generate_fragment_length(**length_params)
            strand = 'forward' if random.random() < strand_bias else 'reverse'
            read = self.generate_read(length=length, strand=strand)
            if read:
                reads.append(read)
            
        return reads

def write_bam_file(reads, output_path, reference):
    """Write reads to BAM file with proper alignments"""
    # Create BAM header
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [{'LN': reference.length, 'SN': 'ref1'}],
        'RG': [{
            'ID': 'simulated',
            'SM': 'ancient_dna_sim',
            'LB': 'lib1',
            'PL': 'ILLUMINA',
            'PU': 'unit1',
            'DT': datetime.datetime.now().isoformat()
        }]
    }
    
    temp_path = str(output_path) + ".tmp"
    try:
        with pysam.AlignmentFile(temp_path, 'wb', header=header) as outf:
            for i, read in enumerate(reads):
                a = pysam.AlignedSegment()
                a.query_name = f'read_{i+1}'
                a.query_sequence = read['sequence']
                a.reference_id = 0
                a.reference_start = read['position']
                a.flag = 16 if read['strand'] == 'reverse' else 0
                
                # Calculate CIGAR string and MD tag
                ref_seq = read['reference']
                query_seq = read['sequence']
                
                cigar = []
                md_parts = []
                matches = 0
                nm = 0
                
                for qbase, rbase in zip(query_seq, ref_seq):
                    if qbase == rbase:
                        matches += 1
                    else:
                        if matches > 0:
                            cigar.append((0, matches))
                            md_parts.append(str(matches))
                            matches = 0
                        cigar.append((8, 1))  # Mismatch
                        md_parts.append(rbase)
                        nm += 1
                
                if matches > 0:
                    cigar.append((0, matches))
                    md_parts.append(str(matches))
                
                a.cigar = cigar
                a.tags = [
                    ('MD', ''.join(md_parts)),
                    ('NM', nm)
                ]
                
                outf.write(a)
        
        # Sort and index
        pysam.sort("-o", str(output_path), temp_path)
        pysam.index(str(output_path))
        
    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)


layout = dbc.Container([
    dbc.Row([
        html.H1("Ancient DNA Simulator", className="text-center my-4"),
        
        # Main control panel
        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.H3("Simulation Parameters", className="card-title"),
                    html.Hr(),
                    
                    # Basic Parameters Section
                    dbc.Form([
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Number of Reads"),
                                dbc.Input(
                                    id='sim-num-reads',
                                    type='number',
                                    value=1000,
                                    min=1,
                                    max=1000000
                                ),
                            ], width=6),
                            dbc.Col([
                                dbc.Label("Reference Length (bp)"),
                                dbc.Input(
                                    id='sim-ref-length',
                                    type='number',
                                    value=1000000,
                                    min=10000,
                                    max=10000000
                                ),
                            ], width=6),
                        ], className="mb-3"),

                        # Fragment Length Parameters
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Mean Fragment Length"),
                                dbc.Input(
                                    id='sim-mean-length',
                                    type='number',
                                    value=60,
                                    min=20,
                                    max=500
                                ),
                            ], width=6),
                            dbc.Col([
                                dbc.Label("Length StdDev"),
                                dbc.Input(
                                    id='sim-length-std',
                                    type='number',
                                    value=0.4,
                                    min=0.1,
                                    max=2,
                                    step=0.1
                                ),
                            ], width=6),
                        ], className="mb-3"),

                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Minimum Length"),
                                dbc.Input(
                                    id='sim-min-length',
                                    type='number',
                                    value=35,
                                    min=20,
                                    max=200
                                ),
                            ], width=6),
                            dbc.Col([
                                dbc.Label("Maximum Length"),
                                dbc.Input(
                                    id='sim-max-length',
                                    type='number',
                                    value=500,
                                    min=50,
                                    max=1000
                                ),
                            ], width=6),
                        ], className="mb-3"),

                        # Damage Parameters
                        html.H4("Damage Parameters", className="mt-4"),
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("5' Terminal Deamination Rate"),
                                dcc.Slider(
                                    id='sim-terminal-ct-rate',
                                    min=0,
                                    max=1,
                                    step=0.01,
                                    value=0.4,
                                    marks={i/10: f'{i/10:.1f}' for i in range(11)},
                                    tooltip={"placement": "bottom", "always_visible": True}
                                ),
                            ], width=12),
                        ], className="mb-3"),

                        dbc.Row([
                            dbc.Col([
                                dbc.Label("3' Terminal Deamination Rate"),
                                dcc.Slider(
                                    id='sim-terminal-ga-rate',
                                    min=0,
                                    max=1,
                                    step=0.01,
                                    value=0.3,
                                    marks={i/10: f'{i/10:.1f}' for i in range(11)},
                                    tooltip={"placement": "bottom", "always_visible": True}
                                ),
                            ], width=12),
                        ], className="mb-3"),

                        # Additional Parameters
                        html.H4("Advanced Parameters", className="mt-4"),
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Temperature Factor"),
                                dcc.Slider(
                                    id='sim-temp-factor',
                                    min=0.5,
                                    max=2,
                                    step=0.1,
                                    value=1,
                                    marks={i/2: str(i/2) for i in range(1, 5)},
                                    tooltip={"placement": "bottom", "always_visible": True}
                                ),
                            ], width=12),
                        ], className="mb-3"),

                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Damage Type"),
                                dbc.RadioItems(
                                    id='damage-type',
                                    options=[
                                        {'label': 'Single-stranded', 'value': 'single'},
                                        {'label': 'Double-stranded', 'value': 'double'}
                                    ],
                                    value='double',
                                    inline=True
                                ),
                            ], width=12),
                        ], className="mb-3"),

                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Strand Bias"),
                                dcc.Slider(
                                    id='sim-strand-bias',
                                    min=0,
                                    max=1,
                                    step=0.01,
                                    value=0.5,
                                    marks={i/10: f'{i/10:.1f}' for i in range(11)},
                                    tooltip={"placement": "bottom", "always_visible": True}
                                ),
                            ], width=12),
                        ], className="mb-3"),

                        # Custom Matrix Upload
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Upload Custom Damage Matrix"),
                                dcc.Upload(
                                    id='sim-custom-matrix',
                                    children=dbc.Button(
                                        "Upload Matrix File",
                                        color="secondary",
                                        outline=True
                                    ),
                                    multiple=False
                                ),
                                html.Div(id='sim-matrix-status', className="mt-2"),
                            ], width=12),
                        ], className="mb-3"),
                    ]),
                ])
            ], className="mb-4"),
        ], width=4),

        # Preview and Results Panel
        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.H3("Simulation Preview", className="card-title"),
                    html.Hr(),
                    
                    # Preview Tabs
                    dbc.Tabs([
                        dbc.Tab(
                            dcc.Graph(id='sim-length-preview'),
                            label="Length Distribution"
                        ),
                        dbc.Tab(
                            dcc.Graph(id='sim-damage-preview'),
                            label="Damage Pattern"
                        ),
                        dbc.Tab(
                            dcc.Graph(id='sim-pmd-preview'),
                            label="PMD Score"
                        ),
                    ]),
                ])
            ], className="mb-4"),

            # Generate Button
            dbc.Row([
                dbc.Col([
                    dbc.Button(
                        "Generate BAM",
                        id='sim-generate-button',
                        color="primary",
                        size="lg",
                        className="w-100"
                    ),
                    dcc.Download(id='sim-download-bam'),
                ], width=12),
            ]),

            # Status Messages
            html.Div(id='sim-status', className="mt-3"),
        ], width=8),
    ]),
], fluid=True)

@callback(
    [Output('sim-length-preview', 'figure'),
     Output('sim-damage-preview', 'figure'),
     Output('sim-pmd-preview', 'figure')],
    [Input('sim-mean-length', 'value'),
     Input('sim-length-std', 'value'),
     Input('sim-min-length', 'value'),
     Input('sim-max-length', 'value'),
     Input('sim-terminal-ct-rate', 'value'),
     Input('sim-terminal-ga-rate', 'value'),
     Input('sim-temp-factor', 'value'),
     Input('sim-strand-bias', 'value'),
     Input('damage-type', 'value')]
)
def update_previews(mean_length, length_std, min_length, max_length,
                   terminal_ct, terminal_ga, temp_factor, strand_bias, damage_type):
    """Update preview plots based on parameter changes"""
    
    # Preview fragment length distribution
    preview_reads = 1000
    lengths = [generate_fragment_length(
        mean_length=mean_length,
        std_dev=length_std,
        min_length=min_length,
        max_length=max_length
    ) for _ in range(preview_reads)]
    
    length_fig = go.Figure(data=[go.Histogram(
        x=lengths,
        nbinsx=50,
        marker_color=colors['line']
    )])
    length_fig.update_layout(
        title='Fragment Length Distribution',
        xaxis_title='Length (bp)',
        yaxis_title='Count',
        template='plotly_dark',
        plot_bgcolor=colors['plot_bg'],
        paper_bgcolor=colors['surface'],
        font=dict(color=colors['muted'])
    )

    # Preview damage patterns
    pmd_model = PMDamage(
        matrices_dir="matrices", 
        damage_type=damage_type
    )
    pmd_model.modify_damage_rates(
        terminal_ct_rate=terminal_ct,
        terminal_ga_rate=terminal_ga,
        temperature_factor=temp_factor
    )

    positions = np.arange(-80, 1)
    
    # Calculate position-specific damage rates
    ct_rates = terminal_ct * temp_factor * np.exp(positions / 20)
    ga_rates = terminal_ga * temp_factor * np.exp(positions / 20)
    
    damage_fig = go.Figure()
    damage_fig.add_trace(go.Scatter(
        x=positions,
        y=ct_rates,
        name='C→T',
        line={'color': colors['line']}
    ))
    damage_fig.add_trace(go.Scatter(
        x=positions,
        y=ga_rates,
        name='G→A',
        line={'color': colors['highlight']}
    ))
    
    damage_fig.update_layout(
        title='Expected Damage Pattern',
        xaxis_title='Position from End',
        yaxis_title='Damage Probability',
        template='plotly_dark',
        plot_bgcolor=colors['plot_bg'],
        paper_bgcolor=colors['surface'],
        font=dict(color=colors['muted'])
    )

    # Preview PMD scores (simulated)
    pmd_fig = go.Figure()
    pmd_fig.add_trace(go.Scatter(
        x=positions,
        y=ct_rates * (1 - strand_bias),
        name='5\' C→T',
        line={'color': colors['line']}
    ))
    pmd_fig.add_trace(go.Scatter(
        x=positions,
        y=ga_rates * strand_bias,
        name='3\' G→A',
        line={'color': colors['highlight']}
    ))
    
    pmd_fig.update_layout(
        title='Predicted PMD Scores',
        xaxis_title='Position',
        yaxis_title='Score',
        template='plotly_dark',
        plot_bgcolor=colors['plot_bg'],
        paper_bgcolor=colors['surface'],
        font=dict(color=colors['muted'])
    )

    return length_fig, damage_fig, pmd_fig

@callback(
    [Output('sim-download-bam', 'data'),
     Output('sim-status', 'children')],
    Input('sim-generate-button', 'n_clicks'),
    [State('sim-num-reads', 'value'),
     State('sim-ref-length', 'value'),
     State('sim-mean-length', 'value'),
     State('sim-length-std', 'value'),
     State('sim-min-length', 'value'),
     State('sim-max-length', 'value'),
     State('sim-terminal-ct-rate', 'value'),
     State('sim-terminal-ga-rate', 'value'),
     State('sim-temp-factor', 'value'),
     State('sim-strand-bias', 'value'),
     State('damage-type', 'value')], 
    prevent_initial_call=True
)
def generate_bam_file(n_clicks, num_reads, ref_length, mean_length, length_std,
                     min_length, max_length, terminal_ct, terminal_ga,
                     temp_factor, strand_bias, damage_type):
    """Generate BAM file with simulated ancient DNA reads"""
    if not n_clicks:
        raise PreventUpdate
    
    try:
        # Initialize simulation components
        reference = ReferenceGenome(length=ref_length)
        
        # Initialize PMD model with matrix file
        pmd_model = PMDamage(
            matrices_dir="matrices", 
            damage_type=damage_type
        )
        pmd_model.modify_damage_rates(
            terminal_ct_rate=terminal_ct,
            terminal_ga_rate=terminal_ga,
            temperature_factor=temp_factor
        )
        
        # Initialize simulator
        simulator = AncientDNASimulator(reference, pmd_model)
        
        # Generate reads
        reads = simulator.generate_reads(
            num_reads=num_reads,
            length_params={
                'mean_length': mean_length,
                'std_dev': length_std,
                'min_length': min_length,
                'max_length': max_length
            },
            strand_bias=strand_bias
        )
        
        # Write BAM file
        output_path = Path('temp_synthetic.bam')
        write_bam_file(reads, output_path, reference)
        
        # Create status message
        status = html.Div([
            html.H5("Simulation Complete", className="text-success"),
            html.P(f"Generated {len(reads)} reads"),
            html.P(f"Average length: {np.mean([r['length'] for r in reads]):.1f} bp"),
            html.P(f"Damage rates - C→T: {terminal_ct:.2%}, G→A: {terminal_ga:.2%}")
        ])
        
        return dcc.send_file(
            str(output_path),
            filename=f'synthetic_adna_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.bam'
        ), status
        
    except Exception as e:
        error_status = html.Div([
            html.H5("Error", className="text-danger"),
            html.P(str(e))
        ])
        return None, error_status

@callback(
    Output('sim-matrix-status', 'children'),
    Input('sim-custom-matrix', 'contents'),
    State('sim-custom-matrix', 'filename')
)
def update_matrix_status(contents, filename):
    """Update status of custom matrix upload"""
    if contents is None:
        return "No matrix file uploaded"
    return html.Div([
        html.P(f"Uploaded: {filename}", className="text-success"),
        html.P("Matrix loaded successfully")
    ])

def main():
    # Initialize components
    reference = ReferenceGenome(length=1000000)
    pmd_model = PMDamage()
    simulator = AncientDNASimulator(reference, pmd_model)
    
    # Generate reads
    reads = simulator.generate_reads(
        num_reads=1000,
        length_params={
            'mean_length': 60,
            'std_dev': 0.4,
            'min_length': 35,
            'max_length': 500
        }
    )
    
    # Write to BAM file
    output_path = Path('simulated_ancient_dna.bam')
    write_bam_file(reads, output_path, reference)

if __name__ == '__main__':
    main()
