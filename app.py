import pysam
import sys
import threading
from Bio import SeqIO
from Bio.Seq import Seq
import os
import numpy as np
import dash
import dash_bootstrap_components as dbc
from dash import dcc, html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash_iconify import DashIconify
from plotly.subplots import make_subplots
import dash_dangerously_set_inner_html
import plotly.graph_objects as go
import plotly.express as px
import plotly.colors
import pandas as pd
import base64
import tempfile
import uuid
import subprocess
from functools import lru_cache
from celery import Celery
import yaml
import zipfile
from ansi2html import Ansi2HTMLConverter
from itertools import islice
import re
import dust_module


from utils.tasks import celery_app, run_mapdamage_task, run_centrifuge_task
from pages import simulation
from utils.files import TempDirectoryManager

###############################################################################################

 ######  ########    ###    ########  ########    ######  ######## ######## ##     ## ########  
##    ##    ##      ## ##   ##     ##    ##      ##    ## ##          ##    ##     ## ##     ## 
##          ##     ##   ##  ##     ##    ##      ##       ##          ##    ##     ## ##     ## 
 ######     ##    ##     ## ########     ##       ######  ######      ##    ##     ## ########  
      ##    ##    ######### ##   ##      ##            ## ##          ##    ##     ## ##        
##    ##    ##    ##     ## ##    ##     ##      ##    ## ##          ##    ##     ## ##        
 ######     ##    ##     ## ##     ##    ##       ######  ########    ##     #######  ##

###############################################################################################

try:
    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f)
except FileNotFoundError:
    config = {}

default_centrifuge_db_path = config.get('centrifuge_db_path', '')

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SOLAR, "/assets/custom.css"], suppress_callback_exceptions=True) # general dash theme, can be changed but colors need to be adjusted accordingly (if you want a pretty dashboard)


# custom temp directory for storing the data (needs better implementation, but works fine)
temp_manager = TempDirectoryManager(max_age_hours=12)
CUSTOM_TEMP_DIR = temp_manager.create_temp_dir()
os.environ["TMPDIR"] = CUSTOM_TEMP_DIR

# Set TMPDIR environment variable to ensure pysam uses the custom temporary directory
os.environ["TMPDIR"] = CUSTOM_TEMP_DIR

# color palette
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



#################################################################

##     ## ######## ##       ########  ######## ########   ######  
##     ## ##       ##       ##     ## ##       ##     ## ##    ## 
##     ## ##       ##       ##     ## ##       ##     ## ##       
######### ######   ##       ########  ######   ########   ######  
##     ## ##       ##       ##        ##       ##   ##         ## 
##     ## ##       ##       ##        ##       ##    ##  ##    ## 
##     ## ######## ######## ##        ######## ##     ##  ######

#################################################################

def rgb_to_hex(rgb):
    """Convert an RGB tuple to a hex string, rounding the values to integers - small unnecessary helper function"""
    return '#{:02x}{:02x}{:02x}'.format(int(round(rgb[0])), int(round(rgb[1])), int(round(rgb[2])))

def prepare_histogram_data(read_lengths):
    if not read_lengths:
        # If read_lengths is empty, return empty arrays
        return np.array([]), np.array([]), np.array([])

    lengths = [length for length, _, _, _, _, _ in read_lengths]
    bins = np.arange(min(lengths), max(lengths) + 1, 1)
    hist, bin_edges = np.histogram(lengths, bins=bins)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    return bin_centers, hist, bins


def save_uploaded_file(content, filename):
    import base64
    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)
    temp_dir = os.path.join(CUSTOM_TEMP_DIR, 'uploads')
    os.makedirs(temp_dir, exist_ok=True)
    file_path = os.path.join(temp_dir, filename)
    with open(file_path, 'wb') as f:
        f.write(decoded)
    return file_path

def complement_base(base):
    base = base.upper()
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return complement.get(base, 'N')

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    reversed_seq = reversed(seq.upper())
    return ''.join([complement.get(base, 'N') for base in reversed_seq])

def parse_misincorporation(file_path):
    df = pd.read_csv(file_path, sep='\t', comment='#')
    return df

def parse_dnacomp(file_path):
    df = pd.read_csv(file_path, sep='\t', comment='#')
    return df

def create_misincorporation_plot(df):
    # Filter data as needed
    df = df[df['End'] == '5p']
    fig = go.Figure()
    mutations = ['C>T', 'G>A', 'A>G', 'T>C']
    for mut in mutations:
        fig.add_trace(go.Scatter(
            x=df['Pos'],
            y=df[mut],
            mode='lines',
            name=mut
        ))
    fig.update_layout(
        title='Misincorporation Pattern',
        xaxis_title='Position',
        yaxis_title='Frequency',
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        font=dict(color=colors['muted'])
    )
    return fig

def create_dnacomp_plot(df):
    fig = go.Figure()
    bases = ['A', 'C', 'G', 'T']
    for base in bases:
        fig.add_trace(go.Scatter(
            x=df['Pos'],
            y=df[base],
            mode='lines',
            name=base
        ))
    fig.update_layout(
        title='DNA Composition',
        xaxis_title='Position',
        yaxis_title='Frequency',
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        font=dict(color=colors['muted'])
    )
    return fig

def display_certainty_metrics(mcmc_df):
    misincorporation_rate = mcmc_df['delta'].iloc[0]
    fragmentation_rate = mcmc_df['lambda'].iloc[0]

    return html.Div([
        html.H4("Certainty Metrics"),
        html.P(f"Misincorporation Rate (delta): {misincorporation_rate:.4f}"),
        html.P(f"Fragmentation Rate (lambda): {fragmentation_rate:.4f}")
    ])

def extract_mismatches(read):
    """
    Extract mismatches from a read using the MD tag and CIGAR string.
    Returns a list of dictionaries with keys 'position', 'ref_base', 'read_base'.
    Positions are adjusted based on the read's orientation.
    """
    try:
        md_tag = read.get_tag('MD')
    except KeyError:
        # MD tag not present
        return []

    # Initialize variables
    mismatches = []
    query_sequence = read.query_sequence
    if query_sequence is None:
        return []
    if read.is_reverse:
        query_sequence = reverse_complement(query_sequence)

    cigar_tuples = read.cigartuples
    md_string = md_tag
    seq_pos = 0  # Position in query_sequence
    ref_pos = 0  # Position in reference

    import re
    md_tokens = re.findall(r'(\d+|\^[A-Z]+|[A-Z])', md_string)

    for token in md_tokens:
        if token.isdigit():
            num_matches = int(token)
            seq_pos += num_matches
            ref_pos += num_matches
        elif token.startswith('^'):
            # Deletion in reference (insertion in read)
            deletion_length = len(token) - 1
            ref_pos += deletion_length
        else:
            # Mismatch
            ref_base = token
            if seq_pos >= len(query_sequence):
                break
            read_base = query_sequence[seq_pos]
            mismatch = {
                'position': seq_pos,
                'ref_base': ref_base,
                'read_base': read_base,
                'strand': '-' if read.is_reverse else '+'
            }
            mismatches.append(mismatch)
            seq_pos += 1
            ref_pos += 1

    return mismatches




################################

########  #### ##       ########
##         ##  ##       ##      
##         ##  ##       ##      
########   ##  ##       ######  
##         ##  ##       ##      
##         ##  ##       ##      
##        #### ######## ########

################################

def process_uploaded_file(content, filename):
    file_ext = filename.split('.')[-1].lower()
    if file_ext not in ['bam', 'sam', 'fasta', 'fq', 'fna', 'fa', 'fastq']:
        raise ValueError(f"Unsupported file format: {file_ext}")

    # Map file extensions to Bio.SeqIO format names
    if file_ext in ['fa', 'fna', 'fasta']:
        file_format = 'fasta'
    elif file_ext in ['fq', 'fastq']:
        file_format = 'fastq'
    elif file_ext in ['bam', 'sam']:
        file_format = file_ext
    else:
        raise ValueError(f"Unsupported file extension: {file_ext}")

    
    content_type, content_string = content.split(',') # content_type is necessary, VScode just doesn't see its reference

    # Decode the base64-encoded string
    decoded = base64.b64decode(content_string)

    unique_filename = f"{uuid.uuid4()}_{filename}"
    temp_file_path = os.path.join(CUSTOM_TEMP_DIR, unique_filename)

    # Write the decoded content to a temporary file
    with open(temp_file_path, "wb") as fp:
        fp.write(decoded)

    return temp_file_path, file_format

def export_selected_reads(
    read_lengths, selected_lengths, selected_cg_content, output_file_path,
    original_file_path, file_format, selected_nm=None, filter_ct=None, exact_ct_changes=None
):
    if file_format in ['bam', 'sam']:
        with pysam.AlignmentFile(original_file_path, "rb" if file_format == 'bam' else "r", check_sq=False) as bamfile:
            header = bamfile.header.to_dict()

        # Retain all references in the header to avoid mismatch errors
        with pysam.AlignmentFile(output_file_path, "wb" if file_format == 'bam' else "w", header=header) as outputfile:
            for length, cg, nm, seq, read, mapq in read_lengths:
                deamination_pattern = get_deamination_pattern(read)
                ct_changes = count_ct_changes(read)

                # set filters based on C>T changes
                if filter_ct is not None:
                    has_ct = 'C>T' in deamination_pattern or 'G>A' in deamination_pattern
                    if filter_ct and not has_ct:
                        continue
                    elif not filter_ct and has_ct:
                        continue

                # set NM filter
                if selected_nm is not None and selected_nm != 'all' and nm != selected_nm:
                    continue

                # set exact C>T changes filter
                if exact_ct_changes is not None and ct_changes != exact_ct_changes:
                    continue

                # set selected lengths filter
                if selected_lengths and length not in selected_lengths:
                    continue  # Skip this read

                # set selected CG content filter (if applicable)
                if selected_cg_content and cg not in selected_cg_content:
                    continue  # Skip this read

                outputfile.write(read)
    else:
        with open(output_file_path, "w") as outputfile:
            for length, cg, _, seq, record in read_lengths:
                # set selected lengths filter
                if selected_lengths and length not in selected_lengths:
                    continue  # Skip this record

                # set selected CG content filter (if applicable)
                if selected_cg_content and cg not in selected_cg_content:
                    continue  # Skip this record

                SeqIO.write(record, outputfile, file_format)



# cached results of file-processing to avoid reloading the data, but not strictly necessary
@lru_cache(maxsize=100)
def load_and_process_file(file_content, filename, selected_nm, filter_ct, exclusively_ct, ct_count_value, ct_checklist_tuple, mapq_range):
    temp_file_path, file_format = process_uploaded_file(file_content, filename)
    deduped_file_path = os.path.join(
        CUSTOM_TEMP_DIR, f"deduped_{uuid.uuid4()}_{filename}"
    )
    read_lengths = remove_duplicates_with_strand_check(
        temp_file_path, deduped_file_path, file_format
    )

    if file_format in ['bam', 'sam']:
        # set MAPQ filter only for BAM/SAM files
        min_mapq, max_mapq = mapq_range
        read_lengths = [
            r for r in read_lengths
            if r[5] is not None and min_mapq <= r[5] <= max_mapq
        ]
        
        if filter_ct is not None:
            # Filter reads that have at least one C>T or G>A change
            read_lengths = [r for r in read_lengths if ('C>T' in get_deamination_pattern(r[4]) or 'G>A' in get_deamination_pattern(r[4])) == filter_ct]

        if exclusively_ct:
            # Keep reads where all mismatches are C>T or G>A
            read_lengths = [r for r in read_lengths if has_only_ct_mismatches(r[4])]

        # First, adjust NM values if subtract_ct is selected
        if 'subtract_ct' in ct_checklist_tuple:
            read_lengths = adjust_nm_for_ct_changes(read_lengths, subtract_ct=True, ct_count_value=ct_count_value)

        # Then, the C>T count filter if specified
        if ct_count_value != 'any':
            read_lengths = [r for r in read_lengths if count_ct_changes(r[4]) == int(ct_count_value)]

        # Finally, the NM filter on the adjusted NM values
        if selected_nm != 'all':
            read_lengths = [r for r in read_lengths if r[2] == selected_nm]

    return read_lengths, file_format, deduped_file_path


def remove_duplicates_with_strand_check(input_file, output_file, file_format, apply_dust=False):
    if file_format in ['bam', 'sam']:
        bamfile = pysam.AlignmentFile(input_file, "rb" if file_format == 'bam' else "r", check_sq=False)
        outputfile = pysam.AlignmentFile(output_file, "wb" if file_format == 'bam' else 'w', header=bamfile.header)
    else:
        outputfile = open(output_file, "w")
    
    unique_reads = set()
    read_lengths = []
    
    if file_format in ['bam', 'sam']:
        for read in bamfile:
            sequence = read.query_sequence
            if apply_dust:
                masked_seq = dust_module.dust_mask(sequence)
                if set(masked_seq) == {'N'}:  # DUST masks low complexity to 'N'
                    continue
           
            reverse_complement = str(Seq(sequence).reverse_complement())
            previous_length = len(unique_reads)
            unique_reads.add(sequence)
            unique_reads.add(reverse_complement)
            if len(unique_reads) > previous_length:
                outputfile.write(read)
                mapq = read.mapping_quality
                nm_tag = read.get_tag('NM') if read.has_tag('NM') else None
                read_lengths.append((len(sequence), calculate_cg_content(sequence), nm_tag, sequence, read, mapq))
            else:
                continue
    else:
        for record in SeqIO.parse(input_file, file_format):
            sequence = str(record.seq)
            if apply_dust:
                masked_seq = dust_module.dust_mask(sequence)
                if set(masked_seq) == {'N'}:  # DUST masks low complexity to 'N'
                    continue
           
            reverse_complement = str(Seq(sequence).reverse_complement())
            previous_length = len(unique_reads)
            unique_reads.add(sequence)
            unique_reads.add(reverse_complement)
            if len(unique_reads) > previous_length:
                SeqIO.write(record, outputfile, file_format)
                read_lengths.append((len(sequence), calculate_cg_content(sequence), None, sequence, record, None))
            else:
                continue
    
    if file_format in ['bam', 'sam']:
        bamfile.close()
    outputfile.close()
    return read_lengths




#####################################################################

########  ##        #######  ######## ######## #### ##    ##  ######  
##     ## ##       ##     ##    ##       ##     ##  ###   ## ##    ## 
##     ## ##       ##     ##    ##       ##     ##  ####  ## ##       
########  ##       ##     ##    ##       ##     ##  ## ## ## ##   ####
##        ##       ##     ##    ##       ##     ##  ##  #### ##    ## 
##        ##       ##     ##    ##       ##     ##  ##   ### ##    ## 
##        ########  #######     ##       ##    #### ##    ##  ######

#####################################################################

def get_mismatches(read):
    mismatches = []

    if read.is_unmapped:
        return mismatches

    read_seq = read.query_sequence.upper()
    qual = read.query_qualities
    if qual is None:
        qual = [0] * len(read_seq)  # Assign zero quality if not available

    ref_seq = read.get_reference_sequence().upper()

    for q_pos, r_pos in read.get_aligned_pairs(matches_only=True):
        read_base = read_seq[q_pos]
        ref_base = ref_seq[r_pos - read.reference_start]
        base_qual = qual[q_pos]
        mapping_qual = read.mapping_quality
        is_reverse = read.is_reverse

        if read_base != ref_base:
            mismatch_info = {
                'read_pos': q_pos,
                'ref_pos': r_pos,
                'read_base': read_base,
                'ref_base': ref_base,
                'base_quality': base_qual,
                'mapping_quality': mapping_qual,
                'is_reverse': is_reverse
            }
            mismatches.append(mismatch_info)

    return mismatches


def calculate_alignment_stats(file_path, file_format):
    if file_format in ['bam', 'sam']:
        bamfile = pysam.AlignmentFile(file_path, "rb" if file_format == 'bam' else 'r', check_sq=False)
        
        total_reads = 0
        mapped_reads = 0
        duplicate_reads = 0
        total_mismatches = 0
        mismatch_counts = {}
        mismatch_details = []

        for read in bamfile:
            total_reads += 1
            if not read.is_unmapped:
                mapped_reads += 1
                nm_tag = read.get_tag('NM') if read.has_tag('NM') else 0
                total_mismatches += nm_tag

                mismatches = get_mismatches(read)
                mismatch_details.extend(mismatches)
                for mismatch in mismatches:
                    key = f"{mismatch['ref_base']}>{mismatch['read_base']}"
                    mismatch_counts[key] = mismatch_counts.get(key, 0) + 1

            if read.is_duplicate:
                duplicate_reads += 1
        
        bamfile.close()

        mapped_percentage = (mapped_reads / total_reads) * 100 if total_reads > 0 else 0
        duplicate_percentage = (duplicate_reads / total_reads) * 100 if total_reads > 0 else 0

        stats = {
            'Total Reads': total_reads,
            'Mapped Reads': mapped_reads,
            'Mapped Percentage': mapped_percentage,
            'Duplicate Reads': duplicate_reads,
            'Duplicate Percentage': duplicate_percentage,
            'Total Mismatches': total_mismatches,
            'Mismatch Counts': mismatch_counts,
            'Mismatch Details': mismatch_details
        }
    else:
        total_reads = sum(1 for _ in SeqIO.parse(file_path, file_format))
        stats = {
            'Total Reads': total_reads,
            'Mapped Reads': 'N/A',
            'Mapped Percentage': 'N/A',
            'Duplicate Reads': 'N/A',
            'Duplicate Percentage': 'N/A',
            'Total Mismatches': 'N/A',
            'Mismatch Counts': 'N/A',
            'Mismatch Details': 'N/A'
        }
    return stats

def filter_mismatches(mismatches, min_base_quality=20, min_mapping_quality=0):
    filtered_mismatches = []
    for mismatch in mismatches:
        base_quality = mismatch.get('base_quality')
        mapping_quality = mismatch.get('mapping_quality')
        
        if base_quality is None or mapping_quality is None:
            continue  # Skip mismatches with missing quality values
        
        if base_quality >= min_base_quality and mapping_quality >= min_mapping_quality:
            filtered_mismatches.append(mismatch)
    return filtered_mismatches

def has_only_ct_mismatches(read):
    if read.is_unmapped:
        return False  # Unmapped reads cannot be considered

    mismatches = []
    for q_pos, r_pos, ref_base in read.get_aligned_pairs(with_seq=True):
        if q_pos is not None and ref_base is not None:
            read_base = read.query_sequence[q_pos].upper()
            ref_base = ref_base.upper()
            if read_base != ref_base:
                mismatches.append({'ref_base': ref_base, 'read_base': read_base})

    if not mismatches:
        return False  

    for mismatch in mismatches:
        ref_base = mismatch['ref_base']
        read_base = mismatch['read_base']
        if read.is_reverse:
            # Reverse strand: G>A corresponds to C>T on reverse strand
            if not (ref_base == 'G' and read_base == 'A'):
                return False
        else:
            # Forward strand
            if not (ref_base == 'C' and read_base == 'T'):
                return False

    return True  # All mismatches are C>T or G>A




def categorize_mismatches(mismatches): # have to check whether this makes sense
    categories = {
        'Transitions': 0,
        'Transversions': 0,
        'C>T/G>A Damage': 0,
        'Other': 0
    }
    for mismatch in mismatches:
        ref_base = mismatch['ref_base']
        read_base = mismatch['read_base']
        mutation = f"{ref_base}>{read_base}"
        if mutation in ['A>G', 'G>A', 'C>T', 'T>C']:
            categories['Transitions'] += 1
            if mutation in ['C>T', 'G>A']:
                categories['C>T/G>A Damage'] += 1
        elif mutation in ['A>C', 'A>T', 'C>A', 'C>G', 'G>C', 'G>T', 'T>A', 'T>G']:
            categories['Transversions'] += 1
        else:
            categories['Other'] += 1
    return categories


def get_mismatch_pattern(read):
    if not isinstance(read, pysam.AlignedSegment):
        return {}
    
    if read.is_unmapped:
        return {}
    
    mismatch_counts = {
        'C>T': 0, 'G>A': 0,
        'A>C': 0, 'A>G': 0, 'A>T': 0,
        'T>A': 0, 'T>C': 0, 'T>G': 0,
        'C>A': 0, 'C>G': 0,
        'G>C': 0, 'G>T': 0
    }

    read_seq = read.query_sequence.upper()
    ref_seq = read.get_reference_sequence().upper()

    if read.is_reverse:
        ref_seq = str(Seq(ref_seq).reverse_complement())
        read_seq = str(Seq(read_seq).reverse_complement())

    for q_pos, r_pos in read.get_aligned_pairs(matches_only=True):
        ref_base = ref_seq[r_pos - read.reference_start]
        read_base = read_seq[q_pos]
        mismatch_key = f"{ref_base}>{read_base}"
        if mismatch_key in mismatch_counts:
            mismatch_counts[mismatch_key] += 1

    return mismatch_counts


def calculate_cg_content(sequence):
    if len(sequence) == 0:
        return 0
    cg_count = sequence.count('C') + sequence.count('G')
    return cg_count / len(sequence)

def create_read_length_histogram(bin_centers, hist, selected_file, read_lengths, bins):
    read_length_fig = go.Figure()

    # Histogram bars
    read_length_fig.add_trace(go.Bar(
        x=bin_centers - 0.5,
        y=hist,
        name=f'{selected_file} Read Count',
        opacity=1,
        hovertemplate='%{x}: %{y}<extra></extra>',
        marker_color=colors['line'],
    ))

    
    mean_cg_content = calculate_mean_cg_content(read_lengths, bins)

    # CG content markers (including zeros)
    read_length_fig.add_trace(go.Scatter(
        x=[x - 0.5 for x in bin_centers],
        y=mean_cg_content,
        name=f'{selected_file} Average CG Content',
        yaxis='y2',
        mode='markers',
        marker=dict(color=colors['highlight'], size=8),
    ))

    # Lines between neighboring bins
    for i in range(len(bin_centers) - 1):
        if mean_cg_content[i] > 0 and mean_cg_content[i + 1] > 0:
            read_length_fig.add_trace(go.Scatter(
                x=[bin_centers[i] - 0.5, bin_centers[i + 1] - 0.5],
                y=[mean_cg_content[i], mean_cg_content[i + 1]],
                mode='lines',
                line=dict(color=colors['highlight'], width=2),
                yaxis='y2',
                showlegend=False
            ))

    read_length_fig.update_layout(
        title=dict(text='Read Length Distribution with Average CG Content', font=dict(color=colors['muted'])),
        uirevision='read_length',
        xaxis=dict(title='Read Length', color=colors['muted'], tickfont=dict(color=colors['muted'])),
        yaxis=dict(
            title=dict(text='Read Count', font=dict(color='#1CA8FF')),
            showline=True,
            linecolor='#1CA8FF',
            tickfont=dict(color=colors['muted']),
            gridcolor=colors['grid'],
            gridwidth=1
        ),
        yaxis2=dict(
            title=dict(text='Average CG Content', font=dict(color='#FF4081')),
            overlaying='y',
            side='right',
            range=[0, 1],
            fixedrange=True,
            showline=True,
            linecolor='#FF4081',
            tickfont=dict(color=colors['muted']),
            gridcolor=colors['grid'],
            gridwidth=1
        ),
        legend=dict(x=0.85, y=1.15, font=dict(color=colors['muted'])),
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        dragmode='select',
        selectdirection='h',
        newselection_line_width=2,
        newselection_line_color=colors['secondary'],
        newselection_line_dash='dashdot'
    )

    return read_length_fig

def create_cg_content_histogram(cg_contents_for_length, selected_lengths):
    cg_content_fig = go.Figure()

    if cg_contents_for_length:
        cg_content_fig.add_trace(go.Histogram(
            x=cg_contents_for_length,
            nbinsx=20,
            marker_color=colors['marker']
        ))
        selected_lengths_str = ', '.join(map(str, selected_lengths)) if selected_lengths else 'All'
        cg_content_fig.update_layout(
            title=f'CG Content Distribution for Selected Read Length(s): {selected_lengths_str}',
            xaxis=dict(title='CG Content', color=colors['muted']),
            yaxis=dict(title='Frequency', color=colors['muted']),
            paper_bgcolor=colors['plot_bg'],
            plot_bgcolor=colors['plot_bg'],
            font=dict(color=colors['muted']),
            dragmode='select',
            newselection_line_width=2,
            newselection_line_color=colors['secondary'],
            newselection_line_dash='dashdot'
        )
    else:
        cg_content_fig.update_layout(
            title='No data available for selected read length(s)',
            paper_bgcolor=colors['plot_bg'],
            plot_bgcolor=colors['plot_bg'],
            font=dict(color=colors['muted'])
        )

    return cg_content_fig

def create_mismatch_type_bar_chart(mismatch_counts):
    mutations = list(mismatch_counts.keys())
    counts = list(mismatch_counts.values())
    
    fig = go.Figure(data=[go.Bar(
        x=mutations,
        y=counts,
        marker_color=colors['highlight']
    )])
    
    fig.update_layout(
        title='Mismatch Types Distribution',
        xaxis=dict(title='Mismatch Type', color=colors['muted']),
        yaxis=dict(title='Count', color=colors['muted']),
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        font=dict(color=colors['muted'])
    )
    return fig

def create_damage_pattern_plot(mismatches):
    positions = [mismatch['read_pos'] for mismatch in mismatches]
    mutations = [f"{mismatch['ref_base']}>{mismatch['read_base']}" for mismatch in mismatches]
    
    df = pd.DataFrame({'Position': positions, 'Mutation': mutations})
    damage_df = df[df['Mutation'].isin(['C>T', 'G>A'])]
    
    fig = px.histogram(damage_df, x='Position', nbins=50, color='Mutation',
                       title='Damage Patterns Along Read Positions')
    fig.update_layout(
        xaxis_title='Read Position',
        yaxis_title='Count',
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        font=dict(color=colors['muted'])
    )
    return fig


def create_detailed_mismatch_pie_chart(stats):
    mismatch_counts = stats.get('Mismatch Counts')
    
    # for everything that has no mismatches 
    if mismatch_counts == 'N/A' or not mismatch_counts:
        return go.Figure().update_layout(
            title_text="No mismatch data available",
            paper_bgcolor=colors['plot_bg'],
            font=dict(color=colors['muted'])
        )
    sorted_keys = sorted(mismatch_counts.keys(), key=lambda x: (x[0], x[2]))  # Sort by base
    
    # The rest of the function (plotting) is way too complex for a pretty visualization
    # Define color ramps for each base using Plotly colors
    color_map = {
        'A': plotly.colors.sequential.Blues,
        'C': plotly.colors.sequential.Oranges,
        'G': plotly.colors.sequential.Greens,
        'T': plotly.colors.sequential.Reds_r,
        'N': plotly.colors.sequential.solar # otherwise not really visible
    }
    
    # Assign colors to each mismatch type
    colors_list = []
    brighter_colors = []
    for key in sorted_keys:
        base = key[0]
        ramp_length = len(color_map[base])
        color_index = sorted_keys.index(key) % ramp_length
        color = color_map[base][color_index]
        colors_list.append(color)
        
        # Convert color to RGB tuple
        if color.startswith('#'):
            rgb_color = plotly.colors.hex_to_rgb(color)
        elif color.startswith('rgb'):
            rgb_color = tuple(map(int, color[4:-1].split(',')))
        else:
            raise ValueError(f"Unsupported color format: {color}")
        
        # Find intermediate color between white and the original color
        brightened_color = plotly.colors.find_intermediate_color((255, 255, 255), rgb_color, 0.6)
        brightened_color_hex = rgb_to_hex(brightened_color)
        brighter_colors.append(brightened_color_hex)
    
    labels = sorted_keys
    values = [mismatch_counts[key] for key in sorted_keys]
    
    # Create custom texttemplate with brighter label colors
    texttemplate = [
        f'<span style="color:{brighter_colors[i]}">%{{label}}: %{{percent}}</span>' 
        for i in range(len(labels))
    ]
    
    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        hole=.4,
        marker=dict(colors=colors_list),
        textinfo='label+percent',
        textposition='outside',
        sort=False,
        texttemplate=texttemplate
    )])
    
    fig.update_layout(
        title_text="Detailed SNP Distribution",
        title_x=0.5,
        title_y=0.9,
        title_font=dict(size=24, color=colors['muted']),
        annotations=[dict(text='SNPs', x=0.5, y=0.5, font_size=20, showarrow=False, font=dict(color=colors['muted']))],
        paper_bgcolor=colors['stats_tab_color'],
        plot_bgcolor=colors['stats_tab_color'],
        font=dict(color=colors['muted']),
        margin=dict(t=160, b=20, l=20, r=20),
        showlegend=False
    )
        
    return fig


def calculate_mean_cg_content(read_lengths, bins):
    mean_cg_contents = []
    for i in range(len(bins) - 1):
        bin_start = bins[i]
        bin_end = bins[i + 1]
        cg_values = [cg for read_length, cg, _, _, _, _ in read_lengths if bin_start <= read_length < bin_end]
        if cg_values:
            mean_cg = np.mean(cg_values)
        else:
            mean_cg = 0
        mean_cg_contents.append(mean_cg)
    return mean_cg_contents

def calculate_confidence_intervals(counts, total_bases, confidence_level=0.95):
    """
    Calculate binomial confidence intervals for mismatch frequencies.
    
    Parameters:
    -----------
    counts : array-like
        Number of mismatches at each position
    total_bases : array-like
        Total number of bases at each position
    confidence_level : float
        Confidence level for interval calculation (default: 0.95)
    
    Returns:
    --------
    tuple : (lower_bounds, upper_bounds)
        Arrays containing lower and upper confidence bounds
    """
    from scipy import stats
    
    frequencies = np.array(counts) / np.maximum(total_bases, 1)
    z_score = stats.norm.ppf((1 + confidence_level) / 2)
    
    # Wilson score interval calculation
    denominator = 1 + z_score**2/np.array(total_bases)
    center = (frequencies + z_score**2/(2*np.array(total_bases))) / denominator
    
    spread = z_score * np.sqrt(
        (frequencies * (1 - frequencies) + z_score**2/(4*np.array(total_bases))) / np.array(total_bases)
    ) / denominator
    
    lower_bounds = np.maximum(0, center - spread)
    upper_bounds = np.minimum(1, center + spread)
    
    return lower_bounds, upper_bounds

def calculate_mismatch_frequency(read_lengths, file_format, sequencing_type):
    if file_format not in ['bam', 'sam']:
        return {}
        
    max_distance = 50  # Standard distance from read ends to analyze
    counts = {
        '5_prime': {
            'C>T_counts': [0]*max_distance, 
            'G>A_counts': [0]*max_distance, 
            'other_counts': [0]*max_distance,
            'C>T': [0]*max_distance, 
            'G>A': [0]*max_distance, 
            'other': [0]*max_distance,
            'total_bases': [0]*max_distance
        },
        '3_prime': {
            'C>T_counts': [0]*max_distance, 
            'G>A_counts': [0]*max_distance, 
            'other_counts': [0]*max_distance,
            'C>T': [0]*max_distance, 
            'G>A': [0]*max_distance, 
            'other': [0]*max_distance,
            'total_bases': [0]*max_distance
        }
    }

    for read_tuple in read_lengths:
        read = read_tuple[4]  # Get read object from tuple
        if read.is_unmapped:
            continue
            
        sequence = read.query_sequence
        if sequence is None:
            continue
            
        seq_length = len(sequence)
        if seq_length == 0:
            continue

        # Get aligned positions and mismatches
        aligned_pairs = read.get_aligned_pairs(with_seq=True)
        if not aligned_pairs:
            continue

        for query_pos, ref_pos, ref_base in aligned_pairs:
            if query_pos is None or ref_pos is None or ref_base is None:
                continue  # Skip insertions/deletions
                
            if query_pos >= seq_length:
                continue

            read_base = sequence[query_pos].upper()
            ref_base = ref_base.upper()
            
            # Calculate positions from both ends
            dist_5p = query_pos
            dist_3p = seq_length - query_pos - 1
            
            if dist_5p < max_distance:
                counts['5_prime']['total_bases'][dist_5p] += 1
                if read_base != ref_base:
                    if read.is_reverse:
                        if ref_base == 'G' and read_base == 'A':
                            counts['5_prime']['G>A_counts'][dist_5p] += 1
                        else:
                            counts['5_prime']['other_counts'][dist_5p] += 1
                    else:
                        if ref_base == 'C' and read_base == 'T':
                            counts['5_prime']['C>T_counts'][dist_5p] += 1
                        else:
                            counts['5_prime']['other_counts'][dist_5p] += 1

            if dist_3p < max_distance:
                counts['3_prime']['total_bases'][dist_3p] += 1
                if read_base != ref_base:
                    if read.is_reverse:
                        if ref_base == 'C' and read_base == 'T':
                            counts['3_prime']['C>T_counts'][dist_3p] += 1
                        else:
                            counts['3_prime']['other_counts'][dist_3p] += 1
                    else:
                        if ref_base == 'G' and read_base == 'A':
                            counts['3_prime']['G>A_counts'][dist_3p] += 1
                        else:
                            counts['3_prime']['other_counts'][dist_3p] += 1

    # Calculate frequencies
    frequencies = {}
    
    if sequencing_type == 'single':
        frequencies = {
            'C>T': [], 'G>A': [], 'other': [],
            'C>T_counts': [], 'G>A_counts': [], 'other_counts': [],
            'total_bases': []
        }
        
        for i in range(max_distance):
            total_5p = counts['5_prime']['total_bases'][i]
            total_3p = counts['3_prime']['total_bases'][i]
            
            frequencies['total_bases'].append(total_5p + total_3p)
            
            for mut_type in ['C>T', 'G>A', 'other']:
                count_5p = counts['5_prime'][f'{mut_type}_counts'][i]
                count_3p = counts['3_prime'][f'{mut_type}_counts'][i]
                
                frequencies[f'{mut_type}_counts'].append(count_5p + count_3p)
                if frequencies['total_bases'][-1] > 0:
                    frequencies[mut_type].append(
                        (count_5p + count_3p) / frequencies['total_bases'][-1]
                    )
                else:
                    frequencies[mut_type].append(0)
    else:
        frequencies = {
            '5_prime': {},
            '3_prime': {}
        }
        
        for end in ['5_prime', '3_prime']:
            frequencies[end] = {
                'C>T': [], 'G>A': [], 'other': [],
                'C>T_counts': [], 'G>A_counts': [], 'other_counts': [],
                'total_bases': counts[end]['total_bases']
            }
            
            for i in range(max_distance):
                total = counts[end]['total_bases'][i]
                
                for mut_type in ['C>T', 'G>A', 'other']:
                    count = counts[end][f'{mut_type}_counts'][i]
                    
                    frequencies[end][f'{mut_type}_counts'].append(count)
                    if total > 0:
                        frequencies[end][mut_type].append(count / total)
                    else:
                        frequencies[end][mut_type].append(0)

    return frequencies


def create_mismatch_frequency_plot(frequencies, sequencing_type):
    if not frequencies:
        fig = go.Figure()
        fig.update_layout(
            title_text="Mismatch Frequency Data Not Available",
            paper_bgcolor=colors['plot_bg'],
            font=dict(color=colors['muted'])
        )
        return fig

    fig = go.Figure()
    
    if sequencing_type == 'single':
        x_values = list(range(len(frequencies['C>T'])))
        
        for mut_type, color, symbol in [
            ('C>T', colors['line'], 'circle'),
            ('G>A', colors['highlight'], 'square'),
            ('other', colors['highlight2'], 'diamond')
        ]:
            # Calculate confidence intervals
            lower_bounds, upper_bounds = calculate_confidence_intervals(
                frequencies[f'{mut_type}_counts'],
                frequencies['total_bases']
            )
            
            # Main line
            fig.add_trace(go.Scatter(
                x=x_values,
                y=frequencies[mut_type],
                name=mut_type,
                mode='lines+markers',
                line=dict(color=color),
                marker=dict(symbol=symbol)
            ))
            
            # Confidence interval
            fig.add_trace(go.Scatter(
                x=x_values + x_values[::-1],
                y=np.concatenate([upper_bounds, lower_bounds[::-1]]),
                fill='toself',
                fillcolor=f'rgba{tuple(list(plotly.colors.hex_to_rgb(color)) + [0.2])}',
                line=dict(width=0),
                showlegend=False,
                name=f'{mut_type} 95% CI'
            ))
            
            # Add sample size indicators
            sizes = frequencies['total_bases']
            normalized_sizes = np.interp(sizes, (min(sizes), max(sizes)), (2, 10))
            
            fig.add_trace(go.Scatter(
                x=x_values,
                y=frequencies[mut_type],
                mode='markers',
                marker=dict(
                    size=normalized_sizes,
                    color=color,
                    symbol=symbol,
                    line=dict(width=1, color='white')
                ),
                showlegend=False
            ))
    else:
        for end, title in [('5_prime', "5' End"), ('3_prime', "3' End")]:
            for mut_type, color, symbol in [
                ('C>T', colors['line'], 'circle'),
                ('G>A', colors['highlight'], 'square'),
                ('other', colors['highlight2'], 'diamond')
            ]:
                name = f"{mut_type} ({title})"
                x_values = list(range(len(frequencies[end][mut_type])))
                
                # Calculate confidence intervals
                lower_bounds, upper_bounds = calculate_confidence_intervals(
                    frequencies[end][f'{mut_type}_counts'],
                    frequencies[end]['total_bases']
                )
                
                # Main line
                fig.add_trace(go.Scatter(
                    x=x_values,
                    y=frequencies[end][mut_type],
                    name=name,
                    mode='lines+markers',
                    line=dict(color=color),
                    marker=dict(symbol=symbol)
                ))
                
                # Confidence interval
                fig.add_trace(go.Scatter(
                    x=x_values + x_values[::-1],
                    y=np.concatenate([upper_bounds, lower_bounds[::-1]]),
                    fill='toself',
                    fillcolor=f'rgba{tuple(list(plotly.colors.hex_to_rgb(color)) + [0.2])}',
                    line=dict(width=0),
                    showlegend=False,
                    name=f'{name} 95% CI'
                ))
                
                # Add sample size indicators
                sizes = frequencies[end]['total_bases']
                normalized_sizes = np.interp(sizes, (min(sizes), max(sizes)), (2, 10))
                
                fig.add_trace(go.Scatter(
                    x=x_values,
                    y=frequencies[end][mut_type],
                    mode='markers',
                    marker=dict(
                        size=normalized_sizes,
                        color=color,
                        symbol=symbol,
                        line=dict(width=1, color='white')
                    ),
                    showlegend=False
                ))

    fig.update_layout(
        title='Mismatch Frequency vs Distance from Read End with 95% Confidence Intervals',
        xaxis_title='Position',
        yaxis_title='Frequency',
        template='plotly_dark',
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        font=dict(color=colors['muted']),
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99,
            bgcolor='rgba(0,0,0,0.5)'
        )
    )

    # Add hover text with sample size information
    for trace in fig.data:
        if trace.mode == 'lines+markers':
            pos = trace.x
            freq = trace.y
            if sequencing_type == 'single':
                n = frequencies['total_bases']
            else:
                end = '5_prime' if "5' End" in trace.name else '3_prime'
                n = frequencies[end]['total_bases']
            
            hover_text = [
                f'Position: {p}<br>'
                f'Frequency: {f:.4f}<br>'
                f'Sample size: {s:,}<br>'
                f'95% CI: [{l:.4f}, {u:.4f}]'
                for p, f, s, l, u in zip(
                    pos, freq, n,
                    lower_bounds, upper_bounds
                )
            ]
            trace.hovertemplate = '%{text}<extra></extra>'
            trace.text = hover_text

    return fig



def display_mapdamage_results(output_dir):
    try:
        # Parse misincorporation data
        misincorporation_df = parse_misincorporation(os.path.join(output_dir, 'misincorporation.txt'))
        # Create misincorporation plot
        misincorporation_fig = create_misincorporation_plot(misincorporation_df)

        # Parse DNA composition data
        dnacomp_df = parse_dnacomp(os.path.join(output_dir, 'dnacomp.txt'))
        dnacomp_fig = create_dnacomp_plot(dnacomp_df)

        # Display certainty metrics if available
        metrics_content = []
        mcmc_path = os.path.join(output_dir, 'Stats_out_MCMC.csv')
        if os.path.exists(mcmc_path):
            mcmc_df = pd.read_csv(mcmc_path)
            metrics_content.append(display_certainty_metrics(mcmc_df))

        return html.Div([
            html.H4("Misincorporation Plot"),
            dcc.Graph(figure=misincorporation_fig),
            html.H4("DNA Composition Plot"),
            dcc.Graph(figure=dnacomp_fig),
            *metrics_content
        ])
    except Exception as e:
        return html.Div(f"Error processing mapDamage2 results: {e}", className="text-danger")


def get_damage_patterns(file_path, file_format, filter_ct=None):
    five_prime_end = []
    three_prime_end = []

    if file_format in ['bam', 'sam']:
        bamfile = pysam.AlignmentFile(file_path, "rb" if file_format == 'bam' else "r", check_sq=False)
        for read in bamfile:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            seq = read.query_sequence
            deamination_pattern = get_deamination_pattern(read)
            has_ct = 'C>T' in deamination_pattern or 'G>A' in deamination_pattern
            if filter_ct is not None:
                if filter_ct and not has_ct:
                    continue
                elif not filter_ct and has_ct:
                    continue

            if read.is_reverse:
                five_prime_end.append(seq[-1])
                three_prime_end.append(seq[0])
            else:
                five_prime_end.append(seq[0])
                three_prime_end.append(seq[-1])
        bamfile.close()
    else:
        for record in SeqIO.parse(file_path, file_format):
            seq = str(record.seq)
            five_prime_end.append(seq[0])
            three_prime_end.append(seq[-1])

    return five_prime_end, three_prime_end

def get_deamination_pattern(read):
    pattern = []
    read_seq = read.query_sequence.upper()
    ref_seq = read.get_reference_sequence().upper()

    if read.is_reverse:
        ref_seq = str(Seq(ref_seq).reverse_complement())
        read_seq = str(Seq(read_seq).reverse_complement())

    for q_pos, r_pos in read.get_aligned_pairs(matches_only=True):
        ref_base = ref_seq[r_pos - read.reference_start]
        read_base = read_seq[q_pos]
        if ref_base.upper() == 'C' and read_base.upper() == 'T':
            pattern.append('C>T')
        elif read.is_reverse and ref_base.upper() == 'G' and read_base.upper() == 'A':
            pattern.append('G>A')

    return pattern

def calculate_overall_cg_content(read_lengths):
    cg_contents = [cg for _, cg, _, _, _, _ in read_lengths]
    return cg_contents

def count_ct_changes(read): # function is older than calculate_mismatch_frequency, a bit redundant
    if isinstance(read, pysam.AlignedSegment):
        pattern = get_deamination_pattern(read)
        return pattern.count('C>T') + pattern.count('G>A')
    else:
        return 0  # For FASTA/FASTQ records, you can't determine C>T changes
        
def adjust_nm_for_ct_changes(read_lengths, subtract_ct, ct_count_value): # for filtering section, for subtracting all C>Ts from NM without removing them from the list
    adjusted_read_lengths = []
    for length, cg, nm, seq, read, mapq in read_lengths:
        ct_changes = count_ct_changes(read)
        
        # C>T count filter first if specified
        if ct_count_value != 'any' and ct_changes != int(ct_count_value):
            continue

        if subtract_ct and isinstance(read, pysam.AlignedSegment):
            if nm is not None:
                # Subtract C>T changes from NM
                adjusted_nm = max(0, nm - ct_changes)
            else:
                adjusted_nm = None
        else:
            adjusted_nm = nm

        adjusted_read_lengths.append((length, cg, adjusted_nm, seq, read, mapq))
    
    return adjusted_read_lengths



def get_read_length_data(read_lengths):
    lengths = [length for length, _, _, _, _, _ in read_lengths]
    length_counts = pd.Series(lengths).value_counts().sort_index()
    text_output = "Read Length Distribution:\n"
    for length, count in length_counts.items():
        text_output += f"Length {length}: {count} reads\n"
    return text_output

def get_overall_cg_content_data(read_lengths):
    cg_contents = [cg for _, cg, _, _, _, _ in read_lengths]
    cg_series = pd.Series(cg_contents)
    text_output = "Overall CG Content Statistics:\n"
    text_output += f"Mean CG Content: {cg_series.mean():.4f}\n"
    text_output += f"Median CG Content: {cg_series.median():.4f}\n"
    text_output += f"Standard Deviation: {cg_series.std():.4f}\n"
    return text_output

def get_read_length_CG_average(read_lengths):
    df = pd.DataFrame(read_lengths, columns=['length', 'cg', 'adjusted_nm', 'seq', 'read', 'mapq'])
    grouped = df.groupby('length')['cg'].mean()
    text_output = "Mean CG Content by Read Length:\n"
    for length, mean_cg in grouped.items():
        text_output += f"Length {length}: Mean CG Content: {mean_cg:.4f}\n"
    return text_output

def get_damage_patterns_data(five_prime_end, three_prime_end):
    five_prime_counts = pd.Series(five_prime_end).value_counts()
    three_prime_counts = pd.Series(three_prime_end).value_counts()
    text_output = "Damage Patterns at Read Ends:\n"
    text_output += "5' End Base Counts:\n"
    for base, count in five_prime_counts.items():
        text_output += f"Base {base}: {count} occurrences\n"
    text_output += "3' End Base Counts:\n"
    for base, count in three_prime_counts.items():
        text_output += f"Base {base}: {count} occurrences\n"
    return text_output

def get_mismatch_frequency_data(frequencies, sequencing_type):
    if not frequencies:
        return "Mismatch frequency data not available for this file format.\n"

    if sequencing_type == 'single':
        text_output = "Mismatch Frequency vs Distance from Read End:\n"
        max_distance = len(frequencies['C>T'])
        text_output += "Position\tC>T Frequency\tG>A Frequency\tOther Frequency\n"
        for i in range(max_distance):
            ct_freq = frequencies['C>T'][i]
            ga_freq = frequencies['G>A'][i]
            other_freq = frequencies['other'][i]
            text_output += f"{i}\t{ct_freq:.4f}\t{ga_freq:.4f}\t{other_freq:.4f}\n"
    elif sequencing_type == 'double':
        text_output = "Mismatch Frequency vs Distance from Read Ends:\n"
        text_output += "\n5' End:\n"
        max_distance_5p = len(frequencies['5_prime']['C>T'])
        text_output += "Position\tC>T Frequency\tG>A Frequency\tOther Frequency\n"
        for i in range(max_distance_5p):
            ct_freq = frequencies['5_prime']['C>T'][i]
            ga_freq = frequencies['5_prime']['G>A'][i]
            other_freq = frequencies['5_prime']['other'][i]
            text_output += f"{i}\t{ct_freq:.4f}\t{ga_freq:.4f}\t{other_freq:.4f}\n"

        text_output += "\n3' End:\n"
        max_distance_3p = len(frequencies['3_prime']['C>T'])
        text_output += "Position\tC>T Frequency\tG>A Frequency\tOther Frequency\n"
        for i in range(max_distance_3p):
            ct_freq = frequencies['3_prime']['C>T'][i]
            ga_freq = frequencies['3_prime']['G>A'][i]
            other_freq = frequencies['3_prime']['other'][i]
            text_output += f"{i}\t{ct_freq:.4f}\t{ga_freq:.4f}\t{other_freq:.4f}\n"
    else:
        text_output = "Unknown sequencing type.\n"

    return text_output


def get_alignment_stats_text(stats):
    text_output = "Alignment Statistics:\n"
    for key, value in stats.items():
        if key not in ['Mismatch Counts', 'Mismatch Details']:
            text_output += f"{key}: {value}\n"
    # mismatch counts added 
    mismatch_counts = stats.get('Mismatch Counts', {})
    if mismatch_counts and mismatch_counts != 'N/A':
        text_output += "\nMismatch Counts:\n"
        for mismatch, count in mismatch_counts.items():
            text_output += f"{mismatch}: {count}\n"
    # mismatch categories added, but idk if this is correct
    if stats.get('Mismatch Details') != 'N/A':
        categories = categorize_mismatches(stats['Mismatch Details'])
        text_output += "\nMismatch Categories:\n"
        for category, count in categories.items():
            text_output += f"{category}: {count}\n"
    return text_output


def create_mapq_histogram(mapq_scores):
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=mapq_scores,
        nbinsx=60,  
        marker_color=colors['marker']
    ))
    fig.update_layout(
        title='Mapping Quality (MAPQ) Score Distribution',
        xaxis=dict(title='MAPQ Score', color=colors['muted']),
        yaxis=dict(title='Frequency', color=colors['muted']),
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        font=dict(color=colors['muted']),
        dragmode='select',
        newselection_line_width=2,
        newselection_line_color=colors['secondary'],
        newselection_line_dash='dashdot'
    )
    return fig

def get_filter_options_text(selected_nm, ct_checklist, ct_count_value, mapq_range_tuple):
    text_output = "Current Filters Applied:\n"
    text_output += f"NM Filter: {selected_nm}\n"
    ct_filters = ', '.join(ct_checklist) if ct_checklist else 'None'
    text_output += f"C>T Checklist Filters: {ct_filters}\n"
    text_output += f"C>T Change Count: {ct_count_value}\n"
    text_output += f"MAPQ Range: {mapq_range_tuple[0]} - {mapq_range_tuple[1]}\n"
    return text_output

def create_damage_pattern_figure(five_prime_end, three_prime_end, selected_file):
    damage_df = pd.DataFrame({
        'Base': five_prime_end + three_prime_end,
        'End': ['5_prime'] * len(five_prime_end) + ['3_prime'] * len(three_prime_end)
    })
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=damage_df[damage_df['End'] == '5_prime']['Base'],
        name=f'{selected_file} 5\' end',
        opacity=0.7,
        marker_color='blue'
    ))
    fig.add_trace(go.Histogram(
        x=damage_df[damage_df['End'] == '3_prime']['Base'],
        name=f'{selected_file} 3\' end',
        opacity=0.7,
        marker_color='red'
    ))
    fig.update_layout(
        barmode='group',
        title='Damage Patterns at Read Ends',
        xaxis=dict(title='Base', color=colors['muted']),
        yaxis=dict(title='Frequency', color=colors['muted']),
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        font=dict(color=colors['muted']),
        legend=dict(x=0.85, y=1.15, font=dict(color=colors['muted']))
    )
    return fig






########################################################

##          ###    ##    ##  #######  ##     ## ######## 
##         ## ##    ##  ##  ##     ## ##     ##    ##    
##        ##   ##    ####   ##     ## ##     ##    ##    
##       ##     ##    ##    ##     ## ##     ##    ##    
##       #########    ##    ##     ## ##     ##    ##    
##       ##     ##    ##    ##     ## ##     ##    ##    
######## ##     ##    ##     #######   #######     ##

#######################################################

settings_offcanvas = dbc.Offcanvas(
    [
        html.H4("Settings", className="mb-3", style={"color": colors['muted']}),

        # General Settings Section
        dbc.Card([
            dbc.CardHeader(
                dbc.Row([
                    dbc.Col(DashIconify(icon="lucide:sun", width=20, color=colors['on_surface']), width="auto"),
                    dbc.Col(html.Span("General Settings", style={"margin-left": "8px"})),
                ], align="center"),
                style={"background-color": colors['surface'], "color": colors['muted']}
            ),
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Theme", style={"color": colors['muted'], "margin-bottom": "4px"}),
                        dcc.Dropdown(
                            id="theme-dropdown",
                            options=[
                                {"label": "Dark Mode", "value": "dark"},
                                {"label": "Light Mode", "value": "light"}
                            ],
                            value="dark",
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
                ], className="g-2"),
            ]),
        ], className="mb-3", style={"background-color": colors['surface']}),

        # Data Processing Section
        dbc.Card([
            dbc.CardHeader(
                dbc.Row([
                    dbc.Col(DashIconify(icon="lucide:database", width=20, color=colors['on_surface']), width="auto"),
                    dbc.Col(html.Span("Data Processing", style={"margin-left": "8px"})),
                ], align="center"),
                style={"background-color": colors['surface'], "color": colors['muted']}
            ),
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Maximum Read Length", style={"color": colors['muted'], "margin-bottom": "4px"}),
                        dcc.Input(
                            id="max-read-length",
                            type="number",
                            value=1000,
                            style={
                                "background-color": colors['surface'],
                                "color": colors['on_surface'],
                                "border": f"1px solid {colors['muted']}",
                                "width": "100%",
                                "margin-bottom": "10px"
                            }
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label("Minimum Quality Score", style={"color": colors['muted'], "margin-bottom": "4px"}),
                        dcc.Input(
                            id="min-quality-score",
                            type="number",
                            value=20,
                            style={
                                "background-color": colors['surface'],
                                "color": colors['on_surface'],
                                "border": f"1px solid {colors['muted']}",
                                "width": "100%",
                                "margin-bottom": "10px"
                                }
                        ),
                    ], width=6),
                ], className="g-2"),
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Centrifuge DB Path", style={"color": colors['muted'], "margin-bottom": "4px"}),
                        dbc.InputGroup([
                            dbc.Input(
                                id="centrifuge-db-path",
                                type="text",
                                value=default_centrifuge_db_path,
                                placeholder="/path/to/centrifuge/db",
                                style={
                                    "background-color": colors['surface'],
                                    "color": colors['on_surface'],
                                    "border": f"1px solid {colors['muted']}",
                                    "width": "100%",
                                }
                            ),
                        ], className="mb-3"),
                    ], width=12),
                ], className="g-2"),
            ]),
        ], className="mb-3", style={"background-color": colors['surface']}),

        # Visualization Preferences Section
        dbc.Card([
            dbc.CardHeader(
                dbc.Row([
                    dbc.Col(DashIconify(icon="lucide:pie-chart", width=20, color=colors['on_surface']), width="auto"),
                    dbc.Col(html.Span("Visualization Preferences", style={"margin-left": "8px"})),
                ], align="center"),
                style={"background-color": colors['surface'], "color": colors['muted']}
            ),
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Color Themes", style={"color": colors['muted'], "margin-bottom": "4px"}),
                        dcc.Dropdown(
                            id="color-theme-dropdown",
                            options=[
                                {"label": "Default", "value": "default"},
                                {"label": "Cool Blues", "value": "cool"},
                                {"label": "Warm Reds", "value": "warm"}
                            ],
                            value="default",
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
                ], className="g-2"),
            ]),
        ], className="mb-3", style={"background-color": colors['surface']}),

        # Performance Settings Section
        dbc.Card([
            dbc.CardHeader(
                dbc.Row([
                    dbc.Col(DashIconify(icon="lucide:cpu", width=20, color=colors['on_surface']), width="auto"),
                    dbc.Col(html.Span("Performance Settings", style={"margin-left": "8px"})),
                ], align="center"),
                style={"background-color": colors['surface'], "color": colors['muted']}
            ),
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Max Threads", style={"color": colors['muted'], "margin-bottom": "4px"}),
                        dcc.Input(
                            id="max-threads",
                            type="number",
                            value=4,
                            style={
                                "background-color": colors['surface'],
                                "color": colors['on_surface'],
                                "border": f"1px solid {colors['muted']}",
                                "width": "100%",
                                "margin-bottom": "10px"
                            }
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label("Cache Size (MB)", style={"color": colors['muted'], "margin-bottom": "4px"}),
                        dcc.Input(
                            id="cache-size",
                            type="number",
                            value=100,
                            style={
                                "background-color": colors['surface'],
                                "color": colors['on_surface'],
                                "border": f"1px solid {colors['muted']}",
                                "width": "100%",
                                "margin-bottom": "10px"
                            }
                        ),
                    ], width=6),
                ], className="g-2"),
            ]),
        ], className="mb-3", style={"background-color": colors['surface']}),

        # Export Settings Section
        dbc.Card([
            dbc.CardHeader(
                dbc.Row([
                    dbc.Col(DashIconify(icon="lucide:file-text", width=20, color=colors['on_surface']), width="auto"),
                    dbc.Col(html.Span("Export Settings", style={"margin-left": "8px"})),
                ], align="center"),
                style={"background-color": colors['surface'], "color": colors['muted']}
            ),
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Default Export Format", style={"color": colors['muted'], "margin-bottom": "4px"}),
                        dcc.Dropdown(
                            id="export-format-dropdown",
                            options=[
                                {"label": "CSV", "value": "csv"},
                                {"label": "TSV", "value": "tsv"},
                                {"label": "BAM", "value": "bam"},
                                {"label": "SAM", "value": "sam"}
                                
                            ],
                            value="csv",
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
                ], className="g-2"),
            ]),
        ], className="mb-3", style={"background-color": colors['surface']}),

        # set and Reset Buttons
        html.Div(
            [
                dbc.Button("Apply Settings", id="apply-settings-button", color="success", className="me-2"),
                dbc.Button("Reset to Default", id="reset-settings-button", color="secondary"),
            ],
            className="d-flex justify-content-end mt-4"
        )
    ],
    id="settings-offcanvas",
    title="Settings",
    placement="end",
    is_open=False,
    style={"width": "450px", "background-color": colors['surface'], "color": colors['muted']}
)

# Console Offcanvas
console_offcanvas = dbc.Offcanvas(
    [
        html.H5("Console Log", className="mb-3", style={"color": colors['muted']}),
        html.Div(id='console-log', className='console-log', style={
            'whiteSpace': 'pre-wrap',
            'height': '100%',
            'overflowY': 'scroll',
            'backgroundColor': colors['background'],
            'padding': '10px',
            'border': f'1px solid {colors["muted"]}',
            'borderRadius': '4px',
            'fontFamily': 'monospace'
        })
    ],
    id="console-offcanvas",
    title="Console Log",
    placement="end",
    is_open=False,
    style={"width": "900px", "backgroundColor": colors['surface']}
)



#old header
def create_header():
    return dbc.Container([
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.Div([
                        DashIconify.DashIconify(icon="emojione:dna", width=48, height=48),
                        html.H1("GRAFT", className="ms-3 mb-0"),
                    ], className="d-flex align-items-center"),
                    html.P("Genomic Read Analysis and Filtering Tool", className="mb-0 mt-1"),
                ], className="header-content")
            ], width=12, className="header-wrapper")
        ], className="header-row")
    ], fluid=True, className="header-container")

layout_file_viewer = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H2("File Viewer", className="card-title mb-4"),
            dbc.Card([
                dbc.CardBody([
                    # File Selection Dropdown
                    dcc.Dropdown(
                        id='viewer-file-selector',
                        options=[],  # Will be populated dynamically
                        value=None,
                        clearable=False,
                        placeholder="Select a file to view",
                        className="mb-3",
                        style={"width": "100%"}
                    ),

                    # Filtering Options Section
                    dbc.Accordion([
                        dbc.AccordionItem([
                            html.Div([
                                html.Label("Mapping Quality (MAPQ) Range:", className="mb-2"),
                                dcc.RangeSlider(
                                    id='viewer-mapq-range',
                                    min=0,
                                    max=255,
                                    step=1,
                                    value=[0, 255],
                                    marks={i: {'label': str(i), 'style': {'color': colors['muted']}} for i in range(0, 256, 25)},
                                    tooltip={"placement": "bottom", "always_visible": True},
                                    updatemode='drag',
                                    className="mb-3",
                                ),
                                dbc.Label("Minimum Base Quality: ", className="mt-3"),
                                dcc.Input(
                                    id='viewer-min-base-quality',
                                    type='number',
                                    value=20,
                                    className="mb-3",
                                    style={
                                        "background-color": colors['surface'],
                                        "color": colors['on_surface'],
                                        "border": f"1px solid {colors['muted']}",
                                        "width": "100%"
                                    }
                                ),
                                html.Label("Mismatch Filter (BAM/SAM only):", className="mb-2"),
                                dcc.Dropdown(
                                    id='viewer-nm-dropdown',
                                    options=[
                                        {'label': 'All Reads', 'value': 'all'},
                                        {'label': 'NM:i:0', 'value': 0},
                                        {'label': 'NM:i:1', 'value': 1},
                                        {'label': 'NM:i:2', 'value': 2},
                                        {'label': 'NM:i:3', 'value': 3},
                                        {'label': 'NM:i:4', 'value': 4},
                                        {'label': 'NM:i:5', 'value': 5}
                                    ],
                                    value='all',
                                    clearable=False,
                                    className="mb-3",
                                ),
                                html.Label("C>T Change Filter:", className="mb-2"),
                                dcc.Checklist(
                                    id='viewer-ct-checklist',
                                    options=[
                                        {'label': ' Show only reads with C>T changes', 'value': 'only_ct'},
                                        {'label': ' Show reads with only C>T changes', 'value': 'exclusively_ct'},
                                        {'label': ' Exclude C>T changed reads', 'value': 'exclude_ct'},
                                        {'label': ' Subtract C>T changes from NM', 'value': 'subtract_ct'}
                                    ],
                                    value=[],
                                    className="mb-3",
                                ),
                                html.Label("C>T Change Count:", className="mb-2"),
                                dcc.Dropdown(
                                    id='viewer-ct-count-dropdown',
                                    options=[
                                        {'label': 'Any', 'value': 'any'},
                                        {'label': 'One', 'value': 1},
                                        {'label': 'Two', 'value': 2},
                                        {'label': 'Three', 'value': 3},
                                        {'label': 'Four', 'value': 4},
                                        {'label': 'Five', 'value': 5}
                                    ],
                                    value='any',
                                    clearable=False,
                                    className="mb-4",
                                ),
                            ])
                        ], title="Filtering Options", style={"background-color": colors['surface'], "color": colors['on_surface']})
                    ], start_collapsed=True, className="mb-4"),

                    # Navigation Controls and Search Options
                    html.Div([
                        dbc.Row([
                            dbc.Col([
                                dbc.InputGroup([
                                    dbc.Input(id='viewer-search-input', placeholder='Search by read name...', size="sm"),
                                    dbc.Button('Search', id='viewer-search-button', n_clicks=0, size="sm")
                                ], size="sm"),
                            ], width=3, className="mb-3 mb-sm-0"),
                            dbc.Col([
                                dbc.InputGroup([
                                    dbc.Input(id='viewer-sequence-search-input', placeholder='Search by sequence...', size="sm"),
                                    dbc.Button('Search', id='viewer-sequence-search-button', n_clicks=0, size="sm")
                                ], size="sm"),
                            ], width=3, className="mb-3 mb-sm-0"),
                            dbc.Col([
                                dbc.InputGroup([
                                    dbc.InputGroupText("Sort by"),
                                    dbc.Select(
                                        id='viewer-sort-field',
                                        options=[
                                            {'label': 'Name', 'value': 'name'},
                                            {'label': 'Length', 'value': 'length'},
                                            {'label': 'Mismatches', 'value': 'mismatches'}
                                        ],
                                        value='name',
                                        size="sm",  # Use bs_size instead of size
                                    ),
                                    dbc.InputGroupText("Order"),
                                    dbc.Select(
                                        id='viewer-sort-order',
                                        options=[
                                            {'label': '↑', 'value': 'asc'},
                                            {'label': '↓', 'value': 'desc'}
                                        ],
                                        value='asc',
                                        size="sm",  # Use bs_size instead of size
                                    ),
                                ], size="sm"),
                            ], width=4, className="mb-3 mb-sm-0"),
                            dbc.Col([
                                dbc.InputGroup([
                                    dbc.Input(id='viewer-page-input', type='number', min=1, placeholder='Page', size="sm"),
                                    dbc.Button('Go', id='viewer-go-button', n_clicks=0, size="sm", className="me-2"),
                                    dbc.Button("←", id='viewer-prev-button', n_clicks=0, disabled=True, size="sm", className="me-2"),
                                    dbc.Button("→", id='viewer-next-button', n_clicks=0, disabled=False, size="sm"),
                                ], size="sm"),
                            ], width=2),
                        ], className="g-2 align-items-center"),
                        html.Div(id='viewer-page-info', className="mt-2 text-muted small"),
                    ], className="mb-4"),


                    # File Viewer Content
                    dcc.Loading(
                        id="loading-viewer",
                        type="circle",
                        children=html.Div(id='file-viewer-content', style={"max-height": "1000px", "overflow": "auto", "padding": "10px", "background-color": colors['surface'], "border-radius": "4px"})
                    ),

                    dcc.Store(id='viewer-page-index', data=0),
                ])
            ], className="mb-4")
        ], width=11),
    ], justify='center'),
], fluid=True, className="px-4 py-3", style={"backgroundColor": colors['background']})




# Define the Navbar separately for reuse 
navbar = dbc.Navbar(
    dbc.Container([
        dcc.Link(
            dbc.Row([
                dbc.Col(DashIconify(icon="mdi:dna", width=40, color=colors['secondary'])),
                dbc.Col(dbc.NavbarBrand("GRAFT", className="ms-2")),
            ],
            align="center",
            className="g-0",
            ),
            href="/",
            style={"textDecoration": "none"},
        ),
        dbc.NavbarToggler(id="navbar-toggler"),
        dbc.Collapse(
            dbc.Nav([
                dbc.NavItem(dcc.Link("Home", href="/", className="nav-link")),
                dbc.NavItem(dcc.Link("Simulation", href="/simulation", className="nav-link")),
                dbc.NavItem(dcc.Link("Convert", href="/convert", className="nav-link")),
                dbc.NavItem(dcc.Link("File Viewer", href="/file-viewer", className="nav-link")),
                dbc.NavItem(dbc.Button("Settings", id="settings-button", color="light", outline=True)),
            ], className="ms-auto", navbar=True),
            id="navbar-collapse",
            navbar=True,
        ),
        console_offcanvas,
    ]),
    color=colors['surface'],
    dark=True,
    className="mb-4",
)



app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    navbar,
    html.Div(id='page-content'),
    # Stores and Download components
    dcc.Store(id='file-store', data={'files': []}),
    dcc.Store(id='selected-lengths', data=[]),
    dcc.Store(id='read-length-selection-store'),
    dcc.Store(id='selected-cg-contents', data=[]),
    dcc.Interval(id='centrifuge-interval', interval=5000, n_intervals=0, disabled=True),
    dcc.Interval(id='console-interval', interval=2000, n_intervals=0),
    dcc.Store(id='centrifuge-task-id'),
    dcc.Store(id='centrifuge-output-path'),
    dcc.Download(id="download-file"),
])



# Converter Page Layout
layout_convert = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H2("Batch File Conversion", className="card-title mb-4"),
            dcc.Upload(
                id='convert-upload-data',
                children=html.Div([
                    DashIconify(icon="mdi:cloud-upload", width=64, color=colors['secondary']),
                    html.Div("Drag and Drop or Click to Select Files", className="mt-2"),
                    html.Div("Supported formats: BAM, SAM, FASTA, FASTQ", className="text-muted"),
                ]),
                multiple=True,
                className="upload-box",
            ),
            html.Div(id='convert-upload-filenames', className='mt-2'),
            html.Label("Select Output Format:", className="mt-3"),
            dcc.Dropdown(
                id='convert-output-format',
                options=[
                    {'label': 'BAM', 'value': 'bam'},
                    {'label': 'SAM', 'value': 'sam'},
                    {'label': 'FASTA', 'value': 'fasta'},
                    {'label': 'FASTQ', 'value': 'fastq'},
                ],
                value='sam',
                clearable=False,
                className="mb-3",
            ),
            dbc.Button("Convert Files", id="convert-button", color="primary", className="mb-3", disabled=True),
            dcc.Download(id="convert-download-file"),
            html.Div(id='convert-status', className='mt-2'),
        ], width=8),
    ], justify='center'),
], fluid=True, className="px-4 py-3", style={"backgroundColor": colors['background']})



layout_main = dbc.Container([
    # Main content
    dbc.Row([
        # Sidebar
        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.H4("File Selection", className="card-title mb-4"),
                    dcc.Dropdown(
                        id='file-selector',
                        options=[],
                        value=None,
                        clearable=False,
                        placeholder="Select a file",
                        className="mb-4",
                    ),
                    html.Label("Sequencing Type:", className="mt-3 mb-2"),
                    dcc.RadioItems(
                        id='sequencing-type',
                        options=[
                            {'label': ' Single-stranded', 'value': 'single'},
                            {'label': ' Double-stranded', 'value': 'double'}
                        ],
                        value='single',
                        labelStyle={'display': 'block'},
                        className="mb-4",
                    ),
                    html.H5("Filtering Options", className="mt-4 mb-3"),
                    html.Label("Mapping Quality (MAPQ) Range:", className="mb-2"),
                    dcc.RangeSlider(
                        id='mapq-range',
                        min=0,
                        max=255,
                        step=1,
                        value=[0, 255],
                        marks={i: {'label': str(i), 'style': {'color': colors['muted']}} for i in range(0, 256, 25)},
                        tooltip={"placement": "bottom", "always_visible": True},
                        updatemode='drag',
                        className="mb-3",
                    ),
                    dbc.Label("Minimum Base Quality: ", className="mt-3"),
                    dcc.Input(
                        id='min-base-quality',
                        type='number',
                        value=20,
                        className="mb-3",
                        style={
                            "background-color": colors['surface'],
                            "color": colors['on_surface'],
                            "border": f"1px solid {colors['muted']}",
                            "width": "100%",
                            "margin-bottom": "10px"
                        }
                    ),
                    html.Label("Mismatch Filter (BAM/SAM only):", className="mb-2"),
                    dcc.Dropdown(
                        id='nm-dropdown',
                        options=[
                            {'label': 'All Reads', 'value': 'all'},
                            {'label': 'NM:i:0', 'value': 0},
                            {'label': 'NM:i:1', 'value': 1},
                            {'label': 'NM:i:2', 'value': 2},
                            {'label': 'NM:i:3', 'value': 3},
                            {'label': 'NM:i:4', 'value': 4},
                            {'label': 'NM:i:5', 'value': 5}
                        ],
                        value='all',
                        clearable=False,
                        className="mb-3",
                    ),
                    html.Label("C>T Change Filter:", className="mb-2"),
                    dcc.Checklist(
                        id='ct-checklist',
                        options=[
                            {'label': ' Show only reads with C>T changes', 'value': 'only_ct'},
                            {'label': ' Show reads with only C>T changes', 'value': 'exclusively_ct'},
                            {'label': ' Exclude C>T changed reads', 'value': 'exclude_ct'},
                            {'label': ' Subtract C>T changes from NM', 'value': 'subtract_ct'}
                        ],
                        value=[],
                        className="mb-3",
                    ),
                    html.Label("C>T Change Count:", className="mb-2"),
                    dcc.Dropdown(
                        id='ct-count-dropdown',
                        options=[
                            {'label': 'Any', 'value': 'any'},
                            {'label': 'One', 'value': 1},
                            {'label': 'Two', 'value': 2},
                            {'label': 'Three', 'value': 3},
                            {'label': 'Four', 'value': 4},
                            {'label': 'Five', 'value': 5}
                        ],
                        value='any',
                        clearable=False,
                        className="mb-4",
                    ),
                    dbc.Button("Run Centrifuge Analysis", id="centrifuge-button", color="primary", className="w-100 mb-3"),
                    dcc.Loading(
                        id="loading-centrifuge",
                        type="circle",
                        children=html.Div(id="centrifuge-output"),
                    ),
                    dbc.Button("Export Selected Reads", id="export-button", color="secondary", className="w-100 mb-3"),
                    dbc.Button("Show Alignment Statistics", id="stats-button", color="info", className="w-100 mb-3"),
                    dbc.Button("Clear Selection", id="clear-selection-button", outline=True, color="danger", className="w-100"),
                ])
            ], className="sticky-top")
        ], width=3, className="mb-4"),

        # Main content area
        dbc.Col([
            dbc.Row([
                # Upload Files Section
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H2("Upload Files", className="card-title mb-4", id="upload"),
                            dcc.Upload(
                                id='upload-data',
                                children=html.Div([
                                    DashIconify(icon="mdi:cloud-upload", width=64, color=colors['secondary']),
                                    html.Div("Drag and Drop or Click to Select Files", className="mt-2"),
                                    html.Div("Supported formats: BAM, SAM, FASTA, FASTQ", className="text-muted"),
                                ]),
                                multiple=True,
                                className="upload-box",
                            ),
                        ])
                    ], className="mb-4")
                ], width=12),
            ]),

            dbc.Row([
                # Analysis Results Section
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H2("Analysis Results", className="card-title mb-4", id="analyze"),
                            dcc.Loading(
                                id="loading-graphs",
                                type="circle",
                                children=dbc.Tabs([
                                    dbc.Tab(dcc.Graph(id='read-length-histogram'), label="Read Length"),
                                    dbc.Tab(dcc.Graph(id='overall-cg-histogram'), label="Overall CG"),
                                    dbc.Tab(dcc.Graph(id='cg-content-histogram'), label="CG Distribution"),
                                    dbc.Tab(dcc.Graph(id='damage-patterns'), label="Damage Patterns"),
                                    dbc.Tab(dcc.Graph(id='mismatch-frequency-plot'), label="Mismatch Frequency"),
                                    dbc.Tab(dcc.Graph(id='mismatch-type-bar-chart'), label="Mismatch Types"),
                                    dbc.Tab(dcc.Graph(id='damage-pattern-plot'), label="Damage Patterns Along Reads"),
                                    dbc.Tab(dcc.Graph(id='mapq-histogram'), label="MAPQ Distribution"),
                                    dbc.Tab(html.Pre(id='data-summary', style={'whiteSpace': 'pre-wrap'}), label="Data Summary"),
                                ]),
                            ),
                        ])
                    ], className="mb-4")
                ], width=12),
            ]),

            dbc.Row([
                # MapDamage2 Analysis Section
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H2("mapDamage2 Analysis", className="card-title mb-4"),
                            html.Label("Select BAM File:", className="mb-2"),
                            dcc.Dropdown(
                                id='mapdamage-bam-selector',
                                options=[],  # Will be populated dynamically
                                value=None,
                                clearable=False,
                                placeholder="Select a BAM file",
                                className="mb-4",
                            ),
                            html.Label("Select Reference Genome (FASTA):", className="mb-2"),
                            dcc.Dropdown(
                                id='mapdamage-ref-selector',
                                options=[],  # Will be populated dynamically
                                value=None,
                                clearable=False,
                                placeholder="Select a reference genome",
                                className="mb-4",
                            ),
                            dbc.Button("Run mapDamage2 Analysis", id="run-mapdamage-button", color="primary"),
                            html.Div(id="mapdamage-status", className="mt-3"),
                            dcc.Loading(
                                id="loading-mapdamage",
                                type="circle",
                                children=[
                                    html.Div(id='mapdamage-output')
                                ]
                            ),
                            dcc.Interval(id='mapdamage-interval', interval=5000, n_intervals=0, disabled=True),
                            dcc.Store(id='mapdamage-task-id'),
                            dcc.Store(id='mapdamage-output-dir'),
                        ])
                    ])
                ], width=12),
            ])
        ], width=9),
    ]),

    # Offcanvas for alignment stats
    dbc.Offcanvas(
        html.Div(id="alignment-stats"),
        id="stats-offcanvas",
        title="Alignment Summary Statistics",
        placement="end",
        backdrop=False,
    ),

    # Settings offcanvas
    settings_offcanvas,

    # # Stores and Download component
    # dcc.Store(id='file-store', data={'files': []}),
    # dcc.Store(id='selected-lengths', data=[]),
    # dcc.Interval(id='centrifuge-interval', interval=5000, n_intervals=0, disabled=True),
    # dcc.Interval(id='console-interval', interval=2000, n_intervals=0),
    # dcc.Store(id='centrifuge-task-id'),
    # dcc.Store(id='centrifuge-output-path'),
    # dcc.Download(id="download-file"),

], fluid=True, className="px-4 py-3", style={"backgroundColor": colors['background']})

app.validation_layout = html.Div([
    app.layout,
    layout_main,
    layout_convert,
    layout_file_viewer
])


# Add custom CSS, have to merge with external assets css
app.index_string = '''
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>GRAFT - Genomic Read Analysis and Filtering Tool</title>
        {%favicon%}
        {%css%}
        <style>
            body {
                background-color: ''' + colors['background'] + ''';
                color: ''' + colors['on_surface'] + ''';
            }
            .card {
                background-color: ''' + colors['surface'] + ''';
                border: none;
                border-radius: 8px;
                box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
            }
            .nav-link {
                color: ''' + colors['muted'] + ''' !important;
            }
            .nav-link:hover {
                color: ''' + colors['on_surface'] + ''' !important;
            }
            .console-log pre {
                background-color: transparent;
                color: inherit;
                font-family: monospace;
                white-space: pre-wrap;
            }
            .upload-box {
                border: 2px dashed ''' + colors['secondary'] + ''';
                border-radius: 8px;
                padding: 40px;
                text-align: center;
                cursor: pointer;
                transition: all 0.3s ease;
            }
            .upload-box:hover {
                background-color: rgba(255, 255, 255, 0.05);
            }
            /* Custom styling for dropdowns */
            .Select-control {
                background-color: ''' + colors['surface'] + ''' !important;
                border-color: ''' + colors['muted'] + ''' !important;
            }
            .Select-menu-outer {
                background-color: ''' + colors['surface'] + ''' !important;
            }
            .Select-option {
                background-color: ''' + colors['surface'] + ''' !important;
                color: ''' + colors['on_surface'] + ''' !important;
            }
            .Select-option:hover {
                background-color: rgba(255, 255, 255, 0.1) !important;
            }
        </style>
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



###################################################################################

 ######     ###    ##       ##       ########     ###     ######  ##    ##  ######  
##    ##   ## ##   ##       ##       ##     ##   ## ##   ##    ## ##   ##  ##    ## 
##        ##   ##  ##       ##       ##     ##  ##   ##  ##       ##  ##   ##       
##       ##     ## ##       ##       ########  ##     ## ##       #####     ######  
##       ######### ##       ##       ##     ## ######### ##       ##  ##         ## 
##    ## ##     ## ##       ##       ##     ## ##     ## ##    ## ##   ##  ##    ## 
 ######  ##     ## ######## ######## ########  ##     ##  ######  ##    ##  ######

##################################################################################


################### FILE VIEWER ##################

@app.callback(
    Output('viewer-ct-checklist', 'value'),
    [Input('viewer-ct-checklist', 'value')]
)
def enforce_mutual_exclusivity(selected_values):
    if len(selected_values) > 1:
        # Keep only the last selected option
        return [selected_values[-1]]
    else:
        return selected_values



@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/convert':
        return layout_convert
    elif pathname == '/simulation':
        return simulation.layout
    elif pathname == '/file-viewer':
        return layout_file_viewer  # New File Viewer layout
    else:
        return layout_main

@app.callback(
    Output('viewer-file-selector', 'options'),
    [Input('file-store', 'data')]
)
def update_viewer_file_selector(store_data):
    files = store_data.get('files', [])
    options = [{'label': f['filename'], 'value': f['filename']} for f in files]
    return options

@app.callback(
    [
        Output('file-viewer-content', 'children'),
        Output('viewer-prev-button', 'disabled'),
        Output('viewer-next-button', 'disabled'),
        Output('viewer-page-info', 'children'),
        Output('viewer-page-index', 'data')
    ],
    [
        Input('viewer-file-selector', 'value'),
        Input('viewer-prev-button', 'n_clicks'),
        Input('viewer-next-button', 'n_clicks'),
        Input('viewer-go-button', 'n_clicks'),
        Input('viewer-search-button', 'n_clicks'),
        Input('viewer-mapq-range', 'value'),
        Input('viewer-min-base-quality', 'value'),
        Input('viewer-nm-dropdown', 'value'),
        Input('viewer-ct-checklist', 'value'),
        Input('viewer-ct-count-dropdown', 'value'),
        Input('viewer-page-input', 'value'),
        Input('viewer-search-input', 'value'),
        Input('viewer-sequence-search-button', 'n_clicks'),
        Input('viewer-sort-field', 'value'),
        Input('viewer-sort-order', 'value'),
    ],
    [
        State('viewer-sequence-search-input', 'value'),
        State('file-store', 'data'),
        State('viewer-page-index', 'data')
    ]
)
def display_file_content(
    selected_file,
    prev_clicks,
    next_clicks,
    go_clicks,
    search_clicks,
    mapq_range,
    min_base_quality,
    selected_nm,
    ct_checklist,
    ct_count_value,
    page_input_value,
    search_input_value,
    sequence_search_clicks,
    sort_field,
    sort_order,
    sequence_search_input_value,
    store_data,
    page_index
):
    if selected_file is None:
        return html.Div("No file selected.", className="text-muted"), True, True, "", page_index

    # Retrieve the file content
    file_data = next((f for f in store_data['files'] if f['filename'] == selected_file), None)
    if file_data is None:
        return html.Div("File not found.", className="text-muted"), True, True, "", page_index

    temp_file_path, file_format = process_uploaded_file(file_data['content'], file_data['filename'])

    ct_checklist_tuple = tuple(ct_checklist)

    # Prepare C>T filters
    filter_ct = None
    exclusively_ct = False
    subtract_ct = False
    if 'only_ct' in ct_checklist_tuple:
        filter_ct = True
    elif 'exclude_ct' in ct_checklist_tuple:
        filter_ct = False
    elif 'exclusively_ct' in ct_checklist_tuple:
        exclusively_ct = True
    if 'subtract_ct' in ct_checklist_tuple:
        subtract_ct = True
    
    mapq_range_tuple = tuple(mapq_range)

    # load_and_process_file
    read_lengths, file_format, deduped_file_path = load_and_process_file(
        file_data['content'], file_data['filename'], selected_nm, filter_ct, exclusively_ct, ct_count_value,
        ct_checklist_tuple, mapq_range_tuple
    )


    # set additional filters based on base quality
    if file_format in ['bam', 'sam']:
        filtered_reads = []
        for length, cg, nm, seq, read, mapq in read_lengths:
            if read.mapping_quality < mapq_range[0] or read.mapping_quality > mapq_range[1]:
                continue
            if min_base_quality is not None:
                base_qualities = read.query_qualities
                if base_qualities is not None and min(base_qualities) < min_base_quality:
                    continue
            filtered_reads.append((length, cg, nm, seq, read, mapq))
        read_lengths = filtered_reads

    # Implement pagination and search functionality
    records_per_page = 20
    ctx = dash.callback_context
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # Initialize matching_indices to all reads
    matching_indices = list(range(len(read_lengths)))

    if triggered_id == 'viewer-prev-button':
        page_index = max(0, page_index - 1)
    elif triggered_id == 'viewer-next-button':
        page_index = page_index + 1
    elif triggered_id == 'viewer-go-button' and page_input_value:
        page_index = min(max(0, page_input_value - 1), (len(read_lengths) - 1) // records_per_page)
    elif triggered_id == 'viewer-search-button' and search_input_value:
        # Partial name matching
        search_name = search_input_value
        matching_indices = [i for i, (_, _, _, _, read, _) in enumerate(read_lengths)
                            if search_name in (read.query_name if isinstance(read, pysam.AlignedSegment) else read.id)]
        if not matching_indices:
            return html.Div(f"No reads found containing '{search_name}' in their names.", className="text-danger"), True, True, "", page_index
        page_index = 0  # Reset to first page after search
    elif triggered_id == 'viewer-sequence-search-button' and sequence_search_input_value:
        # Sequence content search
        search_seq = sequence_search_input_value.upper()
        matching_indices = [i for i, (_, _, _, _, read, _) in enumerate(read_lengths)
                            if search_seq in (read.query_sequence.upper() if isinstance(read, pysam.AlignedSegment) else str(read.seq).upper())]
        if not matching_indices:
            return html.Div(f"No reads found containing the sequence '{search_seq}'.", className="text-danger"), True, True, "", page_index
        page_index = 0  # Reset to first page after search
    else:
        page_index = 0

    # Adjust the reads to include only matching_indices
    filtered_read_lengths = [read_lengths[i] for i in matching_indices]

    # Implement sorting before pagination
    if sort_field and sort_order:
        if sort_field == 'name':
            # Sort by read name
            key_func = lambda x: x[4].query_name if isinstance(x[4], pysam.AlignedSegment) else x[4].id
        elif sort_field == 'length':
            # Sort by read length
            key_func = lambda x: x[0]  # length is at index 0
        elif sort_field == 'mismatches':
            # Sort by number of mismatches (NM tag)
            key_func = lambda x: x[2]  # nm is at index 2
        else:
            key_func = None

        if key_func:
            reverse = sort_order == 'desc'
            filtered_read_lengths.sort(key=key_func, reverse=reverse)

    # Update total_pages based on filtered reads
    total_pages = max(1, (len(filtered_read_lengths) + records_per_page - 1) // records_per_page)

    # Ensure page_index is within valid range
    page_index = min(page_index, total_pages - 1)

    start_index = page_index * records_per_page
    end_index = start_index + records_per_page
    page_reads = filtered_read_lengths[start_index:end_index]

    # Prepare content based on file format
    sequences = []
    if file_format in ['fasta', 'fa', 'fna', 'fastq', 'fq']:
        for _, _, _, _, record, _ in page_reads:
            sequences.append(html.Div([
                html.H5(record.id, className='sequence-header'),
                html.Pre(str(record.seq), className='sequence-content'),
                # For FASTQ, include qualities
                html.Pre(str(record.letter_annotations.get("phred_quality", "")), className='sequence-quality')
            ], className='sequence-block'))
    elif file_format in ['bam', 'sam']:
        for _, _, _, _, read, _ in page_reads:
            # Retrieve reference sequence if possible
            ref_seq = get_reference_sequence(read)
            # Highlight mismatches
            read_seq_html = highlight_mismatches(read)
            # Retrieve metadata
            metadata = read.get_tags()
            metadata_info = ', '.join([f"{tag}:{value}" for tag, value in metadata])

            sequences.append(html.Div([
                html.Div([
                    html.Span(read.query_name, className='sequence-header'),
                    html.Span(f" (Flag: {read.flag}, MapQ: {read.mapping_quality})", className='alignment-info'),
                ], className='header-row'),
                html.Div([
                    html.Span('Metadata:', className='sequence-label'),
                    html.Span(metadata_info, className='metadata-info'),
                ]),
                html.Div([
                    html.Span('CIGAR:', className='sequence-label'),
                    html.Span(read.cigarstring),
                    html.Span(' | '),
                    html.Span('NM:', className='sequence-label'),
                    html.Span(read.get_tag('NM')),
                    html.Span(' | '),
                    html.Span('MD:', className='sequence-label'),
                    html.Span(read.get_tag('MD')),
                ], className='metadata-info'),
                html.Div([
                    html.Span('Ref:', className='sequence-label'),
                    html.Pre(ref_seq if ref_seq else 'N/A', className='reference-sequence'),
                ]),
                html.Div([
                    html.Span('Read:', className='sequence-label'),
                    read_seq_html,  # This will include highlighted mismatches
                ]),
            ], className='sequence-block'))
    else:
        return html.Div("Unsupported file format.", className="text-danger"), True, True, "", page_index

    prev_disabled = page_index == 0
    next_disabled = page_index >= total_pages - 1
    page_info = f"Page {page_index + 1} of {total_pages}"

    return html.Div(sequences), prev_disabled, next_disabled, page_info, page_index




def get_reference_sequence(read): # quite unneccessary, was initially implemented to use if fasta refseq was available
    if read.is_unmapped or read.reference_id < 0:
        return None
    # need access to the reference genome sequence
    # Assuming you have a dictionary or function that maps reference names to sequences
    # For example, reference_sequences = {'chr1': 'ATGCT...'}
    ref_name = read.reference_name
    ref_start = read.reference_start
    ref_end = read.reference_end
    ref_seq = reference_sequences.get(ref_name)
    if ref_seq:
        return ref_seq[ref_start:ref_end]
    else:
        return None

def highlight_mismatches(read):
    read_seq = read.query_sequence
    aligned_pairs = read.get_aligned_pairs(with_seq=True)
    alignment = []

    for query_pos, ref_pos, ref_base in aligned_pairs:
        if query_pos is not None:
            read_base = read_seq[query_pos]

            if ref_pos is not None and ref_base is not None:
                # Compare read base and reference base
                read_base_upper = read_base.upper()
                ref_base_upper = ref_base.upper()

                if read_base_upper != ref_base_upper:
                    # Mismatch found
                    if read.is_reverse:
                        # Read is mapped to reverse strand
                        if ref_base_upper == 'G' and read_base_upper == 'A':
                            # G>A change on reverse strand (corresponds to C>T on forward strand)
                            alignment.append(f'<span style="color: orange; font-weight: bold;">{read_base}</span>')
                        else:
                            # Other mismatches
                            alignment.append(f'<span style="color: red; font-weight: bold;">{read_base}</span>')
                    else:
                        # Read is mapped to forward strand
                        if ref_base_upper == 'C' and read_base_upper == 'T':
                            # C>T change on forward strand
                            alignment.append(f'<span style="color: orange; font-weight: bold;">{read_base}</span>')
                        else:
                            # Other mismatches
                            alignment.append(f'<span style="color: red; font-weight: bold;">{read_base}</span>')
                else:
                    # Match
                    alignment.append(f'<span style="color: lightblue;">{read_base}</span>')
            else:
                # Reference base not available (e.g., soft clipping)
                alignment.append(f'<span style="color: gray;">{read_base}</span>')
        else:
            # Insertion or deletion relative to reference
            alignment.append(f'<span style="color: purple; font-style: italic;">-</span>')

    highlighted_seq = ''.join(alignment)
    return dash_dangerously_set_inner_html.DangerouslySetInnerHTML(
        f'<pre style="margin: 0; font-size: 14px;">{highlighted_seq}</pre>'
    )

reference_sequences = {}

@app.callback(
    Output('reference-store', 'data'),
    [Input('file-store', 'data')]
)
def load_reference_sequences(file_store_data): # unused for now, not necessary
    global reference_sequences
    for f in file_store_data.get('files', []):
        if f['format'] in ['fasta', 'fa', 'fna']:
            temp_file_path, _ = process_uploaded_file(f['content'], f['filename'])
            # Parse the reference genome and store sequences
            for record in SeqIO.parse(temp_file_path, 'fasta'):
                reference_sequences[record.id] = str(record.seq)
    return {'status': 'Reference genome loaded'}

    



###################################


##### Console #####################
    
@app.callback(
    Output("console-offcanvas", "is_open"),
    [Input("console-button", "n_clicks")],
    [State("console-offcanvas", "is_open")],
)
def toggle_console(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open

@app.callback(
    Output('console-log', 'children'),
    [Input('console-interval', 'n_intervals')]
)
def update_console_log(n):
    log_file_path = os.path.join(os.path.dirname(__file__), 'celery.log')
    try:
        with open(log_file_path, 'r') as f:
            lines = f.readlines()
            last_lines = lines[-200:]  # Adjust if needed
            log_text = ''.join(last_lines)
            
            # Add classes for different log levels
            log_text = log_text.replace('[ERROR]', '<span class="log-error">[ERROR]</span>')
            log_text = log_text.replace('[WARNING]', '<span class="log-warning">[WARNING]</span>')
            log_text = log_text.replace('[INFO]', '<span class="log-info">[INFO]</span>')
            log_text = log_text.replace('[DEBUG]', '<span class="log-debug">[DEBUG]</span>')
            
            # Convert ANSI to HTML
            conv = Ansi2HTMLConverter(inline=True)
            html_content = conv.convert(log_text, full=False)
            
            # Wrap the HTML content in a div with appropriate styling
            styled_html = f'''
            <div style="color: #f8f9fa; background-color: {colors['background']};">
                {html_content}
            </div>
            '''
            return dash_dangerously_set_inner_html.DangerouslySetInnerHTML(styled_html)
    except Exception as e:
        return f'Error reading log file: {e}'

#######################################################

##### Converter ######################################

    

@app.callback(
    Output('convert-button', 'disabled'),
    Input('convert-upload-data', 'contents')
)
def toggle_convert_button(contents):
    return contents is None


@app.callback(
    [Output('convert-download-file', 'data'),
     Output('convert-status', 'children')],
    [Input('convert-button', 'n_clicks')],
    [State('convert-upload-data', 'contents'),
     State('convert-upload-data', 'filename'),
     State('convert-output-format', 'value')],
    prevent_initial_call=True
)
def convert_files(n_clicks, contents_list, filenames, output_format):
    if not n_clicks:
        return None, ''
    if contents_list is None:
        return None, 'No files uploaded.'

    converted_files = []
    status_messages = []

    for content, filename in zip(contents_list, filenames):
        try:
            
            temp_input_path, input_format = process_uploaded_file(content, filename)

            
            supported_conversions = {
                ('bam', 'sam'),
                ('sam', 'bam'),
                ('bam', 'fasta'),
                ('sam', 'fasta'),
                ('bam', 'fastq'),
                ('sam', 'fastq'),
                ('fastq', 'fasta'),
                ('fasta', 'fastq'),
                ('fastq', 'bam'),
                ('fastq', 'sam'),
                ('fasta', 'bam'),
                ('fasta', 'sam'),
            }

            input_format = input_format.lower()
            output_format = output_format.lower()

            if (input_format, output_format) not in supported_conversions:
                status_messages.append(f'Conversion from {input_format.upper()} to {output_format.upper()} is not supported for {filename}.')
                continue

            # temp file path, reused as in main analysis tab
            output_filename = f"{os.path.splitext(filename)[0]}_converted.{output_format}"
            temp_output_path = os.path.join(CUSTOM_TEMP_DIR, output_filename)

            # conversion
            if input_format in ['bam', 'sam'] and output_format in ['bam', 'sam']:
                # sam/bam
                with pysam.AlignmentFile(temp_input_path, "rb" if input_format == 'bam' else 'r') as infile:
                    with pysam.AlignmentFile(temp_output_path, "wb" if output_format == 'bam' else 'w', header=infile.header) as outfile:
                        for read in infile:
                            outfile.write(read)
            elif input_format in ['bam', 'sam'] and output_format in ['fasta', 'fastq']:
                # bam/sam to fasta/q
                with pysam.AlignmentFile(temp_input_path, "rb" if input_format == 'bam' else 'r') as infile:
                    with open(temp_output_path, 'w') as outfile:
                        for read in infile:
                            if read.is_unmapped:
                                continue
                            seq = read.query_sequence
                            qual = read.query_qualities
                            read_name = read.query_name
                            if output_format == 'fasta':
                                outfile.write(f'>{read_name}\n{seq}\n')
                            elif output_format == 'fastq':
                                if qual is None:
                                    qual_str = 'I' * len(seq)  # Assign dummy quality
                                else:
                                    qual_str = ''.join([chr(q + 33) for q in qual])
                                outfile.write(f'@{read_name}\n{seq}\n+\n{qual_str}\n')
            elif input_format == 'fastq' and output_format in ['bam', 'sam']:
                # fastq to bam/sam
                header = {
                    'HD': {'VN': '1.0'},
                    'SQ': [{'LN': 1, 'SN': 'chrUn'}]  # Dummy reference
                }
                with pysam.AlignmentFile(temp_output_path, "wb" if output_format == 'bam' else 'w', header=header) as outfile:
                    for record in SeqIO.parse(temp_input_path, 'fastq'):
                        # new unmapped read
                        a = pysam.AlignedSegment()
                        a.query_name = record.id
                        a.query_sequence = str(record.seq)
                        a.flag = 4  # Unmapped
                        a.mapping_quality = 0
                        a.query_qualities = pysam.qualitystring_to_array(record.letter_annotations["phred_quality"])
                        outfile.write(a)
            elif input_format == 'fasta' and output_format == 'fastq':
                # fasta to fastq
                with open(temp_output_path, 'w') as outfile:
                    for record in SeqIO.parse(temp_input_path, 'fasta'):
                        record.letter_annotations["phred_quality"] = [40] * len(record.seq)  # Assign dummy quality
                        SeqIO.write(record, outfile, 'fastq')
            elif input_format == 'fasta' and output_format in ['bam', 'sam']:
                # fasta to bam/sam
                header = {
                    'HD': {'VN': '1.0'},
                    'SQ': [{'LN': 1, 'SN': 'chrUn'}]  # Dummy reference
                }
                with pysam.AlignmentFile(temp_output_path, "wb" if output_format == 'bam' else 'w', header=header) as outfile:
                    for record in SeqIO.parse(temp_input_path, 'fasta'):
                        # new unmapped read again
                        a = pysam.AlignedSegment()
                        a.query_name = record.id
                        a.query_sequence = str(record.seq)
                        a.flag = 4  # Unmapped
                        a.mapping_quality = 0
                        # Assign dummy quality scores
                        a.query_qualities = [40] * len(record.seq)
                        outfile.write(a)
            elif input_format in ['fastq', 'fasta'] and output_format in ['fastq', 'fasta']:
                # fasta/fastq
                records = list(SeqIO.parse(temp_input_path, input_format))
                if input_format == 'fasta' and output_format == 'fastq':
                    for record in records:
                        record.letter_annotations["phred_quality"] = [40] * len(record.seq)  # Assign dummy quality
                SeqIO.write(records, temp_output_path, output_format)
            else:
                status_messages.append(f'Conversion from {input_format.upper()} to {output_format.upper()} is not implemented for {filename}.')
                continue

            converted_files.append(temp_output_path)
            status_messages.append(f'Successfully converted {filename}.')

        except Exception as e:
            status_messages.append(f'Error converting {filename}: {str(e)}')

    if not converted_files:
        return None, html.Div([
            html.H5("Conversion Status:"),
            html.Ul([html.Li(msg) for msg in status_messages])
        ])

    if len(converted_files) == 1:
        # download directly if only one file
        return dcc.send_file(converted_files[0]), html.Div([
            html.H5("Conversion Status:"),
            html.Ul([html.Li(msg) for msg in status_messages])
        ])
    else:
        # Create a ZIP archive if more than one file
        zip_filename = f'converted_files_{uuid.uuid4()}.zip'
        zip_path = os.path.join(CUSTOM_TEMP_DIR, zip_filename)

        with zipfile.ZipFile(zip_path, 'w') as zipf:
            for file_path in converted_files:
                zipf.write(file_path, os.path.basename(file_path))

        return dcc.send_file(zip_path), html.Div([
            html.H5("Conversion Status:"),
            html.Ul([html.Li(msg) for msg in status_messages])
        ])



@app.callback(
    Output('convert-upload-filenames', 'children'),
    Input('convert-upload-data', 'filename')
)
def update_filenames(filenames):
    if filenames:
        file_list = html.Ul([html.Li(f) for f in filenames])
        return html.Div([
            html.Strong("Uploaded files:"),
            file_list
        ])
    else:
        return ''

#######################################################


######### Mapdamage Section ###########################


@app.callback(
    [Output('mapdamage-output', 'children'),
     Output('mapdamage-interval', 'disabled')],
    [Input('mapdamage-interval', 'n_intervals')],
    [State('mapdamage-task-id', 'data'),
     State('mapdamage-output-dir', 'data')]
)
def update_mapdamage_output(n_intervals, task_id, output_dir):
    if not task_id:
        raise dash.exceptions.PreventUpdate

    task_result = run_mapdamage_task.AsyncResult(task_id)
    state = task_result.state

    if state == 'PENDING':
        return html.Div("mapDamage2 analysis is pending...", className="text-info"), False
    elif state == 'STARTED':
        return html.Div("mapDamage2 analysis is running...", className="text-info"), False
    elif state == 'SUCCESS':
        # Disable the interval
        # Process and display the results
        return display_mapdamage_results(output_dir), True
    elif state == 'FAILURE':
        error = str(task_result.info)
        return html.Div(f"Error running mapDamage2: {error}", className="text-danger"), True
    else:
        return html.Div(f"mapDamage2 analysis state: {state}", className="text-info"), False



@app.callback(
    [Output('mapdamage-status', 'children'),
     Output('mapdamage-task-id', 'data')],
    [Input('run-mapdamage-button', 'n_clicks')],
    [State('mapdamage-bam-selector', 'value'),
     State('mapdamage-ref-selector', 'value'),
     State('file-store', 'data')]
)
def run_mapdamage(n_clicks, selected_bam, selected_ref, store_data):
    if not n_clicks:
        return dash.no_update, None

    if selected_bam is None:
        return html.Div("Please select a BAM file.", className="text-danger"), None

    if selected_ref is None:
        return html.Div("Please select a reference genome.", className="text-danger"), None

    # Retrieve BAM file path
    bam_file_data = next((f for f in store_data['files'] if f['filename'] == selected_bam), None)
    if bam_file_data is None:
        return html.Div("Selected BAM file not found.", className="text-danger"), None
    bam_path = process_uploaded_file(bam_file_data['content'], bam_file_data['filename'])[0]

    # Retrieve Reference genome path
    ref_file_data = next((f for f in store_data['files'] if f['filename'] == selected_ref), None)
    if ref_file_data is None:
        return html.Div("Selected reference genome not found.", className="text-danger"), None
    ref_path = process_uploaded_file(ref_file_data['content'], ref_file_data['filename'])[0]

    # Output directory
    output_dir = os.path.join(CUSTOM_TEMP_DIR, f"mapdamage_results_{uuid.uuid4()}")
    os.makedirs(output_dir, exist_ok=True)

    # Start the mapDamage2 task
    task = run_mapdamage_task.delay(bam_path, ref_path, output_dir)

    return html.Div(f"mapDamage2 analysis started. Task ID: {task.id}", className="text-info"), task.id


#######################################################



######## Settings and stats menu ###############################


# Callback to toggle settings offcanvas
@app.callback(
    Output("settings-offcanvas", "is_open"),
    [Input("settings-button", "n_clicks"),
     Input('apply-settings-button', 'n_clicks')],
    [State("settings-offcanvas", "is_open"),
     State('centrifuge-db-path', 'value')],
    prevent_initial_call=True
)
def toggle_settings(settings_click, apply_click, is_open, centrifuge_db_path):
    ctx = dash.callback_context

    if not ctx.triggered:
        return is_open
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
        if button_id == 'settings-button':
            return not is_open  # Toggle the offcanvas
        elif button_id == 'apply-settings-button':
            # Save settings
            config['centrifuge_db_path'] = centrifuge_db_path
            with open('config.yaml', 'w') as f:
                yaml.dump(config, f)
            return False  # Close the offcanvas
    return is_open


@app.callback(
    Output('alignment-stats', 'children'),
    [Input('file-selector', 'value'),
     Input('file-store', 'data')]
)
def update_alignment_stats(selected_file, store_data):
    if selected_file is None:
        return html.Div("No file selected.", className="text-muted")

    file_data = next((f for f in store_data['files'] if f['filename'] == selected_file), None)
    if file_data is None:
        return html.Div("No data available.", className="text-muted")

    temp_file_path, file_format = process_uploaded_file(file_data['content'], file_data['filename'])

    stats = calculate_alignment_stats(temp_file_path, file_format)

    stats_text = get_alignment_stats_text(stats)

    # Only include mismatch information if it's available
    mismatch_info = []
    if isinstance(stats.get('Mismatch Counts'), dict):  # This ensures it's only for SAM/BAM files
        mismatch_info = [
            html.Li(f"C>T Changes: {stats['Mismatch Counts'].get('C>T', 0)}"),
            html.Li(f"G>A Changes: {stats['Mismatch Counts'].get('G>A', 0)}"),
            html.Li(f"Other Mismatches: {sum(stats['Mismatch Counts'].values()) - stats['Mismatch Counts'].get('C>T', 0) - stats['Mismatch Counts'].get('G>A', 0)}")
        ]

    mismatch_pie_chart = create_detailed_mismatch_pie_chart(stats)

    return html.Div([
        html.Ul([
            html.Li(f"Total Reads: {stats['Total Reads']}"),
            html.Li(f"Mapped Reads: {stats['Mapped Reads']} ({stats['Mapped Percentage']:.2f}%)" if stats['Mapped Reads'] != 'N/A' else "Mapped Reads: N/A"),
            html.Li(f"Duplicate Reads: {stats['Duplicate Reads']} ({stats['Duplicate Percentage']:.2f}%)" if stats['Duplicate Reads'] != 'N/A' else "Duplicate Reads: N/A"),
            html.Li(f"Total Mismatches: {stats['Total Mismatches']}")
        ] + mismatch_info),
        dcc.Graph(figure=mismatch_pie_chart)
    ])

@app.callback(
    Output('stats-offcanvas', 'is_open'),
    [Input('stats-button', 'n_clicks')],
    [State('stats-offcanvas', 'is_open')]
)
def toggle_offcanvas(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open

#############################################################



####### CENTRIFUGE IMPLEMENTATION ###########################

def parse_and_display_centrifuge_results(report_path):
    try:
        df = pd.read_csv(report_path, sep='\t')
    except Exception as e:
        return html.Div(f"Error reading Centrifuge report: {e}", className="text-danger")
    
    # Process the DataFrame as needed
    top_n = 10
    df_top = df.head(top_n)

    fig = go.Figure(data=[go.Bar(
        x=df_top['name'],
        y=df_top['numUniqueReads'],
        text=df_top['taxID'],
        hoverinfo='text',
    )])

    fig.update_layout(
        title='Top Taxonomic Classifications',
        xaxis_title='Taxon',
        yaxis_title='Number of Unique Reads',
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        font=dict(color=colors['muted'])
    )

    return dcc.Graph(figure=fig)


# @celery_app.task
# def run_centrifuge_task(file_path, output_path, report_path, db_path):
#     centrifuge_cmd = [
#         "centrifuge",
#         "-x", db_path,
#         "-U", file_path,
#         "-S", output_path,
#         "--report-file", report_path
#     ]
#     try:
#         subprocess.run(centrifuge_cmd, check=True)
#     except subprocess.CalledProcessError as e:
#         raise e


@app.callback(
    [Output('centrifuge-output', 'children'),
     Output('centrifuge-interval', 'disabled'),
     Output('centrifuge-task-id', 'data'),
     Output('centrifuge-output-path', 'data')],
    [Input('centrifuge-button', 'n_clicks')],
    [State('file-selector', 'value'),
     State('file-store', 'data'),
     State('centrifuge-db-path', 'value')],
    prevent_initial_call=True
)
def run_centrifuge(n_clicks, selected_file, store_data, db_path):
    if selected_file is None:
        return html.Div("No file selected for analysis.", className="text-danger"), True, None, None

    if not db_path:
        return html.Div("Centrifuge database path not specified. Please set it in the settings.", className="text-danger"), True, None, None

    
    file_data = next((f for f in store_data['files'] if f['filename'] == selected_file), None)
    if file_data is None:
        return html.Div("File data not found.", className="text-danger"), True, None, None

    temp_file_path, file_format = process_uploaded_file(file_data['content'], file_data['filename'])

    # file must be in fastq (centrifuge only supports fastq; might add automatic conversion)
    if file_format not in ['fastq', 'fq']:
        return html.Div("Centrifuge requires FASTQ files. Please upload a FASTQ file.", className="text-danger"), True, None, None

    
    centrifuge_output_path = os.path.join(CUSTOM_TEMP_DIR, f"centrifuge_output_{uuid.uuid4()}.txt")
    centrifuge_report_path = os.path.join(CUSTOM_TEMP_DIR, f"centrifuge_report_{uuid.uuid4()}.txt")

    
    task = run_centrifuge_task.delay(temp_file_path, centrifuge_output_path, centrifuge_report_path, db_path)

    # return message
    return html.Div(f"Centrifuge analysis started. Task ID: {task.id}", className="text-success"), False, task.id, centrifuge_report_path



def update_centrifuge_output(n_intervals, task_id, centrifuge_report_path):
    if not task_id:
        raise dash.exceptions.PreventUpdate

    task_result = run_centrifuge_task.AsyncResult(task_id)
    if task_result.state == 'PENDING':
        return html.Div("Centrifuge analysis is pending...", className="text-info"), False
    elif task_result.state == 'STARTED':
        return html.Div("Centrifuge analysis is running...", className="text-info"), False
    elif task_result.state == 'SUCCESS':
        # Disable the interval
        return parse_and_display_centrifuge_results(centrifuge_report_path), True
    elif task_result.state == 'FAILURE':
        return html.Div(f"Error running Centrifuge: {task_result.info}", className="text-danger"), True
    else:
        return html.Div(f"Centrifuge analysis state: {task_result.state}", className="text-info"), False
    


#####################################################################

############### Main Callbacks ######################################


@app.callback(
    Output('ct-checklist', 'value'),
    Input('ct-checklist', 'value')
)
def update_ct_checklist(ct_value):
    if len(ct_value) > 1:
        # Keep only the last selected option
        return [ct_value[-1]]
    return ct_value


@app.callback(
    [
        Output('file-store', 'data'),
        Output('file-selector', 'options'),
        Output('file-selector', 'value'),
        Output('mapdamage-bam-selector', 'options'),
        Output('mapdamage-ref-selector', 'options')
    ],
    [
        Input('upload-data', 'contents'),
        Input('upload-data', 'filename')
    ],
    [
        State('file-store', 'data'),
        State('file-selector', 'value')
    ]
)
def update_file_store(contents, filenames, store_data, current_selection):
    if store_data is None:
        store_data = {'files': []}

    if contents is None:
        # Use existing options from store_data
        options = [{'label': f['filename'], 'value': f['filename']} for f in store_data['files']]
        bam_options = [{'label': f['filename'], 'value': f['filename']} for f in store_data['files'] if f['format'] == 'bam']
        ref_options = [{'label': f['filename'], 'value': f['filename']} for f in store_data['files'] if f['format'] in ['fasta', 'fa', 'fna']]
        return store_data, options, current_selection, bam_options, ref_options

    files = store_data.get('files', [])
    new_file_added = False
    for content, filename in zip(contents, filenames):
        if filename not in [f['filename'] for f in files]:
            file_format = filename.split('.')[-1].lower()
            if file_format not in ['bam', 'sam', 'fasta', 'fa', 'fna', 'fq', 'fastq']:
                continue
            files.append({'filename': filename, 'content': content, 'format': file_format})
            new_file_added = True

    store_data['files'] = files
    options = [{'label': f['filename'], 'value': f['filename']} for f in files]

    # Populate BAM file options for mapDamage2
    bam_options = [{'label': f['filename'], 'value': f['filename']} for f in files if f['format'] == 'bam']

    # Populate FASTA file options for reference genome selector
    ref_options = [{'label': f['filename'], 'value': f['filename']} for f in files if f['format'] in ['fasta', 'fa', 'fna']]

    # Automatically select the first file if no file is currently selected
    if current_selection is None and new_file_added:
        current_selection = files[0]['filename']

    return store_data, options, current_selection, bam_options, ref_options



# @app.callback(
#     Output('reference-store', 'data'),
#     [Input('upload-reference', 'contents')],
#     [State('upload-reference', 'filename')]
# )
# def store_reference_genome(content, filename):
#     if content is None:
#         return dash.no_update
#     file_format = filename.split('.')[-1].lower()
#     if file_format not in ['fasta', 'fa', 'fna']:
#         return dash.no_update  # Only accept FASTA files
#     # Save the reference genome
#     ref_path = save_uploaded_file(content, filename)
#     return {'filename': filename, 'path': ref_path}


# @app.callback(
#     [
#         Output('cg-content-histogram', 'figure'),
#         Output('selected-lengths', 'data')
#     ],
#     [
#         Input('read-length-histogram', 'selectedData'),
#     ],
#     [
#         State('read-length-histogram', 'figure'),
#         State('file-store', 'data'),
#         State('file-selector', 'value'),
#         State('nm-dropdown', 'value'),
#         State('ct-checklist', 'value'),
#         State('ct-count-dropdown', 'value'),
#         State('mapq-range', 'value'),
#     ]
# )
# def update_selection(
#     read_length_selectedData, figure_state, store_data, selected_file,
#     selected_nm, ct_checklist, ct_count_value, mapq_range
# ):
#     if selected_file is None or read_length_selectedData is None:
#         raise dash.exceptions.PreventUpdate

#     ct_checklist_tuple = tuple(ct_checklist)
#     mapq_range_tuple = tuple(mapq_range)

#     # Prepare C>T filters
#     filter_ct = None
#     exclusively_ct = False
#     subtract_ct = False
#     if 'only_ct' in ct_checklist_tuple:
#         filter_ct = True
#     elif 'exclude_ct' in ct_checklist_tuple:
#         filter_ct = False
#     elif 'exclusively_ct' in ct_checklist_tuple:
#         exclusively_ct = True
#     if 'subtract_ct' in ct_checklist_tuple:
#         subtract_ct = True

#     # Extract selected lengths
#     selected_x_values = [point['x'] for point in read_length_selectedData['points']]
#     selected_lengths = [int(x) for x in selected_x_values]

#     # If no selection is made, consider all lengths
#     if not selected_lengths:
#         file_data = next((f for f in store_data['files'] if f['filename'] == selected_file), None)
#         if file_data is None:
#             raise dash.exceptions.PreventUpdate
#         read_lengths, _, _ = load_and_process_file(
#             file_data['content'], file_data['filename'], selected_nm,
#             filter_ct, exclusively_ct, ct_count_value, ct_checklist_tuple, mapq_range_tuple
#         )
#         selected_lengths = list(set(length for length, _, _, _, _, _ in read_lengths))

#     # Retrieve read_lengths
#     file_data = next((f for f in store_data['files'] if f['filename'] == selected_file), None)
#     if file_data is None:
#         raise dash.exceptions.PreventUpdate

#     # Process the file to get read_lengths
#     read_lengths, _, _ = load_and_process_file(
#         file_data['content'], file_data['filename'], selected_nm,
#         filter_ct, exclusively_ct, ct_count_value, ct_checklist_tuple, mapq_range_tuple
#     )

#     # Filter read_lengths based on selected_lengths
#     cg_contents_for_length = [cg for length, cg, _, _, _, _ in read_lengths if length in selected_lengths]

#     # Create CG Content Histogram
#     cg_content_fig = create_cg_content_histogram(cg_contents_for_length, selected_lengths)

#     # Store selected_lengths in the figure's custom data
#     figure_state['custom_data'] = {'selected_lengths': selected_lengths}

#     return cg_content_fig, figure_state


@app.callback(
    [
        Output('read-length-histogram', 'figure'),
        Output('overall-cg-histogram', 'figure'),
        Output('cg-content-histogram', 'figure'),
        Output('damage-patterns', 'figure'),
        Output('selected-lengths', 'data'),
        Output('mismatch-frequency-plot', 'figure'),
        Output('data-summary', 'children'),
        Output('mapq-histogram', 'figure'),
        Output('mismatch-type-bar-chart', 'figure'),
        Output('damage-pattern-plot', 'figure')
    ],
    [
        Input('file-selector', 'value'),
        Input('file-store', 'data'),
        Input('mapq-range', 'value'),
        Input('nm-dropdown', 'value'),
        Input('ct-checklist', 'value'),
        Input('ct-count-dropdown', 'value'),
        Input('min-base-quality', 'value'),
        Input('read-length-histogram', 'selectedData'),
        Input('read-length-histogram', 'clickData'),
        Input('clear-selection-button', 'n_clicks'),
        Input('sequencing-type', 'value')
    ],
    [
        State('selected-lengths', 'data'),
        State('read-length-histogram', 'figure')
    ]
)


def update_histograms(
    selected_file, 
    store_data, 
    mapq_range, 
    selected_nm, 
    ct_checklist, 
    ct_count_value,
    min_base_quality, 
    read_length_selectedData, 
    read_length_clickData, 
    clear_selection_nclicks, 
    sequencing_type,
    current_selected_lengths, 
    current_fig
):
    # Check if no file is selected
    if selected_file is None:
        empty_figure = go.Figure()
        empty_text = ""
        return (
            empty_figure, empty_figure, empty_figure, empty_figure, [],
            empty_figure, empty_text, empty_figure, empty_figure, empty_figure
        )


    # Initialize selected_lengths
    selected_lengths = []

    # Handle clear selection
    if clear_selection_nclicks:
        selected_lengths = []
    else:
        # Handle selection from selectedData or clickData
        selected_lengths = current_selected_lengths or []
        
        if read_length_selectedData and read_length_selectedData["points"]:
            selected_x_values = [point['x'] for point in read_length_selectedData['points']]
            # Update selected_lengths by intersecting current selection
            if selected_lengths:
                selected_lengths = list(set(selected_lengths).intersection(set(selected_x_values)))
            else:
                selected_lengths = selected_x_values
        elif read_length_clickData and read_length_clickData["points"]:
            # Handle clickData (single bar click)
            clicked_length = int(read_length_clickData['points'][0]['x'])
            if clicked_length not in selected_lengths:
                selected_lengths.append(clicked_length)

    # Convert ct_checklist to tuple
    ct_checklist_tuple = tuple(ct_checklist)

    # Prepare C>T filters
    filter_ct = None
    exclusively_ct = False
    subtract_ct = False
    if 'only_ct' in ct_checklist_tuple:
        filter_ct = True
    elif 'exclude_ct' in ct_checklist_tuple:
        filter_ct = False
    elif 'exclusively_ct' in ct_checklist_tuple:
        exclusively_ct = True
    if 'subtract_ct' in ct_checklist_tuple:
        subtract_ct = True

    # Convert mapq_range to tuple
    mapq_range_tuple = tuple(mapq_range)

    file_data = next((f for f in store_data['files'] if f['filename'] == selected_file), None)
    if file_data is None:
        empty_figure = go.Figure()
        empty_text = ""
        return (
            empty_figure, empty_figure, empty_figure, empty_figure, [],
            empty_figure, empty_figure, empty_text, empty_figure, empty_figure, empty_figure
        )

    # Initialize default empty figures for mismatch-related plots
    empty_mismatch_fig = go.Figure().update_layout(
        title='Mismatch Frequency Not Available for This File Format',
        xaxis=dict(title='', color=colors['muted']),
        yaxis=dict(title='', color=colors['muted']),
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        font=dict(color=colors['muted'])
    )
    mismatch_freq_fig = empty_mismatch_fig
    old_mismatch_freq_fig = empty_mismatch_fig

    # Adjust the call to load_and_process_file
    read_lengths, file_format, deduped_file_path = load_and_process_file(
        file_data['content'], file_data['filename'], selected_nm, filter_ct,
        exclusively_ct, ct_count_value, ct_checklist_tuple, mapq_range_tuple
    )

    # Handle mismatch frequency based on file format
    if file_format in ['bam', 'sam']:
        frequencies = calculate_mismatch_frequency(read_lengths, file_format, sequencing_type)
        if frequencies:
            mismatch_freq_fig = create_mismatch_frequency_plot(frequencies, sequencing_type)
    else:
        frequencies = None
        frequencies_old = None

    bin_centers, hist, bins = prepare_histogram_data(read_lengths)

    # Create read length histogram with uirevision
    read_length_fig = create_read_length_histogram(bin_centers, hist, selected_file, read_lengths, bins)

    overall_cg_content = calculate_overall_cg_content(read_lengths)
    overall_cg_fig = go.Figure()
    overall_cg_fig.add_trace(go.Histogram(
        x=overall_cg_content,
        nbinsx=20,
        marker_color=colors['marker']
    ))
    overall_cg_fig.update_layout(
        title='Overall CG Content Distribution',
        xaxis=dict(title='CG Content', color=colors['muted']),
        yaxis=dict(title='Frequency', color=colors['muted']),
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        font=dict(color=colors['muted'])
    )

    # Generate CG Content Histogram based on selected_lengths
    cg_contents_for_length = [
        cg for length, cg, _, _, _, _ in read_lengths if length in selected_lengths
    ]
    cg_content_fig = create_cg_content_histogram(cg_contents_for_length, selected_lengths)

    # Generate Damage Patterns Figure
    five_prime_end, three_prime_end = get_damage_patterns(deduped_file_path, file_format, filter_ct=filter_ct)
    damage_fig = create_damage_pattern_figure(five_prime_end, three_prime_end, selected_file)

    mapq_scores = [mapq for _, _, _, _, _, mapq in read_lengths if mapq is not None]
    # Create MAPQ histogram
    mapq_fig = create_mapq_histogram(mapq_scores)

    # Calculate alignment stats
    temp_file_path, file_format = process_uploaded_file(file_data['content'], file_data['filename'])
    stats = calculate_alignment_stats(temp_file_path, file_format)

    # Initialize empty figures for mismatch type and damage pattern
    mismatch_type_fig = go.Figure()
    damage_pattern_fig = go.Figure()

    # Set filters to mismatches
    if stats.get('Mismatch Details') != 'N/A':
        filtered_mismatches = filter_mismatches(
            stats['Mismatch Details'],
            min_base_quality=min_base_quality,
            min_mapping_quality=mapq_range[0]
        )
        # Recalculate mismatch counts
        filtered_mismatch_counts = {}
        for mismatch in filtered_mismatches:
            key = f"{mismatch['ref_base']}>{mismatch['read_base']}"
            filtered_mismatch_counts[key] = filtered_mismatch_counts.get(key, 0) + 1

        # Plotting
        mismatch_type_fig = create_mismatch_type_bar_chart(filtered_mismatch_counts)
        damage_pattern_fig = create_damage_pattern_plot(filtered_mismatches)

    # Data summary after generating all plots
    read_length_text = get_read_length_data(read_lengths)
    cg_content_text = get_overall_cg_content_data(read_lengths)
    damage_patterns_text = get_damage_patterns_data(five_prime_end, three_prime_end)
    mismatch_frequency_text = get_mismatch_frequency_data(frequencies, sequencing_type)
    stats_text = get_alignment_stats_text(stats)
    read_length_cg_average_text = get_read_length_CG_average(read_lengths)
    filter_options_text = get_filter_options_text(selected_nm, ct_checklist, ct_count_value, mapq_range_tuple)
    
    selected_cg_contents = [cg for length, cg, _, _, _, _ in read_lengths if length in selected_lengths]
    if selected_cg_contents:
        cg_series = pd.Series(selected_cg_contents)
        cg_distribution_text = (
            "CG Content for Selected Reads:\n"
            f"Mean: {cg_series.mean():.4f}\n"
            f"Median: {cg_series.median():.4f}\n"
            f"Std Dev: {cg_series.std():.4f}\n"
        )
    else:
        cg_distribution_text = "No CG content data for selected reads.\n"

    # Combined data summary text
    data_summary_text = (
        filter_options_text + "\n" +
        read_length_text + "\n" +
        read_length_cg_average_text + "\n" +
        cg_content_text + "\n" +
        cg_distribution_text + "\n" +
        damage_patterns_text + "\n" +
        mismatch_frequency_text + "\n" +
        stats_text
    )

    # Return updated figures and data
    return (
        read_length_fig, overall_cg_fig, cg_content_fig, damage_fig,
        selected_lengths, mismatch_freq_fig, data_summary_text, mapq_fig,
        mismatch_type_fig, damage_pattern_fig
    )







    


@app.callback(
    Output("download-file", "data"),
    [Input("export-button", "n_clicks")],
    [
        State('file-selector', 'value'),
        State('file-store', 'data'),
        State('nm-dropdown', 'value'),
        State('ct-checklist', 'value'),
        State('ct-count-dropdown', 'value'),
        State('mapq-range', 'value'),
        State('selected-lengths', 'data'),
    ],
    prevent_initial_call=True
)
def export_reads(
    n_clicks, selected_file, store_data, selected_nm, ct_checklist, ct_count_value, mapq_range, selected_lengths
):
    if n_clicks is None or n_clicks == 0:
        raise dash.exceptions.PreventUpdate

    if selected_file is None:
        return None

    if not selected_lengths:
        # No reads are selected; prevent export
        print("No read lengths selected. Export aborted.")
        return None

    ct_checklist_tuple = tuple(ct_checklist)
    mapq_range_tuple = tuple(mapq_range)

    # Prepare C>T filters
    filter_ct = None
    exclusively_ct = False
    subtract_ct = False
    if 'only_ct' in ct_checklist_tuple:
        filter_ct = True
    elif 'exclude_ct' in ct_checklist_tuple:
        filter_ct = False
    elif 'exclusively_ct' in ct_checklist_tuple:
        exclusively_ct = True
    if 'subtract_ct' in ct_checklist_tuple:
        subtract_ct = True

    exact_ct_changes = None if ct_count_value == 'any' else int(ct_count_value)

    file_data = next((f for f in store_data['files'] if f['filename'] == selected_file), None)
    if file_data is None:
        return None

    # Re-process the file with all filters applied
    read_lengths, file_format, deduped_file_path = load_and_process_file(
        file_data['content'], file_data['filename'], selected_nm, filter_ct,
        exclusively_ct, ct_count_value, ct_checklist_tuple, mapq_range_tuple
    )

    # set any additional adjustments if necessary
    if subtract_ct:
        read_lengths = adjust_nm_for_ct_changes(read_lengths, subtract_ct, ct_count_value)

    selected_cg_content = []  # Update if CG content selection is implemented

    output_file_path = os.path.join(CUSTOM_TEMP_DIR, f"selected_{uuid.uuid4()}_{selected_file}")

    # Export the selected reads with filters applied
    export_selected_reads(
        read_lengths, selected_lengths, selected_cg_content, output_file_path,
        deduped_file_path, file_format, selected_nm, filter_ct, exact_ct_changes
    )

    return dcc.send_file(output_file_path)














##########################################################################################################

##     ##    ###    #### ##    ##    ######## ##     ## ##    ##  ######  ######## ####  #######  ##    ## 
###   ###   ## ##    ##  ###   ##    ##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## 
#### ####  ##   ##   ##  ####  ##    ##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## 
## ### ## ##     ##  ##  ## ## ##    ######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ## 
##     ## #########  ##  ##  ####    ##       ##     ## ##  #### ##          ##     ##  ##     ## ##  #### 
##     ## ##     ##  ##  ##   ###    ##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### 
##     ## ##     ## #### ##    ##    ##        #######  ##    ##  ######     ##    ####  #######  ##    ##

##########################################################################################################

if __name__ == '__main__':
    # Define the path for your log file
    log_file_path = os.path.join(os.path.dirname(__file__), 'celery.log')

    def log_to_file(process, log_file_path):
        with open(log_file_path, 'w') as log_file:
            while True:
                output = process.stdout.readline()
                if process.poll() is not None and output == b'':
                    break
                if output:
                    log_file.write(output)
                    log_file.flush()

    # Launch the Celery worker
    if os.environ.get('WERKZEUG_RUN_MAIN') == 'true':
        worker_process = subprocess.Popen(
            [
                sys.executable, '-m', 'celery', '-A', 'tasks.celery_app', 'worker',
                '--loglevel=info'
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
            universal_newlines=True
        )

        # Start a thread to log output to a file
        log_thread = threading.Thread(target=log_to_file, args=(worker_process, log_file_path))
        log_thread.daemon = True
        log_thread.start()

        print(f"Started Celery worker with PID: {worker_process.pid}")

    else:
        worker_process = None

    try:
        app.run_server(debug=True)
    finally:
        if worker_process:
            print("Terminating Celery worker...")
            worker_process.terminate()
            worker_process.wait()
