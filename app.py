import pysam
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

#### START SETUP  ###################################################################################################
#####################################################################################################################

try:
    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f)
except FileNotFoundError:
    config = {}

default_centrifuge_db_path = config.get('centrifuge_db_path', '')

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SOLAR, "/assets/custom.css"], suppress_callback_exceptions=True) # general dash theme, can be changed but colors need to be adjusted accordingly (if you want a pretty dashboard)
celery_app = Celery('tasks', broker='redis://localhost:6379/0')

# custom temp directory for storing the data (needs better implementation, but works fine)
CUSTOM_TEMP_DIR = tempfile.mkdtemp()
os.makedirs(CUSTOM_TEMP_DIR, exist_ok=True)

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






#### HELPER FUNCTIONS ###############################################################################################
#####################################################################################################################

def rgb_to_hex(rgb):
    """Convert an RGB tuple to a hex string, rounding the values to integers - small unnecessary helper function"""
    return '#{:02x}{:02x}{:02x}'.format(int(round(rgb[0])), int(round(rgb[1])), int(round(rgb[2])))

def prepare_histogram_data(read_lengths):
    lengths = [length for length, _, _, _, _, _ in read_lengths]
    bins = np.arange(min(lengths), max(lengths) + 1, 1)
    hist, bin_edges = np.histogram(lengths, bins=bins)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    return bin_centers, hist, bins






#### FILE FUNCTIONS #################################################################################################
#####################################################################################################################

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


# cached results of file-processing to avoid reloading the data, but not strictly necessary
@lru_cache(maxsize=100)
def load_and_process_file(file_content, filename, selected_nm, filter_ct, ct_count_value, ct_checklist_tuple, mapq_range):
    temp_file_path, file_format = process_uploaded_file(file_content, filename)
    deduped_file_path = os.path.join(
        CUSTOM_TEMP_DIR, f"deduped_{uuid.uuid4()}_{filename}"
    )
    read_lengths = remove_duplicates_with_strand_check(
        temp_file_path, deduped_file_path, file_format
    )

    if file_format in ['bam', 'sam']:
        # Apply MAPQ filter only for BAM/SAM files
        min_mapq, max_mapq = mapq_range
        read_lengths = [
            r for r in read_lengths
            if r[5] is not None and min_mapq <= r[5] <= max_mapq
        ]
        
        if filter_ct is not None:
            read_lengths = [r for r in read_lengths if ('C>T' in get_deamination_pattern(r[4]) or 'G>A' in get_deamination_pattern(r[4])) == filter_ct]

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


def remove_duplicates_with_strand_check(input_file, output_file, file_format):
    if file_format in ['bam', 'sam']:
        bamfile = pysam.AlignmentFile(input_file, "rb" if file_format == 'bam' else "r", check_sq=False)
        outputfile = pysam.AlignmentFile(output_file, "wb" if file_format == 'bam' else 'w', header=bamfile.header)
    else:
        outputfile = open(output_file, "w")

    unique_reads = set()
    read_lengths = []

    if file_format in ['bam', 'sam']:
        for read in bamfile:
            # Remove the check for read.is_unmapped
            sequence = read.query_sequence
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
        # Unchanged for FASTA/FASTQ
        for record in SeqIO.parse(input_file, file_format):
            sequence = str(record.seq)
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





#### PLOT FUNCTIONS #################################################################################################
#####################################################################################################################

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
        marker_color=colors['line']
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
        'T': plotly.colors.sequential.Reds_r # otherwise not really visible
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


def calculate_mismatch_frequency(read_lengths, file_format):
    if file_format not in ['bam', 'sam']:
        return None  # Mismatch frequency is not applicable

    # Initialize arrays to store mismatch counts at each distance
    max_distance = 30  # Define maximum distance to consider, can be changed
    counts = {'C>T': [0] * max_distance, 'G>A': [0] * max_distance, 'other': [0] * max_distance}
    total_reads = len(read_lengths)

    if total_reads == 0:
        return None  # No reads to process

    # Iterate over each read to count mismatches by position
    for _, _, _, _, read, _ in read_lengths:
        mismatches = get_mismatch_pattern(read)
        if not mismatches:
            continue  
        for mismatch_type, count in mismatches.items():
            if mismatch_type in counts:
                for position in range(min(len(counts[mismatch_type]), count)):
                    counts[mismatch_type][position] += 1
            else:
                for position in range(min(len(counts['other']), count)):
                    counts['other'][position] += 1

    if sum(sum(v) for v in counts.values()) == 0:
        return None  # No mismatches found

    # Normalize counts to calculate frequency
    frequencies = {key: [c / total_reads for c in counts[key]] for key in counts}
    return frequencies


def create_mismatch_frequency_plot(frequencies):
    if not frequencies:
        fig = go.Figure()
        fig.update_layout(
            title_text="Mismatch Frequency Data Not Available",
            paper_bgcolor=colors['plot_bg'],
            font=dict(color=colors['muted'])
        )
        return fig
    
    fig = go.Figure()

    # Define the distance range
    x_values = list(range(len(frequencies['C>T'])))

    # traces for each type, not really consistent with the rest of the styling theme
    fig.add_trace(go.Scatter(
        x=x_values,
        y=frequencies['C>T'],
        mode='lines+markers',
        name='C>T',
        line=dict(color=colors['line'], width=2, dash='solid'),
        marker=dict(size=6, color=colors['line'], symbol='circle')
    ))

    fig.add_trace(go.Scatter(
        x=x_values,
        y=frequencies['G>A'],
        mode='lines+markers',
        name='G>A',
        line=dict(color=colors['highlight'], width=2, dash='dash'),
        marker=dict(size=6, color=colors['highlight'], symbol='square')
    ))

    fig.add_trace(go.Scatter(
        x=x_values,
        y=frequencies['other'],
        mode='lines+markers',
        name='Other',
        line=dict(color=colors['highlight2'], width=2, dash='dot'),
        marker=dict(size=6, color=colors['highlight2'], symbol='diamond')
    ))

    
    fig.update_layout(
        title=dict(
            text='Mismatch Frequency vs Distance from Read End',
            font=dict(size=20, color=colors['muted']),
            x=0.5,
            y=0.95
        ),
        xaxis=dict(
            title='Distance from Read End (bp)',
            color=colors['muted'],
            tickfont=dict(size=12, color=colors['muted']),
            gridcolor=colors['grid'],
            showline=True,
            linewidth=1,
            linecolor=colors['primary'],
            zeroline=False,
        ),
        yaxis=dict(
            title='Mismatch Frequency',
            color=colors['primary'],
            tickfont=dict(size=12, color=colors['muted']),
            gridcolor=colors['grid'],
            showline=True,
            linewidth=1,
            linecolor=colors['primary'],
            zeroline=False,
            rangemode='tozero',
        ),
        legend=dict(
            x=0.85,
            y=1.15,
            bgcolor='rgba(0, 0, 0, 0.5)',
            bordercolor='rgba(255, 255, 255, 0.3)',
            borderwidth=1,
            font=dict(color=colors['muted'])
        ),
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        dragmode='select',
        selectdirection='h',
        newselection_line_width=2,
        newselection_line_color=colors['secondary'],
        newselection_line_dash='dashdot',
        margin=dict(t=50, b=50, l=50, r=50)
    )

    
    fig.update_traces(hovertemplate='Distance: %{x} bp<br>Frequency: %{y:.4f}<extra></extra>')

    return fig

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

def export_selected_reads(read_lengths, selected_lengths, selected_cg_content, output_file_path, original_file_path, file_format, selected_nm=None, filter_ct=None, exact_ct_changes=None):
    if file_format in ['bam', 'sam']:
        with pysam.AlignmentFile(original_file_path, "rb" if file_format == 'bam' else "r", check_sq=False) as bamfile:
            header = bamfile.header.to_dict()

        # Retain all references in the header to avoid mismatch errors
        with pysam.AlignmentFile(output_file_path, "wb" if file_format == 'bam' else "w", header=header) as outputfile:
            for length, cg, nm, seq, read, mapq in read_lengths:
                deamination_pattern = get_deamination_pattern(read)
                ct_changes = count_ct_changes(read)

                if filter_ct is not None:
                    has_ct = 'C>T' in deamination_pattern or 'G>A' in deamination_pattern
                    if filter_ct and not has_ct:
                        continue
                    elif not filter_ct and has_ct:
                        continue

                if selected_nm is not None and selected_nm != 'all' and nm != selected_nm:
                    continue

                if exact_ct_changes is not None and ct_changes != exact_ct_changes:
                    continue

                if (selected_lengths and length not in selected_lengths) or (selected_cg_content and cg not in selected_cg_content):
                    continue

                outputfile.write(read)
    else:
        with open(output_file_path, "w") as outputfile:
            for length, cg, _, seq, record in read_lengths:
                if (selected_lengths and length not in selected_lengths) or (selected_cg_content and cg not in selected_cg_content):
                    continue
                SeqIO.write(record, outputfile, file_format)

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

def get_mismatch_frequency_data(frequencies):
    if not frequencies:
        return "Mismatch frequency data not available for this file format.\n"
    text_output = "Mismatch Frequency vs Distance from Read End:\n"
    max_distance = len(frequencies['C>T'])
    text_output += "Position\tC>T Frequency\tG>A Frequency\tOther Frequency\n"
    for i in range(max_distance):
        ct_freq = frequencies['C>T'][i]
        ga_freq = frequencies['G>A'][i]
        other_freq = frequencies['other'][i]
        text_output += f"{i}\t{ct_freq:.4f}\t{ga_freq:.4f}\t{other_freq:.4f}\n"
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








#### DASHBOARD LAYOUT ###############################################################################################
#####################################################################################################################

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

        # Apply and Reset Buttons
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


# Define the Navbar separately for reuse 
navbar = dbc.Navbar(
    dbc.Container([
        html.A(
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
                dbc.NavItem(dbc.NavLink("Home", href="/")),
                dbc.NavItem(dbc.NavLink("Convert", href="/convert")),
                dbc.NavItem(dbc.Button("Settings", id="settings-button", color="light", outline=True)),
            ], className="ms-auto", navbar=True),  # Align to the right
            id="navbar-collapse",
            navbar=True,
        ),
    ]),
    color=colors['surface'],
    dark=True,
    className="mb-4",
)

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    navbar,
    html.Div(id='page-content'),
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
                            {'label': ' Show only C>T changes', 'value': 'only_ct'},
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
                    # dbc.Button("Apply Filters (?)", id="apply-filters-button", color="primary", className="w-100 mb-3"),
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
            ], className="mb-4"),

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
            ]),
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

    # Stores and Download component
    dcc.Store(id='file-store', data={'files': []}),
    dcc.Store(id='selected-lengths', data=[]),
    dcc.Interval(id='centrifuge-interval', interval=5000, n_intervals=0, disabled=True),
    dcc.Store(id='centrifuge-task-id'),
    dcc.Store(id='centrifuge-output-path'),
    dcc.Download(id="download-file"),

], fluid=True, className="px-4 py-3", style={"backgroundColor": colors['background']})

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





##### CALLBACK FUNCTIONS ############################################################################################
#####################################################################################################################

@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])

def display_page(pathname):
    if pathname == '/convert':
        return layout_convert
    else:
        return layout_main
    

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

### CENTRIFUGE IMPLEMENTATION

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


@celery_app.task
def run_centrifuge_task(file_path, output_path, report_path, db_path):
    centrifuge_cmd = [
        "centrifuge",
        "-x", db_path,
        "-U", file_path,
        "-S", output_path,
        "--report-file", report_path
    ]
    try:
        subprocess.run(centrifuge_cmd, check=True)
    except subprocess.CalledProcessError as e:
        raise e


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



def update_centrifuge_output(n_intervals, task_id, centrifuge_report_path): # visual error, n_intervals is being accessed
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
    


######


@app.callback(
    Output('stats-offcanvas', 'is_open'),
    [Input('stats-button', 'n_clicks')],
    [State('stats-offcanvas', 'is_open')]
)
def toggle_offcanvas(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open

@app.callback(
    Output('ct-checklist', 'value'),
    Input('ct-checklist', 'value')
)
def update_ct_checklist(ct_value):
    if len(ct_value) > 1:
        return ct_value[-1:]
    return ct_value

@app.callback(
    [Output('file-store', 'data'),
     Output('file-selector', 'options'),
     Output('file-selector', 'value')],
    [Input('upload-data', 'contents'),
     Input('upload-data', 'filename')],
    [State('file-store', 'data'),
     State('file-selector', 'value')]
)

def update_file_store(contents, filenames, store_data, current_selection):
    if contents is None:
        return store_data, [], current_selection

    files = store_data['files']
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
    
    # Automatically select the first file if no file is currently selected
    if current_selection is None and new_file_added:
        current_selection = files[0]['filename']

    return store_data, options, current_selection


@app.callback(
    [Output('read-length-histogram', 'figure'),
     Output('overall-cg-histogram', 'figure'),
     Output('cg-content-histogram', 'figure'),
     Output('damage-patterns', 'figure'),
     Output('selected-lengths', 'data'),
     Output('mismatch-frequency-plot', 'figure'),
     Output('data-summary', 'children'),
     Output('mapq-histogram', 'figure'),
     Output('mismatch-type-bar-chart', 'figure'),
     Output('damage-pattern-plot', 'figure')],
    [Input('file-selector', 'value'),
     Input('file-store', 'data'),
     Input('mapq-range', 'value'),
     Input('nm-dropdown', 'value'),
     Input('ct-checklist', 'value'),
     Input('ct-count-dropdown', 'value'),
     Input('min-base-quality', 'value'),
     Input('read-length-histogram', 'selectedData'),
     Input('read-length-histogram', 'clickData'),
     Input('clear-selection-button', 'n_clicks')],
    [State('selected-lengths', 'data'),
     State('read-length-histogram', 'figure')]
)
def update_histograms(selected_file, store_data, mapq_range, selected_nm, ct_checklist, ct_count_value,
                      min_base_quality,
                      read_length_selectedData, read_length_clickData, clear_selection_nclicks,
                      current_selected_lengths, current_fig):

    # Check if no file is selected
    if selected_file is None:
        empty_figure = go.Figure()
        empty_text = ""
        return empty_figure, empty_figure, empty_figure, empty_figure, [], empty_figure, empty_text, empty_figure, empty_figure, empty_figure

    # If clear-selection button is clicked, reset selection
    if clear_selection_nclicks:
        current_selected_lengths = []
        read_length_selectedData = None
        read_length_clickData = None

    # Convert ct_checklist to tuple
    ct_checklist_tuple = tuple(ct_checklist)

    # Prepare C>T filters
    filter_ct = None
    subtract_ct = False
    if 'only_ct' in ct_checklist_tuple:
        filter_ct = True
    elif 'exclude_ct' in ct_checklist_tuple:
        filter_ct = False
    if 'subtract_ct' in ct_checklist_tuple:
        subtract_ct = True

    # Convert mapq_range to tuple
    mapq_range_tuple = tuple(mapq_range)

    file_data = next((f for f in store_data['files'] if f['filename'] == selected_file), None)
    if file_data is None:
        empty_figure = go.Figure()
        empty_text = ""
        return empty_figure, empty_figure, empty_figure, empty_figure, [], empty_figure, empty_text, empty_figure, empty_figure, empty_figure

    # Load and process file content
    read_lengths, file_format, deduped_file_path = load_and_process_file(
        file_data['content'], file_data['filename'], selected_nm, filter_ct, ct_count_value, ct_checklist_tuple, mapq_range_tuple
    )

    # Handle mismatch frequency based on file format
    if file_format in ['bam', 'sam']:
        frequencies = calculate_mismatch_frequency(read_lengths, file_format)
        mismatch_freq_fig = create_mismatch_frequency_plot(frequencies)
    else:
        frequencies = None
        mismatch_freq_fig = go.Figure().update_layout(
            title='Mismatch Frequency Not Available for This File Format',
            xaxis=dict(title='', color=colors['muted']),
            yaxis=dict(title='', color=colors['muted']),
            paper_bgcolor=colors['plot_bg'],
            plot_bgcolor=colors['plot_bg'],
            font=dict(color=colors['muted'])
        )

    
    bin_centers, hist, bins = prepare_histogram_data(read_lengths)

    
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

    # selection data for CG content hist
    selected_lengths = current_selected_lengths
    if read_length_selectedData:
        selected_x_values = [point['x'] + 0.5 for point in read_length_selectedData['points']]
        selected_lengths = [int(x) for x in selected_x_values]

    if read_length_clickData:
        clicked_length = int(read_length_clickData['points'][0]['x'] + 0.5)
        selected_lengths = [clicked_length]

    cg_contents_for_length = [cg for length, cg, _, _, _, _ in read_lengths if length in selected_lengths]
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

    # Apply filters to mismatches
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

        # plotting
        mismatch_type_fig = create_mismatch_type_bar_chart(filtered_mismatch_counts)

        
        damage_pattern_fig = create_damage_pattern_plot(filtered_mismatches)
    else:
        mismatch_type_fig = go.Figure()
        damage_pattern_fig = go.Figure()

    # data summary after generating all plots
    read_length_text = get_read_length_data(read_lengths)
    cg_content_text = get_overall_cg_content_data(read_lengths)
    damage_patterns_text = get_damage_patterns_data(five_prime_end, three_prime_end)
    mismatch_frequency_text = get_mismatch_frequency_data(frequencies)
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

    # combined data
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

    # main return
    return (
        read_length_fig, overall_cg_fig, cg_content_fig, damage_fig,
        selected_lengths, mismatch_freq_fig, data_summary_text, mapq_fig,
        mismatch_type_fig, damage_pattern_fig
    )



    


@app.callback(
    Output("download-file", "data"),
    [Input("export-button", "n_clicks"),
     State('read-length-histogram', 'selectedData'),
     State('cg-content-histogram', 'clickData'),
     State('file-selector', 'value'),
     State('file-store', 'data'),
     State('nm-dropdown', 'value'),
     State('ct-checklist', 'value'),
     State('ct-count-dropdown', 'value')],
    prevent_initial_call=True
)
def export_reads(n_clicks, read_length_selectedData, cg_clickData, selected_file, store_data, selected_nm, ct_checklist, ct_count_value):
    if selected_file is None:
        return None

    filter_ct = None
    subtract_ct = False
    if 'only_ct' in ct_checklist:
        filter_ct = True
    elif 'exclude_ct' in ct_checklist:
        filter_ct = False
    elif 'subtract_ct' in ct_checklist:
        subtract_ct = True

    exact_ct_changes = None if ct_count_value == 'any' else int(ct_count_value)

    file_data = next((f for f in store_data['files'] if f['filename'] == selected_file), None)
    if file_data is None:
        return None

    temp_file_path, file_format = process_uploaded_file(file_data['content'], file_data['filename'])

    deduped_file_path = os.path.join(CUSTOM_TEMP_DIR, f"deduped_{uuid.uuid4()}_{selected_file}")
    read_lengths = remove_duplicates_with_strand_check(temp_file_path, deduped_file_path, file_format)

    if subtract_ct:
        read_lengths = adjust_nm_for_ct_changes(read_lengths, subtract_ct)

    selected_lengths = []
    selected_cg_content = []

    if read_length_selectedData is not None:
        selected_x_values = [point['x'] + 0.5 for point in read_length_selectedData['points']]
        selected_lengths = [int(x) for x in selected_x_values]

    if cg_clickData is not None:
        selected_cg_content.append(cg_clickData['points'][0]['x'])

    output_file_path = os.path.join(CUSTOM_TEMP_DIR, f"selected_{uuid.uuid4()}_{selected_file}")
    export_selected_reads(read_lengths, selected_lengths, selected_cg_content, output_file_path, deduped_file_path, file_format, selected_nm, filter_ct, exact_ct_changes)

    return dcc.send_file(output_file_path)

if __name__ == '__main__':
    app.run_server(debug=True)