import plotly.graph_objects as go
import numpy as np
import pandas as pd
from config import colors
from .utils import apply_graft_theme


def create_read_length_histogram(processed_reads_data_list_of_dicts, selected_lengths=None):
    fig = go.Figure()
    if not processed_reads_data_list_of_dicts:
        fig.update_layout(title="No data available for Read Length Distribution")
        return apply_graft_theme(fig)

    lengths = [read_dict['length'] for read_dict in processed_reads_data_list_of_dicts]
    if not lengths:
        fig.update_layout(title="No valid read lengths found")
        return apply_graft_theme(fig)

    min_len = min(lengths) if lengths else 0
    max_len = max(lengths) if lengths else 1
    bin_size = 1
    bins = np.arange(min_len - 0.5 * bin_size, max_len + 1.5 * bin_size, bin_size)
    if len(bins) < 2:
        bins = np.array([min_len - 0.5, min_len + 0.5, min_len + 1.5]) if lengths else np.array([0.5, 1.5])

    hist, bin_edges = np.histogram(lengths, bins=bins)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    df = pd.DataFrame(processed_reads_data_list_of_dicts)
    if df.empty or 'length' not in df.columns or 'cg' not in df.columns:
        # If no CG data, create a series of NaNs for all bin_centers to ensure no line is drawn
        mean_cg_per_bin = pd.Series(data=[np.nan]*len(bin_centers), index=bin_centers)
    else:
        bin_labels = bin_centers if len(np.unique(bin_centers)) == len(bin_centers) else False
        df['bin'] = pd.cut(df['length'], bins=bin_edges, labels=bin_labels, right=False, include_lowest=True)
        
        # Group by bin and calculate mean CG
        grouped_cg = df.groupby('bin', observed=True)['cg'].mean()
        
        # Reindex to include all bin_centers. For bins with no data, this will introduce NaNs.
        mean_cg_per_bin = grouped_cg.reindex(bin_centers)
        # No .fillna(0). If a bin_center had no reads, mean_cg_per_bin[bin_center] will be NaN.

    fig.add_trace(go.Bar(
        x=bin_centers,
        y=hist,
        name='Read Count',
        marker_color=colors['line'],
        opacity=0.8,
        hoverinfo='x+y'
    ))

    fig.add_trace(go.Scatter(
        x=mean_cg_per_bin.index,
        y=mean_cg_per_bin.values, # NaNs for gaps
        name='Average CG Content',
        mode='lines+markers',
        yaxis='y2',
        line=dict(color=colors['highlight'], width=2),
        marker=dict(color=colors['highlight'], size=6, symbol='circle-open'),
        hoverinfo='x+y',
        connectgaps=False # leave NaN values open
    ))

    apply_graft_theme(fig)
    fig.update_layout(
        title='Read Length Distribution and Average CG Content',
        xaxis_title='Read Length (bp)',
        yaxis=dict(
            title='Read Count',
            titlefont_color=colors['line'],
            tickfont_color=colors['muted']
        ),
        yaxis2=dict(
            title='Average CG Content',
            titlefont_color=colors['highlight'],
            tickfont_color=colors['muted'],
            overlaying='y',
            side='right',
            range=[0, 1],
            showgrid=False,
        ),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        barmode='overlay',
        clickmode='event+select',
        dragmode='select',
    )

    if selected_lengths:
        selected_indices = [i for i, center in enumerate(bin_centers) if int(center) in selected_lengths]
        colors_bar = [colors['secondary'] if i in selected_indices else colors['line'] for i in range(len(bin_centers))]
        opacities_bar = [1.0 if i in selected_indices else 0.7 for i in range(len(bin_centers))]
        fig.data[0].marker.color = colors_bar
        fig.data[0].marker.opacity = opacities_bar

    return fig


def create_overall_cg_histogram(processed_reads_data_list_of_dicts):
    """Creates the histogram for overall CG content distribution."""
    fig = go.Figure()
    if not processed_reads_data_list_of_dicts:
        fig.update_layout(title="No data available for Overall CG Content")
        return apply_graft_theme(fig)


    cg_contents = [read_dict['cg'] for read_dict in processed_reads_data_list_of_dicts if read_dict.get('cg') is not None]
    if not cg_contents:
         fig.update_layout(title="No valid CG content found")
         return apply_graft_theme(fig)

    fig.add_trace(go.Histogram(
        x=cg_contents,
        name='CG Content',
        nbinsx=30,
        marker_color=colors['marker'],
        opacity=0.8
    ))

    apply_graft_theme(fig)
    fig.update_layout(
        title='Overall CG Content Distribution',
        xaxis_title='CG Content', # Keep as proportion, format as %
        yaxis_title='Frequency',
        xaxis_tickformat='.0%', # Format x-axis as percentage
        bargap=0.1
    )
    return fig


def create_cg_content_histogram_selected(processed_reads_data_list_of_dicts, selected_lengths):
    """Creates CG content histogram specifically for selected read lengths."""
    fig = go.Figure()
    if not processed_reads_data_list_of_dicts or not selected_lengths:
        len_str = "No selection" if not selected_lengths else f"{len(selected_lengths)} lengths"
        fig.update_layout(title=f"CG Content for Selected Reads ({len_str}) - No Data")
        return apply_graft_theme(fig)

    selected_lengths_set = set(selected_lengths)

    cg_contents_selected = [
        read_dict['cg'] for read_dict in processed_reads_data_list_of_dicts
        if read_dict.get('length') in selected_lengths_set and read_dict.get('cg') is not None
    ]

    if not cg_contents_selected:
        len_str = ', '.join(map(str, sorted(list(selected_lengths_set))[:5]))
        if len(selected_lengths_set) > 5: len_str += "..."
        fig.update_layout(title=f"CG Content for Length(s) {len_str} - No Reads Found")
        return apply_graft_theme(fig)

    if len(selected_lengths_set) == 1:
         title = f"CG Content Distribution for Read Length: {list(selected_lengths_set)[0]} bp"
    elif len(selected_lengths_set) < 10:
        len_str = ', '.join(map(str, sorted(list(selected_lengths_set))))
        title = f"CG Content Distribution for Lengths: {len_str}"
    else:
         title = f"CG Content Distribution for {len(selected_lengths_set)} Selected Lengths"

    fig.add_trace(go.Histogram(
        x=cg_contents_selected,
        name='CG Content',
        nbinsx=25,
        marker_color=colors['highlight2'],
        opacity=0.8
    ))

    apply_graft_theme(fig)
    fig.update_layout(
        title=title,
        xaxis_title='CG Content',
        yaxis_title='Frequency',
        xaxis_tickformat='.0%',
        bargap=0.1
    )
    return fig


def create_mapq_histogram(mapq_scores):
    """Creates the histogram for Mapping Quality (MAPQ) scores."""
    fig = go.Figure()
    if mapq_scores is None or len(mapq_scores) == 0:
        fig.update_layout(title="MAPQ Distribution (No Data Available)")
        return apply_graft_theme(fig)

    fig.add_trace(go.Histogram(
        x=mapq_scores,
        name='MAPQ',

        xbins=dict(start=0, end=max(mapq_scores) + 1 if mapq_scores else 1, size=5), # Bins of size 5
        marker_color=colors['secondary'],
        opacity=0.8
    ))

    apply_graft_theme(fig)
    fig.update_layout(
        title='Mapping Quality (MAPQ) Score Distribution',
        xaxis_title='MAPQ Score',
        yaxis_title='Frequency',
        bargap=0.1
    )
    return fig
