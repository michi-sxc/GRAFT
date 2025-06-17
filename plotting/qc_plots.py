import plotly.graph_objects as go
from config import colors
from .utils import apply_graft_theme

def create_mismatch_type_bar_chart(mismatch_counts_by_type):
    """Creates a bar chart showing the counts of each specific mismatch type."""
    fig = go.Figure()
    if not mismatch_counts_by_type:
        fig.update_layout(title="Mismatch Type Distribution (No Data)")
        return apply_graft_theme(fig)

    # Sort by ref base then mismatch type for consistent order
    sorted_items = sorted(mismatch_counts_by_type.items(), key=lambda item: (item[0][0], item[0][2]))
    mutations = [item[0] for item in sorted_items]
    counts = [item[1] for item in sorted_items]

    # Filter out zero counts to avoid clutter, unless all are zero
    if any(c > 0 for c in counts):
        mutations_nz = [m for m, c in zip(mutations, counts) if c > 0]
        counts_nz = [c for c in counts if c > 0]
    else:
        mutations_nz = mutations
        counts_nz = counts


    # Assign colors based on reference base (my implementation could be changed)
    color_discrete_map = {
        'A': colors['line'], 'C': colors['highlight'],
        'G': colors['secondary'], 'T': colors['marker']
    }
    bar_colors = [color_discrete_map.get(mut[0], colors['muted']) for mut in mutations_nz]

    fig.add_trace(go.Bar(
        x=mutations_nz,
        y=counts_nz,
        text=counts_nz, # Count on bar
        textposition='auto',
        marker_color=bar_colors,
        hoverinfo='x+y'
    ))

    apply_graft_theme(fig)
    fig.update_layout(
        title='Mismatch Type Distribution',
        xaxis_title='Mismatch Type (Reference > Read)',
        yaxis_title='Count',
        xaxis_tickangle=-45, # Angle ticks for readability
    )
    return fig


def create_mismatch_pie_chart(stats_data):
    """Creates a detailed SNP distribution pie chart from alignment stats."""
    fig = go.Figure()
    mismatch_counts = get_mismatch_counts_by_type(stats_data.get('Mismatch Details', []))

    # Using dummy stats_data for standalone example if not provided
    if stats_data.get('File Format', 'N/A') == 'N/A': # Simplified check
        stats_data['File Format'] = 'bam/sam'


    if not mismatch_counts or stats_data.get('File Format', 'N/A') != 'bam/sam':
        fig.update_layout(title_text="SNP Distribution (No BAM/SAM Data)")
        return apply_graft_theme(fig)

    labels = [m for m, c in mismatch_counts.items() if c > 0]
    values = [c for c in mismatch_counts.values() if c > 0]

    if not labels:
        fig.update_layout(title_text="SNP Distribution (No Mismatches Found)")
        return apply_graft_theme(fig)

    ref_base_colors = {
        'A': '#1f77b4', 'C': '#ff7f0e', 'G': '#2ca02c', 'T': '#d62728'
    }
    pie_colors = [ref_base_colors.get(label[0], '#8c564b') for label in labels]

    fig.add_trace(go.Pie(
        labels=labels,
        values=values,
        hole=.4,
        marker=dict(colors=pie_colors, line=dict(color=colors['surface'], width=1)),
        textinfo='label+percent',
        insidetextorientation='radial',
        hoverinfo='label+value+percent',
        sort=False,
        pull=[0.02] * len(labels),
        # Adjust domain: pie uses bottom 90% of plot area's height,
        # leaving 10% space at the top of the plot area. Fix is pretty ugly.
        domain=dict(x=[0, 0.65], y=[0, 0.90]) 
    ))

    apply_graft_theme(fig) # Apply theme first

    # Then, explicitly set layout properties to override/ensure desired state
    fig.update_layout(
        title=dict(
            text="",
            y=0.97, 
            x=0.5,   
            xanchor='center',
            yanchor='top',
            font=dict(
                size=18,
                color=colors['on_surface']
            ),
            pad=dict(t=0, b=15) # Padding below the title text (b=bottom)
        ),
        showlegend=True,
        legend=dict(
            font=dict(color=colors['on_surface'], size=11),
            bgcolor='rgba(0,0,0,0)',
            bordercolor='rgba(0,0,0,0)',
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=0.68,
            tracegroupgap=5,
            itemsizing='constant'
        ),

        margin=dict(t=80, b=20, l=20, r=20), 
        paper_bgcolor=colors['stats_tab_color'],
        plot_bgcolor=colors['stats_tab_color'],
        height=400, 
    )
    fig.update_traces(textfont_color=colors['on_surface'], textfont_size=10)

    return fig

# Helper function needed by create_mismatch_pie_chart
def get_mismatch_counts_by_type(mismatch_details):
    """Counts occurrences of each specific mismatch type (A>C, A>G, etc.)."""
    counts = {}
    bases = ['A', 'C', 'G', 'T']
    for ref in bases:
        for read_b in bases:
            if ref != read_b:
                counts[f"{ref}>{read_b}"] = 0

    for mismatch in mismatch_details:
        ref = mismatch.get('ref_base', 'N').upper()
        read_b = mismatch.get('read_base', 'N').upper()
        key = f"{ref}>{read_b}"
        if key in counts:
            counts[key] += 1
    return counts
