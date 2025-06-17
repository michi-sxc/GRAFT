import plotly.graph_objects as go
import plotly.colors
import numpy as np
import pandas as pd
from config import colors
from .utils import apply_graft_theme, calculate_confidence_intervals # Use local CI calc

def create_damage_at_ends_figure(five_prime_pct, three_prime_pct):
    """Creates the bar chart for base frequencies at read ends."""
    fig = go.Figure()

    bases = ['A', 'C', 'G', 'T', 'N'] # Include N category
    five_freqs = [five_prime_pct.get(b, 0) for b in bases]
    three_freqs = [three_prime_pct.get(b, 0) for b in bases]

    # Filter out bases with 0 frequency in both ends to avoid clutter
    bases_to_plot = [b for i, b in enumerate(bases) if five_freqs[i] > 0 or three_freqs[i] > 0]
    five_freqs_plot = [five_prime_pct.get(b, 0) for b in bases_to_plot]
    three_freqs_plot = [three_prime_pct.get(b, 0) for b in bases_to_plot]

    if not bases_to_plot:
         fig.update_layout(title="Damage Patterns at Read Ends (No Data)")
         return apply_graft_theme(fig)

    # 5' end frequencies
    fig.add_trace(go.Bar(
        x=bases_to_plot,
        y=five_freqs_plot,
        name="5' Terminus", # Use terminus instead of end
        marker_color=colors['highlight'], # Use highlight color
        text=[f'{freq:.1f}%' for freq in five_freqs_plot],
        textposition='auto',
        hoverinfo='x+y+name'
    ))

    # 3' end frequencies
    fig.add_trace(go.Bar(
        x=bases_to_plot,
        y=three_freqs_plot,
        name="3' Terminus",
        marker_color=colors['line'], # Use line color
        text=[f'{freq:.1f}%' for freq in three_freqs_plot],
        textposition='auto',
        hoverinfo='x+y+name'
    ))

    apply_graft_theme(fig)
    fig.update_layout(
        title='Base Composition at Read Termini',
        xaxis_title='Base',
        yaxis_title='Frequency (%)',
        yaxis_range=[0, max(max(five_freqs_plot, default=0), max(three_freqs_plot, default=0)) * 1.1 + 5], # Dynamic range with padding
        barmode='group',
        bargap=0.15, # Standard gap
        bargroupgap=0.1,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig


def create_mismatch_frequency_plot(frequencies, sequencing_type):
    """Plots mismatch frequency vs distance from read end with confidence intervals."""
    fig = go.Figure()
    if not frequencies:
        fig.update_layout(title="Mismatch Frequency Data Not Available")
        return apply_graft_theme(fig)

    # Define plot properties for different mutation types
    plot_defs = {
        'C>T_CpG': {'color': colors['highlight'], 'symbol': 'circle', 'name': 'C>T (CpG)'},
        'C>T_nonCpG': {'color': colors['line'], 'symbol': 'square', 'name': 'C>T (non-CpG)'},
        'other': {'color': colors['highlight2'], 'symbol': 'diamond', 'name': 'Other Mismatch'}
    }

    max_freq_observed = 0.0 # To set y-axis range

    if sequencing_type == 'single':
        x_values = np.arange(len(frequencies.get('C>T_CpG', []))) # Positions 0 to N-1

        for mut_type, props in plot_defs.items():
            if mut_type not in frequencies: continue # Skip if data missing

            freq_data = frequencies[mut_type]
            
            counts_data = frequencies.get(f'{mut_type}_mismatches', []) 

            if mut_type == 'C>T_CpG':
                site_counts = frequencies.get('C_in_CpG_sites', [])   
            elif mut_type == 'C>T_nonCpG':
                site_counts = frequencies.get('C_not_in_CpG_sites', [])
            else: # 'other'
                site_counts = frequencies.get('total_considered_sites', []) 

            if len(freq_data) == 0 or len(counts_data) == 0 or len(site_counts) == 0: continue

            max_freq_observed = max(max_freq_observed, np.max(freq_data) if len(freq_data) > 0 else 0)

            # Calculate Confidence Intervals
            lower_ci, upper_ci = calculate_confidence_intervals(counts_data, site_counts)

            # Hover text
            hover_texts = [
                f"Pos: {p}<br>Type: {props['name']}<br>Freq: {f:.3%}<br>Count: {c}<br>Sites: {s}<br>95% CI: [{l:.3%}, {u:.3%}]"
                for p, f, c, s, l, u in zip(x_values, freq_data, counts_data, site_counts, lower_ci, upper_ci)
            ]

            # Add CI Area Trace (first, so it's behind lines)
            fig.add_trace(go.Scatter(
                x=np.concatenate([x_values, x_values[::-1]]), # x_values forward then backward
                y=np.concatenate([upper_ci, lower_ci[::-1]]), # upper bounds forward, lower bounds backward
                fill='toself',
                fillcolor=f'rgba({plotly.colors.hex_to_rgb(props["color"])[0]}, {plotly.colors.hex_to_rgb(props["color"])[1]}, {plotly.colors.hex_to_rgb(props["color"])[2]}, 0.2)', # Transparent fill
                line=dict(width=0), # No border line for fill area
                name=f"{props['name']} 95% CI",
                hoverinfo='none',
                showlegend=False
            ))

            # Add Main Line Trace
            fig.add_trace(go.Scatter(
                x=x_values,
                y=freq_data,
                name=props['name'],
                mode='lines+markers',
                line=dict(color=props['color'], width=2),
                marker=dict(symbol=props['symbol'], size=6, color=props['color']),
                text=hover_texts, # hover text
                hoverinfo='text' # Show only text on hover
            ))

    else: # Double-stranded
        for end_label, end_title in [('5_prime', "5' End"), ('3_prime', "3' End")]:
            if end_label not in frequencies: continue # Check if end_label (e.g. '5_prime') exists as key
            
            # frequencies[end_label] must be dictionary before trying to access sub-keys
            if not isinstance(frequencies[end_label], dict):
                #logger.warning(f"Expected a dictionary for frequencies['{end_label}'], but got {type(frequencies[end_label])}. Skipping.")
                continue

            x_values = np.arange(len(frequencies[end_label].get('C>T_CpG', []))) # Get length from a known frequency array

            for mut_type, props in plot_defs.items():
                if mut_type not in frequencies[end_label]: continue

                freq_data = frequencies[end_label][mut_type]
                counts_data = frequencies[end_label].get(f'{mut_type}_counts', []) # MATCHES KEY FROM STATS

                site_counts_arr = [] # Use a more generic name
                if mut_type == 'C>T_CpG':
                    site_counts_arr = frequencies[end_label].get('CpG_C_sites', [])       # MATCHES KEY
                elif mut_type == 'C>T_nonCpG':
                    site_counts_arr = frequencies[end_label].get('nonCpG_C_sites', [])   # MATCHES KEY
                else: # 'other'
                    site_counts_arr = frequencies[end_label].get('total_bases', [])      # MATCHES KEY
                
                # Renamed to avoid conflict if site_counts was already defined
                if len(freq_data) == 0 or len(counts_data) == 0 or len(site_counts_arr) == 0: continue

                max_freq_observed = max(max_freq_observed, np.max(freq_data) if len(freq_data) > 0 else 0)

                lower_ci, upper_ci = calculate_confidence_intervals(counts_data, site_counts_arr) # Use site_counts_arr

                trace_name = f"{props['name']} ({end_title})"
                hover_texts = [
                    f"Pos: {p} ({end_title})<br>Type: {props['name']}<br>Freq: {f:.3%}<br>Count: {c}<br>Sites: {s}<br>95% CI: [{l:.3%}, {u:.3%}]"
                    for p, f, c, s, l, u in zip(x_values, freq_data, counts_data, site_counts_arr, lower_ci, upper_ci) # Use site_counts_arr
                ]

                # Add CI Area Trace
                fig.add_trace(go.Scatter(
                    x=np.concatenate([x_values, x_values[::-1]]),
                    y=np.concatenate([upper_ci, lower_ci[::-1]]),
                    fill='toself',
                    fillcolor=f'rgba({plotly.colors.hex_to_rgb(props["color"])[0]}, {plotly.colors.hex_to_rgb(props["color"])[1]}, {plotly.colors.hex_to_rgb(props["color"])[2]}, 0.15)', # More transparent for overlap
                    line=dict(width=0),
                    name=f"{trace_name} 95% CI",
                    hoverinfo='none',
                    showlegend=False,
                    legendgroup=end_label # Group traces by end for legend interaction
                ))

                # Add Main Line Trace
                fig.add_trace(go.Scatter(
                    x=x_values,
                    y=freq_data,
                    name=trace_name,
                    mode='lines+markers',
                    line=dict(color=props['color'], width=2, dash=('dash' if end_label=='3_prime' else 'solid')), # Dashed line for 3' end
                    marker=dict(symbol=props['symbol'], size=6, color=props['color']),
                    text=hover_texts,
                    hoverinfo='text',
                    legendgroup=end_label
                ))

    apply_graft_theme(fig)
    fig.update_layout(
        title='Mismatch Frequency vs Distance from Read End',
        xaxis_title='Distance from Read End (bp)',
        yaxis_title='Mismatch Frequency',
        yaxis_range=[0, max(0.01, max_freq_observed * 1.1)], # Set range based on observed max, minimum 1%
        yaxis_tickformat='.1%', # Format y-axis as percentage with 1 decimal
        hovermode='closest', # Show hover for closest point
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig


def create_damage_pattern_heatmap(mismatch_details):
    """Creates a heatmap showing mismatch types along the read."""
    fig = go.Figure()
    if not mismatch_details:
        fig.update_layout(title="Damage Pattern Heatmap (No Mismatch Data)")
        return apply_graft_theme(fig)

    # Create DataFrame for easier manipulation
    df = pd.DataFrame(mismatch_details)
    if df.empty:
        fig.update_layout(title="Damage Pattern Heatmap (No Mismatch Data)")
        return apply_graft_theme(fig)

    # --- Determine max position ---
    # Consider query_alignment_length if available and relevant, otherwise use read_pos
    max_pos = df['read_pos'].max() if not df.empty else 0
    if max_pos == 0: # Handle case where all mismatches are at pos 0
         max_pos = 1

    # --- Define Damage Categories ---
    categories = ['C>T', 'G>A', 'Transition (Other)', 'Transversion', 'Other/N']
    category_map = {cat: i for i, cat in enumerate(categories)}
    heatmap_counts = np.zeros((len(categories), max_pos + 1), dtype=int)
    total_bases_at_pos = np.zeros(max_pos + 1, dtype=int) # Approx total bases considered at position

    # --- Populate Heatmap Data ---
    # Estimate total bases (tricky without original alignment, approximate by max count)
    # A better way requires passing the full read data
    df_agg = df.groupby('read_pos').size()
    for pos, count in df_agg.items():
         if pos <= max_pos:
             total_bases_at_pos[pos] = count # Approximation: total mismatches at pos

    # Count mismatches per category per position
    for _, row in df.iterrows():
        pos = row.get('read_pos', -1)
        if pos < 0 or pos > max_pos: continue

        ref = row.get('ref_base', 'N').upper()
        read = row.get('read_base', 'N').upper()
        is_rev = row.get('is_reverse', False)

        cat_key = 'Other/N' # Default
        if ref != 'N' and read != 'N' and ref != read:
            mut = f"{ref}>{read}"
            is_transition = (ref in 'AG' and read in 'AG') or (ref in 'CT' and read in 'CT')

            if is_rev: # Biological change on reverse strand
                if mut == 'G>A': cat_key = 'C>T' # Treat G>A on rev as C>T damage
                elif mut == 'C>T': cat_key = 'Transition (Other)' # Treat C>T on rev as G>A (non-damage transition)
                elif is_transition: cat_key = 'Transition (Other)'
                else: cat_key = 'Transversion'
            else: # Forward strand
                if mut == 'C>T': cat_key = 'C>T'
                elif mut == 'G>A': cat_key = 'G>A' # Keep G>A separate
                elif is_transition: cat_key = 'Transition (Other)'
                else: cat_key = 'Transversion'

        if cat_key in category_map:
             heatmap_counts[category_map[cat_key], pos] += 1


    # Calculate Frequencies (handle division by zero)
    # Note: Denominator choice is important. Could be used as total mismatches at pos,
    # or total bases sequenced at pos? Using total mismatches here.
    total_mismatches_at_pos_safe = np.maximum(1, total_bases_at_pos) # Avoid div zero
    heatmap_freq = heatmap_counts / total_mismatches_at_pos_safe[:, np.newaxis].T # Divide each column by total for that pos


    # --- Create Heatmap ---
    hover_texts = []
    for r_idx, cat in enumerate(categories):
        row_texts = []
        for c_idx in range(max_pos + 1):
            count = heatmap_counts[r_idx, c_idx]
            total = total_bases_at_pos[c_idx]
            freq = heatmap_freq[r_idx, c_idx]
            row_texts.append(f"Pos: {c_idx}<br>Type: {cat}<br>Count: {count}<br>Total Mismatches at Pos: {total}<br>Freq: {freq:.2%}")
        hover_texts.append(row_texts)


    fig = go.Figure(data=go.Heatmap(
        z=heatmap_freq,
        x=np.arange(max_pos + 1),
        y=categories,
        colorscale='Viridis', #'YlOrRd' might be better for damage focus
        colorbar=dict(title='Frequency', tickformat='.1%'),
        text=hover_texts,
        hoverongaps=False,
        hoverinfo='text'
    ))

    apply_graft_theme(fig)
    fig.update_layout(
        title='Mismatch Type Frequency Along Read Position',
        xaxis_title='Position in Read (bp)',
        yaxis_title='Mismatch Type',
        yaxis=dict(autorange='reversed'), # Put C>T at top if desired
        xaxis=dict(dtick=10 if max_pos > 50 else 5) # Adjust tick frequency based on length
    )
    return fig
