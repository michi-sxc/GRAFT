import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
from scipy import stats 
from config import colors 

# --- Theming ---

def apply_graft_theme(fig):
    """Applies the standard GRAFT dark theme to a Plotly figure."""
    fig.update_layout(
        template='plotly_dark', # Base dark theme
        paper_bgcolor=colors['plot_bg'],
        plot_bgcolor=colors['plot_bg'],
        font=dict(color=colors['muted']),
        # Default modebar buttons
        modebar=dict(
            add=['downloadSVG'],
            remove=['lasso2d', 'select2d']
        )
    )
    # Customize axes colors maybe?
    fig.update_xaxes(gridcolor=colors['grid'], linecolor=colors['muted'], zerolinecolor=colors['grid'])
    fig.update_yaxes(gridcolor=colors['grid'], linecolor=colors['muted'], zerolinecolor=colors['grid'])
    return fig

# --- Plot Export ---

def create_light_theme_figure(fig_dark):
    """Creates a light-themed copy of a Plotly figure for export."""
    if not fig_dark: return go.Figure() # Handle empty input

    # use copy to avoid modifying the original dark figure
    fig_light = go.Figure(fig_dark)

    fig_light.update_layout(
        template='plotly_white', # better for publication
        paper_bgcolor='white',
        plot_bgcolor='white',
        font=dict(color='black'), 
    )
    # Update axes for light theme
    fig_light.update_xaxes(gridcolor='lightgrey', linecolor='black', zerolinecolor='lightgrey', color='black', title_font_color='black', tickfont_color='black')
    fig_light.update_yaxes(gridcolor='lightgrey', linecolor='black', zerolinecolor='lightgrey', color='black', title_font_color='black', tickfont_color='black')

    # Specific trace color adjustments 
    # Adjust according to layout specs
    color_map_dark_to_light = {
        colors['line']: '#1f77b4', # Blue
        colors['marker']: '#ff7f0e', # Orange
        colors['highlight']: '#d62728', # Red
        colors['highlight2']: '#ffbb78', # Light orange
        colors['secondary']: '#2ca02c', # Green
        
    }

    for trace in fig_light.data:
        # Update marker colors
        if hasattr(trace, 'marker') and trace.marker:
            original_color = trace.marker.color
            if isinstance(original_color, str) and original_color in color_map_dark_to_light:
                trace.marker.color = color_map_dark_to_light[original_color]
            # Handle marker outlines, etc. if necessary
            if trace.marker.line:
                 trace.marker.line.color = 'black'

        # Update line colors
        if hasattr(trace, 'line') and trace.line:
             original_color = trace.line.color
             if isinstance(original_color, str) and original_color in color_map_dark_to_light:
                trace.line.color = color_map_dark_to_light[original_color]

        # Update fill colors (e.g., for confidence intervals)
        # Note: This is trickier as fillcolor might be rgba

    return fig_light


def export_plot_to_file(fig_dark, output_path, width=1200, height=800, scale=2, format='png'):
    """Exports a Plotly figure to a file using a light theme."""
    if not fig_dark:
        print(f"Warning: Attempted to export an empty figure to {output_path}.")
        return

    print(f"Exporting plot to {output_path} (Format: {format}, Size: {width}x{height}, Scale: {scale})")
    fig_light = create_light_theme_figure(fig_dark)

    try:
        pio.write_image(fig_light, output_path, format=format, width=width, height=height, scale=scale)
        print(f"Successfully exported plot: {output_path}")
    except ValueError as ve:
         # Often indicates kaleido/orca is not installed or found
         print(f"Export Error (ValueError): {ve}")
         print("Please ensure 'kaleido' is installed (`pip install kaleido`) for static image export.")
    except Exception as e:
        print(f"Failed to export plot {output_path}: {e}")


# --- Confidence Interval Calculation --- (Moved from plotting to utils)
def calculate_confidence_intervals(counts, total_bases, confidence_level=0.95):
    """
    Calculate binomial confidence intervals (Wilson score interval) for frequencies.
    Handles zero denominators safely.
    """
    counts_arr = np.array(counts)
    total_bases_arr = np.array(total_bases)

    # Ensure non-negative counts and totals >= counts
    counts_arr = np.maximum(0, counts_arr)
    total_bases_arr = np.maximum(counts_arr, total_bases_arr) # Total must be at least count
    total_bases_arr_safe = np.maximum(1, total_bases_arr) # Avoid division by zero in formulas

    alpha = 1 - confidence_level
    z_score = stats.norm.ppf(1 - alpha / 2)

    # Wilson score interval calculation
    p_hat = counts_arr / total_bases_arr_safe # Observed proportion
    n = total_bases_arr # Use original totals here for accuracy where possible

    factor = z_score**2 / n
    denominator = 1 + factor

    center_adjusted = p_hat + factor / 2
    term_under_sqrt = p_hat * (1 - p_hat) / n + factor / (4 * n)
    # Handle potential negative values under sqrt due to floating point errors
    term_under_sqrt[term_under_sqrt < 0] = 0
    interval_width = z_score * np.sqrt(term_under_sqrt)

    lower = (center_adjusted - interval_width) / denominator
    upper = (center_adjusted + interval_width) / denominator

    # Correct bounds for 0 counts or 100% counts on safe totals
    lower[counts_arr == 0] = 0
    upper[counts_arr == total_bases_arr] = 1

    # Ensure bounds are within [0, 1]
    lower_bounds = np.maximum(0, lower)
    upper_bounds = np.minimum(1, upper)

    # Handle cases where total_bases was originally 0 - interval should be undefined or [0,0] or [0,1]
    lower_bounds[total_bases_arr == 0] = 0 # Or np.nan
    upper_bounds[total_bases_arr == 0] = 0 # Or np.nan / 1.0 depending on convention

    return lower_bounds, upper_bounds
