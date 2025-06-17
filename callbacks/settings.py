import dash
from dash import Input, Output, State, callback_context
import logging

from app import app 
logger = logging.getLogger(__name__)

# --- Toggle Settings Offcanvas ---
@app.callback(
    Output("settings-offcanvas", "is_open"),
    Input("settings-button", "n_clicks"),
    Input("apply-settings-button", "n_clicks"), # Close on apply
    State("settings-offcanvas", "is_open"),
    prevent_initial_call=True,
)
def toggle_settings_offcanvas(settings_click, apply_click, is_open):
    """Toggles the settings offcanvas visibility."""
    triggered_id = callback_context.triggered_id
    if triggered_id == "settings-button":
        return not is_open
    elif triggered_id == "apply-settings-button":
        return False # Close the offcanvas when Apply is clicked
    return is_open


# --- Apply Settings (Example - currently placeholder) ---
@app.callback(
    Output("apply-settings-button", "disabled"), # Can disable button until changes are made
    # Inputs for all settings components
    Input("plot-format-select", "value"),
    Input("plot-width-input", "value"),
    Input("plot-height-input", "value"),
    Input("plot-scale-input", "value"),
    # Add inputs for other settings like performance toggles etc.
    # State("show-performance-stats", "value"),
    prevent_initial_call=True
)
def settings_changed(plot_format, width, height, scale): # Add other settings args
    """
    Handles applying settings. Currently just logs and enables/disables button.
    Actual application logic might happen elsewhere or trigger other callbacks.
    """
    # In a real scenario, you might compare current values to stored/default values
    # to determine if the Apply button should be enabled.
    logger.info("Settings changed - Apply button active (placeholder logic).")
    # Here you could save settings to config.yaml or a user session if needed
    # For example:
    # current_config = {
    #     'export_settings': {'format': plot_format, 'width': width, 'height': height, 'scale': scale},
    #     # 'performance_settings': {'show_stats': show_perf_stats}
    # }
    # try:
    #     with open('config.yaml', 'w') as f:
    #         yaml.dump(current_config, f, default_flow_style=False)
    #     logger.info("Settings saved to config.yaml")
    # except Exception as e:
    #     logger.error(f"Failed to save settings to config.yaml: {e}")

    return False # Keep apply button enabled once settings are potentially changed

# --- Reset Settings ---
@app.callback(
    [Output("plot-format-select", "value"),
     Output("plot-width-input", "value"),
     Output("plot-height-input", "value"),
     Output("plot-scale-input", "value")],
    # Add outputs for other settings components
    Input("reset-settings-button", "n_clicks"),
    prevent_initial_call=True
)
def reset_settings_to_default(n_clicks):
    """Resets settings components to their default values."""
    if not n_clicks:
        raise dash.exceptions.PreventUpdate

    # Define default values
    defaults = {
        'format': 'png',
        'width': 1200,
        'height': 800,
        'scale': 2,
        # Add other defaults
    }
    logger.info("Resetting settings to default values.")
    return defaults['format'], defaults['width'], defaults['height'], defaults['scale'] # Return defaults in order