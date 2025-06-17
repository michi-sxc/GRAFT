import uuid
import dash
from dash import Input, Output, State
import dash_bootstrap_components as dbc
import logging


from app import app 
from processing.read_operations import FileMerger 

logger = logging.getLogger(__name__)

# --- Open Merge Modal ---
@app.callback(
    Output('merge-modal', 'is_open'),
    Input('open-merge-modal-button', 'n_clicks'),
    State('merge-modal', 'is_open'),
    prevent_initial_call=True
)
def open_merge_modal(n_clicks, is_open):
    """Opens the merge file selection modal."""
    if n_clicks:
        return not is_open
    return is_open

# --- Populate Merge Checklist ---
@app.callback(
    [Output('merge-file-checklist', 'options'),
     Output('merge-file-checklist', 'value'),
     Output('merge-validation-alert', 'children'),
     Output('merge-validation-alert', 'is_open'),
     Output('execute-merge-button', 'disabled')],
    Input('merge-modal', 'is_open'), # Update when modal opens
    State('file-store', 'data')
)
def populate_merge_checklist(is_open, store_data):
    """Populates the checklist in the merge modal with compatible files."""
    if not is_open or not store_data or not store_data.get('files'):
        return [], [], "", False, True # No options, hidden alert, disabled button

    files = store_data['files']
    options = [{'label': f['filename'], 'value': f['filename']} for f in files]

    # Basic validation - check if there are at least 2 files
    if len(files) < 2:
        alert_msg = "You need at least two uploaded files to perform a merge."
        return options, [], alert_msg, True, True # Show alert, disable button
    else:
        # Further validation could happen here (e.g., pre-checking types)
        # But the main validation happens on button click
        return options, [], "", False, False # Enable button by default if >1 file


# --- Merge Execution and Download Trigger ---
# This callback now primarily handles triggering the merge and storing results.
# A separate callback handles the actual download.
@app.callback(
    [Output('merge-status', 'children'),
     Output('merge-progress', 'value'),
     Output('merged-file-data-store', 'data'), # Store result for download callback
     Output('merge-modal', 'is_open', allow_duplicate=True)], # Close modal on success/error
    Input('execute-merge-button', 'n_clicks'),
    State('merge-file-checklist', 'value'),
    State('file-store', 'data'),
    State('merge-modal', 'is_open'), # Keep track of modal state
    prevent_initial_call=True,
    # background=True, # Candidate for background callback
    # progress=[Output('merge-progress', 'value'), Output('merge-status', 'children')] # If using background
)
# def execute_merge(set_progress, n_clicks, selected_filenames, store_data, is_open_state): # If using background
def execute_merge(n_clicks, selected_filenames, store_data, is_open_state):
    """Validates selection and initiates the merge process."""
    if not n_clicks or not selected_filenames:
        return "Please select files to merge.", 0, None, is_open_state # Stay in modal

    if len(selected_filenames) < 2:
         # This validation should ideally disable the button earlier
         return dbc.Alert("Please select at least two files.", color="warning"), 0, None, is_open_state

    files_to_merge_data = [f for f in store_data['files'] if f['filename'] in selected_filenames]

    if not files_to_merge_data:
         return dbc.Alert("Selected files not found in store.", color="danger"), 0, None, is_open_state

    merger = FileMerger() # Create instance

    # --- Validation ---
    is_valid, message = merger.validate_files_for_merge(files_to_merge_data)
    if not is_valid:
        return dbc.Alert(message, color="danger"), 0, None, is_open_state # Show error in modal

    # --- Progress Update Function ---
    # If using background callbacks, use set_progress directly
    # Otherwise, this function is called synchronously within merge_files
    def progress_updater(value, msg):
         print(f"Merge Progress: {value}% - {msg}") # Log progress
         # If not using background callback's set_progress, this won't update UI live
         # Need other mechanisms (like polling or websockets) for live sync progress

    # --- Execute Merge ---
    try:
        logger.info(f"Starting merge for: {selected_filenames}")
        # merged_b64_content = merger.merge_files(progress_callback=set_progress) # If background
        merged_b64_content = merger.merge_files(progress_callback=progress_updater) # If synchronous

        output_filename = f"merged_{merger.output_format}_{uuid.uuid4().hex[:8]}.{merger.output_format}"

        # Store results for the download callback
        merged_data_for_store = {
            'filename': output_filename,
            'content_b64': merged_b64_content, # Store encoded content
            'format': merger.output_format
        }
        success_msg = dbc.Alert(f"Merge successful! File '{output_filename}' is ready for download.", color="success")
        logger.info(f"Merge successful. Output: {output_filename}")
        return success_msg, 100, merged_data_for_store, False # Close modal on success

    except Exception as e:
        logger.error(f"Merge execution failed: {e}", exc_info=True)
        error_msg = dbc.Alert(f"Merge Error: {e}", color="danger")
        # Don't close modal on error, show message inside
        return error_msg, 0, None, True


# --- Download Callback for Merged File ---
@app.callback(
    Output('download-merged-file', 'data'),
    Input('merged-file-data-store', 'data'), # Triggered when merge completes successfully
    prevent_initial_call=True
)
def download_merged_file_callback(merged_data):
    """Sends the merged file content to the browser for download."""
    if not merged_data or not merged_data.get('content_b64'):
        raise dash.exceptions.PreventUpdate

    filename = merged_data.get('filename', 'merged_file.out')
    content_b64 = merged_data.get('content_b64')

    # Decode base64 content before sending
    # dcc.send_bytes might be better here if available and handling large files
    try:
        # decoded_content = base64.b64decode(content_b64)
        # return dcc.send_bytes(decoded_content, filename)

        # Alternative: send base64 string directly using dictionary format
         return dict(content=content_b64, filename=filename, base64=True)

    except Exception as e:
        logger.error(f"Error preparing merged file for download: {e}")
        # How to signal error back to user from download callback? Difficult.
        # Maybe update a status message elsewhere.
        return None