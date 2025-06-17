import dash
from dash import Input, Output, State, html, callback_context, DiskcacheManager
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import numpy as np
import plotly.graph_objects as go
import pysam
import diskcache
from Bio import SeqIO

from app import app
from processing.read_operations import load_reads_from_file
from processing.filters import filter_read_tuples, adjust_nm_for_ct, deduplicate_reads
from processing.stats import (
    calculate_alignment_stats, get_damage_read_end, calculate_mismatch_frequency_vs_end,
    get_mismatch_counts_by_type, categorize_mismatches
)
from plotting.histograms import (
    create_read_length_histogram, create_overall_cg_histogram,
    create_cg_content_histogram_selected, create_mapq_histogram
)
from plotting.damage_plots import create_damage_at_ends_figure, create_mismatch_frequency_plot, create_damage_pattern_heatmap
from plotting.qc_plots import create_mismatch_type_bar_chart, create_mismatch_pie_chart
import logging

logger = logging.getLogger(__name__)


if not hasattr(app, 'background_callback_manager'): # Check if already set
    try:
        cache = diskcache.Cache("./graft_cache") # Ensure this directory is writable
        app.background_callback_manager = DiskcacheManager(cache)
        logger.info("DiskcacheManager for background callbacks initialized in callbacks/main.py.")
    except Exception as e:
        logger.error(f"Failed to initialize DiskcacheManager in callbacks/main.py: {e}")
        # Fallback or raise error if background callbacks are critical
        app.background_callback_manager = None


# --- Persistent Settings Store ---

# Store the IDs of components whose values need to be restored
# This helps manage the potentially long list of outputs
PERSISTENT_COMPONENT_IDS_VALUES = {
    'file-selector': 'value',
    'mapq-range': 'value',
    'min-base-quality': 'value',
    'soft-clip-handling': 'value',
    'nm-dropdown': 'value',
    'ct-checklist': 'value',
    'ct-count-dropdown': 'value',
    'sequencing-type': 'value',
    'apply-sequence-deduplication-switch': 'value',
    'selected-lengths': 'data', # This is a store, will be restored to memory
}
# For file selector options
PERSISTENT_FILE_SELECTOR_OPTIONS_ID = 'file-selector'


@app.callback(
    Output('file-selector', 'value', allow_duplicate=True),
    Input('url', 'pathname'),
    State('persistent-settings-store', 'data'),
    prevent_initial_call=True
)
def handle_navigation_file_selector(pathname, stored_settings):
    """Clear or restore file selection based on navigation."""
    if pathname != '/':
        # Not on main page, file-selector might not exist.
        # Do not attempt to update it.
        # logger.debug("handle_navigation_file_selector: Not on main page, preventing update.")
        raise dash.exceptions.PreventUpdate # MODIFIED
    else:
        # On main page, file-selector should exist.
        # Restore from persistent storage if available
        # logger.debug("handle_navigation_file_selector: On main page, attempting to restore file selector value.")
        if stored_settings and 'file-selector' in stored_settings:
            return stored_settings['file-selector']
        return None


@app.callback(
    [Output(comp_id, prop) for comp_id, prop in PERSISTENT_COMPONENT_IDS_VALUES.items()] +
    [Output(PERSISTENT_FILE_SELECTOR_OPTIONS_ID, 'options')],
    Input('main-page-loaded-signal', 'data'), # MODIFIED Input
    State('persistent-settings-store', 'data'),
    State('persistent-file-metadata-store', 'data'),
    State('url', 'pathname'), # Keep for an additional check
    prevent_initial_call=True # ADDED: Crucial for signal-based trigger
)
def load_persistent_settings(main_page_signal_data, stored_settings, stored_file_metadata, current_pathname):
    # This callback now runs only after display_page has loaded the main layout and sent a signal.
    # prevent_initial_call=True ensures it doesn't run when main_page_signal_data is initially None.

    if not main_page_signal_data or current_pathname != '/':
        # logger.debug(
        #    f"load_persistent_settings: Invalid trigger. Signal: {main_page_signal_data}, Pathname: {current_pathname}. Preventing update."
        # )
        raise dash.exceptions.PreventUpdate

    # logger.debug(f"load_persistent_settings: Loading persistent settings for main page. Pathname: {current_pathname}")

    if stored_settings is None and stored_file_metadata is None:
        # logger.debug("load_persistent_settings: No stored settings or file metadata found. Preventing update.")
        raise dash.exceptions.PreventUpdate

    output_values = []
    
    # Restore component values
    for comp_id, prop in PERSISTENT_COMPONENT_IDS_VALUES.items():
        if stored_settings and comp_id in stored_settings:
            output_values.append(stored_settings[comp_id])
        else:
            if comp_id == 'selected-lengths':
                output_values.append([])
            elif comp_id == 'mapq-range':
                output_values.append([0, 255])
            elif comp_id == 'ct-checklist':
                 output_values.append([])
            else:
                default_map = {
                    'file-selector': None, 
                    'min-base-quality': 20,
                    'soft-clip-handling': 'show', 
                    'nm-dropdown': 'all',
                    'ct-count-dropdown': 'any', 
                    'sequencing-type': 'single',
                    'apply-sequence-deduplication-switch': False
                }
                output_values.append(default_map.get(comp_id, dash.no_update))

    # Restore file selector options
    file_selector_options = []
    if stored_file_metadata:
        for meta in stored_file_metadata:
            file_selector_options.append({'label': meta['filename'], 'value': meta['filename']})
    
    output_values.append(file_selector_options)
    
    # logger.debug(f"Loading from persistent-settings-store: {stored_settings}")
    # logger.debug(f"Generated output values for components: {output_values[:-1]}")
    # logger.debug(f"Generated file selector options: {file_selector_options}")

    return tuple(output_values)


# --- File Upload and Selection ---
@app.callback(
    [Output('file-store', 'data'), # Memory store for current session's file contents
     Output('file-selector', 'options', allow_duplicate=True),
     Output('file-selector', 'value', allow_duplicate=True),
     Output('upload-status-messages', 'children'),
     Output('persistent-file-metadata-store', 'data', allow_duplicate=True)], # ADDED allow_duplicate
    Input('upload-data', 'contents'),
    State('upload-data', 'filename'),
    State('file-store', 'data'), # Current memory contents
    State('persistent-file-metadata-store', 'data'), # Previously persisted metadata
    # State('file-selector', 'options'), # Not strictly needed, can reconstruct from persisted_metadata
    State('file-selector', 'value'),
    State('url', 'pathname'),
    prevent_initial_call=True
)
def update_file_store_and_selector(
    list_of_contents, list_of_filenames,
    current_memory_file_store_data,
    persisted_metadata_store_data,
    current_selected_file_value,
    current_pathname
):
    if current_pathname != '/':
        raise PreventUpdate

    ctx = dash.callback_context # To check if actually triggered by upload
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None

    memory_files_list = (current_memory_file_store_data or {}).get('files', [])
    
    # Persisted metadata (for file selector options)
    # This represents all files "known" to the application from previous uploads
    persisted_metadata_list = persisted_metadata_store_data or []
    persisted_metadata_dict = {meta['filename']: meta for meta in persisted_metadata_list}

    messages = []
    new_files_added_to_memory_this_upload = False
    made_changes_to_persistent_metadata = False # Flag to see if we need to update the store

    if triggered_id == 'upload-data' and list_of_contents:
        existing_filenames_in_memory = {f['filename'] for f in memory_files_list}

        for content, filename in zip(list_of_contents, list_of_filenames):
            file_format_parts = filename.split('.')
            file_ext = file_format_parts[-1].lower()
            true_format = file_format_parts[-2].lower() if len(file_format_parts) > 1 and file_ext in ['gz', 'bz2'] else file_ext

            if true_format not in ['bam', 'sam', 'fasta', 'fa', 'fna', 'fq', 'fastq']:
                messages.append(dbc.Alert(f"Unsupported file format for '{filename}'. Skipping.", color="danger", duration=4000))
                continue
            
            # Handle re-upload of a file already in memory for this session
            if filename in existing_filenames_in_memory:
                messages.append(dbc.Alert(f"File '{filename}' re-uploaded to current session. Content updated.", color="info", duration=3000))
                for f_data in memory_files_list:
                    if f_data['filename'] == filename:
                        f_data['content'] = content # Update content
                        f_data['format'] = true_format # Ensure format is also up-to-date
                        break
                # Ensure it's in persisted metadata if somehow it wasn't (e.g. error in previous logic)
                if filename not in persisted_metadata_dict:
                    persisted_metadata_dict[filename] = {'filename': filename, 'format': true_format}
                    made_changes_to_persistent_metadata = True
                continue # Next file in upload batch

            # New file for this session's memory store
            new_file_data_for_memory = {'filename': filename, 'content': content, 'format': true_format}
            memory_files_list.append(new_file_data_for_memory)
            new_files_added_to_memory_this_upload = True
            messages.append(dbc.Alert(f"File '{filename}' uploaded to session.", color="success", duration=3000))

            # Add to persisted_metadata_dict if not already there
            if filename not in persisted_metadata_dict:
                persisted_metadata_dict[filename] = {'filename': filename, 'format': true_format}
                made_changes_to_persistent_metadata = True
    
    # Options for the dropdown are always derived from the full list of known files
    final_options = [{'label': meta['filename'], 'value': meta['filename']} for meta in persisted_metadata_dict.values()]
    
    new_selected_value = current_selected_file_value
    if new_files_added_to_memory_this_upload and list_of_filenames:
        new_selected_value = list_of_filenames[-1] # Select the last uploaded file from the batch

    if new_selected_value is None and final_options:
        new_selected_value = final_options[0]['value']
    
    if new_selected_value and not any(opt['value'] == new_selected_value for opt in final_options):
        new_selected_value = final_options[0]['value'] if final_options else None

    updated_memory_file_store = {'files': memory_files_list}
    
    # Only update persistent-file-metadata-store if there were actual changes to its content
    output_persisted_metadata = list(persisted_metadata_dict.values()) if made_changes_to_persistent_metadata else dash.no_update
    if made_changes_to_persistent_metadata:
         logger.info(f"Updating persistent-file-metadata-store.") # Data can be large
    # else:
    #      logger.debug("No changes to persistent_metadata_dict, so not updating persistent-file-metadata-store.")

    return updated_memory_file_store, final_options, new_selected_value, messages, output_persisted_metadata


# --- Data Loading and Processing Store ---
# This callback acts as the central point for loading and filtering data based on UI controls.
# It stores the results in intermediate stores, triggering plot updates.
@app.callback(
    [Output('processed-data-store', 'data'),
     Output('alignment-stats-store', 'data'),
     Output('mismatch-frequencies-store', 'data'),
     Output('main-loader-status', 'children'), # For text status
     Output('main-progress-bar', 'value'),    # For progress bar value
     Output('main-progress-bar', 'style')],   # To show/hide progress bar
    [Input('file-selector', 'value'),
     Input('mapq-range', 'value'),
     Input('min-base-quality', 'value'),
     Input('soft-clip-handling', 'value'),
     Input('nm-dropdown', 'value'),
     Input('ct-checklist', 'value'),
     Input('ct-count-dropdown', 'value'),
     Input('sequencing-type', 'value'),
     Input('apply-sequence-deduplication-switch', 'value')],
    State('file-store', 'data'),
    State('url', 'pathname'),  # Add this to check current page
    background=True,
    manager=app.background_callback_manager if hasattr(app, 'background_callback_manager') else None,
    progress=[
        Output('main-loader-status', 'children'),
        Output('main-progress-bar', 'value')
    ],
    prevent_initial_call=True,
)
def load_and_process_data(
    set_progress,
    selected_file, mapq_range, min_base_quality, soft_clip_option,
    selected_nm, ct_checklist, ct_count_value,
    sequencing_type,
    apply_sequence_deduplication_value,
    file_store_data,
    current_pathname
):
    # Check if we're on the main page
    if current_pathname != '/':
        logger.debug(f"load_and_process_data: Not on main page (pathname: {current_pathname}). Preventing update.")
        raise PreventUpdate

    # Default empty/error structures for main outputs
    default_processed_data = {'reads': [], 'format': None, 'error': "Initialization", 'filename': None, 'filters_applied': {}}
    default_stats_data = {'stats': None, 'error': "Initialization"}
    default_mismatch_freq_data = {'frequencies': None, 'format': None, 'error': "Initialization"}

    # Initial state for progress outputs (progress bar hidden initially by its style in layout)
    # The style output in the main callback will make it visible when processing starts
    initial_progress_message = "Waiting for file selection or filter change..."
    
    # --- Function to update progress via set_progress ---
    # Total steps estimated for progress calculation
    total_major_steps = 8 # 1.Load, 2.DefineFilters, 3.NMAdjust, 4.Filter, 5.Dedupe, 6.MismatchFreq, 7.Serialize, 8.AlignStats

    def _update_progress_ui(step_num, message):
        if set_progress: # Check if set_progress is available (it should be for background callbacks)
            percentage = int((step_num / total_major_steps) * 100)
            set_progress((html.Div(message, className="text-info small"), percentage))

    _update_progress_ui(0, "Initializing data processing...")

    if not selected_file or not file_store_data:
        logger.info("No file selected or file store empty in load_and_process_data.")
        # For background, return tuple matching all outputs (main + style for progress bar)
        return default_processed_data, default_stats_data, default_mismatch_freq_data, html.Div(initial_progress_message), 0, {'display': 'none'}


    file_data = next((f for f in file_store_data['files'] if f['filename'] == selected_file), None)
    if not file_data:
        logger.warning(f"File content for '{selected_file}' not found in session memory. Please re-upload to process.")
        # More descriptive error message for the user
        error_msg_for_ui = f"File content for '{selected_file}' not in current session. Please re-upload the file to generate plots."
        
        # Ensure filters_applied is part of the error state, even if empty
        # This helps data_summary_tab if it expects this key
        current_filters = { # Capture current filter settings to display them even on error
            'mapq_range': tuple(mapq_range),
            'min_base_quality': int(min_base_quality) if min_base_quality is not None else 0,
            'soft_clip_option': soft_clip_option,
            'filter_ct': True if 'only_ct' in ct_checklist else False if 'exclude_ct' in ct_checklist else None,
            'exclusively_ct': 'exclusively_ct' in ct_checklist,
            'ct_count_value': ct_count_value if ct_count_value == 'any' else int(ct_count_value),
            'selected_nm': selected_nm if selected_nm == 'all' else int(selected_nm),
            'subtract_ct_from_nm': 'subtract_ct' in ct_checklist,
            'apply_sequence_deduplication': apply_sequence_deduplication_value
        }

        processed_error = {
            'reads': [], 'format': None, 'error': error_msg_for_ui, 'filename': selected_file,
            'filters_applied': current_filters # Store current filter settings
        }
        stats_error = {'stats': None, 'error': error_msg_for_ui}
        freq_error = {'frequencies': None, 'format': None, 'error': error_msg_for_ui}
        
        # Update UI status and hide progress bar
        return processed_error, stats_error, freq_error, html.Div(error_msg_for_ui, className="text-warning small fw-bold"), 0, {'display': 'none'}

    logger.info(f"Processing file: {selected_file} with filters: MAPQ={mapq_range}, BaseQ>={min_base_quality}, SoftClip={soft_clip_option}, NM={selected_nm}, CT_Check={ct_checklist}, CT_Count={ct_count_value}, SeqType={sequencing_type}, SeqDedupe={apply_sequence_deduplication_value}")
    _update_progress_ui(0, f"Starting processing for {selected_file}...")


    try:
        # --- 1. Load Raw Reads ---
        _update_progress_ui(1, f"Loading raw reads from {selected_file}...")
        raw_reads_data_tuples, file_format_inferred, _ = load_reads_from_file(
            file_data['content'],
            file_data['filename'],
            soft_clip_option=soft_clip_option
        )
        logger.debug(f"Loaded {len(raw_reads_data_tuples)} raw reads. Format: {file_format_inferred}")

        # --- 2. Define Filters Dictionary ---
        _update_progress_ui(2, "Preparing filters...")
        filter_ct_flag = None
        if 'only_ct' in ct_checklist: filter_ct_flag = True
        elif 'exclude_ct' in ct_checklist: filter_ct_flag = False

        filters = {
            'mapq_range': tuple(mapq_range),
            'min_base_quality': int(min_base_quality) if min_base_quality is not None else 0,
            'soft_clip_option': soft_clip_option,
            'filter_ct': filter_ct_flag,
            'exclusively_ct': 'exclusively_ct' in ct_checklist,
            'ct_count_value': ct_count_value if ct_count_value == 'any' else int(ct_count_value),
            'selected_nm': selected_nm if selected_nm == 'all' else int(selected_nm),
            'subtract_ct_from_nm': 'subtract_ct' in ct_checklist,
            'apply_sequence_deduplication': apply_sequence_deduplication_value
        }
        logger.info(f"Filters being applied in load_and_process_data: {filters}")


        # --- 3. Adjust NM based on C>T ---
        _update_progress_ui(3, "Adjusting NM values for C>T changes...")
        nm_adjusted_reads_tuples = adjust_nm_for_ct(raw_reads_data_tuples, filters['subtract_ct_from_nm'])
        logger.debug(f"Count after NM adjustment: {len(nm_adjusted_reads_tuples)}")

        # --- 4. Apply Filters ---
        _update_progress_ui(4, "Applying general read filters...")
        filtered_reads_tuples = filter_read_tuples(nm_adjusted_reads_tuples, file_format_inferred, filters)
        logger.debug(f"Count after general filters: {len(filtered_reads_tuples)}")

        # --- 5. Deduplicate ---
        _update_progress_ui(5, "Applying sequence deduplication (if selected)...")
        if filters.get('apply_sequence_deduplication', False):
            final_reads_tuples = deduplicate_reads(filtered_reads_tuples)
            logger.info(f"Sequence deduplication applied. Original: {len(filtered_reads_tuples)}, Final: {len(final_reads_tuples)}")
        else:
            final_reads_tuples = filtered_reads_tuples
            logger.info("Sequence deduplication skipped.")
        logger.info(f"Final processed reads (tuples) count for {selected_file}: {len(final_reads_tuples)}")

        if not final_reads_tuples:
            no_reads_msg = "No reads remained after filtering and deduplication."
            logger.warning(no_reads_msg)
            processed_data_output = {**default_processed_data, 'error': no_reads_msg, 'filename': selected_file, 'filters_applied': filters}
            stats_data_output = {**default_stats_data, 'error': no_reads_msg}
            mismatch_frequencies_output = {**default_mismatch_freq_data, 'error': no_reads_msg, 'format': file_format_inferred}
            _update_progress_ui(total_major_steps, "Processing complete: No reads passed filters.")
            return processed_data_output, stats_data_output, mismatch_frequencies_output, html.Div(no_reads_msg, className="text-warning small"), 100, {'display': 'none'}

        # --- 6. Calculate Mismatch Frequencies ---
        _update_progress_ui(6, "Calculating mismatch frequencies...")
        mismatch_frequencies_output = {'frequencies': None, 'format': file_format_inferred, 'error': None}
        if file_format_inferred in ['bam', 'sam']: # Already checked final_reads_tuples and tuple structure implicitly
            try:
                logger.info(f"PRE-CALC: Attempting mismatch frequency calculation for {len(final_reads_tuples)} BAM/SAM reads. Current mismatch_frequencies_output['frequencies'] is {mismatch_frequencies_output['frequencies']}")
                frequencies_calculated = calculate_mismatch_frequency_vs_end(final_reads_tuples, sequencing_type)
                # Detailed logging for frequencies_calculated (as added in previous step)
                if frequencies_calculated is not None:
                    logger.info(f"POST-CALC: calculate_mismatch_frequency_vs_end returned non-None. Type: {type(frequencies_calculated)}")
                    if isinstance(frequencies_calculated, dict) and sequencing_type == 'single':
                        logger.info(f"POST-CALC (single): 'C>T_CpG' (first 5 freqs): {frequencies_calculated.get('C>T_CpG', [])[:5]}, 'C>T_CpG_counts' (first 5): {frequencies_calculated.get('C>T_CpG_counts', [])[:5]}")
                    elif isinstance(frequencies_calculated, dict) and sequencing_type == 'double':
                         logger.info(f"POST-CALC (double): '5_prime' data (first 5 C>T_CpG freqs): {frequencies_calculated.get('5_prime', {}).get('C>T_CpG', [])[:5]}")
                else:
                    logger.warning("POST-CALC: calculate_mismatch_frequency_vs_end returned None!")

                mismatch_frequencies_output['frequencies'] = frequencies_calculated
                logger.info(f"POST-ASSIGN: mismatch_frequencies_output['frequencies'] is now set. Type: {type(mismatch_frequencies_output['frequencies'])}. Is None: {mismatch_frequencies_output['frequencies'] is None}")
                logger.info("Mismatch frequency calculation successful (message from load_and_process_data).")
            except ValueError as ve:
                logger.error(f"ValueError calculating mismatch frequencies for {selected_file}: {ve}", exc_info=True)
                mismatch_frequencies_output['error'] = f"Data structure error during frequency calc: {ve}"
            except Exception as freq_e:
                logger.error(f"General error calculating mismatch frequencies for {selected_file}: {freq_e}", exc_info=True)
                mismatch_frequencies_output['error'] = str(freq_e)
        else: # Not BAM/SAM
            mismatch_frequencies_output['error'] = "Mismatch frequency calculation not applicable for this file type."
        logger.info(f"FINAL-PRE-RETURN (mismatch_freq): mismatch_frequencies_output to be stored: error='{mismatch_frequencies_output.get('error')}', frequencies_is_none={mismatch_frequencies_output.get('frequencies') is None}")


        # --- 7. Serialize Reads for processed-data-store ---
        _update_progress_ui(7, "Serializing reads for display...")
        serializable_reads = []
        for length, cg, nm_val, seq_str, read_obj, mapq_val in final_reads_tuples:
            record_info = {
                'length': length, 'cg': cg, 'nm': nm_val, 'seq': seq_str, 'mapq': mapq_val,
                'original_nm': read_obj.get_tag('NM') if isinstance(read_obj, pysam.AlignedSegment) and read_obj.has_tag('NM') else nm_val,
            }
            if isinstance(read_obj, pysam.AlignedSegment):
                record_info.update({
                    'name': read_obj.query_name, 'flag': read_obj.flag,
                    'is_mapped': not read_obj.is_unmapped, 'is_reverse': read_obj.is_reverse,
                    'cigarstring': read_obj.cigarstring,
                })
            elif isinstance(read_obj, SeqIO.SeqRecord):
                 record_info.update({
                     'name': read_obj.id, 'description': read_obj.description, 'is_reverse': False
                 })
            serializable_reads.append(record_info)
        processed_data_output = {
            'reads': serializable_reads, 'format': file_format_inferred,
            'filename': selected_file, 'filters_applied': filters, 'error': None
        }

        # --- 8. Calculate Alignment Stats ---
        _update_progress_ui(total_major_steps, "Calculating alignment statistics...") # Last step before completion
        alignment_stats = calculate_alignment_stats(final_reads_tuples) # Uses original objects in final_reads_tuples
        stats_data_output = {'stats': alignment_stats, 'error': None}

        _update_progress_ui(total_major_steps, f"Processing for {selected_file} complete!")
        # Final return for main outputs, and hide progress bar
        return processed_data_output, stats_data_output, mismatch_frequencies_output, html.Div(f"Processing for {selected_file} complete!", className="text-success small"), 100, {'display': 'none'}

    except Exception as e:
        logger.error(f"Error processing file {selected_file} in load_and_process_data (background): {e}", exc_info=True)
        error_msg = f"Error processing {selected_file}: {str(e)}"
        processed_error_state = {**default_processed_data, 'error': error_msg, 'filename': selected_file}
        stats_error_state = {**default_stats_data, 'error': error_msg}
        freq_error_state = {**default_mismatch_freq_data, 'error': error_msg}
        return processed_error_state, stats_error_state, freq_error_state, html.Div(error_msg, className="text-danger small"), 0, {'display': 'none'}
    
# --- Plotting Callbacks ---
# These callbacks are triggered by the 'processed-data-store'
@app.callback(
    Output('read-length-histogram', 'figure'),
    Input('processed-data-store', 'data'),
    Input('selected-lengths', 'data')
)
def update_read_length_plot(processed_data, selected_lengths):
    if processed_data is None: # Explicit check for None first
        return go.Figure().update_layout(title="Read Length: Waiting for data...")
    # Now we know processed_data is not None, so .get() is safe
    if processed_data.get('error'):
        return go.Figure().update_layout(title=f"Read Length: {processed_data.get('error', 'Error loading data')}")
    if not processed_data.get('reads'): # Check if reads list is empty
        return go.Figure().update_layout(title="Read Length: No reads to display")
    
    # If you reach here, processed_data is valid and has reads
    return create_read_length_histogram(processed_data['reads'], selected_lengths)

@app.callback(
    Output('overall-cg-histogram', 'figure'),
    Input('processed-data-store', 'data')
)
def update_overall_cg_plot(processed_data):
    if processed_data is None:
        return go.Figure().update_layout(title="Overall CG: Waiting for data...")
    if processed_data.get('error'):
        return go.Figure().update_layout(title=f"Overall CG: {processed_data.get('error', 'Error loading data')}")
    if not processed_data.get('reads'):
        return go.Figure().update_layout(title="Overall CG: No reads to display")
    return create_overall_cg_histogram(processed_data['reads'])

@app.callback(
    Output('cg-content-histogram', 'figure'),
    Input('processed-data-store', 'data'),
    Input('selected-lengths', 'data') # Triggered by selection changes
)
def update_selected_cg_plot(processed_data, selected_lengths):
    if processed_data is None:
        return go.Figure().update_layout(title="CG Content (Selected): Waiting for data...")
    if processed_data.get('error'):
        return go.Figure().update_layout(title=f"CG Content (Selected): {processed_data.get('error', 'Error loading data')}")
    if not processed_data.get('reads'):
        return go.Figure().update_layout(title="CG Content (Selected): No reads to display")
    return create_cg_content_histogram_selected(processed_data['reads'], selected_lengths)

@app.callback(
    Output('damage-patterns', 'figure'),
    Input('processed-data-store', 'data')
)
def update_damage_ends_plot(processed_data):
    if processed_data is None:
        return go.Figure().update_layout(title="End Patterns: Waiting for data...")
    if processed_data.get('error'):
        return go.Figure().update_layout(title=f"End Patterns: {processed_data.get('error', 'Error loading data')}")
    if not processed_data.get('reads'):
        return go.Figure().update_layout(title="End Patterns: No reads to display")
    if not ('seq' in processed_data['reads'][0] and 'is_reverse' in processed_data['reads'][0]):
        logger.warning("Damage patterns plot cannot be generated: 'seq' or 'is_reverse' missing from processed read data.")
        return go.Figure().update_layout(title="End Patterns (Sequence data missing)")
    five_prime_pct, three_prime_pct = get_damage_read_end(processed_data['reads'])
    return create_damage_at_ends_figure(five_prime_pct, three_prime_pct)

@app.callback(
    Output('mismatch-frequency-plot', 'figure'),
    # Input('processed-data-store', 'data'), # Not strictly needed if mismatch_frequencies_data has all info
    Input('sequencing-type', 'value'),
    Input('mismatch-frequencies-store', 'data')
)
def update_mismatch_freq_plot(sequencing_type, mismatch_frequencies_data): 
    if mismatch_frequencies_data is None:
        return go.Figure().update_layout(title="Damage Freq: Waiting for calculation...")

    # Centralized error display
    if mismatch_frequencies_data.get('error'):
        # This will now correctly show "ValueError in frequency calculation: too many values to unpack (expected 6)"
        # OR "Mismatch frequency calculation not applicable for this file type."
        return go.Figure().update_layout(title=f"Damage Freq Error: {mismatch_frequencies_data.get('error')}")

    # If no error, but format is still not right (should be caught by error above, but as a safeguard)
    if mismatch_frequencies_data.get('format') not in ['bam', 'sam']:
         return go.Figure().update_layout(title="Damage Frequency (Requires BAM/SAM; current file is not or error occurred)")

    # If frequencies themselves are None or empty (but no error was reported)
    if mismatch_frequencies_data.get('frequencies') is None:
        return go.Figure().update_layout(title="Damage Freq: Calculation successful but no frequency data generated.")

    return create_mismatch_frequency_plot(mismatch_frequencies_data['frequencies'], sequencing_type)

@app.callback(
    Output('mismatch-type-bar-chart', 'figure'),
    Input('alignment-stats-store', 'data')
)
def update_mismatch_type_plot(stats_data):
    if stats_data is None:
        return go.Figure().update_layout(title="Mismatch Types: Waiting for data...")
    if stats_data.get('error'):
        return go.Figure().update_layout(title=f"Mismatch Types: {stats_data.get('error', 'Error loading data')}")
    if not stats_data.get('stats'):
        return go.Figure().update_layout(title="Mismatch Types: No statistics available")
    if stats_data['stats'].get('File Format') != 'bam/sam':
        return go.Figure().update_layout(title="Mismatch Types (Requires BAM/SAM)")
    mismatch_counts = get_mismatch_counts_by_type(stats_data['stats'].get('Mismatch Details', []))
    return create_mismatch_type_bar_chart(mismatch_counts)


@app.callback(
    Output('damage-pattern-plot', 'figure'),
    Input('alignment-stats-store', 'data')
)
def update_damage_heatmap_plot(stats_data):
    if stats_data is None:
        return go.Figure().update_layout(title="Damage Heatmap: Waiting for data...")
    if stats_data.get('error'):
        return go.Figure().update_layout(title=f"Damage Heatmap: {stats_data.get('error', 'Error loading data')}")
    if not stats_data.get('stats'):
        return go.Figure().update_layout(title="Damage Heatmap: No statistics available")
    if stats_data['stats'].get('File Format') != 'bam/sam':
        return go.Figure().update_layout(title="Damage Heatmap (Requires BAM/SAM)")
    return create_damage_pattern_heatmap(stats_data['stats'].get('Mismatch Details', []))


@app.callback(
    Output('mapq-histogram', 'figure'),
    Input('alignment-stats-store', 'data')
)
def update_mapq_dist_plot(stats_data):
    if stats_data is None:
        return go.Figure().update_layout(title="MAPQ Distribution: Waiting for data...")
    if stats_data.get('error'):
        return go.Figure().update_layout(title=f"MAPQ Distribution: {stats_data.get('error', 'Error loading data')}")
    if not stats_data.get('stats'):
        return go.Figure().update_layout(title="MAPQ Distribution: No statistics available")
    if stats_data['stats'].get('File Format') != 'bam/sam':
        return go.Figure().update_layout(title="MAPQ Distribution (Requires BAM/SAM)")
    return create_mapq_histogram(stats_data['stats'].get('MAPQ Scores', []))


# --- Selection Handling ---
@app.callback(
    Output('selected-lengths', 'data', allow_duplicate=True),
    Input('read-length-histogram', 'selectedData'),
    Input('read-length-histogram', 'clickData'),
    Input('clear-selection-button', 'n_clicks'),
    Input('file-selector', 'value'),
    State('selected-lengths', 'data'),
    State('url', 'pathname'),  
    prevent_initial_call=True,
)
def update_selected_lengths(selected_data, click_data, n_clicks_clear, selected_file, current_selection, current_pathname):
    """Manages the list of selected read lengths based on plot interactions."""
    
    # Check if we're on the main page
    if current_pathname != '/':
        raise PreventUpdate
    
    triggered_id = callback_context.triggered_id
    current_selection = current_selection or []

    if triggered_id == 'clear-selection-button' or triggered_id == 'file-selector':
        return []  # Clear selection

    new_selection = set(current_selection)

    if triggered_id == 'read-length-histogram' and selected_data and selected_data.get('points'):
        # Box select - replace current selection with the points in the box
        selected_x_values = set()
        for point in selected_data['points']:
            # selectedData often gives range for bars, use 'x' which is the bin center
            if 'x' in point:
                 # Assuming bin centers represent the length
                 selected_x_values.add(int(round(point['x'])))
        logger.debug(f"Box selection update: {selected_x_values}")
        # If shift is held, add to selection? Dash doesn't easily exposes modifiers.
        # Default: Overwrite selection with box select.
        return sorted(list(selected_x_values))

    elif triggered_id == 'read-length-histogram' and click_data and click_data.get('points'):
        # Single click - toggle selection of the clicked bar
        clicked_length = int(round(click_data['points'][0]['x']))
        if clicked_length in new_selection:
            new_selection.remove(clicked_length)
            logger.debug(f"Click deselection: {clicked_length}")
        else:
            new_selection.add(clicked_length)
            logger.debug(f"Click selection: {clicked_length}")
        return sorted(list(new_selection))

    # If no relevant input triggered, return the current state
    # This prevents selectedData=None after clearing selection from resetting it incorrectly
    return current_selection


# --- Data Summary Tab ---
@app.callback(
    Output('data-summary', 'children'),
    Input('processed-data-store', 'data'),
    Input('alignment-stats-store', 'data'),
    Input('selected-lengths', 'data'),
    Input('sequencing-type', 'value'),
    Input('mismatch-frequencies-store', 'data'), 
)
def update_data_summary_tab(processed_data, stats_data, selected_lengths, sequencing_type, mismatch_frequencies_data):
    if processed_data is None:
        return "Waiting for data..."
    if processed_data.get('error'):
        return f"Error: {processed_data.get('error', 'Error loading data')}"
    if not processed_data.get('reads'):
        return "No reads to display"

    reads_list_of_dicts = processed_data.get('reads', [])
    filters_applied = processed_data.get('filters_applied', {})
    stats = stats_data.get('stats') if stats_data else None

    # --- Generate Text Sections ---
    summary = f"--- Summary for: {processed_data.get('filename', 'N/A')} ---\n\n"

    # Applied Filters
    summary += "**Applied Filters:**\n"
    summary += f"- MAPQ Range: {filters_applied.get('mapq_range', 'N/A')}\n"
    summary += f"- Min Base Quality: {filters_applied.get('min_base_quality', 'N/A')}\n"
    summary += f"- Soft Clipping: {filters_applied.get('soft_clip_option', 'N/A')}\n"
    summary += f"- NM Filter: {filters_applied.get('selected_nm', 'N/A')}\n"
    ct_filt = filters_applied.get('filter_ct')
    ct_excl = filters_applied.get('exclusively_ct')
    ct_sub = 'subtract_ct' in filters_applied.get('ct_checklist', []) # Reconstruct from filters dict if possible
    ct_filt_str = "Require C>T/G>A" if ct_filt is True else "Exclude C>T/G>A" if ct_filt is False else "None"
    if ct_excl: ct_filt_str = "Require ONLY C>T/G>A"
    summary += f"- C>T Logic: {ct_filt_str}" + (" (Subtracted from NM)" if ct_sub else "") + "\n"
    summary += f"- C>T Count: {filters_applied.get('ct_count_value', 'N/A')}\n\n"

    # Basic Counts
    summary += "**Read Counts (After Filtering & Deduplication):**\n"
    summary += f"- Total Reads Displayed: {len(reads_list_of_dicts)}\n\n"

    # Read Length Info
    if reads_list_of_dicts:
        lengths = [r_dict['length'] for r_dict in reads_list_of_dicts]
        summary += "**Read Length Distribution:**\n"
        summary += f"- Min Length: {min(lengths)}\n"
        summary += f"- Max Length: {max(lengths)}\n"
        summary += f"- Average Length: {np.mean(lengths):.2f}\n"
        summary += f"- Median Length: {np.median(lengths)}\n\n"
    else:
        summary += "**Read Length Distribution:** No reads passed filters.\n\n"

    # CG Content Info
    if reads_list_of_dicts:
        # CORRECTED: Use 'cg' key
        cg_contents = [r_dict['cg'] for r_dict in reads_list_of_dicts if r_dict.get('cg') is not None]
        if cg_contents:
             summary += "**Overall CG Content:**\n"
             summary += f"- Average CG: {np.mean(cg_contents):.3f}\n"
             summary += f"- Median CG: {np.median(cg_contents):.3f}\n"
             summary += f"- Std Dev CG: {np.std(cg_contents):.3f}\n\n"

             # CG for selected lengths
             if selected_lengths:
                 # CORRECTED: Use 'cg' key
                 cg_selected = [r_dict['cg'] for r_dict in reads_list_of_dicts if r_dict.get('length') in selected_lengths and r_dict.get('cg') is not None]
                 if cg_selected:
                     summary += f"**CG Content for {len(selected_lengths)} Selected Length(s):**\n"
                     summary += f"- Average CG (Selected): {np.mean(cg_selected):.3f}\n"
                     summary += f"- Median CG (Selected): {np.median(cg_selected):.3f}\n\n"
                 else:
                      summary += "**CG Content for Selected Length(s):** No reads with selected lengths found with CG data.\n\n" # More specific
             else:
                 summary += "**CG Content for Selected Length(s):** No lengths selected.\n\n"
        else:
            summary += "**Overall CG Content:** No CG data available in processed reads.\n\n" # More specific
    else:
        summary += "**Overall CG Content:** No reads passed filters.\n\n"

    # Alignment Stats Summary
    summary += "**Alignment Statistics (Based on Filtered Reads):**\n"
    if stats and stats_data.get('error') is None:
         summary += f"- Total Reads Analyzed: {stats.get('Total Reads', 'N/A')}\n"
         if stats.get('File Format') == 'bam/sam':
             summary += f"- Mapped Reads: {stats.get('Mapped Reads', 'N/A')} ({stats.get('Mapped Percentage', 0.0):.2f}%)\n"
             summary += f"- Duplicate Reads: {stats.get('Duplicate Reads', 'N/A')} ({stats.get('Duplicate Percentage', 0.0):.2f}%)\n"
             summary += f"- Total Mismatches (NM Sum): {stats.get('Total Mismatches (NM)', 'N/A')}\n" # This is from NM tags
             summary += f"- Avg Mismatches (NM) per Mapped Read: {stats.get('Avg Mismatches per Read', 0.0):.3f}\n"

             # Mismatch Categories from detailed mismatch list
             mismatch_details = stats.get('Mismatch Details', [])
             if mismatch_details:
                 categories = categorize_mismatches(mismatch_details)
                 counts_by_type = get_mismatch_counts_by_type(mismatch_details) # This gives specific A>C etc.
                 summary += "\n**Mismatch Breakdown (from MD/CIGAR analysis):**\n" # Clarify source
                 # CORRECTED KEY:
                 summary += f"- Total Distinct Mismatches Found: {categories.get('Total Analyzed Mismatches', 0)}\n"
                 summary += f"  - Transitions: {categories.get('Transitions', 0)}\n"
                 summary += f"    - C>T (Fwd Strand Ref): {categories.get('C>T Damage (Fwd)', 0)}\n"
                 summary += f"    - G>A (Fwd Strand Ref): {categories.get('G>A Damage (Rev)', 0)}\n" # Note: G>A (Rev) can mean C>T on original template
                 # To get "Other Transitions":
                 other_transitions = categories.get('Transitions', 0) - categories.get('C>T Damage (Fwd)', 0) - categories.get('G>A Damage (Rev)', 0)
                 summary += f"    - Other Transitions: {other_transitions}\n"
                 summary += f"  - Transversions: {categories.get('Transversions', 0)}\n"
                 if categories.get('Other Mismatches', 0) > 0:
                     summary += f"  - Other/Uncategorized: {categories.get('Other Mismatches', 0)}\n"
             summary += "\n"
         else:
             summary += "- (Detailed alignment statistics not applicable for FASTA/FASTQ)\n\n"
    else:
        summary += f"- Error loading stats: {stats_data.get('error', 'Unknown error') if stats_data else 'Stats data unavailable'}\n\n" # Handle if stats_data itself is None


    # Damage Patterns (Only if BAM/SAM)
    if stats and stats.get('File Format') == 'bam/sam':
        summary += "**Damage Patterns:**\n"
        five_pct, three_pct = get_damage_read_end(reads_list_of_dicts) # Recalculate based on final reads
        summary += "- 5' Terminus Composition:\n"
        for base, pct in sorted(five_pct.items()):
            if pct > 0.1: summary += f"    {base}: {pct:.1f}%\n"
        summary += "- 3' Terminus Composition:\n"
        for base, pct in sorted(three_pct.items()):
             if pct > 0.1: summary += f"    {base}: {pct:.1f}%\n"

      
    # Mismatch Frequency (Brief summary from Store)
        summary += "\n**Mismatch Frequency (First 10bp from Store):**\n"
        if mismatch_frequencies_data and not mismatch_frequencies_data.get('error') and mismatch_frequencies_data.get('frequencies'):
            freqs = mismatch_frequencies_data['frequencies']
            # Check if freqs is not None and has content before proceeding
            if freqs:
                if sequencing_type == 'single':
                    cpg_5 = freqs.get('C>T_CpG', [])[:10]
                    noncpg_5 = freqs.get('C>T_nonCpG', [])[:10]
                    other_5 = freqs.get('other', [])[:10]
                    if len(cpg_5) > 0 or len(noncpg_5) > 0 or len(other_5) > 0:
                        summary += f"- C>T (CpG, Avg): {np.mean(cpg_5) if len(cpg_5)>0 else 0:.2%}\n"
                        summary += f"- C>T (nonCpG, Avg): {np.mean(noncpg_5) if len(noncpg_5)>0 else 0:.2%}\n"
                        summary += f"- Other (Avg): {np.mean(other_5) if len(other_5)>0 else 0:.2%}\n"
                    else:
                        summary += "- No frequency data available for single-end analysis in summary.\n"
                else: # Paired-end
                    cpg_5 = freqs.get('5_prime', {}).get('C>T_CpG', [])[:10]
                    noncpg_5 = freqs.get('5_prime', {}).get('C>T_nonCpG', [])[:10]
                    other_5 = freqs.get('5_prime', {}).get('other', [])[:10]
                    cpg_3 = freqs.get('3_prime', {}).get('C>T_CpG', [])[:10]
                    noncpg_3 = freqs.get('3_prime', {}).get('C>T_nonCpG', [])[:10]
                    other_3 = freqs.get('3_prime', {}).get('other', [])[:10]
                    if any(len(arr)>0 for arr in [cpg_5, noncpg_5, other_5, cpg_3, noncpg_3, other_3]):
                        summary += f"- 5' C>T (CpG, Avg): {np.mean(cpg_5) if len(cpg_5)>0 else 0:.2%}\n"
                        summary += f"- 5' C>T (nonCpG, Avg): {np.mean(noncpg_5) if len(noncpg_5)>0 else 0:.2%}\n"
                        summary += f"- 5' Other (Avg): {np.mean(other_5) if len(other_5)>0 else 0:.2%}\n"
                        summary += f"- 3' C>T (CpG, Avg): {np.mean(cpg_3) if len(cpg_3)>0 else 0:.2%}\n"
                        summary += f"- 3' C>T (nonCpG, Avg): {np.mean(noncpg_3) if len(noncpg_3)>0 else 0:.2%}\n" # Typo fix: nonc_3 -> noncpg_3
                        summary += f"- 3' Other (Avg): {np.mean(other_3) if len(other_3)>0 else 0:.2%}\n"
                    else:
                        summary += "- No frequency data available for paired-end analysis in summary.\n"
            else: # freqs is None or empty
                summary += "- No frequency data generated or available in store for summary.\n"
        elif mismatch_frequencies_data and mismatch_frequencies_data.get('error'):
            summary += f"- Error in stored frequency data: {mismatch_frequencies_data.get('error')}\n"
        else:
            summary += "- Mismatch frequency data not available in store for summary.\n"
        summary += "\n"

    


    return summary


# --- Alignment Stats Offcanvas ---
@app.callback(
    Output('stats-offcanvas', 'is_open'),
    Input('stats-button', 'n_clicks'),
    State('stats-offcanvas', 'is_open'),
    prevent_initial_call=True,
)
def toggle_stats_offcanvas(n_clicks, is_open):
    """Toggles the visibility of the alignment statistics offcanvas."""
    if n_clicks:
        return not is_open
    return is_open

# Replace the existing update_alignment_stats_offcanvas callback with this:
@app.callback(
    [Output('alignment-stats-content', 'children'),
     Output('total-reads-stat', 'children'),
     Output('mapped-reads-stat', 'children'),
     Output('duplicate-reads-stat', 'children'),
     Output('avg-mismatches-stat', 'children'),
     Output('mismatch-details-content', 'children'),
     Output('stats-pie-chart', 'figure')],
    Input('alignment-stats-store', 'data'),
    Input('stats-offcanvas', 'is_open')
)
def update_alignment_stats_offcanvas(stats_data, is_open):
    """Populates the alignment statistics offcanvas content with new design."""
    if not is_open:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update
    
    if stats_data is None:
        empty_msg = html.Div("Alignment statistics will load when a file is processed.", 
                           className="text-center text-muted p-5")
        return empty_msg, "—", "—", "—", "—", empty_msg, go.Figure()
    
    if stats_data.get('error'):
        error_alert = dbc.Alert(f"Error loading stats: {stats_data['error']}", 
                              color="danger", className="m-3")
        return error_alert, "—", "—", "—", "—", error_alert, go.Figure()

    stats = stats_data.get('stats')
    if not stats:
        empty_msg = html.Div("No alignment statistics available.", 
                           className="text-center text-muted p-5")
        return empty_msg, "—", "—", "—", "—", empty_msg, go.Figure()

    # Extract stats values
    total_reads = stats.get('Total Reads', 0)
    mapped_reads = stats.get('Mapped Reads', 0)
    mapped_pct = stats.get('Mapped Percentage', 0.0)
    duplicate_reads = stats.get('Duplicate Reads', 0)
    duplicate_pct = stats.get('Duplicate Percentage', 0.0)
    avg_mismatches = stats.get('Avg Mismatches per Mapped Read', 0.0)

    # Format stat values for cards
    total_reads_display = f"{total_reads:,}"
    
    if stats.get('File Format') == 'bam/sam':
        mapped_display = f"{mapped_reads:,} ({mapped_pct:.1f}%)"
        duplicate_display = f"{duplicate_reads:,} ({duplicate_pct:.1f}%)"
        avg_mm_display = f"{avg_mismatches:.2f}"
    else:
        mapped_display = "N/A"
        duplicate_display = "N/A"
        avg_mm_display = "N/A"

    # Create mismatch details content
    mismatch_details_content = []
    if stats.get('File Format') == 'bam/sam':
        mismatch_details = stats.get('Mismatch Details', [])
        if mismatch_details:
            categories = categorize_mismatches(mismatch_details)
            
            # Create progress bars for mismatch types
            total_analyzed = categories.get('Total Analyzed Mismatches', 0)
            if total_analyzed > 0:
                mismatch_types = [
                    ('C>T Damage', categories.get('C>T Damage (Fwd)', 0), 'primary'),
                    ('G>A Damage', categories.get('G>A Damage (Rev)', 0), 'info'),
                    ('Other Transitions', 
                     categories.get('Transitions', 0) - 
                     categories.get('C>T Damage (Fwd)', 0) - 
                     categories.get('G>A Damage (Rev)', 0), 'success'),
                    ('Transversions', categories.get('Transversions', 0), 'warning'),
                    ('Other', categories.get('Other Mismatches', 0), 'secondary')
                ]
                
                for name, count, color in mismatch_types:
                    if count > 0:
                        pct = (count / total_analyzed) * 100
                        mismatch_details_content.append(
                            html.Div([
                                html.Div([
                                    html.Span(name, className="small"),
                                    html.Span(f"{count:,} ({pct:.1f}%)", 
                                            className="small float-end text-muted")
                                ], className="d-flex justify-content-between mb-1"),
                                dbc.Progress(
                                    value=pct,
                                    color=color,
                                    style={"height": "8px"},
                                    className="mb-3"
                                )
                            ])
                        )
                
                mismatch_details_content.append(
                    html.Hr(className="my-2")
                )
                mismatch_details_content.append(
                    html.P(f"Total analyzed: {total_analyzed:,} mismatches", 
                          className="small text-muted text-center mb-0")
                )
        else:
            mismatch_details_content = html.P("No mismatch details available", 
                                             className="text-muted text-center")
    else:
        mismatch_details_content = html.P("Mismatch analysis requires BAM/SAM format", 
                                         className="text-muted text-center")

    # Create the pie chart
    pie_chart_fig = create_mismatch_pie_chart(stats)
    # Update pie chart styling for dark theme
    pie_chart_fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        margin=dict(t=40, b=40, l=40, r=40),
        height=350,
        font=dict(color='#E0E0E0'),
        title_font_color='#E0E0E0'
    )

    # Main content placeholder (empty since we're using individual outputs)
    main_content = html.Div()

    return (main_content, total_reads_display, mapped_display, duplicate_display, 
            avg_mm_display, mismatch_details_content, pie_chart_fig)


# Add this callback to handle the sequencing type button group
@app.callback(
    [Output('seq-single', 'active'),
     Output('seq-double', 'active'),
     Output('sequencing-type', 'value', allow_duplicate=True)],
    [Input('seq-single', 'n_clicks'),
     Input('seq-double', 'n_clicks')],
    [State('seq-single', 'active'),
     State('seq-double', 'active')],
    prevent_initial_call=True
)
def update_sequencing_type_buttons(single_clicks, double_clicks, single_active, double_active):
    """Handle the sequencing type button group clicks."""
    ctx = dash.callback_context
    if not ctx.triggered:
        return single_active, double_active, 'single' if single_active else 'double'
    
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if button_id == 'seq-single':
        return True, False, 'single'
    elif button_id == 'seq-double':
        return False, True, 'double'
    
    return single_active, double_active, 'single' if single_active else 'double'

# --- Callback to handle mutual exclusivity of C>T checklist ---
@app.callback(
    Output('ct-checklist', 'value', allow_duplicate=True),
    Input('ct-checklist', 'value'),
    prevent_initial_call=True
)
def update_ct_checklist_exclusive(selected_values):
    """Ensure only one of 'only_ct', 'exclusively_ct', 'exclude_ct' is selected."""
    if not selected_values:
        return []

    exclusive_options = {'only_ct', 'exclusively_ct', 'exclude_ct'}
    last_triggered_id = dash.callback_context.triggered[0]['prop_id'].split('.')[0]
    last_value = selected_values[-1] # The one just clicked

    # If the last clicked item is one of the exclusive ones
    if last_value in exclusive_options:
        # Keep only the last clicked exclusive option and any non-exclusive ones (like 'subtract_ct')
        new_value = [val for val in selected_values if val not in exclusive_options or val == last_value]
        return new_value
    else:
        # If the last clicked was non-exclusive (e.g., 'subtract_ct'), just return the current selection
        return selected_values

@app.callback(
    Output('persistent-settings-store', 'data'), # Only this output now
    [Input('file-selector', 'value'),
     Input('mapq-range', 'value'),
     Input('min-base-quality', 'value'),
     Input('soft-clip-handling', 'value'),
     Input('nm-dropdown', 'value'),
     Input('ct-checklist', 'value'),
     Input('ct-count-dropdown', 'value'),
     Input('sequencing-type', 'value'),
     Input('apply-sequence-deduplication-switch', 'value'),
     Input('selected-lengths', 'data')],
    State('url', 'pathname'),
    prevent_initial_call=True
)
def save_persistent_settings(
    selected_file, mapq_range, min_base_quality, soft_clip_option,
    selected_nm, ct_checklist, ct_count_value, sequencing_type,
    apply_sequence_deduplication_value, selected_lengths_data,
    current_pathname
):
    if current_pathname != '/':
        # logger.debug(f"save_persistent_settings: Not on main page (pathname: {current_pathname}). Preventing update.")
        raise dash.exceptions.PreventUpdate

    settings_to_save = {
        'file-selector': selected_file,
        'mapq-range': mapq_range,
        'min-base-quality': min_base_quality,
        'soft-clip-handling': soft_clip_option,
        'nm-dropdown': selected_nm,
        'ct-checklist': ct_checklist,
        'ct-count-dropdown': ct_count_value,
        'sequencing-type': sequencing_type,
        'apply-sequence-deduplication-switch': apply_sequence_deduplication_value,
        'selected-lengths': selected_lengths_data,
    }
    
    # logger.debug(f"Saving to persistent-settings-store: {settings_to_save}")
    return settings_to_save # Only return settings, not metadata