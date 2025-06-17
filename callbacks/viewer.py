import re
import dash
from dash import Input, Output, State, html
import pysam
from Bio import SeqIO
import math

from app import app 
from config import colors, RECORDS_PER_PAGE_VIEWER # Import config
from processing.file_parser import parse_uploaded_file
from processing.read_operations import load_reads_from_file # Reuse loading logic
from processing.filters import filter_read_tuples, adjust_nm_for_ct # Reuse filtering logic
from processing.stats import get_mismatches_from_read # For highlighting
from utils.helpers import complement_base, reverse_complement # For highlighting
import logging

logger = logging.getLogger(__name__)

# --- Viewer File Selection ---
@app.callback(
    Output('viewer-file-selector', 'options'),
    Input('file-store', 'data') # Triggered when files are uploaded/changed
)
def update_viewer_file_options(store_data):
    """Populates the file selector dropdown in the viewer tab."""
    if not store_data or not store_data.get('files'):
        return []
    return [{'label': f['filename'], 'value': f['filename']} for f in store_data['files']]


# --- Viewer Data Loading and Filtering ---
# This is similar to the main page but specific to the viewer's filters
@app.callback(
    Output('viewer-processed-data-store', 'data'),
    [Input('viewer-file-selector', 'value'),
     Input('viewer-mapq-range', 'value'),
     Input('viewer-min-base-quality', 'value'),
     Input('viewer-nm-dropdown', 'value'),
     Input('viewer-ct-checklist', 'value'),
     Input('viewer-ct-count-dropdown', 'value'),
     Input('viewer-softclip-display', 'value')],
    State('file-store', 'data'),
    prevent_initial_call=True,
)
def load_and_filter_viewer_data(
    selected_file, mapq_range, min_base_quality, selected_nm,
    ct_checklist, ct_count_value, softclip_display_filter_option, file_store_data
):
    if not selected_file or not file_store_data:
        return {'reads': [], 'format': None, 'filename': None, 'error': "No file selected", 'filters_applied': {}}

    file_data = next((f for f in file_store_data['files'] if f['filename'] == selected_file), None)
    if not file_data:
         return {'reads': [], 'format': None, 'filename': selected_file, 'error': "File data not found", 'filters_applied': {}}

    logger.info(f"Viewer: Processing {selected_file} for viewer display...")

    try:
        # --- 1. Load Raw Reads (tuples with objects) ---
        # Determine soft_clip_option for initial loading (affects length calculation if 'exclude_regions')
        # For viewer filtering, 'exclude' on softclip_display_filter_option means 'exclude_all'
        initial_load_soft_clip_opt = 'exclude_regions' if softclip_display_filter_option == 'exclude_regions_length_effect_only' else 'show' # Or 'show' generally
        
        raw_reads_data_tuples, file_format_inferred, _, file_header = load_reads_from_file(
            file_data['content'],
            file_data['filename'],
            soft_clip_option=initial_load_soft_clip_opt # Viewer mainly cares about full seq for display
        )

        # --- 2. Define Filters for Viewer ---
        filter_ct_flag_viewer = None
        if 'only_ct' in ct_checklist: filter_ct_flag_viewer = True
        elif 'exclude_ct' in ct_checklist: filter_ct_flag_viewer = False

        viewer_filters = {
            'mapq_range': tuple(mapq_range),
            'min_base_quality': int(min_base_quality) if min_base_quality is not None else 0,
            # soft_clip_option for filtering reads (not just for length calculation)
            'soft_clip_option': 'exclude_all' if softclip_display_filter_option == 'exclude' else 'show',
            'filter_ct': filter_ct_flag_viewer,
            'exclusively_ct': 'exclusively_ct' in ct_checklist,
            'ct_count_value': ct_count_value if ct_count_value == 'any' else int(ct_count_value),
            'selected_nm': selected_nm if selected_nm == 'all' else int(selected_nm),
            'subtract_ct_from_nm': 'subtract_ct' in ct_checklist
        }

        # --- 3. Adjust NM (operates on tuples with objects) ---
        nm_adjusted_reads_tuples = adjust_nm_for_ct(raw_reads_data_tuples, viewer_filters['subtract_ct_from_nm'])

        # --- 4. Apply Filters (operates on tuples with objects) ---
        filtered_reads_tuples = filter_read_tuples(nm_adjusted_reads_tuples, file_format_inferred, viewer_filters)

        # --- 5. Deduplication (Optional for viewer - might be slow and change read order) ---
        # For now, let's skip deduplication for the viewer to maintain original order after filtering,
        # unless specific use case requires it.
        # final_reads_tuples = deduplicate_reads(filtered_reads_tuples)
        final_reads_tuples = filtered_reads_tuples

        logger.info(f"Viewer: Prepared {len(final_reads_tuples)} reads for display after filtering.")

        # --- 6. Serialize Reads for viewer-processed-data-store ---
        # The viewer needs more info than the main page plots for detailed rendering.
        serializable_viewer_reads = []
        
        # Store the header for potential export use (same as main page)
        file_header = None
        if file_format_inferred in ['bam', 'sam'] and final_reads_tuples:
            first_read_obj = final_reads_tuples[0][4]
            if isinstance(first_read_obj, pysam.AlignedSegment):
                try:
                    file_header = first_read_obj.header.to_dict() if hasattr(first_read_obj, 'header') else None
                except:
                    logger.warning("Could not extract header from read object in viewer")
        
        for length, cg, nm_val, seq_str, read_obj, mapq_val in final_reads_tuples:
            record_info = {
                'length': length,
                'cg': cg,
                'nm': nm_val,
                'seq': seq_str,
                'mapq': mapq_val,
            }
            if isinstance(read_obj, pysam.AlignedSegment):
                # Store comprehensive mapping information (same as main page)
                mismatches_list = get_mismatches_from_read(read_obj)

                record_info.update({
                    'name': read_obj.query_name,
                    'flag': read_obj.flag,
                    'is_mapped': not read_obj.is_unmapped,
                    'is_reverse': read_obj.is_reverse,
                    'cigarstring': read_obj.cigarstring,
                    'md_tag': read_obj.get_tag('MD') if read_obj.has_tag('MD') else None,
                    'original_nm_tag': read_obj.get_tag('NM') if read_obj.has_tag('NM') else None,
                    'query_qualities': list(read_obj.query_qualities) if read_obj.query_qualities is not None else None,
                    'mismatches': mismatches_list,
                    # Additional mapping information for complete reconstruction
                    'reference_id': read_obj.reference_id,
                    'reference_start': read_obj.reference_start,
                    'reference_end': read_obj.reference_end,
                    'next_reference_id': read_obj.next_reference_id,
                    'next_reference_start': read_obj.next_reference_start,
                    'template_length': read_obj.template_length,
                    'mapping_quality': read_obj.mapping_quality,
                    'tags': dict(read_obj.get_tags()) if hasattr(read_obj, 'get_tags') else {},
                    'reference_name': read_obj.reference_name if hasattr(read_obj, 'reference_name') else None,
                    'next_reference_name': read_obj.next_reference_name if hasattr(read_obj, 'next_reference_name') else None,
                })
            elif isinstance(read_obj, SeqIO.SeqRecord):
                 record_info.update({
                     'name': read_obj.id,
                     'description': read_obj.description,
                     'is_reverse': False,
                     'query_qualities': list(read_obj.letter_annotations.get("phred_quality", [])) if read_obj.letter_annotations.get("phred_quality") else None,
                     'mismatches': []
                 })
            serializable_viewer_reads.append(record_info)

        return {
            'reads': serializable_viewer_reads,
            'format': file_format_inferred,
            'filename': selected_file,
            'filters_applied': viewer_filters,
            'file_header': file_header,  # Store header for export
            'error': None
        }

    except Exception as e:
        logger.error(f"Viewer: Error processing {selected_file}: {e}", exc_info=True)
        return {'reads': [], 'format': None, 'filename': selected_file, 'error': f"Error processing for viewer: {str(e)}", 'filters_applied': {}}


# --- Viewer Content Display ---
@app.callback(
    [Output('file-viewer-content', 'children'),
     Output('viewer-pagination', 'max_value'),
     Output('viewer-page-info', 'children'),
     Output('viewer-pagination', 'active_page')],
    [Input('viewer-processed-data-store', 'data'),
     Input('viewer-pagination', 'active_page'),
     Input('viewer-search-button', 'n_clicks'),
     Input('viewer-sequence-search-button', 'n_clicks'),
     Input('viewer-sort-field', 'value'),
     Input('viewer-sort-order', 'value'),
     Input('viewer-softclip-display', 'value')],
    [State('viewer-search-input', 'value'),
     State('viewer-sequence-search-input', 'value')],
    prevent_initial_call=True,
)
def display_viewer_content(
    processed_data_viewer, active_page, search_clicks, seq_search_clicks,
    sort_field, sort_order, softclip_display_style,
    search_term, seq_search_term
):
    if not processed_data_viewer or processed_data_viewer.get('error'):
         return html.Div(f"Error: {processed_data_viewer.get('error', 'No data loaded for viewer')}", className="text-danger m-3"), 1, "Page 0 of 0 (0 Records)", 1

    # processed_data_viewer['reads'] is now a list of dictionaries
    all_filtered_reads_dicts = processed_data_viewer.get('reads', [])
    file_format = processed_data_viewer.get('format')

    if not all_filtered_reads_dicts:
         return html.Div("No reads match the current filters for viewer.", className="text-muted m-3"), 1, "Page 0 of 0 (0 Records)", 1

    active_page = active_page or 1

    # --- Apply Search ---
    reads_to_display_dicts = all_filtered_reads_dicts # Start with all filtered reads
    triggered_id = dash.callback_context.triggered_id

    if triggered_id == 'viewer-search-button' and search_term:
        term = search_term.lower()
        reads_to_display_dicts = [
            r_dict for r_dict in reads_to_display_dicts if term in r_dict.get('name', '').lower()
        ]
        active_page = 1
    elif triggered_id == 'viewer-sequence-search-button' and seq_search_term:
        term = seq_search_term.upper()
        reads_to_display_dicts = [
             r_dict for r_dict in reads_to_display_dicts if term in r_dict.get('seq', '').upper()
        ]
        active_page = 1

    # --- Apply Sorting ---
    if sort_field and reads_to_display_dicts:
         reverse = sort_order == 'desc'
         if sort_field == 'name':
             key_func = lambda x: x.get('name', '')
         elif sort_field == 'length':
             key_func = lambda x: x.get('length', 0)
         elif sort_field == 'mismatches': # Using 'original_nm_tag' for consistency with typical NM usage
             key_func = lambda x: x.get('original_nm_tag') if x.get('original_nm_tag') is not None else float('inf')
         elif sort_field == 'mapq':
             key_func = lambda x: x.get('mapq') if x.get('mapq') is not None else -1
         else:
             key_func = None

         if key_func:
            try:
                reads_to_display_dicts.sort(key=key_func, reverse=reverse)
            except Exception as e:
                 logger.warning(f"Sorting failed in viewer: {e}")

    # --- Pagination ---
    total_records = len(reads_to_display_dicts)
    max_page = max(1, math.ceil(total_records / RECORDS_PER_PAGE_VIEWER))
    active_page = min(active_page, max_page)
    start_index = (active_page - 1) * RECORDS_PER_PAGE_VIEWER
    end_index = start_index + RECORDS_PER_PAGE_VIEWER
    page_reads_dicts = reads_to_display_dicts[start_index:end_index] # This is now a list of dicts

    page_info = f"Page {active_page} of {max_page} ({total_records} Records)"

    # --- Generate HTML Content ---
    content_divs = []
    if file_format in ['fasta', 'fastq']:
        # CORRECTED: Iterate through dictionaries
        for read_dict in page_reads_dicts:
            read_id = read_dict.get('name', 'N/A')
            seq_str = read_dict.get('seq', '')
            qualities = read_dict.get('query_qualities') # This is now a list of ints

            header_html = html.Div(read_id, className='viewer-read-header')
            sequence_html = html.Pre(seq_str, className='viewer-read-sequence')
            quality_html = None
            if file_format == 'fastq' and qualities is not None:
                # Convert Phred int list back to ASCII string for display
                qual_ascii_str = "".join([chr(q + 33) for q in qualities if q is not None])
                quality_html = html.Pre(f"+\n{qual_ascii_str}", className='viewer-read-quality')

            block_content = [header_html, sequence_html]
            if quality_html:
                block_content.append(quality_html)
            content_divs.append(html.Div(block_content, className='viewer-read-block'))

    elif file_format in ['bam', 'sam']:
        # CORRECTED: Iterate through dictionaries
        for read_dict in page_reads_dicts:
            # Access all necessary fields from the read_dict
            read_name = read_dict.get('name', 'N/A')
            flag_val = read_dict.get('flag', 'N/A')
            mapq_val = read_dict.get('mapq', 'N/A')
            read_len = read_dict.get('length', 'N/A') # The length calculated during initial load
            cigar_str = read_dict.get('cigarstring', 'N/A')
            nm_original_val = read_dict.get('original_nm_tag', 'N/A') # Show original NM from tag
            md_tag_val = read_dict.get('md_tag', 'N/A')

            header_line = html.Div([
                 html.Span(f"{read_name}", className="fw-bold me-3"),
                 html.Span(f"FLAGS: {flag_val}", className="me-2 small text-muted"),
                 html.Span(f"MAPQ: {mapq_val}", className="me-2 small text-muted"),
                 html.Span(f"LEN: {read_len}", className="me-2 small text-muted"),
             ], className="viewer-read-header")

            details_line = html.Div([
                  html.Span(f"CIGAR: {cigar_str}", className="me-2 small text-muted"),
                  html.Span(f"NM (tag): {nm_original_val}", className="me-2 small text-muted"),
                  html.Span(f"MD: {md_tag_val}", className="me-2 small text-muted"),
             ], className="viewer-read-details")

            # generate_highlighted_sequence_html now takes the read_dict
            seq_html_content = generate_highlighted_sequence_html(read_dict, softclip_display_style)

            content_divs.append(html.Div([header_line, details_line, seq_html_content], className='viewer-read-block'))
    else:
        return html.Div(f"Unsupported format for viewer: {file_format}", className="text-danger m-3"), 1, "Page 0 of 0 (0 Records)", 1

    if not content_divs:
        # This case handles if page_reads_dicts is empty after pagination logic (e.g. page out of bounds, though active_page should prevent this)
        # or if search/sort results in no reads for the current page view.
        no_reads_msg = "No reads to display for this page."
        if search_term or seq_search_term: # If a search was active
            no_reads_msg = "No reads found matching your search criteria for this page."
        content_divs = html.Div(no_reads_msg, className="text-muted m-3 text-center")


    return content_divs, max_page, page_info, active_page


# --- Helper for Sequence Highlighting ---
def generate_highlighted_sequence_html(read_dict, softclip_display_style):
    """
    Generates HTML for a read sequence with mismatch/softclip highlighting.
    Now takes a serialized read_dict.
    """
    if not read_dict or not read_dict.get('seq'):
        return html.Pre("Sequence Unavailable")

    read_seq = read_dict['seq']
    read_len = len(read_seq)
    query_qualities = read_dict.get('query_qualities', []) # Default to empty list
    
    # Mismatches are pre-calculated and stored in read_dict
    mismatches_list = read_dict.get('mismatches', [])
    # Create a map for quick lookup: read_pos -> (ref_base, read_base_mismatch, is_reverse_mismatch)
    mismatch_map = {
        m['read_pos']: (m['ref_base'], m['read_base'], m.get('is_reverse', read_dict.get('is_reverse', False)))
        for m in mismatches_list
    }

    # Identify soft-clipped positions from stored CIGAR string
    soft_clipped_indices = set()
    cigarstring = read_dict.get('cigarstring')
    if cigarstring: # Only parse if CIGAR is available (BAM/SAM)
        try:
            # Minimal CIGAR parsing for soft clips (pysam is not available here directly on read_obj)
            # A more robust CIGAR parser might be needed if complex CIGARs are common
            # For simplicity, we'll assume simple CIGARs or that pysam.AlignedSegment.cigartuples would be stored if needed.
            # Let's simulate parsing 'S' operations from the string.
            # This is a simplification; a proper CIGAR parser is complex.
            # For now, let's assume the soft_clip_option 'exclude' has filtered these reads if needed.
            # If 'highlight' is chosen, we need this.
            current_q_pos_cigar = 0
            operations = re.findall(r'(\d+)([MIDNSHPX=])', cigarstring)
            for length, op_char in operations:
                length = int(length)
                if op_char == 'S': # Soft clip
                    soft_clipped_indices.update(range(current_q_pos_cigar, current_q_pos_cigar + length))
                if op_char in ['M', 'I', 'S', '=', 'X']: # Consumes query
                    current_q_pos_cigar += length
        except Exception as e:
            logger.warning(f"Error parsing CIGAR '{cigarstring}' for soft clips in viewer: {e}")


    html_elements = []
    for i, base_char in enumerate(read_seq):
        is_softclipped_viewer = i in soft_clipped_indices

        if softclip_display_style == 'hide' and is_softclipped_viewer:
            continue
        # 'exclude' filter should be applied by load_and_filter_viewer_data based on 'soft_clip_option'

        style = {"padding": "0 1px"}
        title = f"Pos: {i}, Base: {base_char}"
        if query_qualities and i < len(query_qualities) and query_qualities[i] is not None:
            title += f", Qual: {query_qualities[i]}"

        base_class = "base-match"

        if is_softclipped_viewer and softclip_display_style == 'highlight':
            style.update({"backgroundColor": colors['muted'], "color": colors['surface'], "borderRadius": "2px"})
            base_class = "base-softclip"
            title += " (Soft Clip)"
        elif i in mismatch_map:
            ref_base_mm, read_base_mm, is_reverse_mm = mismatch_map[i]
            title += f", Ref: {ref_base_mm}"
            is_damage_type = False
            # Check for C>T or G>A type damage based on the stored mismatch info
            # read_dict['is_reverse'] refers to the overall read orientation
            is_read_reverse_oriented = read_dict.get('is_reverse', False)
            if is_read_reverse_oriented: # If read is mapped to reverse strand
                if ref_base_mm == 'G' and read_base_mm == 'A': is_damage_type = True # C>T on original template
            else: # Read mapped to forward strand
                if ref_base_mm == 'C' and read_base_mm == 'T': is_damage_type = True

            if is_damage_type:
                style.update({"color": "black", "backgroundColor": colors['highlight2'], "borderRadius": "2px", "fontWeight": "bold"})
                base_class = "base-damage-ct"
                title += " (C>T/G>A type)"
            else:
                style.update({"color": "white", "backgroundColor": colors['highlight'], "borderRadius": "2px"})
                base_class = "base-mismatch"
                title += " (Mismatch)"

        html_elements.append(html.Span(base_char, style=style, title=title, className=f"seq-base {base_class}"))

    return html.Pre(html_elements, style={"whiteSpace": "pre-wrap", "wordBreak": "break-all", "lineHeight": "1.4"})


# --- Viewer C>T Checklist Mutual Exclusivity ---
# Reusing the same logic as the main page
@app.callback(
    Output('viewer-ct-checklist', 'value'),
    Input('viewer-ct-checklist', 'value'),
    prevent_initial_call=True
)
def update_viewer_ct_checklist_exclusive(selected_values):
    """Ensure only one of 'only_ct', 'exclusively_ct', 'exclude_ct' is selected in viewer."""
    if not selected_values: return []
    exclusive_options = {'only_ct', 'exclusively_ct', 'exclude_ct'}
    last_value = selected_values[-1]
    if last_value in exclusive_options:
        return [val for val in selected_values if val not in exclusive_options or val == last_value]
    return selected_values

# --- Add Store Component to Viewer Layout ---
# Make sure the viewer layout in layouts/viewer.py includes:
# dcc.Store(id='viewer-processed-data-store')
