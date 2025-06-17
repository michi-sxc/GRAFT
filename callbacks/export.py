import dash
from dash import Input, Output, State, dcc
import os
import uuid
import zipfile
import shutil
import logging
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from app import app 
from config import CUSTOM_TEMP_DIR
from processing.read_operations import export_selected_reads
from plotting.utils import export_plot_to_file 

logger = logging.getLogger(__name__)

def export_reads_from_dicts(read_dicts_list, selected_lengths, output_file_path, output_format, file_header=None):
    """
    Export reads from serialized dictionaries while preserving original mapping information.
    This reconstructs the complete SAM/BAM format with proper reference context.
    """
    logger.info(f"Exporting {len(read_dicts_list)} read dictionaries with lengths in {selected_lengths} to {output_file_path} (format: {output_format})")
    count = 0
    selected_lengths_set = set(selected_lengths)

    # Filter reads by selected lengths
    filtered_reads = [r for r in read_dicts_list if r.get('length') in selected_lengths_set]
    
    if not filtered_reads:
        logger.warning("No reads match the selected lengths for export.")
        return

    try:
        if output_format in ['bam', 'sam']:
            # Use provided header or create a minimal one
            if file_header:
                header_dict = file_header
            else:
                # Fallback to minimal header
                header_dict = {
                    'HD': {'VN': '1.6', 'SO': 'unsorted'},
                    'SQ': []
                }
                logger.warning("No header provided for BAM/SAM export, using minimal header")
            
            with pysam.AlignmentFile(output_file_path, "wb" if output_format == 'bam' else "w", header=header_dict) as outfile:
                for read_dict in filtered_reads:
                    # Fully reconstruct AlignedSegment from stored information
                    read = pysam.AlignedSegment()
                    
                    # Basic read information
                    read.query_name = read_dict.get('name', f'read_{count}')
                    read.query_sequence = read_dict.get('seq', '')
                    
                    # Reconstruct all mapping information exactly as it was
                    read.flag = read_dict.get('flag', 4)  # Default to unmapped if missing
                    read.reference_id = read_dict.get('reference_id', -1)
                    read.reference_start = read_dict.get('reference_start', -1)
                    read.mapping_quality = read_dict.get('mapping_quality', 0)
                    read.cigarstring = read_dict.get('cigarstring', '*')
                    read.next_reference_id = read_dict.get('next_reference_id', -1)
                    read.next_reference_start = read_dict.get('next_reference_start', -1)
                    read.template_length = read_dict.get('template_length', 0)
                    
                    # Add quality scores
                    qualities = read_dict.get('query_qualities')
                    if qualities:
                        read.query_qualities = qualities
                    else:
                        # Only add dummy qualities if the read has sequence
                        seq_len = len(read_dict.get('seq', ''))
                        if seq_len > 0:
                            read.query_qualities = [30] * seq_len
                    
                    # Restore all original tags
                    tags = read_dict.get('tags', {})
                    for tag_name, tag_value in tags.items():
                        try:
                            read.set_tag(tag_name, tag_value)
                        except Exception as e:
                            logger.warning(f"Could not set tag {tag_name}={tag_value}: {e}")
                    
                    outfile.write(read)
                    count += 1

        elif output_format in ['fasta', 'fastq']:
            with open(output_file_path, "w") as outfile:
                for read_dict in filtered_reads:
                    # Create SeqRecord from dictionary
                    seq_str = read_dict.get('seq', '')
                    read_id = read_dict.get('name', f'read_{count}')
                    
                    record = SeqRecord(
                        Seq(seq_str),
                        id=read_id,
                        description=read_dict.get('description', '')
                    )
                    
                    # Add qualities for FASTQ
                    if output_format == 'fastq':
                        qualities = read_dict.get('query_qualities')
                        if qualities:
                            record.letter_annotations["phred_quality"] = qualities
                        else:
                            # Assign default quality scores if missing
                            record.letter_annotations["phred_quality"] = [30] * len(seq_str)
                    
                    SeqIO.write(record, outfile, output_format)
                    count += 1
        else:
            raise ValueError(f"Unsupported output format for export: {output_format}")

        logger.info(f"Successfully exported {count} selected reads.")

    except Exception as e:
        logger.error(f"Error during read export to {output_file_path}: {e}", exc_info=True)
        # Clean up partially written file
        if os.path.exists(output_file_path):
            try: 
                os.remove(output_file_path)
            except OSError: 
                pass
        raise

# --- Export Selected Reads ---
@app.callback(
    Output("download-selected-reads", "data"), # Use a unique download component ID
    Input("export-button", "n_clicks"),
    State('processed-data-store', 'data'), # Get the currently filtered reads
    State('selected-lengths', 'data'), # Get the lengths selected in the plot
    prevent_initial_call=True
)
def trigger_export_selected_reads(n_clicks, processed_data, selected_lengths):
    """Exports reads that match the current filters AND selected lengths."""
    if n_clicks is None: # If button hasn't been clicked, don't proceed
        raise dash.exceptions.PreventUpdate

    if not processed_data or processed_data.get('error') or not selected_lengths:
        logger.warning("Export triggered but no data/selection available or error occurred.")
        raise dash.exceptions.PreventUpdate
    
    reads_to_export = processed_data.get('reads', [])
    file_format = processed_data.get('format')
    original_filename = processed_data.get('filename', 'data')

    if not reads_to_export:
         logger.warning("No reads passed filters for export.")
         raise dash.exceptions.PreventUpdate

    # Determine output format (usually same as input, unless specified otherwise)
    # For simplicity, assume same format as input for now
    output_format = file_format
    if output_format not in ['bam', 'sam', 'fasta', 'fastq']:
         logger.error(f"Cannot export in unknown format: {output_format}")
         # Provide user feedback?
         raise dash.exceptions.PreventUpdate

    # Define output path
    base_name = os.path.splitext(original_filename)[0]
    output_filename = f"{base_name}_selected_L{len(selected_lengths)}_{uuid.uuid4().hex[:6]}.{output_format}"
    output_path = os.path.join(CUSTOM_TEMP_DIR, output_filename)

    try:
        # Use the new function that works with dictionaries and preserves mapping info
        file_header = processed_data.get('file_header')
        export_reads_from_dicts(
            read_dicts_list=reads_to_export,
            selected_lengths=selected_lengths,
            output_file_path=output_path,
            output_format=output_format,
            file_header=file_header
        )

        if os.path.exists(output_path):
             return dcc.send_file(output_path, filename=output_filename)
        else:
             logger.error(f"Exported file not found at expected path: {output_path}")
             # Provide feedback to user?
             return None

    except Exception as e:
        logger.error(f"Error during selected reads export: {e}", exc_info=True)
        # Provide feedback
        return None # Prevent download on error


# --- Export All Plots ---
@app.callback(
    Output("download-plots-zip", "data"), # Unique download component ID
    Input("export-plots-button", "n_clicks"),
    [State('read-length-histogram', 'figure'),
     State('overall-cg-histogram', 'figure'),
     State('cg-content-histogram', 'figure'),
     State('damage-patterns', 'figure'),
     State('mismatch-frequency-plot', 'figure'),
     State('mismatch-type-bar-chart', 'figure'),
     State('damage-pattern-plot', 'figure'),
     State('mapq-histogram', 'figure'),
     # Settings from settings offcanvas
     State('plot-width-input', 'value'),
     State('plot-height-input', 'value'),
     State('plot-scale-input', 'value'),
     State('plot-format-select', 'value'),
     # Need filename to name the zip file meaningfully
     State('processed-data-store', 'data')],
    prevent_initial_call=True
)
def trigger_export_all_plots(n_clicks, fig_len, fig_cg_overall, fig_cg_sel, fig_damage_ends,
                             fig_mismatch_freq, fig_mismatch_types, fig_damage_heatmap, fig_mapq,
                             export_width, export_height, export_scale, export_format, processed_data):
    """Exports all currently displayed plots as a ZIP archive."""
    if n_clicks is None: # No click event
        raise dash.exceptions.PreventUpdate

    base_filename = processed_data.get('filename', 'graft_analysis') if processed_data else 'graft_analysis'
    base_filename = os.path.splitext(base_filename)[0] # Remove extension

    figures_to_export = {
        f'{base_filename}_read_length': fig_len,
        f'{base_filename}_cg_overall': fig_cg_overall,
        f'{base_filename}_cg_selected': fig_cg_sel,
        f'{base_filename}_damage_ends': fig_damage_ends,
        f'{base_filename}_mismatch_frequency': fig_mismatch_freq,
        f'{base_filename}_mismatch_types': fig_mismatch_types,
        f'{base_filename}_damage_heatmap': fig_damage_heatmap,
        f'{base_filename}_mapq': fig_mapq,
    }

    # Temporary directory for the individual plot files
    temp_plot_dir = os.path.join(CUSTOM_TEMP_DIR, f"plots_export_{uuid.uuid4().hex[:8]}")
    os.makedirs(temp_plot_dir, exist_ok=True)
    logger.info(f"Created temporary plot export directory: {temp_plot_dir}")

    exported_files = []
    try:
        for name, fig_data in figures_to_export.items():
            if fig_data: # Check if figure has data
                 # Check if figure layout indicates no data
                 if fig_data.get('layout', {}).get('title', {}).get('text', '').endswith('(No Data)'):
                      logger.debug(f"Skipping export for '{name}' as it contains no data.")
                      continue

                 output_plot_path = os.path.join(temp_plot_dir, f"{name}.{export_format}")
                 try:
                    export_plot_to_file(
                        fig_data, output_plot_path,
                        width=export_width, height=export_height,
                        scale=export_scale, format=export_format
                    )
                    if os.path.exists(output_plot_path):
                        exported_files.append(output_plot_path)
                 except Exception as plot_err:
                      logger.error(f"Failed to export plot '{name}': {plot_err}", exc_info=True)
            else:
                logger.debug(f"Skipping export for '{name}' as figure data is empty.")


        if not exported_files:
            logger.warning("No plots were successfully exported.")
            # Provide user feedback? How from download callback?
            return None

        # Create ZIP file
        zip_filename = f"{base_filename}_plots_{export_format}.zip"
        zip_path = os.path.join(CUSTOM_TEMP_DIR, zip_filename)

        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file_path in exported_files:
                zipf.write(file_path, arcname=os.path.basename(file_path))
        logger.info(f"Created plot ZIP archive: {zip_path}")

        return dcc.send_file(zip_path, filename=zip_filename)

    except Exception as e:
        logger.error(f"Error during plot ZIP creation: {e}", exc_info=True)
        return None
    finally:
    
        if os.path.exists(temp_plot_dir):
            try:
                shutil.rmtree(temp_plot_dir)
                logger.info(f"Cleaned up temporary plot directory: {temp_plot_dir}")
            except Exception as clean_err:
                logger.error(f"Error cleaning up plot export directory {temp_plot_dir}: {clean_err}")
        # Zip file cleanup might be handled by dcc.send_file depending on context
