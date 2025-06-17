import dash
from dash import Input, Output, State, html, dcc
import os
import uuid
import zipfile
import logging
import pysam
from Bio import SeqIO

from app import app 
from config import CUSTOM_TEMP_DIR
from processing.file_parser import parse_uploaded_file 

logger = logging.getLogger(__name__)

# --- Enable/Disable Convert Button ---
@app.callback(
    Output('convert-button', 'disabled'),
    Input('convert-upload-data', 'contents')
)
def toggle_convert_button(list_of_contents):
    """Enable convert button only if files are uploaded."""
    return not list_of_contents

# --- Update Filename Display ---
@app.callback(
    Output('convert-upload-filenames', 'children'),
    Input('convert-upload-data', 'filename')
)
def update_converter_filenames_display(list_of_filenames):
    """Show the names of the files uploaded for conversion."""
    if not list_of_filenames:
        return ""
    if len(list_of_filenames) == 1:
        return f"File selected: {list_of_filenames[0]}"
    else:
         return f"{len(list_of_filenames)} files selected: " + ", ".join(list_of_filenames)


# --- File Conversion Logic ---
@app.callback(
    [Output('convert-download-file', 'data'),
     Output('convert-status', 'children')],
    Input('convert-button', 'n_clicks'),
    [State('convert-upload-data', 'contents'),
     State('convert-upload-data', 'filename'),
     State('convert-output-format', 'value'),
     # Add states for new advanced options (using dummy IDs for now, will match layout)
     State('converter-bam-filter-checklist', 'value'), # For BAM -> FASTA/Q
     State('converter-min-mapq', 'value'),             # For BAM -> FASTA/Q
     State('converter-dummy-quality-phred', 'value'), # For BAM/FASTA -> FASTQ (phred base for quality)
     State('converter-header-source', 'value'),       # For FASTA/Q -> BAM/SAM
     State('converter-header-include-checklist', 'value'), # For BAM <-> SAM
     State('converter-bam-compression', 'value'),          # For BAM output
     State('converter-dummy-quality-char-code', 'value') # For FASTA -> FASTQ (actual quality score)
     ], 
    prevent_initial_call=True,
    suppress_callback_exceptions=True
)
def perform_conversion(n_clicks, list_of_contents, list_of_filenames, output_format_value,
                       bam_filters, min_mapq, dummy_phred_offset,
                       header_source, header_include, bam_compression,
                       fasta_to_fastq_dummy_qual_score):
    if not n_clicks or not list_of_contents:
        raise dash.exceptions.PreventUpdate

    # Initialize defaults for options that might not be present in the DOM if their section is hidden
    # This is important because Dash sends None if a State's component isn't in the layout when callback fires
    bam_filters = bam_filters or [] # Default to empty list if None
    min_mapq = min_mapq if min_mapq is not None else 0
    dummy_phred_offset = dummy_phred_offset if dummy_phred_offset is not None else 33
    header_source = header_source or 'minimal_unaligned'
    header_include = header_include or ["include_header"] # Default to including if checkbox not there
    bam_compression = bam_compression if bam_compression is not None else 6
    fasta_to_fastq_dummy_qual_score = fasta_to_fastq_dummy_qual_score if fasta_to_fastq_dummy_qual_score is not None else 40


    converted_files_paths = []
    status_messages = []
    files_to_zip = {}

    for content, filename in zip(list_of_contents, list_of_filenames):
        input_temp_path = None
        output_temp_path = None
        try:
            logger.info(f"Converting {filename} to {output_format_value.upper()} with options: bam_filters={bam_filters}, min_mapq={min_mapq}, header_source='{header_source}', header_include={header_include}, bam_comp={bam_compression}, fasta_q_dummy_qual={fasta_to_fastq_dummy_qual_score}")
            input_temp_path, input_format = parse_uploaded_file(content, filename)

            base_name = os.path.splitext(filename)[0]
            if base_name.lower().endswith(('.gz', '.bz2')):
                 base_name = os.path.splitext(base_name)[0]
            output_filename = f"{base_name}_converted.{output_format_value}"
            output_temp_path = os.path.join(CUSTOM_TEMP_DIR, output_filename)

            if input_format == output_format_value:
                 status_messages.append(html.Div(f"Skipping '{filename}': Input and output formats are the same ({input_format}).", className="text-warning"))
                 continue

            # --- BAM/SAM <-> BAM/SAM ---
            if input_format in ['bam', 'sam'] and output_format_value in ['bam', 'sam']:
                read_mode = "rb" if input_format == 'bam' else "r"
                # Adjust write_mode for compression level if output is BAM
                write_mode = f"wb{bam_compression}" if output_format_value == 'bam' else "w"
                
                pysam_header = None
                if "include_header" in header_include:
                    with pysam.AlignmentFile(input_temp_path, read_mode, check_sq=False) as infile_header_check:
                        pysam_header = infile_header_check.header
                
                with pysam.AlignmentFile(input_temp_path, read_mode, check_sq=False) as infile:
                    with pysam.AlignmentFile(output_temp_path, write_mode, header=pysam_header) as outfile:
                        for read in infile:
                            outfile.write(read)

            # --- BAM/SAM -> FASTA/FASTQ ---
            elif input_format in ['bam', 'sam'] and output_format_value in ['fasta', 'fastq']:
                read_mode = "rb" if input_format == 'bam' else "r"
                with pysam.AlignmentFile(input_temp_path, read_mode, check_sq=False) as infile, \
                     open(output_temp_path, 'w') as outfile:
                    count = 0
                    filtered_out_count = 0
                    for read in infile:
                        # Apply filters
                        if 'exclude_unmapped' in bam_filters and read.is_unmapped:
                            filtered_out_count +=1; continue
                        if 'exclude_secondary' in bam_filters and read.is_secondary:
                            filtered_out_count +=1; continue
                        if 'exclude_supplementary' in bam_filters and read.is_supplementary:
                            filtered_out_count +=1; continue
                        if not read.is_unmapped and read.mapping_quality < min_mapq:
                            filtered_out_count +=1; continue
                        
                        seq = read.query_sequence
                        qual = read.query_qualities
                        read_name = read.query_name
                        if not seq: continue

                        if output_format_value == 'fasta':
                            outfile.write(f'>{read_name}\n{seq}\n')
                        elif output_format_value == 'fastq':
                            if qual is None: # No quality scores in BAM read
                                qual_str = chr(fasta_to_fastq_dummy_qual_score + dummy_phred_offset) * len(seq)
                            else:
                                # Qualities should be within valid ASCII Phred range (0-93 for Phred+33)
                                # dummy_phred_offset is 33 or 64
                                max_qual_phred = 93 if dummy_phred_offset == 33 else 60 # Approx for Phred+64
                                qual_str = ''.join([chr(min(max(q, 0), max_qual_phred) + dummy_phred_offset) for q in qual])
                            outfile.write(f'@{read_name}\n{seq}\n+\n{qual_str}\n')
                        count += 1
                    logger.info(f"Converted {count} reads (filtered out {filtered_out_count}) from {filename}.")


            # --- FASTA/FASTQ -> SAM/BAM ---
            elif input_format in ['fasta', 'fastq'] and output_format_value in ['bam', 'sam']:
                 write_mode = f"wb{bam_compression}" if output_format_value == 'bam' else "w"
                 header_dict = {'HD': {'VN': '1.6', 'SO': 'unsorted'}} # Keep 'unknown' or 'unsorted'
                 if header_source == 'minimal_unaligned': # Default
                     pass # header_dict is already minimal
                 # elif header_source == 'from_ref': # Not implemented yet
                 #     status_messages.append(html.Div(f"Header generation from reference not yet implemented for '{filename}'. Using minimal header.", className="text-warning"))

                 with pysam.AlignmentFile(output_temp_path, write_mode, header=header_dict) as outfile:
                    count = 0
                    for record in SeqIO.parse(input_temp_path, input_format):
                        a = pysam.AlignedSegment()
                        a.query_name = record.id
                        a.query_sequence = str(record.seq)
                        a.flag = 4  # Mark as unmapped
                        a.reference_id = -1
                        a.reference_start = -1
                        a.mapping_quality = 0
                        a.cigarstring = '*'

                        if input_format == 'fastq' and record.letter_annotations.get("phred_quality"):
                             a.query_qualities = record.letter_annotations["phred_quality"]
                        else: # FASTA or FASTQ missing quality
                             # Use fasta_to_fastq_dummy_qual_score for consistency
                             a.query_qualities = [fasta_to_fastq_dummy_qual_score] * len(record.seq)
                        outfile.write(a)
                        count += 1
                    logger.info(f"Converted {count} records from {filename} to unmapped {output_format_value.upper()}.")

            # --- FASTA <-> FASTQ ---
            elif input_format == 'fasta' and output_format_value == 'fastq':
                 with open(output_temp_path, 'w') as outfile:
                    count = 0
                    for record in SeqIO.parse(input_temp_path, 'fasta'):
                        # Use fasta_to_fastq_dummy_qual_score and dummy_phred_offset
                        record.letter_annotations["phred_quality"] = [fasta_to_fastq_dummy_qual_score] * len(record.seq)
                        SeqIO.write(record, outfile, 'fastq') # SeqIO.write handles Phred offset correctly for 'fastq' if qualities are integers
                        count += 1
                    logger.info(f"Converted {count} records from FASTA to FASTQ (dummy quality Phred score: {fasta_to_fastq_dummy_qual_score}).")
            elif input_format == 'fastq' and output_format_value == 'fasta':
                  with open(output_temp_path, 'w') as outfile:
                      count = SeqIO.convert(input_temp_path, 'fastq', outfile, 'fasta')
                  logger.info(f"Converted {count} records from FASTQ to FASTA.")

            else:
                 status_messages.append(html.Div(f"Conversion from {input_format.upper()} to {output_format_value.upper()} is not supported for '{filename}'.", className="text-danger"))
                 continue

            files_to_zip[output_filename] = output_temp_path
            status_messages.append(html.Div(f"Successfully converted '{filename}' to '{output_filename}'.", className="text-success"))

        except Exception as e:
            logger.error(f"Error converting {filename}: {e}", exc_info=True)
            status_messages.append(html.Div(f"Error converting '{filename}': {e}", className="text-danger"))
        finally:
            if input_temp_path and os.path.exists(input_temp_path):
                try: os.remove(input_temp_path)
                except OSError: logger.warning(f"Could not remove input temp file: {input_temp_path}")

    if not files_to_zip:
        return None, html.Div(status_messages)

    if len(files_to_zip) == 1:
        # Download single file directly
        output_filename, output_path = list(files_to_zip.items())[0]
        # Note: dcc.send_file handles cleanup if path is given directly
        # For safety, we might handle cleanup manually after download trigger
        return dcc.send_file(output_path, filename=output_filename), html.Div(status_messages)
    else:
        # Create a ZIP archive for multiple files
        zip_filename = f"converted_files_{uuid.uuid4()}.zip"
        zip_path = os.path.join(CUSTOM_TEMP_DIR, zip_filename)
        try:
            with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
                for output_filename, file_path in files_to_zip.items():
                    if os.path.exists(file_path):
                        zipf.write(file_path, arcname=output_filename)
                        logger.debug(f"Added {output_filename} to zip.")
                    else:
                         logger.warning(f"Converted file not found for zipping: {file_path}")

            # Send the zip file
            return dcc.send_file(zip_path, filename=zip_filename), html.Div(status_messages)

        except Exception as e:
            logger.error(f"Error creating zip file: {e}")
            status_messages.append(html.Div(f"Error creating ZIP file: {e}", className="text-danger"))
            return None, html.Div(status_messages)
        finally:
             # Clean up individual converted files AFTER zipping
             for file_path in files_to_zip.values():
                 if os.path.exists(file_path):
                     try: os.remove(file_path)
                     except OSError: logger.warning(f"Could not remove output temp file: {file_path}")
             # Clean up zip file after sending? dcc.send_file might handle this if path is local.



# --- Dynamically display advanced conversion options ---


@app.callback(
    [Output('adv-opts-bam-to-seq-div', 'style'),
     Output('adv-opts-seq-to-bam-div', 'style'),
     Output('adv-opts-bam-sam-div', 'style'),
     Output('adv-opts-fasta-to-fastq-div', 'style'),
     Output('converter-adv-options-initial-message', 'style')],
    Input('convert-output-format', 'value'),
    Input('convert-upload-data', 'filename')
)
def update_advanced_converter_options(output_format_value, list_of_filenames):
    # Default styles (all hidden initially)
    bam_to_seq_style = {'display': 'none', 'margin-top': '1rem'}
    seq_to_bam_style = {'display': 'none', 'margin-top': '1rem'}
    bam_sam_style = {'display': 'none', 'margin-top': '1rem'}
    fasta_to_fastq_style = {'display': 'none', 'margin-top': '1rem'}
    initial_message_style = {'display': 'block'} 

    if not output_format_value:
        return bam_to_seq_style, seq_to_bam_style, bam_sam_style, fasta_to_fastq_style, initial_message_style

    input_format_inferred = None

    if list_of_filenames:
        first_filename = list_of_filenames[0].lower()
        if any(ext in first_filename for ext in ['.bam', '.sam']):
            input_format_inferred = 'bam_sam'
        elif any(ext in first_filename for ext in ['.fa', '.fna', '.fasta']):
            input_format_inferred = 'fasta'
        elif any(ext in first_filename for ext in ['.fq', '.fastq']):
            input_format_inferred = 'fastq'

    # Determine which section to show
    show_section = False
    if output_format_value in ['fasta', 'fastq'] and input_format_inferred == 'bam_sam':
        bam_to_seq_style = {'display': 'block', 'margin-top': '1rem'}
        show_section = True
    elif output_format_value in ['bam', 'sam'] and input_format_inferred in ['fasta', 'fastq']:
        seq_to_bam_style = {'display': 'block', 'margin-top': '1rem'}
        show_section = True
    elif output_format_value in ['bam', 'sam'] and input_format_inferred == 'bam_sam':
        bam_sam_style = {'display': 'block', 'margin-top': '1rem'}
        show_section = True
    elif output_format_value == 'fastq' and input_format_inferred == 'fasta':
        fasta_to_fastq_style = {'display': 'block', 'margin-top': '1rem'}
        show_section = True
    # Note: No specific advanced options for FASTQ -> FASTA in this setup, so it will show the initial message or nothing if other conditions met.

    if show_section:
        initial_message_style = {'display': 'none'} # Hide initial message if a section is shown

    return bam_to_seq_style, seq_to_bam_style, bam_sam_style, fasta_to_fastq_style, initial_message_style