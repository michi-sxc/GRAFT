import pysam
import os
import base64
import uuid
import logging
from Bio import SeqIO
from config import CUSTOM_TEMP_DIR
from .file_parser import parse_uploaded_file 
from .utils import calculate_cg_content

logger = logging.getLogger(__name__)

# --- Initial Read Loading ---

def load_reads_from_file(file_content, filename, soft_clip_option='show'):
    """
    Loads all reads from a file into a list of tuples.
    Handles different file formats and soft-clipping interpretation for length.
    Input: Base64 content, filename, soft_clip_option ('show', 'exclude_regions', 'exclude_all').
    Returns: List of tuples [(length, cg, nm, seq, read_obj, mapq)], file_format, temp_file_path, header_dict
    """
    logger.info(f"Loading reads from {filename} with soft_clip_option='{soft_clip_option}'")
    temp_file_path, file_format = parse_uploaded_file(file_content, filename)
    reads_data = []
    header_dict = None  # Store the file header

    try:
        if file_format in ['bam', 'sam']:
            # check_sq=False could be used
            with pysam.AlignmentFile(temp_file_path, "rb" if file_format == 'bam' else "r", check_sq=False) as infile:
                # Store the header for later use in exports
                header_dict = infile.header.to_dict() if infile.header else None
                
                for read in infile:
                    # --- Determine Sequence and Length based on Soft Clipping Option ---
                    sequence = read.query_sequence if read.query_sequence is not None else ""
                    read_length = 0
                    is_unmapped = read.is_unmapped

                    if is_unmapped:
                         # For unmapped reads, length is always the sequence length
                         read_length = len(sequence)
                    else:
                        # Apply soft clip handling for length calculation
                        if soft_clip_option == 'exclude_regions':
                             # Length of the aligned portion (excluding soft clips)
                             read_length = read.query_alignment_length
                        elif soft_clip_option == 'exclude_all':
                             read_length = read.query_length # with soft clips
                        else: # 'show' or default
                             read_length = read.query_length 

                    # --- Extract Other Information ---
                    cg = calculate_cg_content(sequence)
                    nm_tag = None
                    mapq = None
                    if not is_unmapped: # Get tags only for mapped reads
                         try:
                             nm_tag = read.get_tag('NM') if read.has_tag('NM') else None
                         except KeyError: pass # Ignore if tag doesn't exist
                         except ValueError: pass # Ignore if tag has wrong type
                         mapq = read.mapping_quality

                    reads_data.append((read_length, cg, nm_tag, sequence, read, mapq))

        elif file_format in ['fasta', 'fastq']: # use inferred format
            with open(temp_file_path, "rt") as handle:
                for record in SeqIO.parse(handle, file_format): # use inferred format
                    sequence = str(record.seq)
                    read_length = len(sequence)
                    cg = calculate_cg_content(sequence)
                    # FASTA/Q: nm=None, read_obj=record, mapq=None
                    reads_data.append((read_length, cg, None, sequence, record, None)) # Ensures 6 elements
        else:
             logger.error(f"Unsupported file format encountered during loading: {file_format}")
             # Should not happen if parse_uploaded_file is correct

        logger.info(f"Successfully loaded {len(reads_data)} reads from {filename}.")
        return reads_data, file_format, temp_file_path, header_dict

    except Exception as e:
        logger.error(f"Error loading reads from {temp_file_path} ({filename}): {e}", exc_info=True)
        # Attempt to clean up the temp file in case of error
        if os.path.exists(temp_file_path):
            try:
                os.remove(temp_file_path)
            except OSError:
                logger.warning(f"Could not remove temp file on error: {temp_file_path}")
        raise # Re-raise the exception to be handled by the callback


# --- Read Export ---

def export_selected_reads(
    read_tuple_list, # Input is the list of tuples already filtered (NM, C>T etc)
    selected_lengths, # List of integer lengths to keep
    output_file_path,
    output_format
):
    """
    Writes selected reads (based on length) to a new file.
    Input:
        read_tuple_list: List of (length, cg, nm, seq, read_obj, mapq) tuples.
        selected_lengths: List or set of read lengths to export.
        output_file_path: Path to write the output file.
        output_format: 'bam', 'sam', 'fasta', or 'fastq'.
    """
    logger.info(f"Exporting reads with lengths in {selected_lengths} to {output_file_path} (format: {output_format})")
    count = 0
    selected_lengths_set = set(selected_lengths) # Faster lookups

    # Determine header for BAM/SAM export if possible
    header = None
    first_read_obj = read_tuple_list[0][4] if read_tuple_list else None
    if output_format in ['bam', 'sam'] and isinstance(first_read_obj, pysam.AlignedSegment):
         # Attempt to get header from the first read object's file (might not always work)
        try:
            # This assumes the read object still holds a reference to its original file handle's header
            if hasattr(first_read_obj, 'header'):
                header = first_read_obj.header
            else: # Create a minimal default header
                logger.warning("Could not retrieve header from read object for export. Using minimal header.")
                header = pysam.AlignmentHeader.from_dict({}) # Empty header initially
        except Exception as e:
             logger.warning(f"Error getting header for export: {e}. Using minimal header.")
             header = pysam.AlignmentHeader.from_dict({})


    try:
        if output_format in ['bam', 'sam']:
            with pysam.AlignmentFile(output_file_path, "wb" if output_format == 'bam' else "w", header=header) as outfile:
                for length, _, _, _, read_obj, _ in read_tuple_list:
                    if length in selected_lengths_set and isinstance(read_obj, pysam.AlignedSegment):
                        outfile.write(read_obj)
                        count += 1
        elif output_format in ['fasta', 'fastq']:
             with open(output_file_path, "w") as outfile:
                 for length, _, _, _, record, _ in read_tuple_list:
                     if length in selected_lengths_set and isinstance(record, SeqIO.SeqRecord):
                          # Ensure FASTQ has quality if needed
                          if output_format == 'fastq' and not record.letter_annotations.get("phred_quality"):
                              logger.debug(f"Assigning dummy quality scores to {record.id} for FASTQ export.")
                              record.letter_annotations["phred_quality"] = [40] * len(record.seq) # Assign default quality
                          SeqIO.write(record, outfile, output_format)
                          count += 1
        else:
             raise ValueError(f"Unsupported output format for export: {output_format}")

        logger.info(f"Successfully exported {count} selected reads.")

    except Exception as e:
         logger.error(f"Error during read export to {output_file_path}: {e}", exc_info=True)
         # Clean up partially written file?
         if os.path.exists(output_file_path):
             try: os.remove(output_file_path)
             except OSError: pass
         raise


# --- File Merging ---

class FileMerger:
    """Handles merging of compatible sequence files."""
    def __init__(self):
        self.progress = 0
        self.status_message = ""
        self.current_files = []
        self.output_format = None
        self.temp_files_to_clean = []

    def _cleanup_temp_files(self):
        for temp_file in self.temp_files_to_clean:
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                    logger.info(f"Cleaned up merge temp file: {temp_file}")
                if os.path.exists(temp_file + '.bai'): # Clean index too
                     os.remove(temp_file + '.bai')
                     logger.info(f"Cleaned up merge temp index file: {temp_file + '.bai'}")
            except Exception as e:
                logger.error(f"Error removing merge temporary file {temp_file}: {e}")
        self.temp_files_to_clean = []


    def validate_files_for_merge(self, files_data):
        """Checks if a list of file data dicts can be merged."""
        if not files_data or len(files_data) < 2:
            return False, "Please select at least two files to merge."

        base_format = files_data[0]['format']
        if base_format not in ['bam', 'sam', 'fasta', 'fastq']:
             return False, f"Unsupported format for merging: {base_format}"

        compatible_formats = {
            'bam': ['bam', 'sam'], # Allow merging SAM into BAM
            'sam': ['sam', 'bam'], # Allow merging BAM into SAM
            'fastq': ['fastq'],
            'fasta': ['fasta']
        }
        # target format
        self.output_format = 'bam' if base_format in ['bam', 'sam'] else base_format

        allowed_set = compatible_formats.get(self.output_format, [self.output_format])

        for file_info in files_data:
            if file_info['format'] not in allowed_set:
                return False, f"Incompatible file formats for merging: Cannot merge '{file_info['format']}' with target '{self.output_format}'. Found in '{file_info['filename']}'."

        self.current_files = files_data
        return True, f"Files validated for merging into format: {self.output_format.upper()}"


    def merge_files(self, progress_callback=None):
        """
        Merges the validated files. Uses background task principles ideally.
        Returns base64 encoded content of the merged file.
        """
        if not self.current_files or self.output_format is None:
            raise RuntimeError("Files not validated before calling merge_files.")

        self._cleanup_temp_files() # Clean up from previous runs, settings are kept however
        merged_content_b64 = None
        output_temp_path = None

        try:
            if self.output_format in ['bam', 'sam']:
                output_temp_path = self._merge_bam_sam_files(progress_callback)
            elif self.output_format in ['fasta', 'fastq']:
                 output_temp_path = self._merge_fasta_fastq_files(progress_callback)
            else:
                 raise ValueError(f"Unexpected output format for merge: {self.output_format}")

            # Read the final merged file and encode it
            if output_temp_path and os.path.exists(output_temp_path):
                 with open(output_temp_path, 'rb') as f:
                     merged_content_b64 = base64.b64encode(f.read()).decode('utf-8')
                 logger.info(f"Successfully read and encoded merged file: {output_temp_path}")
            else:
                 raise IOError("Merged output file was not created or found.")

        except Exception as e:
             logger.error(f"Error during file merging: {e}", exc_info=True)
             self.status_message = f"Error: {e}"
             raise # Re-raise for the callback to handle
        finally:
            # cleanup happens even if encoding fails
            self._cleanup_temp_files() # Includes the final output temp path

        return merged_content_b64


    def _update_progress(self, value, message, callback_func):
        self.progress = value
        self.status_message = message
        if callback_func:
            try:
                callback_func(value, message)
            except Exception as e:
                 logger.warning(f"Progress callback function failed: {e}")


    def _merge_bam_sam_files(self, progress_callback):
        """Handles BAM/SAM merging logic."""
        input_paths = []
        total_reads_estimate = 0 # Estimate for progress
        self._update_progress(0, "Preparing files...", progress_callback)

        # 1. Decode and prepare input files, estimate total reads
        for i, file_data in enumerate(self.current_files):
             try:
                 temp_path, _ = parse_uploaded_file(file_data['content'], file_data['filename'])
                 self.temp_files_to_clean.append(temp_path)
                 input_paths.append(temp_path)
                 # Estimate reads (can be inaccurate but better than nothing for progress)
                 # total_reads_estimate += count_bam_reads(temp_path) # This is slow
                 logger.info(f"Prepared input file for merge: {temp_path}")
                 self._update_progress(int(5 * (i+1)/len(self.current_files)), f"Prepared {file_data['filename']}", progress_callback)
             except Exception as e:
                 raise IOError(f"Failed to prepare input file {file_data['filename']} for merging: {e}")

        # 2. Use pysam merge (simpler and often faster than manual iteration)
        #    pysam.merge requires file paths, not handles.
        #    Using a temporary output file path.
        merged_bam_path = os.path.join(CUSTOM_TEMP_DIR, f"merged_{uuid.uuid4()}.bam")
        self.temp_files_to_clean.append(merged_bam_path) # Must get cleaned

        self._update_progress(10, "Starting merge process...", progress_callback)
        try:
            pysam.merge('-c', '-p', '-f', '-o', merged_bam_path, *input_paths)

            self._update_progress(80, "Merge command completed, finalizing...", progress_callback)
        except pysam.SamtoolsError as e:
             logger.error(f"pysam.merge failed: {e}")
             raise RuntimeError(f"Samtools merge operation failed: {e.stderr}")

        # 3. Sort and Index the merged file (maybe good practice?)
        sorted_bam_path = os.path.join(CUSTOM_TEMP_DIR, f"sorted_merged_{uuid.uuid4()}.bam")
        self.temp_files_to_clean.append(sorted_bam_path)
        self.temp_files_to_clean.append(sorted_bam_path + '.bai') # Add index to cleanup list

        try:
             self._update_progress(85, "Sorting merged file...", progress_callback)
             pysam.sort("-o", sorted_bam_path, merged_bam_path)
             self._update_progress(95, "Indexing sorted file...", progress_callback)
             pysam.index(sorted_bam_path)
             self._update_progress(100, "Merge, sort, and index complete.", progress_callback)
             return sorted_bam_path # Return path of final sorted, indexed file
        except pysam.SamtoolsError as e:
             logger.error(f"pysam sort/index failed: {e}. Returning unsorted merged file.")
             # If sort/index fails, still return merged (but maybe unsorted) file
             self._update_progress(100, "Merge complete (sort/index failed).", progress_callback)
             return merged_bam_path # Return potentially unsorted merged file path
        except Exception as e:
             logger.error(f"Unexpected error during sort/index: {e}")
             self._update_progress(100, "Merge complete (error during sort/index).", progress_callback)
             return merged_bam_path


    def _merge_fasta_fastq_files(self, progress_callback):
        """Handles FASTA/FASTQ merging logic."""
        merged_records = []
        processed_files = 0
        self._update_progress(0, "Preparing files...", progress_callback)

        output_temp_path = os.path.join(CUSTOM_TEMP_DIR, f"merged_{uuid.uuid4()}.{self.output_format}")
        self.temp_files_to_clean.append(output_temp_path)

        used_ids = set() # Track IDs to avoid duplicates

        for i, file_data in enumerate(self.current_files):
            temp_file_path, input_format = parse_uploaded_file(file_data['content'], file_data['filename'])
            self.temp_files_to_clean.append(temp_file_path)
            logger.info(f"Processing {file_data['filename']} for FASTA/Q merge.")
            record_count = 0
            try:
                with open(temp_file_path, "rt") as handle:
                    for record in SeqIO.parse(handle, input_format):
                        original_id = record.id
                        counter = 1
                        # Unique IDs within merged file
                        while record.id in used_ids:
                            record.id = f"{original_id}_dup{counter}"
                            record.description = record.id
                            counter += 1

                        used_ids.add(record.id)
                        merged_records.append(record)
                        record_count += 1

                processed_files += 1
                self._update_progress(int(80 * processed_files / len(self.current_files)), f"Processed {record_count} records from {file_data['filename']}", progress_callback)

            except Exception as e:
                 raise IOError(f"Failed to parse input file {file_data['filename']} for merging: {e}")

        # Merged record to output file
        self._update_progress(90, f"Writing {len(merged_records)} merged records...", progress_callback)
        try:
            with open(output_temp_path, "w") as outfile:
                SeqIO.write(merged_records, outfile, self.output_format)
            self._update_progress(100, "Merge complete.", progress_callback)
            return output_temp_path
        except Exception as e:
            raise IOError(f"Failed to write merged output file {output_temp_path}: {e}")
