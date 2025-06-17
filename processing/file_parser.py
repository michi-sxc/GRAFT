import base64
import os
import uuid
import pysam
from Bio import SeqIO
from config import CUSTOM_TEMP_DIR
import logging
import gzip
import bz2

logger = logging.getLogger(__name__)

def parse_uploaded_file(content, filename):
    """
    Decodes uploaded file content and saves it to a temporary file.
    Determines file format based on extension. Handles basic compression.
    Returns the temporary file path and detected format.
    """
    if not content or not filename:
        raise ValueError("File content or filename is missing.")

    file_ext_parts = filename.lower().split('.')
    main_ext = file_ext_parts[-1]
    second_ext = file_ext_parts[-2] if len(file_ext_parts) > 1 else None

    compression = None
    if main_ext == 'gz':
        compression = 'gzip'
        file_format_ext = second_ext
        output_filename_base = filename[:-3] # Remove .gz
    elif main_ext == 'bz2':
        compression = 'bzip2'
        file_format_ext = second_ext
        output_filename_base = filename[:-4] # Remove .bz2
    else:
        file_format_ext = main_ext
        output_filename_base = filename

    # Determine file format
    if file_format_ext in ['fa', 'fna', 'fasta']:
        file_format = 'fasta'
    elif file_format_ext in ['fq', 'fastq']:
        file_format = 'fastq'
    elif file_format_ext == 'bam':
        file_format = 'bam'
    elif file_format_ext == 'sam':
        file_format = 'sam'
    else:
        raise ValueError(f"Unsupported file format extension: '{file_format_ext}' in filename '{filename}'")

    logger.info(f"Parsing file: {filename}, detected format: {file_format}, compression: {compression}")

    try:
        content_type, content_string = content.split(',', 1) # Split only once
        decoded = base64.b64decode(content_string)
    except Exception as e:
        logger.error(f"Error decoding base64 content for {filename}: {e}")
        raise ValueError(f"Could not decode file content for {filename}. Is it correctly base64 encoded?")

    # Create a unique temporary file path
    unique_filename = f"{uuid.uuid4()}_{output_filename_base}" # Use base name for clarity
    temp_file_path = os.path.join(CUSTOM_TEMP_DIR, unique_filename)

    try:
        # Handle decompression if necessary
        if compression == 'gzip':
            with gzip.open(temp_file_path, "wb") as fp:
                fp.write(decoded)
            logger.info(f"Decompressed gzip file to: {temp_file_path}")
        elif compression == 'bzip2':
             with bz2.open(temp_file_path, "wb") as fp:
                fp.write(decoded)
             logger.info(f"Decompressed bzip2 file to: {temp_file_path}")
        else:
            # Write directly if no compression
            with open(temp_file_path, "wb") as fp:
                fp.write(decoded)
            logger.info(f"Saved uploaded file to: {temp_file_path}")

    except Exception as e:
        logger.error(f"Error writing temporary file {temp_file_path}: {e}")
        # Attempt cleanup if write fails
        if os.path.exists(temp_file_path):
            os.remove(temp_file_path)
        raise IOError(f"Could not write temporary file for {filename}.")

    return temp_file_path, file_format

def count_bam_reads(bam_path):
    """Count total reads in a BAM/SAM file safely."""
    count = 0
    try:
        # Use view for potentially unindexed files, faster for simple counting usually
        # Check_sq=False is important for potentially incomplete headers
        with pysam.AlignmentFile(bam_path, "rb" if bam_path.endswith(".bam") else "r", check_sq=False) as bamfile:
            # pysam's count_reads() might be faster if indexed and complete
            # return bamfile.count_reads()
            for _ in bamfile:
                count += 1
        return count
    except ValueError as e:
        logger.warning(f"ValueError counting reads in {bam_path} (likely missing index or issue): {e}. Iterating manually.")
        # Fallback to pure iteration (already done above, but as a safeguard)
        count = 0
        with pysam.AlignmentFile(bam_path, "rb" if bam_path.endswith(".bam") else "r", check_sq=False) as bamfile:
            for _ in bamfile:
                count += 1
        return count
    except Exception as e:
        logger.error(f"Error counting reads in {bam_path}: {e}")
        return 0 # Return 0 on error

def count_fasta_fastq_reads(file_path, file_format):
    """Count reads in FASTA/FASTQ files."""
    count = 0
    try:
        # SeqIO.parse is memory efficient as it yields records
        with open(file_path, "rt") as handle: # Ensure text mode
            for _ in SeqIO.parse(handle, file_format):
                count += 1
        return count
    except Exception as e:
        logger.error(f"Error counting reads in {file_path} (format: {file_format}): {e}")
        return 0