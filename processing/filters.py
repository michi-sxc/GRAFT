import pysam

from utils.helpers import reverse_complement
from .stats import count_ct_changes, get_mismatches_from_read # Relative import
import logging
import dust_module # deprecated, but kept for compatibility
from Bio.Seq import Seq

logger = logging.getLogger(__name__)

def apply_filters_to_read(read_tuple, file_format, filters):
    """
    Checks if a single read passes the specified filters.
    Input: read_tuple (length, cg, nm, seq, read_obj, mapq), file_format, filters (dict)
    Returns: True if the read passes, False otherwise.
    """
    length, cg, nm, seq, read_obj, mapq = read_tuple
    min_base_quality = filters.get('min_base_quality', 0)
    mapq_range = filters.get('mapq_range', (0, 255))
    soft_clip_option = filters.get('soft_clip_option', 'show')
    filter_ct = filters.get('filter_ct', None) # True=only_ct, False=exclude_ct, None=no filter
    exclusively_ct = filters.get('exclusively_ct', False)
    ct_count_value = filters.get('ct_count_value', 'any') # 'any' or integer
    selected_nm = filters.get('selected_nm', 'all') # 'all' or integer
    apply_dust = filters.get('apply_dust', False) # Add DUST filter option

    # --- DUST Filter ---
    if apply_dust and seq:
        try:
            # dust_mask returns the sequence or raises an error/returns specific value, but not in GitHub repo (deprecated)
            masked_seq = dust_module.dust_mask(seq)
            # Simple check: if the entire sequence is masked (all 'N'), filter it out
            if set(masked_seq) == {'N'}:
                return False
        except Exception as e:
             logger.warning(f"DUST module failed for read {getattr(read_obj, 'query_name', 'N/A')}: {e}")
             # pass

    # --- BAM/SAM Specific Filters ---
    if file_format in ['bam', 'sam']:
        if not isinstance(read_obj, pysam.AlignedSegment):
             logger.warning("Expected AlignedSegment for BAM/SAM format but got different type.")
             return True # Or False, depending on desired behavior for malformed data

        # --- Soft Clipping Filter ---
        has_soft_clipping = 'S' in (read_obj.cigarstring or "")
        if soft_clip_option == 'exclude_all' and has_soft_clipping:
            return False
        # Note: 'exclude_regions' affects length calculation BEFORE filtering, handled in load_reads

        # --- MAPQ Filter ---
        # Handle unmapped reads based on MAPQ filter. If min MAPQ > 0, unmapped reads (MAPQ=0 or 255) fail.
        # If MAPQ is None (shouldn't happen for AlignedSegment but defensive), treat as failing quality checks.
        current_mapq = read_obj.mapping_quality if read_obj.mapping_quality is not None else -1
        if current_mapq < mapq_range[0] or current_mapq > mapq_range[1]:
            return False

        # --- Base Quality Filter ---
        if min_base_quality > 0:
            base_qualities = read_obj.query_qualities
            if base_qualities is None:
                return False # No qualities, fails filter > 0
            # Find minimum quality efficiently
            min_q = min(base_qualities) if base_qualities else 0
            if min_q < min_base_quality:
                return False

        # --- C>T Change Presence Filter (filter_ct) ---
        if filter_ct is not None:
            # Check if C>T or G>A exists
            num_ct = count_ct_changes(read_obj)
            has_ct = num_ct > 0
            if filter_ct and not has_ct: # Require C>T, but none found
                return False
            if not filter_ct and has_ct: # Exclude C>T, but found
                return False

        # --- Exclusively C>T Filter ---
        if exclusively_ct:
            if not has_only_ct_mismatches(read_obj):
                return False

         # --- C>T Change Count Filter ---
        if ct_count_value != 'any':
             # Calculate only if needed
             num_ct = count_ct_changes(read_obj)
             if num_ct != int(ct_count_value):
                 return False

        # --- NM Tag Filter ---
        # This filter needs the "potentially adjusted" NM value if subtract_ct is True
        # Handle adjustment only after loading, so this filter applies to the original  or "adjusted" NM based on when it's called.
        # 'nm' in read_tuple is the value to filter on.
        if selected_nm != 'all':
            if nm is None or nm != selected_nm: # nm value passed in the tuple
                return False

    # --- FASTA/FASTQ Filters (not implemented) ---
    # Will do later

    # If we haven't returned False yet, the read passes
    return True


def has_only_ct_mismatches(read):
    """Checks if all mismatches in a read are C>T (forward) or G>A (reverse)."""
    if not isinstance(read, pysam.AlignedSegment) or read.is_unmapped:
        return False

    mismatches = get_mismatches_from_read(read)
    if not mismatches:
        return False # No mismatches means it doesn't have ONLY C>T

    for mismatch in mismatches:
        ref_base = mismatch['ref_base']
        read_base = mismatch['read_base']
        is_reverse = mismatch['is_reverse']

        is_damage = False
        if is_reverse:
            if ref_base == 'G' and read_base == 'A':
                is_damage = True
        else:
            if ref_base == 'C' and read_base == 'T':
                is_damage = True

        if not is_damage:
            return False # Found a mismatch that isn't C>T damage

    return True # All mismatches found were C>T damage


def adjust_nm_for_ct(read_tuple_list, subtract_ct):
    """
    Adjusts the NM value in a list of read tuples by subtracting C>T changes.
    Input: List of (length, cg, nm, seq, read_obj, mapq) tuples.
    Returns: New list with potentially adjusted 'nm' values.
    """
    if not subtract_ct:
        return read_tuple_list # No adjustment needed

    adjusted_list = []
    for length, cg, nm, seq, read_obj, mapq in read_tuple_list:
        adjusted_nm = nm
        if isinstance(read_obj, pysam.AlignedSegment) and nm is not None:
            ct_changes = count_ct_changes(read_obj)
            adjusted_nm = max(0, nm - ct_changes) # Ensure NM doesn't go below 0

        adjusted_list.append((length, cg, adjusted_nm, seq, read_obj, mapq))

    return adjusted_list


def filter_read_tuples(read_tuple_list, file_format, filters):
    """
    Applies a dictionary of filters to a list of read tuples.
    Returns the filtered list.
    """
    filtered_reads = [
        read_tuple for read_tuple in read_tuple_list
        if apply_filters_to_read(read_tuple, file_format, filters)
    ]
    return filtered_reads


def deduplicate_reads(read_tuple_list):
    unique_seqs = set()
    deduped_reads = []
    removed_count = 0

    for read_data in read_tuple_list:
        # read_data is (length, cg, nm, seq_str, read_obj, mapq_val)
        seq = read_data[3] # Sequence string is at index 3
        if not seq:  # Skip reads with no sequence
            deduped_reads.append(read_data) # Keep reads with no sequence? Or filter them earlier?
            continue

        seq_upper = seq.upper()
        # Calculate reverse complement
        try:
            rev_comp_seq_upper = reverse_complement(seq_upper)
        except Exception as e:
            logger.warning(f"Could not generate reverse complement for deduplication of '{seq_upper[:30]}...': {e}. Keeping read.")
            deduped_reads.append(read_data)
            continue

        # Check if either the sequence or its reverse complement is already seen
        if seq_upper not in unique_seqs and rev_comp_seq_upper not in unique_seqs:
            deduped_reads.append(read_data)
            unique_seqs.add(seq_upper)
            # No need to add rev_comp_seq_upper separately as check is symmetric
            # unique_seqs.add(rev_comp_seq_upper) # redundant with the above check
        else:
            removed_count +=1
            # Optional: Log some info about removed duplicates if count is very high
            if removed_count % 1000 == 0:
                logger.debug(f"Deduplication: Removed {removed_count} reads so far. Example duplicate seq: {seq_upper[:30]}...")


    logger.info(f"Deduplication: Input {len(read_tuple_list)} reads -> {len(deduped_reads)} unique reads ({removed_count} removed).")
    return deduped_reads