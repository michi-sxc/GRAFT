import pysam
import numpy as np
from utils.helpers import complement_base
from config import MAX_MISMATCH_DISTANCE
import logging

logger = logging.getLogger(__name__)

# --- Mismatch Extraction ---
def get_mismatches_from_read(read):
    """
    Extract mismatches from a PySAM AlignedSegment.
    Prioritizes get_aligned_pairs(with_seq=True) as it's robust.
    MD tag parsing can be added as an alternative or for cross-validation
    but get_aligned_pairs is generally sufficient and simpler.
    """
    mismatches = []
    if not isinstance(read, pysam.AlignedSegment) or read.is_unmapped:
        return mismatches

    try:
        query_sequence = read.query_sequence
        query_qualities = read.query_qualities

        if not query_sequence:
            logger.debug(f"Read {read.query_name} has no query sequence. Skipping mismatch extraction.")
            return mismatches

        # Using get_aligned_pairs is usually the most reliable way
        # to get mismatches when reference sequence is available via the BAM.
        aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=True)
        if not aligned_pairs: # Should not happen if read is mapped and has sequence
            logger.debug(f"Read {read.query_name} yielded no aligned pairs. Skipping mismatch extraction.")
            return mismatches

        for q_pos, r_pos, ref_base_from_bam in aligned_pairs:
            # q_pos: position in query (read)
            # r_pos: position in reference
            # ref_base_from_bam: reference base at r_pos from BAM's perspective (e.g., via CIGAR and MD)

            if q_pos is not None and r_pos is not None and ref_base_from_bam is not None:
                # We have an alignment of a query base to a reference base
                if q_pos >= len(query_sequence): # Safety check
                    logger.warning(f"Query position {q_pos} out of bounds for read {read.query_name} (len {len(query_sequence)}).")
                    continue

                read_base = query_sequence[q_pos].upper()
                ref_base = ref_base_from_bam.upper()

                if read_base != ref_base and ref_base != 'N' and read_base != '.': # Standard mismatch
                    base_qual = 0
                    if query_qualities and q_pos < len(query_qualities):
                        base_qual = query_qualities[q_pos]

                    mismatches.append({
                        'read_pos': q_pos,       # 0-based position in the read sequence
                        'ref_pos': r_pos,        # 0-based position in the reference sequence
                        'read_base': read_base,
                        'ref_base': ref_base,
                        'base_quality': base_qual,
                        'mapping_quality': read.mapping_quality,
                        'is_reverse': read.is_reverse
                    })
    except Exception as e:
        logger.warning(f"Error extracting mismatches for read {read.query_name}: {e}", exc_info=False) # Keep exc_info False for less verbose logs unless debugging this specifically
    return mismatches


# --- Alignment Statistics ---
def calculate_alignment_stats(processed_reads_data_tuples):
    """
    Calculates alignment statistics from processed read data tuples.
    Input: List of tuples [(length, cg, nm, seq, read_obj, mapq)],
           where read_obj is a pysam.AlignedSegment or Bio.SeqRecord.
    """
    total_reads = len(processed_reads_data_tuples)
    mapped_reads_count = 0
    duplicate_reads_count = 0
    total_nm_sum = 0
    mismatch_details_list = []
    mapq_scores_list = []

    if total_reads == 0:
        return {
            'Total Reads': 0, 'Mapped Reads': 0, 'Mapped Percentage': 0.0,
            'Duplicate Reads': 0, 'Duplicate Percentage': 0.0,
            'Total Mismatches (NM)': 0, 'Avg Mismatches per Mapped Read': 0.0,
            'MAPQ Scores': [], 'Mismatch Details': [], 'File Format': 'N/A'
        }

    # Infer format from the first read object
    first_read_obj = processed_reads_data_tuples[0][4] if processed_reads_data_tuples else None
    file_format_inferred = 'bam/sam' if isinstance(first_read_obj, pysam.AlignedSegment) else 'fasta/fastq'

    if file_format_inferred == 'bam/sam':
        for length, cg, nm, seq, read_obj, mapq in processed_reads_data_tuples:
            if not isinstance(read_obj, pysam.AlignedSegment): continue # Should not happen

            if not read_obj.is_unmapped:
                mapped_reads_count += 1
                if mapq is not None:
                    mapq_scores_list.append(mapq)
                if nm is not None: # nm is the original NM tag value from the tuple
                    total_nm_sum += nm
                mismatch_details_list.extend(get_mismatches_from_read(read_obj))

            if read_obj.is_duplicate or (read_obj.flag & 1024): # BAM_FDUP for PCR/optical duplicates
                duplicate_reads_count += 1
    else: # FASTA/FASTQ
        # These stats are not applicable or directly derivable
        mapped_reads_count = 'N/A'
        duplicate_reads_count = 'N/A'
        total_nm_sum = 'N/A'

    mapped_percentage = (mapped_reads_count / total_reads) * 100 if isinstance(mapped_reads_count, int) and total_reads > 0 else 0.0
    duplicate_percentage = (duplicate_reads_count / total_reads) * 100 if isinstance(duplicate_reads_count, int) and total_reads > 0 else 0.0
    avg_mismatches_per_mapped = (total_nm_sum / mapped_reads_count) if isinstance(total_nm_sum, int) and isinstance(mapped_reads_count, int) and mapped_reads_count > 0 else 0.0

    return {
        'Total Reads': total_reads,
        'Mapped Reads': mapped_reads_count,
        'Mapped Percentage': mapped_percentage,
        'Duplicate Reads': duplicate_reads_count,
        'Duplicate Percentage': duplicate_percentage,
        'Total Mismatches (NM)': total_nm_sum, # Sum of original NM tags from input
        'Avg Mismatches per Mapped Read': avg_mismatches_per_mapped,
        'MAPQ Scores': mapq_scores_list,
        'Mismatch Details': mismatch_details_list, # List of mismatch dicts from get_mismatches_from_read
        'File Format': file_format_inferred
    }


# --- Mismatch Categorization ---
def categorize_mismatches(mismatch_details_list):
    """Categorizes mismatches from a list of mismatch detail dictionaries."""
    categories = {
        'Transitions': 0, 'Transversions': 0,
        'C>T Damage (Fwd)': 0, 'G>A Damage (Rev)': 0, # More specific naming
        'Other Mismatches': 0,
        'Total Analyzed Mismatches': 0 
    }
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}

    for mismatch in mismatch_details_list:
        ref = mismatch.get('ref_base', 'N').upper()
        read_b = mismatch.get('read_base', 'N').upper()
        # is_rev = mismatch.get('is_reverse', False) # Not strictly needed for categorization here

        if ref == 'N' or read_b == 'N' or ref == read_b:
            continue
        categories['Total Analyzed Mismatches'] += 1

        is_transition = (ref in purines and read_b in purines) or \
                        (ref in pyrimidines and read_b in pyrimidines)
        is_transversion = (ref in purines and read_b in pyrimidines) or \
                          (ref in pyrimidines and read_b in purines)

        if is_transition:
            categories['Transitions'] += 1
            if ref == 'C' and read_b == 'T':
                categories['C>T Damage (Fwd)'] += 1
            elif ref == 'G' and read_b == 'A':
                categories['G>A Damage (Rev)'] += 1 # This is G>A on the ref strand, which is C>T on the original '-' DNA strand if the read mapped to '+'
        elif is_transversion:
            categories['Transversions'] += 1
        else:
            categories['Other Mismatches'] += 1
            logger.debug(f"Mismatch {ref}>{read_b} categorized as 'Other'.")
    return categories


def get_mismatch_counts_by_type(mismatch_details_list):
    """Counts occurrences of each specific mismatch type (A>C, A>G, etc.)."""
    counts = {}
    valid_bases = ['A', 'C', 'G', 'T']
    for r_base in valid_bases:
        for q_base in valid_bases:
            if r_base != q_base:
                counts[f"{r_base}>{q_base}"] = 0

    for mismatch in mismatch_details_list:
        ref = mismatch.get('ref_base', 'N').upper()
        read_b = mismatch.get('read_base', 'N').upper()
        key = f"{ref}>{read_b}"
        if key in counts: # Only count valid base changes
            counts[key] += 1
    return counts


# --- Damage Pattern Calculation ---
def get_damage_read_end(processed_reads_list_of_dicts):
    """
    Calculates base frequency at 5' and 3' ends.
    Input: List of *serialized* read dictionaries.
           Requires 'seq' and 'is_reverse' keys in the dictionaries.
    """
    five_prime_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    three_prime_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    total_reads_analyzed = 0

    # Check if the necessary keys are present in the first record (if any)
    if not processed_reads_list_of_dicts or \
       ('seq' not in processed_reads_list_of_dicts[0]) or \
       ('is_reverse' not in processed_reads_list_of_dicts[0]): # is_reverse might not be there for FASTA/Q
        logger.warning("get_damage_read_end: 'seq' or 'is_reverse' missing from serialized read data. Cannot calculate end patterns.")
        # Return zeroed percentages
        return {b: 0.0 for b in five_prime_counts}, {b: 0.0 for b in three_prime_counts}

    for read_dict in processed_reads_list_of_dicts:
        seq = read_dict.get('seq')
        # is_reverse is crucial for BAM/SAM. For FASTA/Q, it defaults to False (or isn't present).
        is_reverse = read_dict.get('is_reverse', False)

        if not seq or len(seq) == 0:
            continue

        base_5p_read = seq[0].upper()        # 5'-most base as sequenced
        base_3p_read = seq[-1].upper()      # 3'-most base as sequenced

        # Determine biological 5' and 3' ends
        if is_reverse: # This was a read from the reverse strand of the original molecule
            # Biological 5' end of original molecule is complement of 3'-most base of sequenced read
            actual_bio_5p_base = complement_base(base_3p_read)
            # Biological 3' end of original molecule is complement of 5'-most base of sequenced read
            actual_bio_3p_base = complement_base(base_5p_read)
        else: # Forward strand or FASTA/FASTQ
            actual_bio_5p_base = base_5p_read
            actual_bio_3p_base = base_3p_read

        five_prime_counts[actual_bio_5p_base] = five_prime_counts.get(actual_bio_5p_base, 0) + 1
        three_prime_counts[actual_bio_3p_base] = three_prime_counts.get(actual_bio_3p_base, 0) + 1
        total_reads_analyzed += 1

    five_prime_pct = {b: (c / total_reads_analyzed) * 100 if total_reads_analyzed > 0 else 0
                      for b, c in five_prime_counts.items()}
    three_prime_pct = {b: (c / total_reads_analyzed) * 100 if total_reads_analyzed > 0 else 0
                       for b, c in three_prime_counts.items()}

    return five_prime_pct, three_prime_pct


def calculate_mismatch_frequency_vs_end(processed_reads_data_tuples, sequencing_type='single'):
    """
    Calculates C>T (CpG/non-CpG) and other mismatch frequencies vs distance from read ends.
    Input: List of tuples [(length, cg, nm, seq, read_obj, mapq)], where read_obj is a pysam.AlignedSegment.
           This function is called *before* serialization.
    """
    non_duplicate_reads_tuples = []
    bam_duplicates_skipped = 0
    for read_tuple in processed_reads_data_tuples:
        read_obj = read_tuple[4] # read_obj is the 5th element
        if isinstance(read_obj, pysam.AlignedSegment):
            if read_obj.is_duplicate: # Checks for BAM flag 0x400
                bam_duplicates_skipped += 1
                continue
        non_duplicate_reads_tuples.append(read_tuple)
    
    if bam_duplicates_skipped > 0:
        logger.info(f"Mismatch Freq Calc: Skipped {bam_duplicates_skipped} reads marked as BAM duplicates.")
    
    if not non_duplicate_reads_tuples:
        logger.warning("Mismatch Freq Calc: No non-duplicate reads available for calculation.")
        # Return an empty structure or None, consistent with how plotting handles no data
        # For single-stranded:
        if sequencing_type == 'single':
            empty_freq_single = {
                'total_bases': np.array([]), # Or np.zeros(MAX_MISMATCH_DISTANCE) if plot expects full length
                'C>T_CpG_counts': np.array([]), 'C>T_CpG': np.array([]),
                'C_in_CpG_sites': np.array([]),
                'C>T_nonCpG_counts': np.array([]), 'C>T_nonCpG': np.array([]),
                'C_not_in_CpG_sites': np.array([]),
                'other_counts': np.array([]), 'other': np.array([]),
                'G_sites': np.array([]), 'A_sites': np.array([]), 'T_sites': np.array([])
            }
            return empty_freq_single
        # For double-stranded:
        else:
            empty_freq_double = {
                '5_prime': {
                    'total_bases': np.array([]),
                    'C>T_CpG_counts': np.array([]), 'C>T_CpG': np.array([]), 'C_in_CpG_sites': np.array([]),
                    'C>T_nonCpG_counts': np.array([]), 'C>T_nonCpG': np.array([]), 'C_not_in_CpG_sites': np.array([]),
                    'other_counts': np.array([]), 'other': np.array([]),
                    'G_sites': np.array([]), 'A_sites': np.array([]), 'T_sites': np.array([])
                },
                '3_prime': {
                    'total_bases': np.array([]),
                    'C>T_CpG_counts': np.array([]), 'C>T_CpG': np.array([]), 'C_in_CpG_sites': np.array([]),
                    'C>T_nonCpG_counts': np.array([]), 'C>T_nonCpG': np.array([]), 'C_not_in_CpG_sites': np.array([]),
                    'other_counts': np.array([]), 'other': np.array([]),
                    'G_sites': np.array([]), 'A_sites': np.array([]), 'T_sites': np.array([])
                }
            }
            return empty_freq_double

    max_dist = MAX_MISMATCH_DISTANCE
    ends_map = ['5_prime', '3_prime']
    mut_categories = ['C>T_CpG', 'C>T_nonCpG', 'other']
    # Denominators: Count of C's in CpG, C's not in CpG, and all other bases (A,T,G)
    # These are reference-based sites.
    site_categories = ['C_in_CpG', 'C_not_in_CpG', 'G_base', 'A_base', 'T_base']

    # Initialize counts structure
    counts = {
        end: {
            **{f'{mut}_mismatches': np.zeros(max_dist, dtype=int) for mut in mut_categories},
            **{f'{site}_sites': np.zeros(max_dist, dtype=int) for site in site_categories},
            'total_considered_sites': np.zeros(max_dist, dtype=int) # All A,C,G,T sites at position
        } for end in ends_map
    }

    for length, cg, nm, seq_str, read_obj, mapq in non_duplicate_reads_tuples:
        if not isinstance(read_obj, pysam.AlignedSegment) or read_obj.is_unmapped:
            # This check is somewhat redundant now due to the is_duplicate check above,
            # but good as a safeguard if non-AlignedSegment objects somehow get here.
            continue

        query_seq = read_obj.query_sequence
        if not query_seq: continue
        read_len = len(query_seq)
        if read_len == 0: continue

        try:
            # with_seq=True provides the reference base from the BAM file's perspective
            aligned_pairs = read_obj.get_aligned_pairs(matches_only=False, with_seq=True)
        except Exception as e:
            logger.warning(f"Could not get aligned pairs for {read_obj.query_name}: {e}")
            continue

        # Determine CpG context along the reference part of the alignment
        # This requires looking ahead. We can build a map of ref_pos to its context.
        ref_pos_to_context = {} # Key: ref_pos, Value: {'ref_base': 'C', 'next_ref_base': 'G'} or similar
        # Iterate once to get all ref bases and their next ref base if available
        ref_alignment_info = {} # store ref_base at ref_pos
        for qp, rp, rb in aligned_pairs:
            if rp is not None and rb is not None:
                ref_alignment_info[rp] = rb.upper()

        for rp_current, current_rb_upper in ref_alignment_info.items():
            # Check if current_rb_upper is 'C' and if the *next* reference position has 'G'
            next_rb_upper = ref_alignment_info.get(rp_current + 1)
            if current_rb_upper == 'C' and next_rb_upper == 'G':
                ref_pos_to_context[rp_current] = 'CpG' # This C is in a CpG context
            elif current_rb_upper == 'C':
                 ref_pos_to_context[rp_current] = 'C_non_CpG'
            # Context for G needed if considering G>A on reverse strand
            # (which is C>T on original molecule)
            # A 'G' is part of a CpG if the *previous* base was 'C'
            prev_rb_upper = ref_alignment_info.get(rp_current -1)
            if current_rb_upper == 'G' and prev_rb_upper == 'C':
                 ref_pos_to_context[rp_current] = 'G_in_CpG' # This G is in a CpG context (..CG..)
            elif current_rb_upper == 'G':
                 ref_pos_to_context[rp_current] = 'G_non_CpG'


        # Main iteration for counting sites and mismatches
        for q_idx, (q_pos, ref_pos, ref_base_bam) in enumerate(aligned_pairs):
            if q_pos is None or ref_pos is None or ref_base_bam is None:
                continue
            if q_pos >= read_len: continue

            read_base_query = query_seq[q_pos].upper()
            ref_base_original = ref_base_bam.upper()

            # Distances from read ends
            dist_from_5p = q_pos
            dist_from_3p = read_len - 1 - q_pos

            current_ref_context = ref_pos_to_context.get(ref_pos)

            for end_idx, dist_val in enumerate([dist_from_5p, dist_from_3p]):
                if dist_val < max_dist:
                    end_label = ends_map[end_idx]
                    counts[end_label]['total_considered_sites'][dist_val] += 1

                    # --- Site Counting (based on reference base) ---
                    is_ref_c_in_cpg = (ref_base_original == 'C' and current_ref_context == 'CpG')
                    is_ref_c_not_in_cpg = (ref_base_original == 'C' and current_ref_context == 'C_non_CpG')
                    is_ref_g = (ref_base_original == 'G') # Could further split by G_in_CpG vs G_non_CpG
                    is_ref_a = (ref_base_original == 'A')
                    is_ref_t = (ref_base_original == 'T')

                    if read_obj.is_reverse:
                        # Biological C site (G on ref strand for reverse read)
                        if is_ref_g:
                            if ref_pos_to_context.get(ref_pos) == 'G_in_CpG': # G is the G of a C-G pair
                                counts[end_label]['C_in_CpG_sites'][dist_val] += 1
                            else:
                                counts[end_label]['C_not_in_CpG_sites'][dist_val] += 1
                        elif is_ref_c_in_cpg: counts[end_label]['G_base_sites'][dist_val] +=1 # Bio G
                        elif is_ref_c_not_in_cpg: counts[end_label]['G_base_sites'][dist_val] +=1 # Bio G
                        elif is_ref_a: counts[end_label]['T_base_sites'][dist_val] +=1 # Bio T
                        elif is_ref_t: counts[end_label]['A_base_sites'][dist_val] +=1 # Bio A
                    else: # Forward strand
                        if is_ref_c_in_cpg: counts[end_label]['C_in_CpG_sites'][dist_val] += 1
                        elif is_ref_c_not_in_cpg: counts[end_label]['C_not_in_CpG_sites'][dist_val] += 1
                        elif is_ref_g: counts[end_label]['G_base_sites'][dist_val] += 1
                        elif is_ref_a: counts[end_label]['A_base_sites'][dist_val] += 1
                        elif is_ref_t: counts[end_label]['T_base_sites'][dist_val] += 1

                    # --- Mismatch Counting ---
                    if read_base_query != ref_base_original and ref_base_original != 'N' and read_base_query != 'N':
                        mismatch_category_key = 'other' # Default
                        if read_obj.is_reverse:
                            # G->A on ref for reverse read == C->T on original (-) strand
                            if ref_base_original == 'G' and read_base_query == 'A':
                                mismatch_category_key = 'C>T_CpG' if ref_pos_to_context.get(ref_pos) == 'G_in_CpG' else 'C>T_nonCpG'
                        else: # Forward strand
                            # C->T on ref for forward read == C->T on original (+) strand
                            if ref_base_original == 'C' and read_base_query == 'T':
                                mismatch_category_key = 'C>T_CpG' if current_ref_context == 'CpG' else 'C>T_nonCpG'

                        counts[end_label][f'{mismatch_category_key}_mismatches'][dist_val] += 1

    # --- Calculate Frequencies ---
    output_frequencies = {}
    # Define denominators for each mutation type
    # For C>T_CpG, denominator is C_in_CpG_sites
    # For C>T_nonCpG, denominator is C_not_in_CpG_sites
    # For 'other' mismatches, denominator is 'total_considered_sites' minus C_in_CpG and C_not_in_CpG sites
    # (or more simply, total_considered_sites if we assume 'other' can happen anywhere)

    if sequencing_type == 'single':
        output_frequencies = {'total_considered_sites': np.zeros(max_dist, dtype=int)}
        for mut_cat in mut_categories:
            output_frequencies[f'{mut_cat}_mismatches'] = np.zeros(max_dist, dtype=int)
            output_frequencies[mut_cat] = np.zeros(max_dist, dtype=float) # For the frequency
        for site_cat in site_categories:
            output_frequencies[f'{site_cat}_sites'] = np.zeros(max_dist, dtype=int)

        for i in range(max_dist):
            output_frequencies['total_considered_sites'][i] = counts['5_prime']['total_considered_sites'][i] + counts['3_prime']['total_considered_sites'][i]
            for site_cat in site_categories:
                output_frequencies[f'{site_cat}_sites'][i] = counts['5_prime'][f'{site_cat}_sites'][i] + counts['3_prime'][f'{site_cat}_sites'][i]

            for mut_cat in mut_categories:
                total_mismatches = counts['5_prime'][f'{mut_cat}_mismatches'][i] + counts['3_prime'][f'{mut_cat}_mismatches'][i]
                output_frequencies[f'{mut_cat}_mismatches'][i] = total_mismatches
                
                # Determine denominator
                if mut_cat == 'C>T_CpG':
                    num_sites = output_frequencies['C_in_CpG_sites'][i]
                elif mut_cat == 'C>T_nonCpG':
                    num_sites = output_frequencies['C_not_in_CpG_sites'][i]
                else: # 'other'
                    # Other mismatches can occur at G, A, or T sites, or C sites that didn't result in C>T
                    num_sites = output_frequencies['G_base_sites'][i] + \
                                output_frequencies['A_base_sites'][i] + \
                                output_frequencies['T_base_sites'][i] + \
                                (output_frequencies['C_in_CpG_sites'][i] - (counts['5_prime']['C>T_CpG_mismatches'][i] + counts['3_prime']['C>T_CpG_mismatches'][i])) + \
                                (output_frequencies['C_not_in_CpG_sites'][i] - (counts['5_prime']['C>T_nonCpG_mismatches'][i] + counts['3_prime']['C>T_nonCpG_mismatches'][i]))
                    num_sites = max(0, num_sites) # Must be non-negative; could also use total_considered_sites

                output_frequencies[mut_cat][i] = (total_mismatches / num_sites) if num_sites > 0 else 0.0
    else: # Double-stranded
        output_frequencies = {end_label: {} for end_label in ends_map}
        for end_label in ends_map:
            # Align keys with plotting function's expectations
            output_frequencies[end_label]['total_bases'] = np.zeros(max_dist, dtype=int) # For "other" site counts
            for mut_cat in mut_categories:
                output_frequencies[end_label][f'{mut_cat}_counts'] = np.zeros(max_dist, dtype=int) # CHANGED KEY
                output_frequencies[end_label][mut_cat] = np.zeros(max_dist, dtype=float) # This key is for the frequency array itself

            output_frequencies[end_label]['CpG_C_sites'] = np.zeros(max_dist, dtype=int)      # CHANGED KEY
            output_frequencies[end_label]['nonCpG_C_sites'] = np.zeros(max_dist, dtype=int) # CHANGED KEY
            # Store G, A, T sites for calculating 'other' denominator or if 'total_bases' is used by plotting
            output_frequencies[end_label]['G_sites'] = np.zeros(max_dist, dtype=int)
            output_frequencies[end_label]['A_sites'] = np.zeros(max_dist, dtype=int)
            output_frequencies[end_label]['T_sites'] = np.zeros(max_dist, dtype=int)


            for i in range(max_dist):
                output_frequencies[end_label]['total_bases'][i] = counts[end_label]['total_considered_sites'][i]
                output_frequencies[end_label]['CpG_C_sites'][i] = counts[end_label]['C_in_CpG_sites'][i]
                output_frequencies[end_label]['nonCpG_C_sites'][i] = counts[end_label]['C_not_in_CpG_sites'][i]
                output_frequencies[end_label]['G_sites'][i] = counts[end_label]['G_base_sites'][i]
                output_frequencies[end_label]['A_sites'][i] = counts[end_label]['A_base_sites'][i]
                output_frequencies[end_label]['T_sites'][i] = counts[end_label]['T_base_sites'][i]


                for mut_cat in mut_categories:
                    mismatches_at_pos = counts[end_label][f'{mut_cat}_mismatches'][i]
                    output_frequencies[end_label][f'{mut_cat}_counts'][i] = mismatches_at_pos # CHANGED KEY (for mismatch count)
                    
                    num_sites = 0
                    if mut_cat == 'C>T_CpG':
                        num_sites = output_frequencies[end_label]['CpG_C_sites'][i]
                    elif mut_cat == 'C>T_nonCpG':
                        num_sites = output_frequencies[end_label]['nonCpG_C_sites'][i]
                    else: # 'other'
                         # Denominator for 'other' used by plotting is 'total_bases'
                         num_sites = output_frequencies[end_label]['total_bases'][i] # Use this for consistency with plotting
                                     # This num_sites is for calculating the frequency.
                                     # The CI in plotting will use the stored site_counts arrays.

                    output_frequencies[end_label][mut_cat][i] = (mismatches_at_pos / num_sites) if num_sites > 0 else 0.0
        #logger.info(f"Mismatch Freq Calculation - Raw Counts (double-stranded): {counts}") # Optional: Keep for debugging
        #logger.info(f"Mismatch Freq Calculation - Final Output Frequencies (double-stranded): {output_frequencies}") # Optional: Keep for debugging
        
    return output_frequencies


# --- Helper Functions (expect pysam.AlignedSegment) ---
def get_deamination_pattern(read_obj):
    """Identifies if a read contains C>T or G>A damage (biological)."""
    if not isinstance(read_obj, pysam.AlignedSegment) or read_obj.is_unmapped:
        return []
    pattern = []
    mismatches = get_mismatches_from_read(read_obj)
    for mismatch in mismatches:
        ref = mismatch['ref_base']
        read_b = mismatch['read_base']
        if read_obj.is_reverse: # G->A on reference for a reverse read means C->T on original template
            if ref == 'G' and read_b == 'A': pattern.append('C>T')
        else: # C->T on reference for a forward read
            if ref == 'C' and read_b == 'T': pattern.append('C>T')
    return list(set(pattern))


def count_ct_changes(read_obj):
    """Counts biological C>T changes in a read, respecting strand."""
    count = 0
    if not isinstance(read_obj, pysam.AlignedSegment) or read_obj.is_unmapped:
        return 0
    mismatches = get_mismatches_from_read(read_obj)
    for mismatch in mismatches:
        ref = mismatch['ref_base']
        read_b = mismatch['read_base']
        if read_obj.is_reverse:
            if ref == 'G' and read_b == 'A': count += 1
        else:
            if ref == 'C' and read_b == 'T': count += 1
    return count