# Lower-level processing helpers


def calculate_cg_content(sequence):
    """Calculate CG content of a sequence string."""
    if not sequence or len(sequence) == 0:
        return 0.0
    cg_count = sequence.upper().count('C') + sequence.upper().count('G')
    return cg_count / len(sequence)