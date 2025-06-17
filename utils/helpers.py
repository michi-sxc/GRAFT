# Basic helper functions unrelated to core processing/plotting

def rgb_to_hex(rgb):
    """Convert an RGB tuple to a hex string, rounding the values to integers."""
    return '#{:02x}{:02x}{:02x}'.format(int(round(rgb[0])), int(round(rgb[1])), int(round(rgb[2])))

def complement_base(base):
    """Return the complementary base."""
    base = base.upper()
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return complement_map.get(base, 'N')

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    reversed_seq = reversed(seq.upper())
    return ''.join([complement_map.get(base, 'N') for base in reversed_seq])
