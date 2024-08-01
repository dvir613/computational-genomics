from ...assignment3 import genomical_statistics


def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)


def calculate_pssm(sequence):
    # Placeholder for PSSM calculation
    pass


def calculate_custom_feature(sequence):
    # Placeholder for another custom feature
    pass
