from ...assignment3.genomical_statistics import calc_pssm


def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)


def calculate_pssm(sequence):
    # Placeholder for PSSM calculation
    pass


def calculate_custom_feature(sequence):
    # Placeholder for another custom feature
    pass
