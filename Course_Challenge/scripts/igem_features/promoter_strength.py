from functools import partial
import numpy as np
import pandas as pd
from pathlib import Path
from Course_Challenge.utils.consts import *
from typing import List, Tuple, Union, Optional

energy_matrix = get_energy_matrix_for_rna_polymerase()


def sliding_window_promoter_strength(sequence, window_length=41) -> pd.DataFrame:
    """
    Calculate the promoter zones' strength for subsequences of a given DNA/RNA sequence using a sliding window approach.

    Args:
    sequence (str): The nucleotide sequence for which to calculate the promoter strength.
    window_length (int): The length of the sliding window.
    file_path (str): Path to the file containing the RNA polymerase energy matrix.

    Returns:
    np.ndarray: An array of promoter strength values for each window in the sequence.
    """

    idx_energy = np.zeros(len(sequence) - window_length + 1)

    for idx in range(len(sequence) - window_length + 1):
        subsequence = sequence[idx:idx + window_length]
        strength = 0.0
        for i, base in enumerate(subsequence):
            position = i - (window_length - 1)  # Adjust position based on window length
            if position in energy_matrix.index:
                strength += energy_matrix.loc[position, base]
        idx_energy[idx] = strength

    return idx_energy


def calculate_promoter_strength(sequence: str, energy_matrix: pd.DataFrame, start: int, end: int) -> float:
    """
    Calculate the promoter strength using the RNAP energy matrix.

    Args:
    sequence (str): The nucleotide sequence.
    energy_matrix (pd.DataFrame): DataFrame containing the RNAP energy matrix.
    start (int): Start position of the sequence window.
    end (int): End position of the sequence window.

    Returns:
    float: Total energy for the given sequence window.
    """
    strength = 0.0
    for position in range(start, end + 1):
        base = sequence[position - start]  # Adjust position to match sequence index
        if position in energy_matrix.index:
            strength += energy_matrix.loc[position, base]
    return strength




