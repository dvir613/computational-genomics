# from assignment3.genomical_statistics import calc_pssm
import numpy as np
from seqfold import dg
from tqdm import tqdm

from Course_Challenge.utils.consts import NUCLEOTIDES, load_anti_sd_hybridization_energy
mapping_anti_SD_hybridization_energy = load_anti_sd_hybridization_energy()

def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)


def calculate_folding_energy_window(sequence):
    window_length = 40
    # create a numpy array to record the indices scores
    energy_values = np.zeros(len(sequence) - window_length + 1)
    # energy_values = []
    if len(sequence) < window_length:
        return energy_values  # Return empty list if sequence is too short

    # Calculate energy for each window in the sequence using Seqfold
    for idx in tqdm(range(len(sequence) - window_length + 1)):
        subsequence = sequence[idx:idx + window_length]
        result = dg(subsequence, temp=37.0)  # Use Seqfold's fold function
        energy_values[idx] = result  # Append the Gibbs free energy

    return energy_values


def calculate_folding_energy(sequence):
    result = dg(sequence, temp=37.0)  # Use Seqfold's fold function for the entire sequence
    return result  # Return the Gibbs free energy


def energy_difference(variant_seq, control_seq):
    # Calculate folding energy for both variant and control sequences
    variant_energy = calculate_folding_energy(variant_seq)
    control_energy = calculate_folding_energy(control_seq)

    # Calculate the energy difference
    difference = variant_energy - control_energy
    return difference


def calculate_pssm(sequence):
    # Placeholder for PSSM calculation
    pass


# function that calculates the score of matching to the PSSMs (from PSSM.xlsx) for each index in the sequence
# need to run it once for each PSSM from the excel file
def score_pssm_match(pssm, sequence):
    window_length = pssm.shape[1]
    # create a dictionary to record the indices scores
    idx_pssm_score = np.zeros(len(sequence) - window_length + 1)
    for idx in range(len(sequence)-window_length+1):
        subsequence = sequence[idx:idx + window_length]
        # print(subsequence)
        score = sum(pssm[NUCLEOTIDES.index(nt), position] for position, nt in enumerate(subsequence))
        idx_pssm_score[idx] = score
    return idx_pssm_score


def anti_SD_hybridization_energy(sequence):
    window_length = 6  # the length of SD sequence
    # create a numpy array to record the indices scores
    idx_energy = np.zeros(len(sequence) - window_length + 1)
    for idx in range(len(sequence) - window_length + 1):
        subsequence = sequence[idx:idx + window_length]
        # assign the value according to the anti SD hybridization energies
        idx_energy[idx] = mapping_anti_SD_hybridization_energy.get(subsequence, 0)
    return idx_energy


# this data is provided and there is no need to calculate it
# # function that calculate the number of codons that have changed in the variant regarding the control
# def calculate_codon_changes(control_sequence, variant_sequence):
#     control_codons = [control_sequence[i:i+3] for i in range(0, len(control_sequence), 3)]
#     variant_codons = [variant_sequence[i:i+3] for i in range(0, len(variant_sequence), 3)]
#     codon_changes = sum(1 for c_codon, v_codon in zip(control_codons, variant_codons) if c_codon != v_codon)
#     return codon_changes


def calculate_CAI(sequence):
    pass


def calculate_tAI(sequence):
    pass


def calculate_custom_feature(sequence):
    # Placeholder for another custom feature
    pass
