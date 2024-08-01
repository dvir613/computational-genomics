from ...assignment3.genomical_statistics import calc_pssm
import pandas as pd

# load necessary files:
PSSM_matrices = pd.ExcelFile("PSSM.xlsx")

anti_SD_hybridization_energy = pd.read_excel("asd_hyb.xlsx")  # anti-Shine-Dalgarno hybridization energy
anti_SD_hybridization_energy.columns = ['Sequence', 'Energy']
# Create a dictionary from the DataFrame to map string to value
mapping_anti_SD_hybridization_energy = pd.Series(anti_SD_hybridization_energy['Energy'].values,
                                                 index=anti_SD_hybridization_energy['Sequence']).to_dict()


def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)


def calculate_pssm(sequence):
    # Placeholder for PSSM calculation
    pass


# function that calculates the score of matching to the PSSMs (from PSSM.xlsx) for each index in the sequence
# need to run it once for each PSSM from the excel file
def score_pssm_match(pssm, sequence):
    nucleotides = ['A', 'C', 'G', 'T']
    window_length = pssm.shape[1]
    # create a dictionary to record the indices scores
    idx_pssm_score = {}
    for idx in range(len(sequence)-window_length+1):
        idx_pssm_score[idx] = 0
        subsequence = sequence[idx:idx+window_length]
        for position, nt in enumerate(subsequence):
            idx_pssm_score[idx] += pssm[nucleotides.index(nt), position]
    return idx_pssm_score


def anti_SD_hybridization_energy(sequence):
    window_length = 6  # the length of SD sequence
    # create a dictionary to record the indices scores
    idx_energy = {}
    for idx in range(len(sequence) - window_length + 1):
        subsequence = sequence[idx:idx + window_length]
        # assign the value according to the anti SD hybridization energies
        idx_energy[idx] = mapping_anti_SD_hybridization_energy.get(subsequence, 0)
    return idx_energy


def calculate_custom_feature(sequence):
    # Placeholder for another custom feature
    pass
