import pandas as pd
import numpy as np
import os

NUCLEOTIDES = ['A', 'C', 'G', 'T']


# Function to clean sequences (mainly from the quotation marks)
def clean_sequence(sequence):
    return ''.join([nt for nt in sequence if nt in NUCLEOTIDES])


# Helper function to get the absolute path of the data files
def get_data_file_path(filename):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, filename)


# Lazy loading mechanism for DataFrames
PSSM_matrices = None
mapping_anti_SD_hybridization_energy = None
rnap_energy_mat = None


def load_pssm_matrices():
    global PSSM_matrices
    if PSSM_matrices is None:
        PSSM_matrices = pd.ExcelFile(get_data_file_path("PSSM.xlsx"))
    return PSSM_matrices


def load_anti_sd_hybridization_energy():
    global mapping_anti_SD_hybridization_energy
    if mapping_anti_SD_hybridization_energy is None:
        anti_SD_hybridization_energy = pd.read_excel(get_data_file_path("asd_hyb.xlsx"))
        anti_SD_hybridization_energy.columns = ['Sequence', 'Energy']
        anti_SD_hybridization_energy['Sequence'] = anti_SD_hybridization_energy['Sequence'].apply(clean_sequence)
        mapping_anti_SD_hybridization_energy = pd.Series(
            anti_SD_hybridization_energy['Energy'].values,
            index=anti_SD_hybridization_energy['Sequence']
        ).to_dict()
    return mapping_anti_SD_hybridization_energy


def get_energy_matrix_for_rna_polymerase() -> pd.DataFrame:
    """
    from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002811
        data link: https://doi.org/10.1371/journal.pcbi.1002811.s003
        Energy matrix for RNAP in kT. Inferred from an experiment done in
        TK310 with no supplemental cAMP (and hence, no CRP present in the
        cells). The matrix covers base pairs [-41:-1] where 0 denotes the
        transcription start site. Each row corresponds to a given position;
        each column corresponds to a value for that base pair. The columns
        are ordered [A,C,G,T].
    """
    promotor_strength_data = np.loadtxt(get_data_file_path('Energy matrix for RNAP.txt'))
    return pd.DataFrame(promotor_strength_data, columns=list(NUCLEOTIDES))


# rnap_energy_mat = get_energy_matrix_for_rna_polymeras()