# import pandas as pd
#
#
# NUCLEOTIDES = ['A', 'C', 'G', 'T']
#
#
# # Function to clean sequences (mainly from the quotation marks)
# def clean_sequence(sequence):
#     return ''.join([nt for nt in sequence if nt in NUCLEOTIDES])
#
#
# # load necessary files:
# PSSM_matrices = pd.ExcelFile("PSSM.xlsx")
#
# anti_SD_hybridization_energy = pd.read_excel("asd_hyb.xlsx")  # anti-Shine-Dalgarno hybridization energy
# anti_SD_hybridization_energy.columns = ['Sequence', 'Energy']
# anti_SD_hybridization_energy['Sequence'] = anti_SD_hybridization_energy['Sequence'].apply(clean_sequence)
# # Create a dictionary from the DataFrame to map string to value
# mapping_anti_SD_hybridization_energy = pd.Series(anti_SD_hybridization_energy['Energy'].values,
#                                                  index=anti_SD_hybridization_energy['Sequence']).to_dict()



import pandas as pd
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


