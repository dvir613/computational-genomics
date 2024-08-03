'''
The script should generate the features for all the sequences using the functions from features script
'''

import pandas as pd
from sklearn.preprocessing import LabelEncoder
from features import *  # imports everything (functions and constants)

save_data = True

def preprocess_data(excel_file_path):
    # Load the sequence data
    sequence_df = pd.read_excel(excel_file_path, sheet_name='Variants data', engine='openpyxl')
    # Clean the sequences
    sequence_df['Variant sequence'] = sequence_df['Variant sequence'].apply(clean_sequence)

    sequence_df['GC_Content'] = sequence_df['Variant sequence'].apply(calculate_gc_content)

    # Process each PSSM matrix
    for sheet_name in PSSM_matrices.sheet_names:
        pssm = pd.read_excel(PSSM_matrices, sheet_name=sheet_name, index_col=0).to_numpy()
        pssm_scores = sequence_df['Variant sequence'].apply(lambda seq: score_pssm_match(pssm, seq))
        # Add PSSM scores to the DataFrame
        pssm_feature_df = pd.DataFrame(pssm_scores.tolist(), columns=[f'PSSM_{sheet_name}_{idx}' for idx in
                                                                      range(pssm_scores.iloc[0].size)])
        sequence_df = pd.concat([sequence_df, pssm_feature_df], axis=1)

    # Calculate anti-SD hybridization energies and concatenate them to the DataFrame
    anti_sd_features = sequence_df['Variant sequence'].apply(anti_SD_hybridization_energy)
    # Convert anti-SD features into a sparse DataFrame
    anti_sd_df = pd.DataFrame(anti_sd_features.tolist(), columns=[f'anti_SD_hybridization_energy_{i}' for i in
                                                                  range(anti_sd_features.iloc[0].size)]).fillna(0)
    # Concatenate the anti-SD features to the original DataFrame
    sequence_df = pd.concat([sequence_df, anti_sd_df], axis=1)

    return sequence_df


def save_dataframe_in_chunks(df, filename, chunk_size=100000):
    chunks = [df[i:i+chunk_size] for i in range(0, df.shape[0], chunk_size)]
    for i, chunk in enumerate(chunks):
        mode = 'w' if i == 0 else 'a'
        header = (i == 0)
        chunk.to_csv(filename, mode=mode, index=False, header=header)


data_file_path = "../data/Train_data.xlsx"
df = preprocess_data(data_file_path)

if save_data:
    # Save the DataFrame in chunks to manage memory usage better
    save_dataframe_in_chunks(df, "../data/data_with_features.csv")
