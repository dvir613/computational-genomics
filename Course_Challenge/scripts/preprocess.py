'''
The script should generate the features for all the sequences using the functions from features script
'''

import pandas as pd
from sklearn.preprocessing import LabelEncoder
from features import *  # imports everything (functions and constants)
from typing import Optional, Tuple, Union, List

save_data = True


def generate_features(excel_file_path):
    # Load the sequence data
    features_df = pd.read_excel(excel_file_path, sheet_name='Variants data', engine='openpyxl')
    # Clean the sequences
    features_df['Variant sequence'] = features_df['Variant sequence'].apply(clean_sequence)

    features_df['GC_Content'] = features_df['Variant sequence'].apply(calculate_gc_content)

    # Process each PSSM matrix
    for sheet_name in PSSM_matrices.sheet_names:
        pssm = pd.read_excel(PSSM_matrices, sheet_name=sheet_name, index_col=0).to_numpy()
        pssm_scores = features_df['Variant sequence'].apply(lambda seq: score_pssm_match(pssm, seq))
        # Add PSSM scores to the DataFrame
        pssm_feature_df = pd.DataFrame(pssm_scores.tolist(), columns=[f'PSSM_{sheet_name}_{idx}' for idx in
                                                                      range(pssm_scores.iloc[0].size)])
        features_df = pd.concat([features_df, pssm_feature_df], axis=1)

    # Calculate anti-SD hybridization energies and concatenate them to the DataFrame
    anti_sd_features = features_df['Variant sequence'].apply(anti_SD_hybridization_energy)
    # Convert anti-SD features into a sparse DataFrame
    anti_sd_df = pd.DataFrame(anti_sd_features.tolist(), columns=[f'anti_SD_hybridization_energy_{i}' for i in
                                                                  range(anti_sd_features.iloc[0].size)]).fillna(0)
    # Concatenate the anti-SD features to the original DataFrame
    features_df = pd.concat([features_df, anti_sd_df], axis=1)

    # drop the sequences column
    features_df.drop(columns='Variant sequence', inplace=True)

    return features_df


def remove_zero_variance_features(X: pd.DataFrame):
    zero_variance_cols = X.columns[X.var() == 0]
    return X.drop(columns=zero_variance_cols), zero_variance_cols


def preprocess_data(train_excel_file_path, test_excel_file_path):
    X_features_df = generate_features(train_excel_file_path)
    print('number of features before zero variance remove: ', X_features_df.shape[1])
    X_features_df, zero_variance_features = remove_zero_variance_features(X_features_df)
    print('number of features after zero variance remove: ', X_features_df.shape[1])
    Y_features_df = generate_features(test_excel_file_path)
    # drop from the test data features the same features dropped from the train data
    Y_features_df.drop(columns=zero_variance_features, inplace=True)
    return X_features_df, Y_features_df


def save_dataframe_in_chunks(df, filename, chunk_size=100000):
    chunks = [df[i:i+chunk_size] for i in range(0, df.shape[0], chunk_size)]
    for i, chunk in enumerate(chunks):
        mode = 'w' if i == 0 else 'a'
        header = (i == 0)
        chunk.to_csv(filename, mode=mode, index=False, header=header)


train_data_file_path = "../data/Train_data.xlsx"
test_data_file_path = "../data/Test_data.xlsx"
X_train_processed_features_df, X_test_processed_features_df = preprocess_data(train_data_file_path, test_data_file_path)

if save_data:
    # Save the DataFrame in chunks to manage memory usage better
    save_dataframe_in_chunks(X_train_processed_features_df, "../data/Train_data_with_features.csv")
    print("Train data saved in ../data/Train_data_with_features.csv")
    save_dataframe_in_chunks(X_test_processed_features_df, "../data/Test_data_with_features.csv")
    print("Test data saved in ../data/Test_data_with_features")
