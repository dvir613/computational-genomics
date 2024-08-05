'''
The script should generate the features for all the sequences using the functions from features script
'''

from Course_Challenge.utils.consts import load_pssm_matrices, clean_sequence
from features import *  # imports everything (functions and constants)
from igem_features.nucli_features import *
from igem_features.entropy import *
from igem_features.promoter_strength import *
from igem_features.delta_G.TX_prediction import *

save_data = True

PSSM_matrices = load_pssm_matrices()


def generate_features(excel_file_path):
    # Load the sequence data
    variants_df = pd.read_excel(excel_file_path, sheet_name='Variants data', engine='openpyxl')
    # variants_df = variants_df.iloc[:30, :]
    features_df = pd.read_excel(excel_file_path, sheet_name='Features', engine='openpyxl')
    # features_df = features_df.iloc[:30, :]

    # Clean the sequences
    features_df['Variant sequence'] = variants_df['Variant sequence'].apply(clean_sequence)

    # calculate GC content
    features_df['GC_Content'] = features_df['Variant sequence'].apply(calculate_gc_content)

    # # folding energy window=40
    # folding_window_features = features_df['Variant sequence'].apply(calculate_folding_energy_window)
    # # Convert folding energy into a sparse DataFrame
    # folding_window_df = pd.DataFrame(folding_window_features.tolist(), columns=[f'folding_energy_window_40_{i}' for i in
    #                                                               range(folding_window_features.iloc[0].size)]).fillna(0)
    # # Concatenate the folding energy features to the original DataFrame
    # features_df = pd.concat([features_df, folding_window_df], axis=1)
    #
    # # folding energy full
    # features_df['folding_energy'] = features_df['Variant sequence'].apply(calculate_folding_energy)

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

    # add nucli features
    # Pass the entire column to the extract_nucli_features function
    nucli_features_df = extract_nucli_features(features_df['Variant sequence'])

    # Join the new features with the original DataFrame if needed
    features_df = pd.concat([features_df, nucli_features_df], axis=1)

    # # entropy
    # that's not a good feature because it depends on all the sequences that are passed and not calculated per sequence
    # # Pass the entire column to the extract_nucli_features function
    # entropy_feature_df = entropy(features_df['Variant sequence'])
    #
    # # Join the new features with the original DataFrame if needed
    # features_df = pd.concat([features_df, entropy_feature_df], axis=1)

    # promoter strength
    # Calculate promoter strength and concatenate them to the DataFrame
    promoter_energy_features = features_df['Variant sequence'].apply(sliding_window_promoter_strength)
    # Convert promoter strength features into a sparse DataFrame
    promoter_energy_df = pd.DataFrame(promoter_energy_features.tolist(), columns=[f'promoter_energy_{i}' for i in
                                                                  range(promoter_energy_features.iloc[0].size)]).fillna(0)
    # Concatenate the promoter strength features to the original DataFrame
    features_df = pd.concat([features_df, promoter_energy_df], axis=1)

    # delta G
    # Calculate delta G and concatenate them to the DataFrame
    delta_G_features = process_deltaG_sequences(features_df['Variant sequence'])

    # Concatenate the delta G features to the original DataFrame
    features_df = pd.concat([features_df, delta_G_features], axis=1)

    # drop the sequences column
    features_df.drop(columns='Variant sequence', inplace=True)

    return features_df


def remove_zero_variance_features(X: pd.DataFrame):
    zero_variance_cols = X.columns[X.var() == 0]
    return X.drop(columns=zero_variance_cols), zero_variance_cols


def preprocess_data(train_excel_file_path, test_excel_file_path):
    X_features_df = generate_features(train_excel_file_path)
    # print(X_features_df.head(5))
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
