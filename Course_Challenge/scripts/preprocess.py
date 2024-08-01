'''
The script should generate the features for all the sequences using the functions from features script
'''

import pandas as pd
from sklearn.preprocessing import LabelEncoder


def preprocess_data(file_path):
    df = pd.read_csv(file_path)

    # Example preprocessing: Convert nucleotides to numeric
    df['nucleotides'] = df['nucleotides'].apply(lambda x: [int(i) for i in x])

    # Split features and target
    X = df.drop(columns=['target'])
    y = df['target']

    return X, y
