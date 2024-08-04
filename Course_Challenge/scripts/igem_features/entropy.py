# Features from Tamir's article
# https://www.nature.com/articles/s41598-021-89918-6#Sec8

import pandas as pd
import math

from Course_Challenge.utils.consts import NUCLEOTIDES

def get_prob_dict_for_idx(idx_bases):
    counts = {}
    for char in idx_bases:
        if char in counts:
            counts[char] += 1
        else:
            counts[char] = 1
    prob_d = {c: counts[c] / len(idx_bases) for c in counts}
    return prob_d


def get_prob_dict(seq_lst: list):
    d = {}
    seq_len = len(seq_lst[0])
    for i in range(seq_len):
        bases_per_idx = [seq[i] for seq in seq_lst]
        d[i] = get_prob_dict_for_idx(bases_per_idx)
    return d


def entropy(seq: 'pd.Series[str]') -> pd.DataFrame:
    # calc entropy for seq
    def seq_entropy(seq, prob_dict):
        entropy_val = 0
        for i, v in enumerate(seq):
            p_i = prob_dict[i][v]
            entropy_val += p_i * math.log2(p_i)
        return -entropy_val / 2  # divide by 2 for normalization purposes

    prob_dict = get_prob_dict(seq.to_list())
    # df['entropy'] = df[seq_col].apply(seq_entropy, prob_dict=prob_dict)
    # return df
    return pd.DataFrame(seq.apply(seq_entropy, prob_dict=prob_dict).tolist(), columns=['entropy'])


# def clean_sequence(sequence):
#     return ''.join([nt for nt in sequence if nt in NUCLEOTIDES])
#
#
# excel_file_path = r"C:\Users\Dvir\Desktop\Limudim\computational genomics\Course_Challenge\data\Train_data.xlsx"
# # Load the sequence data
# features_df = pd.read_excel(excel_file_path, sheet_name='Variants data', engine='openpyxl')
# features_df = features_df.iloc[:30, :]
# # Clean the sequences
# features_df['Variant sequence'] = features_df['Variant sequence'].apply(clean_sequence)
#
# # Pass the entire column to the extract_nucli_features function
# nucli_features_df = entropy(features_df['Variant sequence'])
#
# # Join the new features with the original DataFrame if needed
# features_df = pd.concat([features_df, nucli_features_df], axis=1)
#
# # nucli_features_df = features_df['Variant sequence'].apply(extract_nucli_features)
# print(features_df.head(5))