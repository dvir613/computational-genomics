import itertools
import random

import numpy as np


# helper function to generate random DNA sequences
def generate_random_sequence(length):
    return ''.join(random.choices('ATCG', k=length))


# function that receives DNA sequences as input and
# returns the nucleotides pairs that their appearance is statistically significant
def find_significant_nt_pairs(sequences, random_check=100, alpha=0.05):
    statistical_significant_nt_pairs = []
    possible_pairs = [''.join(p) for p in itertools.product('ATCG', repeat=2)]
    # generate dictionary to record the appearances of each pair to calculate pValue
    empirical_p_value = {pair: 0 for pair in possible_pairs}
    # calculate the frequencies of each pair in the given sequences
    frequencies = {pair: 0 for pair in possible_pairs}
    for sequence in sequences:
        for pair in possible_pairs:
            frequencies[pair] += sequence.count(pair)

    # calculate the pairs frequencies for random sequences of the same length
    for check in range(random_check):
        random_sequences = [generate_random_sequence(len(sequence)) for sequence in sequences]
        random_frequencies = {pair: 0 for pair in possible_pairs}
        for sequence in random_sequences:
            for pair in possible_pairs:
                random_frequencies[pair] += sequence.count(pair)
        for pair in possible_pairs:
            if random_frequencies[pair] >= frequencies[pair]:
                empirical_p_value[pair] += 1

    # return the pairs that their appearance is statistically significant (pValue < alpha)
    for pair in possible_pairs:
        if (empirical_p_value[pair] / random_check) <= alpha:
            statistical_significant_nt_pairs.append(pair)

    return statistical_significant_nt_pairs


# # example
# sequences = ['ATCGATATATGCATC', 'ATCCATATATAT', 'ATCCTTATCCTTATATATATAT']
# for i in range(30):
#     nt_sig = find_significant_nt_pairs(sequences)
#     print(nt_sig)


# function that calculates PSSM and returns it as a numpy array
# where the order of the rows is A C G T
def calc_pssm(sequences, windowStart, windowEnd):
    nucleotides = ['A', 'C', 'G', 'T']
    window_length = windowEnd - windowStart
    pssm_array = np.zeros((4, window_length))
    for sequence in sequences:
        for position in range(windowStart, windowEnd):
            nt_index = nucleotides.index(sequence[position])
            pssm_array[nt_index, position-windowStart] += 1
    pssm_array = pssm_array/len(sequences)
    return pssm_array


# # example
# sequence_list = ['ACTGACTG', 'ACTGGCTA', 'AGCTCTAA', 'ATTTGCG']
# pssm = calc_pssm(sequence_list, 0, 5)
# print(pssm)


# function that finds the index in the sequence that matches best the given (calculated in the previous function) PSSM
def find_best_match(pssm, sequence):
    nucleotides = ['A', 'C', 'G', 'T']
    window_length = pssm.shape[1]
    # create a dictionary to record the indices scores
    idx_score = {}
    for idx in range(len(sequence)-window_length+1):
        idx_score[idx] = 0
        subsequence = sequence[idx:idx+window_length]
        for position, nt in enumerate(subsequence):
            idx_score[idx] += pssm[nucleotides.index(nt), position]
    return int(max(idx_score, key=idx_score.get))


# # example
# ind = find_best_match(pssm, 'TCGGTCAACTTGTCATGGATT')
# print(ind)

