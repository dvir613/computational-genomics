from collections import Counter
import pandas as pd


# function that reads a csv file as a dataframe
def read_orf(file):
    df = pd.read_csv(file)
    return df


# function that counts the appearances of the codons along each sequence and returns list of all the counters
def count_codons(sequences):
    codons_counter_list = []
    for seq in sequences:
        # Create a list of codons
        codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
        # since we were asked to use the Counter method
        counter = Counter(codons)
        codons_counter_list.append(counter)

        # counter = {}
        # i = 0
        # while i < (len(seq)-2):
        #     codon = seq[i:i+3]
        #     if codon in counter.keys():
        #         i += 3
        #         continue
        #     else:
        #         counter[codon] = seq.count(codon, i, len(seq))
        #         i += 3
        # codons_counter_list.append(counter)

    return codons_counter_list


# find the highest frequency appearing codon in each row
def find_highest_frequency_codon(codons_frequency):
    max_freq_list = []
    for counter in codons_frequency:
        most_common_codon, max_count = counter.most_common(1)[0]  # returns the most common codon and its count
        max_freq_list.append(max_count)
    return max_freq_list



# check
# orf_df = read_orf('E_cooli_ORF.csv')
# print(orf_df.head(5))
#
# codons_frequency = count_codons(orf_df.orf)
# print(codons_frequency[0])
#
# max_freq = find_highest_frequency_codon(codons_frequency)
# print(max_freq[0])


