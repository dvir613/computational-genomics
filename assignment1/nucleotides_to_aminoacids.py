
# a dictionary between nucleotides codons to the amino acids
nt_2_aa_dict = {"aa": "X", "ab": "Y", "ba": "Z", "bb": "W"}


# function that receives nucleotides sequence as a string and returns amino acids as string
def nt_2_aa(nt_vec):
    res = ""
    # each codon is 2 nucleotides. first check if the input sequence has a valid length
    if len(nt_vec) % 2 != 0:
        res = "Invalid length"
    else:
        nt_vec = nt_vec.lower()
        for i in range(int(len(nt_vec)/2)):
            res += nt_2_aa_dict.get(nt_vec[2*i:2*i+2])

    return res

# to check uncomment these examples
# nt_seq1 = "aaBBabAbbb"
# aa_seq1 = nt_2_aa(nt_seq1)
# print(aa_seq1)
#
# nt_seq2 = "aaBBabbbaba"
# aa_seq2 = nt_2_aa(nt_seq2)
# print(aa_seq2)

