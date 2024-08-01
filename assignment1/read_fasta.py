from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
import pandas as pd


# function to read fasta file and convert it to pandas dataframe
def read_fasta(file):
    header = []
    seqs = []
    for record in SeqIO.parse(file, "fasta"):
        header.append(record.description)
        seqs.append(record.seq)

    df_data = {"header": header, "sequence": seqs}
    df = pd.DataFrame(df_data)
    return df


# function to add a column of the corresponding RNA sequences
def convert_to_rna(df):
    df["rna_sequence"] = df["sequence"].apply(transcribe_dna_to_rna)
    return df


# helper function for transcribing DNA sequence of type Seq to RNA sequnce of type Seq
def transcribe_dna_to_rna(dna_seq):
    rna_seq = dna_seq.transcribe()
    return rna_seq


# function to add a column of the corresponding amino acids
def convert_to_aa(df):
    df["aa_sequence"] = df["rna_sequence"].apply(translate_to_amino_acids)
    return df


# helper function to translate RNA sequence to amino acids
def translate_to_amino_acids(rna_seq):
    return rna_seq.translate()


# function to count the number of codons in the sequnces
def calc_seq_len(df):
    df["sequence_len"] = df["sequence"].apply(lambda x: int(len(x)/3))
    return df


# function ton calculate the GC content of a sequence
def GC_content_calc(seq):
    mutable_seq = MutableSeq(str(seq))
    g_content = mutable_seq.count("G")
    c_content = mutable_seq.count("C")
    return 100 * (g_content + c_content) / len(mutable_seq)  # return the percentage of GC in the sequence


# function to add a column of the GC content
def GC_content_total(df):
    df["gc_content"] = df["sequence"].apply(GC_content_calc)
    return df


# check
# df = read_fasta("yeast_genes.fa")
# print(df.shape)
# print(df.head(3))
# df_w_rna = convert_to_rna(df)
# print(df_w_rna.iloc[0])
# df_w_aa = convert_to_aa(df)
# print(df_w_aa.iloc[0])
# df_w_len = calc_seq_len(df)
# print(df_w_len.iloc[0])
# print(df.head(3))
# s = Seq("ATCGGGTACG")
# gc = GC_content_calc(s)
# print(gc)
# df_gc_content = GC_content_total(df_w_len)
# print(df_gc_content.iloc[0])


