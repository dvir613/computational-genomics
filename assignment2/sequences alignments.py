from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq


# function to calculate pairwise global alignment for both DNA and AA sequences according to specified scoring values
def global_alignment(seq1, seq2, match=+1, mismatch=-1, open_gap=-10, extend_gap=-0.5):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap

    alignments = aligner.align(seq1, seq2)

    return alignments.score


# function to calculate pairwise local alignment for both DNA and AA sequences according to specified scoring values
def local_alignment(seq1, seq2, match=+1, mismatch=-1, open_gap=-10, extend_gap=-0.5):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap

    alignments = aligner.align(seq1, seq2)

    return alignments.score


# function to calculate pairwise alignment score according to a specific BLOSUM matrix and defined gap penalties
def global_alignment_blossom(str1, str2, scoring_matrix='BLOSUM50', open_gap=-10, extend_gap=-0.5):
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load(scoring_matrix)
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap

    alignments = aligner.align(str1, str2)

    return alignments.score


# function that receives pairwise alignment and counts the number of perfect matches
def calc_identity_score(str1, str2):
    identity_count = 0
    # Count identical letters in the same position
    for a, b in zip(str1, str2):  # use zip to iterate over both strings simultaneously
        if a == b:
            identity_count += 1

    # Calculate percentage
    identity_percentage = (identity_count / len(str1)) * 100
    return round(identity_percentage, 2)  # return the result rounded to 2 decimal places


# uncomment the next lines to check the functions
# seqA = Seq("CCCGCCGCTGTACGTACGCTAGCTAAA")
# seqB = Seq("CCCGCGGCCGTACGCTCACTAGCTCAGTAA")
# seqC = Seq("FINTQQAHNRGIFGNGARVAVL")
# seqD = Seq("AHNRGIFGNGARVAVLD")
#
# global_score = global_alignment(seqA, seqB)
# print("Global alignment score = %.1f" % global_score)

# global_AA_score = global_alignment(seqC, seqD)
# print("Global alignment AA score = %.1f" % global_AA_score)
#
# local_score = local_alignment(seqA, seqB)
# print("Local alignment score = %.1f" % local_score)

# local_AA_score = local_alignment(seqC, seqD)
# print("Local alignment AA score = %.1f" % local_AA_score)
#
# global_blossom_AA = global_alignment_blossom(seqC, seqD)
# print("Global 'BLOSUM50' AA alignment score = %.1f" % global_blossom_AA)
#
# str1 = 'CCCGC-C--GCTGTACGTACGCTAGCTA'
# str2 = 'CCCGCACTAGCT---CAT-CG---GTAA'
# identity_score = calc_identity_score(str1, str2)
# print("Identity score = %.2f" % identity_score)

