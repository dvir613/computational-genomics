from datetime import datetime
from itertools import product
import pandas as pd
from typing import List, Tuple
from tqdm import tqdm

from Course_Challenge.utils.consts import NUCLEOTIDES

m2 = list(product(NUCLEOTIDES, repeat=2))
m3 = list(product(NUCLEOTIDES, repeat=3))
m4 = list(product(NUCLEOTIDES, repeat=4))
m5 = list(product(NUCLEOTIDES, repeat=5))
k_gap = 2
k_tuple = 2


def kmers(seq: str, k: int) -> List[str]:
    v = []
    for i in range(len(seq) - k + 1):
        v.append(seq[i:i + k])
    return v


def z_curve(sequences: 'pd.Series[str]') -> pd.DataFrame:
    T = sequences.str.count('T')
    A = sequences.str.count('A')
    C = sequences.str.count('C')
    G = sequences.str.count('G')

    x_ = (A + G) - (C + T)
    y_ = (A + C) - (G + T)
    z_ = (A + T) - (C + G)

    return pd.DataFrame({'z_curve_x': x_, 'z_curve_y': y_, 'z_curve_z': z_})

def cumulative_skew(sequences: 'pd.Series[str]') -> pd.DataFrame:
    T = sequences.str.count('T')
    A = sequences.str.count('A')
    C = sequences.str.count('C')
    G = sequences.str.count('G')

    ATSkew = (A - T) / (A + T)
    GCSkew = (G - C) / (G + C)

    return pd.DataFrame({'at_skew': ATSkew, 'gc_skew': GCSkew})


def get_k_gap_description(nucleotides: Tuple[str], before_gap: int, after_gap: int, k: int, gap: str='_') -> str:
    return f'{"".join(nucleotides[:before_gap])}{k*gap}{"".join(nucleotides[before_gap:before_gap+after_gap])}_count'


def mono_mono_k_gap(_kmers: 'pd.Series[List[str]]', g: int) -> pd.DataFrame:
    def count_matches(V, _gGap):
        _count = 0
        for v in V:
            if v[0] == _gGap[0] and v[-1] == _gGap[1]:
                _count += 1
        return _count

    d = {}
    m = m2
    for i in tqdm(range(1, g + 1)):
        V = _kmers[i + 2]
        for gGap in m:
            key = get_k_gap_description(gGap, 1, 1, i)
            d[key] = V.apply(lambda v: count_matches(v, gGap) / len(v))

    return pd.DataFrame(d)


def mono_di_k_gap(_kmers: 'pd.Series[List[str]]', g: int) -> pd.DataFrame:
    def count_matches(V, _gGap):
        _count = 0
        for v in V:
            if v[0] == _gGap[0] and v[-2] == _gGap[1] and v[-1] == _gGap[2]:
                _count += 1
        return _count

    d = {}
    m = m3
    for i in tqdm(range(1, g + 1)):
        V = _kmers[i + 3]
        for gGap in m:
            key = get_k_gap_description(gGap, 1, 2, i)
            d[key] = V.apply(lambda v: count_matches(v, gGap) / len(v))

    return pd.DataFrame(d)


def di_mono_k_gap(_kmers: 'pd.Series[List[str]]', g: int) -> pd.DataFrame:
    def count_matches(V, _gGap):
        _count = 0
        for v in V:
            if v[0] == _gGap[0] and v[1] == _gGap[1] and v[-1] == _gGap[2]:
                _count += 1
        return _count

    d = {}
    m = m3
    for i in tqdm(range(1, g + 1)):
        V = _kmers[i + 3]
        for gGap in m:
            key = get_k_gap_description(gGap, 2, 1, i)
            d[key] = V.apply(lambda v: count_matches(v, gGap) / len(v))

    return pd.DataFrame(d)


def mono_tri_k_gap(_kmers: 'pd.Series[List[str]]', g: int) -> pd.DataFrame:
    def count_matches(V, _gGap):
        _count = 0
        for v in V:
            if v[0] == _gGap[0] and v[-3] == _gGap[1] and v[-2] == _gGap[2] and v[-1] == _gGap[3]:
                _count += 1
        return _count

    d = {}
    m = m4
    for i in tqdm(range(1, g + 1)):
        V = _kmers[i + 4]
        for gGap in m:
            key = get_k_gap_description(gGap, 1, 3, i)
            d[key] = V.apply(lambda v: count_matches(v, gGap) / len(v))

    return pd.DataFrame(d)


def tri_mono_k_gap(_kmers: 'pd.Series[List[str]]', g: int) -> pd.DataFrame:
    def count_matches(V, _gGap):
        _count = 0
        for v in V:
            if v[0] == _gGap[0] and v[1] == _gGap[1] and v[2] == _gGap[2] and v[-1] == _gGap[3]:
                _count += 1
        return _count

    d = {}
    m = m4
    for i in tqdm(range(1, g + 1)):
        V = _kmers[i + 4]
        for gGap in m:
            key = get_k_gap_description(gGap, 3, 1, i)
            d[key] = V.apply(lambda v: count_matches(v, gGap) / len(v))

    return pd.DataFrame(d)


def di_di_k_gap(_kmers: 'pd.Series[List[str]]', g: int) -> pd.DataFrame:
    def count_matches(V, _gGap):
        _count = 0
        for v in V:
            if v[0] == _gGap[0] and v[1] == _gGap[1] and v[-2] == _gGap[2] and v[-1] == _gGap[3]:
                _count += 1
        return _count

    d = {}
    m = m4
    for i in tqdm(range(1, g + 1)):
        V = _kmers[i + 4]
        for gGap in m:
            key = get_k_gap_description(gGap, 2, 2, i)
            d[key] = V.apply(lambda v: count_matches(v, gGap) / len(v))

    return pd.DataFrame(d)


def di_tri_k_gap(_kmers: 'pd.Series[List[str]]', g: int) -> pd.DataFrame:
    def count_matches(V, _gGap):
        _count = 0
        for v in V:
            if v[0] == _gGap[0] and v[1] == _gGap[1] and v[-3] == _gGap[2] and v[-2] == _gGap[3] and v[-1] == _gGap[4]:
                _count += 1
        return _count

    d = {}
    m = m5
    for i in tqdm(range(1, g + 1)):
        V = _kmers[i + 5]
        for gGap in m:
            key = get_k_gap_description(gGap, 2, 3, i)
            d[key] = V.apply(lambda v: count_matches(v, gGap) / len(v))

    return pd.DataFrame(d)


def tri_di_k_gap(_kmers: 'pd.Series[List[str]]', g: int) -> pd.DataFrame:
    def count_matches(V, _gGap):
        _count = 0
        for v in V:
            if v[0] == _gGap[0] and v[1] == _gGap[1] and v[2] == _gGap[2] and v[-2] == _gGap[3] and v[-1] == _gGap[4]:
                _count += 1
        return _count

    d = {}
    m = m5
    for i in tqdm(range(1, g + 1)):
        V = _kmers[i + 5]
        for gGap in m:
            key = get_k_gap_description(gGap, 3, 2, i)
            d[key] = V.apply(lambda v: count_matches(v, gGap) / len(v))

    return pd.DataFrame(d)


def extract_nucli_features(sequences: 'pd.Series[str]') -> pd.DataFrame:
    d = []
    print("Generating KMERS")
    KMERS = [sequences.apply(lambda sequence: kmers(sequence, i)) for i in tqdm(range(5 + k_gap + 1))]

    print(f'start z_curve, time: {datetime.now()}')
    res = z_curve(sequences)
    d.append(res)

    print(f'start cumulative_skew, time: {datetime.now()}')
    res = cumulative_skew(sequences)
    d.append(res)

    print(f'start mono_mono_k_gap, time: {datetime.now()}')
    res = mono_mono_k_gap(KMERS, k_gap)  # 4*(k)*4 = 32
    d.append(res)

    print(f'start mono_di_k_gap, time: {datetime.now()}')
    res = mono_di_k_gap(KMERS, k_gap)  # 4*k*(4^2) = 128
    d.append(res)

    print(f'start mono_tri_k_gap, time: {datetime.now()}')
    res = mono_tri_k_gap(KMERS, k_gap)  # 4*k*(4^3) = 512
    d.append(res)

    print(f'start di_mono_k_gap, time: {datetime.now()}')
    res = di_mono_k_gap(KMERS, k_gap)  # (4^2)*k*(4)    = 128
    d.append(res)

    print(f'start di_di_k_gap, time: {datetime.now()}')
    res = di_di_k_gap(KMERS, k_gap)  # (4^2)*k*(4^2)  = 512
    d.append(res)

    print(f'start di_tri_k_gap, time: {datetime.now()}')
    res = di_tri_k_gap(KMERS, k_gap)  # (4^2)*k*(4^3)  = 2048
    d.append(res)

    print(f'start tri_mono_k_gap, time: {datetime.now()}')
    res = tri_mono_k_gap(KMERS, k_gap)  # (4^3)*k*(4)    = 512
    d.append(res)

    print(f'start tri_di_k_gap, time: {datetime.now()}')
    res = tri_di_k_gap(KMERS, k_gap)  # (4^3)*k*(4^2)  = 2048
    d.append(res)

    return pd.concat(d, axis=1)  # in total with k=2 -> 5943


#
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
# nucli_features_df = extract_nucli_features(features_df['Variant sequence'])
#
# # Join the new features with the original DataFrame if needed
# features_df = pd.concat([features_df, nucli_features_df], axis=1)
#
# # nucli_features_df = features_df['Variant sequence'].apply(extract_nucli_features)
# print(features_df.head(5))