'''
------------------------------------------------------------------------------
Code related to the paper "Identifying sequential residue patterns in bitter
and umami peptides" by Arghya Dutta, Tristan Bereau, and Thomas A. Vilgis.
Please consider citing the paper if you find our code useful in your research.
The BibTeX citation can be found in the README of the repository.
------------------------------------------------------------------------------
'''

import random
import collections
import itertools
import pathlib
import numpy as np
import pandas as pd
import amino_acid_data

'''
------------------------------------------------------------------------------
Parameters
- We set max_seq_length equal to 420 (the LCM of 1..7) because we have library
peptide units upto max_lib_pep_length=7.
- We used 500 (set using the variable `runs`) train--test splits.
- For singlets (i.e. `max_lib_pep_length`=1), set `sample_size` equal to 1. 
See discussion before the function `chunkstring`.
- Note that `max_lib_pep_length`=7 with 500 train--test splits takes almost 4.5
hours to run on a single core of a MacBook Pro M1 (October, 2022).
------------------------------------------------------------------------------
'''
runs = 500
max_seq_length = 420
sample_size = 5
max_lib_pep_length = 5

# ------------------------------------------------------------------------------
# Import the whole data and the train-test splits
# ------------------------------------------------------------------------------

data = pd.read_csv('../data/all-peptides.csv')
# print(data)
cwd = pathlib.Path.cwd()
ts_split_folder = cwd / '../data/ts-splits/'

train_set_paths = sorted(
    [path for path in list(ts_split_folder.iterdir()) if 'train' in path.stem])
train_sets = [pd.read_csv(path) for path in train_set_paths]

test_set_paths = sorted(
    [path for path in list(ts_split_folder.iterdir()) if 'test' in path.stem])
test_sets = [pd.read_csv(path) for path in test_set_paths]


# ------------------------------------------------------------------------------
# Convert the amino acid sequences of the proteins to a string of
# hydrophobicity and charge as defined in the paper.
# ------------------------------------------------------------------------------

def conv(x):
    if x == 'PZ':
        return ('H')
    elif x == 'NZ':
        return 'P'
    elif x == 'NP':
        return '+'
    elif x == 'NN':
        return '-'
    else:
        print('trouble!')


for df in train_sets:
    df['hydChrg'] = df.sequence.apply(
        lambda x: [conv(amino_acid_data.get_CGed_hyd_chrg(i)) for i in x])
for df in test_sets:
    df['hydChrg'] = df.sequence.apply(
        lambda x: [conv(amino_acid_data.get_CGed_hyd_chrg(i)) for i in x])
data['hydChrg'] = data.sequence.apply(
    lambda x: [conv(amino_acid_data.get_CGed_hyd_chrg(i)) for i in x])

# ------------------------------------------------------------------------------
# Construct the peptide library
# ------------------------------------------------------------------------------

singlets = ['H', 'P', '+', '-']
doublets = list(itertools.product(singlets, singlets))
triplets = list(itertools.product(singlets, singlets, singlets))
quartets = list(itertools.product(singlets, singlets, singlets, singlets))
quintets = list(
    itertools.product(singlets, singlets, singlets, singlets, singlets))
sextets = list(
    itertools.product(singlets, singlets, singlets, singlets, singlets, singlets))
septets = list(
    itertools.product(singlets, singlets, singlets, singlets, singlets, singlets, singlets))

# ------------------------------------------------------------------------------
# Make the long peptides by repeating the units
# ------------------------------------------------------------------------------

single_sequences = {
    'sing_' + ''.join(i): (tuple([i] * max_seq_length))
    for i in singlets
}
double_sequences = {
    'doub_' + ''.join(i): i * (max_seq_length // 2)
    for i in doublets
}
triple_sequences = {
    'trip_' + ''.join(i): i * (max_seq_length // 3)
    for i in triplets
}
quartile_sequences = {
    'quart_' + ''.join(i): i * (max_seq_length // 4)
    for i in quartets
}
quintile_sequences = {
    'quint_' + ''.join(i): i * (max_seq_length // 5)
    for i in quintets
}
sextet_sequences = {
    'sext_' + ''.join(i): i * (max_seq_length // 6)
    for i in sextets
}
septet_sequences = {
    'sept_' + ''.join(i): i * (max_seq_length // 7)
    for i in septets
}

# ------------------------------------------------------------------------------
# Note that there are some identical peptides. Remove them.
# ------------------------------------------------------------------------------

all_lib_seqs = {
    1: single_sequences,
    2: double_sequences,
    3: triple_sequences,
    4: quartile_sequences,
    5: quintile_sequences,
    6: sextet_sequences,
    7: septet_sequences
}


def get_unique_lib_seqs(max_lib_pep_length):

    # start by adding singlets
    all_sequences = {**single_sequences}
    unique_sequences = {**single_sequences}

    # If the length is >1, then duplicate sequences will start appearing. We
    # will check the seqs and add only if it is unique. We will add everything
    # to all_sequences, just to check the length of it.

    for i in range(2, max_lib_pep_length+1):
        all_sequences.update(all_lib_seqs[i])
        sequences = all_lib_seqs[i]
        for item in sequences.items():
            if item[1] not in unique_sequences.values():
                unique_sequences[item[0]] = item[1]
            else:
                continue
    return all_sequences, unique_sequences


all_sequences, unique_sequences = get_unique_lib_seqs(max_lib_pep_length)
# print(len(all_sequences), len(unique_sequences))
# print(unique_sequences)

# ------------------------------------------------------------------------------
# Define the comparison metric
# ------------------------------------------------------------------------------


def compare_seq(seq1, seq2):
    """
    Reference: Schilling, C., Mack, T., Lickfett, S., Sieste, S., Ruggeri, F.
    S., Sneideris, T., Dutta, A., Bereau, T., Naraghi, R., Sinske, D., Knowles,
    T. P. J., Synatschke, C. V., Weil, T., Knöll, B., Sequence‐Optimized 
    Peptide Nanofibers as Growth Stimulators for Regeneration of Peripheral 
    Neurons. Adv. Funct. Mater. 2019, 29, 1809112. 
    """
    score = 0.0
    if (len(seq2) >= len(seq1)):
        for i, _ in enumerate(seq1):
            if seq1[i] == seq2[i]:
                score += 1.0
    else:
        print('The second peptide should be longer.')
    return score / len(seq1)

'''
------------------------------------------------------------------------------
- We want to calculate the relative occurrences of the AA types (+, -, H, P) at
each lattice position of the first n (=`sample_size`) best-predicted sequences
from the library.
- At each lattice position, we will calculate the fraction of times one of the
four types occur in those n best library sequences (best means the library
sequence have the highest overlap with the peptides with the desired target
property, bitter or umami).
- Note that if we only consider singlets and choose `sample_size=5`, then we'll
get each AA type with a 25% probability of finding in the 5 best sequences,
because there are only 4 library sequences available. For that, we need to
check only the best one, i.e. `sample_size=1`.
- There are 12 unique library sequences if we consider singlets and doublets.
So, we used `sample_size=5` when max_lib_pep_length>1.
------------------------------------------------------------------------------
'''

def chunkstring(string, length):
    return (string[0 + i:length + i] for i in range(0, len(string), length))


def get_AA_count(df, metric, sample_size):
    """
    Input: 
    - A dataframe with sequence names and values of a metric (e.g. bitterness)
    - Name of the metric
    - How many highest-ranking samples to consider while making the output
    Output: 
    - A dataframe with histogram of AA types at each sequence positions
    """
#     print(df)
    _sample_best_candidates = list(
        df.sort_values(by=metric, axis=0,
                       ascending=False).sequence)[:sample_size]

    sample_best_candidates = [
        list(chunkstring(x.split('_')[1], 1)) for x in _sample_best_candidates
    ]

    # print(sample_best_candidates)
    lat_0 = []
    lat_1 = []
    lat_2 = []
    lat_3 = []
    lat_4 = []
    lat_5 = []
    lat_6 = []

    counts = {}

    for i in sample_best_candidates:
        lat_0.append(i[0])
        counts[0] = dict(collections.Counter(lat_0))

    for i in sample_best_candidates:
        if len(i) > 1:
            lat_1.append(i[1])
            counts[1] = dict(collections.Counter(lat_1))

    for i in sample_best_candidates:
        if len(i) > 2:
            lat_2.append(i[2])
            counts[2] = dict(collections.Counter(lat_2))

    for i in sample_best_candidates:
        if len(i) > 3:
            lat_3.append(i[3])
            counts[3] = dict(collections.Counter(lat_3))

    for i in sample_best_candidates:
        if len(i) > 4:
            lat_4.append(i[4])
            counts[4] = dict(collections.Counter(lat_4))

    for i in sample_best_candidates:
        if len(i) > 5:
            lat_5.append(i[5])
            counts[5] = dict(collections.Counter(lat_5))

    for i in sample_best_candidates:
        if len(i) > 6:
            lat_6.append(i[6])
            counts[6] = dict(collections.Counter(lat_6))
#     print(counts)
    return counts


def get_AA_count_df(df, metric, sample_size):
    counts = get_AA_count(df, metric, sample_size)
    val = pd.DataFrame(counts).T
    val = val.apply(lambda x: x / x.sum(), axis=1)
    return val

# print(get_AA_count_df(temp1, 'bitterness', 5))


# ------------------------------------------------------------------------------
# Find the average sequence
# ------------------------------------------------------------------------------


def get_av_seq(runs, metric,
               train_sets,
               sample_size,
               max_lib_pep_length=max_lib_pep_length):

    # For error bars in av_ovlp bar plots, we need to calculate the mean
    # overlap and the standard error of the mean.

    # First make the keys for the dictionary. They should be like of the form
    # (lattice_position, type)
    position_type = []
    for i in range(max_lib_pep_length):
        position_type.append((i, 'H'))
        position_type.append((i, 'P'))
        position_type.append((i, '+'))
        position_type.append((i, '-'))

#     print(position_type)

#     Then make the dictionary
    store = {}
    for i in position_type:
        store[i] = []

    # Now find average patterns from the train sets

    for run in range(runs):
        seq_data = train_sets[
            run].loc[:, ['sequence', 'bitter', 'umami', 'hydChrg']]

        foo_dict = {}
        # print(seq_data.shape)
        for row in seq_data.itertuples():
            foo_dict[row.Index] = {}
            for key, value in unique_sequences.items():
                score = 0.0
                score = compare_seq(row.hydChrg, list(value))
                foo_dict[row.Index][key] = score

        full = pd.merge(seq_data,
                        pd.DataFrame.from_dict(foo_dict).T,
                        left_index=True,
                        right_index=True)

#         print(full)

        temp1 = (full.drop(['hydChrg', 'umami'],
                           axis=1).groupby('bitter').mean().T)

        temp1.reset_index(inplace=True)
        temp1.columns = ['sequence', 'av_ovlp_wth_umami', 'av_ovlp_wth_bitter']
        temp1['bitterness'] = temp1.av_ovlp_wth_bitter
        temp1['umaminess'] = temp1.av_ovlp_wth_umami

        val = get_AA_count_df(temp1, metric, sample_size)

#         All AA types may not be present as columns in val. Add them

        for typ in ['H', 'P', '+', '-']:
            if typ not in val.columns:
                val[typ] = np.NAN

#         Just make the NANs 0s, they are probabilities after all.
        val.fillna(0, inplace=True)

#         print(val)

        for i in position_type:
            store[i].append(val.loc[i])

    means = pd.DataFrame()
    sems = pd.DataFrame()
    for key in store.keys():
        means.at[key] = pd.Series(store[key]).mean()
        sems.at[key] = pd.Series(store[key]).sem()


#     print(means, sems)

    return means[singlets], sems[singlets]


# ------------------------------------------------------------------------------
# Now let's analyze the bitters!
# ------------------------------------------------------------------------------
av_bitter, sem_bitter = get_av_seq(
    runs=runs, metric='bitterness', train_sets=train_sets,
    sample_size=sample_size,
    max_lib_pep_length=max_lib_pep_length)

# print(av_bitter)
# print(sem_bitter)

# ------------------------------------------------------------------------------
# Save the data for plotting
# ------------------------------------------------------------------------------

av_bitter.to_csv(
    f'../data/for-statistics/av_bitter_max_lib_pep_length_{max_lib_pep_length}.csv',
    index=None)
sem_bitter.to_csv(
    f'../data/for-statistics/sem_bitter_max_lib_pep_length_{max_lib_pep_length}.csv',
    index=None)

# ------------------------------------------------------------------------------
# Extract the best predicted bitter sequence
# ------------------------------------------------------------------------------

best_predicted_bitter = []
for i in av_bitter.index:
    best_predicted_bitter.append(
        av_bitter.loc[i].sort_values(ascending=False).index[0])
# print(best_predicted_bitter)

# ------------------------------------------------------------------------------
# Set the baseline bitter peptide as HHHHH...
# ------------------------------------------------------------------------------

baseline_bitter_peptide = []
for i in range(max_seq_length):
    baseline_bitter_peptide.append('H')
# baseline_bitter_peptide

# ------------------------------------------------------------------------------
# Get test statistics
# ------------------------------------------------------------------------------
def get_test_statistics_bitter(test_sets, bitter_pred):

    bitter_pred = bitter_pred * 50
    perf = pd.DataFrame(columns=[
        'pred_ovlp_bitters_mean', 'pred_ovlp_bitters_sem',
        'all_H_ovlp_bitters_mean', 'all_H_ovlp_bitters_sem',
        'rand_ovlp_bitters_mean', 'rand_ovlp_bitters_sem'
    ])

    for run, test_data in zip(range(runs), test_sets):
        test_data = test_sets[run]
        bitter_sample = test_data.loc[
            test_data.bitter == 1, ['sequence', 'bitter', 'umami', 'hydChrg']]
        # print(bitter_sample.shape)
        # Make a random peptide
        rand_pep = []
        for i in range(len(bitter_pred)):
            rand_pep.append(random.choice(['H', 'P', '+', '-']))

        perf.at[run,
                'pred_ovlp_bitters_mean'] = bitter_sample['hydChrg'].apply(
                    lambda x: compare_seq(x, bitter_pred)).mean()
        perf.at[run,
                'all_H_ovlp_bitters_mean'] = bitter_sample['hydChrg'].apply(
                    lambda x: compare_seq(x, baseline_bitter_peptide)).mean()
        perf.at[run,
                'rand_ovlp_bitters_mean'] = bitter_sample['hydChrg'].apply(
                    lambda x: compare_seq(x, rand_pep)).mean()

        perf.at[run, 'pred_ovlp_bitters_sem'] = bitter_sample['hydChrg'].apply(
            lambda x: compare_seq(x, bitter_pred)).sem()
        perf.at[run,
                'all_H_ovlp_bitters_sem'] = bitter_sample['hydChrg'].apply(
                    lambda x: compare_seq(x, baseline_bitter_peptide)).sem()
        perf.at[run, 'rand_ovlp_bitters_sem'] = bitter_sample['hydChrg'].apply(
            lambda x: compare_seq(x, rand_pep)).sem()


        # print(perf)

    return perf

# ------------------------------------------------------------------------------
# Save the result with a label denoting the maximum library peptide length
# ------------------------------------------------------------------------------

if (len(best_predicted_bitter) == max_lib_pep_length):
    perf_bitter = get_test_statistics_bitter(test_sets, best_predicted_bitter)
else:
    perf_bitter = pd.DataFrame()
    print('Something wrong?!')
perf_bitter.to_csv(
    f'../data/for-statistics/perf_bitter_max_lib_pep_length_{max_lib_pep_length}.csv', index=None)

# ------------------------------------------------------------------------------
# Now let's analyze the Umamis!
# ------------------------------------------------------------------------------
av_umami, sem_umami = get_av_seq(
    runs=runs, metric='umaminess', train_sets=train_sets,
    sample_size=sample_size,
    max_lib_pep_length=max_lib_pep_length)

# print(av_umami)
# print(sem_umami)

# ------------------------------------------------------------------------------
# Save the data for plotting
# ------------------------------------------------------------------------------

av_umami.to_csv(
    f'../data/for-statistics/av_umami_max_lib_pep_length_{max_lib_pep_length}.csv',
    index=None)
sem_umami.to_csv(
    f'../data/for-statistics/sem_umami_max_lib_pep_length_{max_lib_pep_length}.csv',
    index=None)

# ------------------------------------------------------------------------------
# Extract the best predicted umami sequence
# ------------------------------------------------------------------------------

best_predicted_umami = []
for i in av_umami.index:
    best_predicted_umami.append(
        av_umami.loc[i].sort_values(ascending=False).index[0])
# print(best_predicted_umami)

# ------------------------------------------------------------------------------
# Set the baseline umami peptide as -----...
# ------------------------------------------------------------------------------

baseline_umami_peptide = []
for i in range(max_seq_length):
    baseline_umami_peptide.append('-')
# baseline_umami_peptide

# ------------------------------------------------------------------------------
# Get test statistics
# ------------------------------------------------------------------------------


def get_test_statistics_umami(test_sets, umami_pred):

    umami_pred = umami_pred * 50
    perf = pd.DataFrame(columns=[
        'pred_ovlp_umamis_mean', 'pred_ovlp_umamis_sem',
        'all_H_ovlp_umamis_mean', 'all_H_ovlp_umamis_sem',
        'rand_ovlp_umamis_mean', 'rand_ovlp_umamis_sem'
    ])

    # analyse on test data
    for run, test_data in zip(range(runs), test_sets):
        test_data = test_sets[run]
        umami_sample = test_data.loc[
            test_data.umami == 1, ['sequence', 'bitter', 'umami', 'hydChrg']]
        # print(umami_sample.shape)
        # Make a random peptide
        rand_pep = []
        for i in range(len(umami_pred)):
            rand_pep.append(random.choice(['H', 'P', '+', '-']))

        perf.at[run,
                'pred_ovlp_umamis_mean'] = umami_sample['hydChrg'].apply(
                    lambda x: compare_seq(x, umami_pred)).mean()
        perf.at[run,
                'all_neg_ovlp_umamis_mean'] = umami_sample['hydChrg'].apply(
                    lambda x: compare_seq(x, baseline_umami_peptide)).mean()
        perf.at[run,
                'rand_ovlp_umamis_mean'] = umami_sample['hydChrg'].apply(
                    lambda x: compare_seq(x, rand_pep)).mean()

        perf.at[run, 'pred_ovlp_umamis_sem'] = umami_sample['hydChrg'].apply(
            lambda x: compare_seq(x, umami_pred)).sem()
        perf.at[run,
                'all_neg_ovlp_umamis_sem'] = umami_sample['hydChrg'].apply(
                    lambda x: compare_seq(x, baseline_umami_peptide)).sem()
        perf.at[run, 'rand_ovlp_umamis_sem'] = umami_sample['hydChrg'].apply(
            lambda x: compare_seq(x, rand_pep)).sem()

        # print(perf)

    return perf

# ------------------------------------------------------------------------------
# Save the result with a label denoting the maximum library peptide length
# ------------------------------------------------------------------------------


if (len(best_predicted_umami) == max_lib_pep_length):
    perf_umami = get_test_statistics_umami(test_sets, best_predicted_umami)
else:
    perf_umami = pd.DataFrame()
    print('Something wrong?!')
perf_umami.to_csv(
    f'../data/for-statistics/perf_umami_max_lib_pep_length_{max_lib_pep_length}.csv', index=None)
