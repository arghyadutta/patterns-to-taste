from sklearn.model_selection import train_test_split
import pandas as pd

# ------------------------------------------------------------------------------
# Generate stratified train--test splits
# ------------------------------------------------------------------------------
# We generated the train--test splits using the following code. The 500 pairs
# we used are in the `data/ts-splits` folder. 
# NOTE: It will overwrite the ts-splits that we used for the paper. For this
# reason, I commented out the code; uncomment it before you run.

# data = pd.read_csv('../data/all-peptides.csv')
# runs = 3 # set to 500 in the paper
# test_size = 0.2
# for run in range(runs):
#     train_data, test_data = train_test_split(data,
#                                              test_size=test_size,
#                                              stratify=data.bitter)
#     zero_prefixed_num = str(run).zfill(3)
#     train_data.to_csv(
#         f'../data/ts-splits/train_{zero_prefixed_num}.csv',
#         index=None)
#     test_data.to_csv(
#         f'../data/ts-splits/test_{zero_prefixed_num}.csv',
#         index=None)