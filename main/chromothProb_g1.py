" Seaborn - Plot Bar Chart; Author: Samuel Winnall "

################################### Notes #####################################


##############################################################################

import shutil
import errno
import sv_gen
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import csv
import sys


# load data
r_0_TSV = '../output/sumstats/prob_trends/unconnected_G1/lambda_0.tsv'
L0 = pd.read_csv(r_0_TSV, sep="\t")

r_5_TSV = '../output/sumstats/prob_trends/unconnected_G1/lambda_5.tsv'
L5 = pd.read_csv(r_5_TSV, sep="\t")

r_10_TSV = '../output/sumstats/prob_trends/unconnected_G1/lambda_10.tsv'
L10 = pd.read_csv(r_10_TSV, sep="\t")

r_50_TSV = '../output/sumstats/prob_trends/unconnected_G1/lambda_50.tsv'
L50 = pd.read_csv(r_50_TSV, sep="\t")

r_100_TSV = '../output/sumstats/prob_trends/unconnected_G1/lambda_100.tsv'
L100 = pd.read_csv(r_100_TSV, sep="\t")

# create dataframe
prob_df = pd.DataFrame(list(zip(L0, L5, L10, L50, L100)), columns=['0', '5', '10', '50', '100'])

# plot figure
sns_plot = sns.barplot(data=prob_df)
plt.savefig("../output/sumstats/prob_trends/unconnected_G1/chromothProb_g1.png")
