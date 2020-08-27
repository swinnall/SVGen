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
r_1_TSV = '../output/sumstats/prob_trends/nCycles/1cycle.tsv'
L1 = pd.read_csv(r_1_TSV, sep="\t")

r_2_TSV = '../output/sumstats/prob_trends/nCycles/2cycle.tsv'
L2 = pd.read_csv(r_2_TSV, sep="\t")


# create dataframe
prob_df = pd.DataFrame(list(zip(L1, L2)), columns=['1 Division', '2 Divisions'])

# plot figure
sns_plot = sns.barplot(data=prob_df)
plt.savefig("../output/sumstats/prob_trends/nCycles/chromothProb_nCycles.png")
