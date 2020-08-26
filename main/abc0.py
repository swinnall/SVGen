" Seaborn - Plot Bar Chart; Author: Samuel Winnall "

################################### Notes #####################################


##############################################################################

import shutil
import errno
import sv_gen
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys


# load
tips = sns.load_dataset("tips")


# create dataframe
prob_df = pd.DataFrame(dsbData, columns=['Neutral', 'Biased'])


# plot figure
sns_plot = sns.barplot(data=prob_df)
plt.savefig("../output/sumstats/chromothProb.png")
