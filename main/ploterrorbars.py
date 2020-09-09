" Seaborn - Plot Bar Chart; Author: Samuel Winnall "

################################### Notes #####################################


##############################################################################

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import csv

def calcStats(p_df,i):

    # number of points
    n = len(p_df[i])

    # statistics
    q    = [0 for i in range(n)]
    mean = [0 for i in range(n)]
    var  = [0 for i in range(n)]
    sd   = [0 for i in range(n)]

    # calculations
    for j in range(n):
        # q data
        q[j] = 1 - p_df[i][j]
        # calculate mean, Np
        mean[j] = N*p_df[i][j]
        # calculate variance, Npq
        var[j] = N*p_df[i][j]*q[j]
        # calculate sd
        sd[j] = var[j]**0.5

    # set error bar positions
    yerr = [0 for i in range(n)]
    for i in range(n):
        yerr[i] = sd[i]/2

    return mean, yerr


################################################################################

# number of samples
N = 10000

# probability data
pSel = [0.0035, 0.0193, 0.0205, 0.0340] # Selection experiment
pDSB = [0, 0.025, 0.0232, 0.0056, 0.004]
pG1  = [0.0004, 0.0174, 0.0464, 0.1086, 0.1104] # G1 experimnet
pCyc = [0.0019, 0.0177]

p_df = [pSel, pDSB, pG1, pCyc]

# store data
y0 = []; y1 = []; y2 = []; y3 = []
means_df = []
error_df = []
for i in range(len(p_df)):

    # pass probabilities to statistics function
    mean, yerr = calcStats(p_df,i)

    # store output in memory
    means_df.append(mean)
    error_df.append(yerr)


# plot figure
fig, axs = plt.subplots(2, 2)
axs[0, 0].errorbar(['Random', 'Biased', 'Fixed 1', 'Fixed 2'], means_df[0], error_df[0], marker='s', ls='none', ecolor='k')
axs[0, 1].errorbar(['0', '5', '10', '50', '100'], means_df[1], error_df[1], color='orange', marker='s', ls='none', ecolor='k')
axs[1, 0].errorbar(['0', '5', '10', '50', '100'], means_df[2], error_df[2], color='g', marker='s', ls='none', ecolor='k')
axs[1, 1].errorbar(['1', '2'], means_df[3], error_df[3], color='r', marker='s', ls='none', ecolor='k')

#axs[0,0].set(xlabel='Chromosome Selection Type')
#axs[0,1].set(xlabel='Number of DSBs per Cycle')
#axs[1,0].set(xlabel='Number of Unrepaired Breaks Entering S Phase per Cycle')
#axs[1,1].set(xlabel='Number of Cell Divisions')

axs[0,0].set(ylabel='Count')
axs[1,0].set(ylabel='Count')



#plt.figure()
#plt.errorbar(x=columns, y=mean, yerr=yerr, marker='s', ls='none')
#plt.errorbar(x=columns, y=numdf, fmt='none', yerror=[ylow, yupp], ecolor='k', elinewidth=2)

plt.savefig("../output/sumstats/prob_trends/fourPanels.png")
