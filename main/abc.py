" Approximate Bayesian Computation; Author: Samuel Winnall "

################################# Overview ###################################

# 1. run breakpoint_model

# 2. check summary stats against chromothripsis criteria

# 3. calculate distance from summary statistics

# 4. store distance and parameters

# 5. [ insert new parameters into original model script ]

##############################################################################

import shutil
import errno
import sv_gen
import matplotlib.pyplot as plt
import pandas as pd


def copy(src, dest):
    try:
        shutil.copytree(src, dest)
    except OSError as e:
        # If the error was caused because the source wasn't a directory
        if e.errno == errno.ENOTDIR:
            shutil.copy(src, dest)
        else:
            print('Directory not copied. Error: %s' % e)



def readCN():

    # read SVGen copy number info
    r_cn_TSV = '../output/0' + str(dest) + '/cn_data.tsv'
    cn_df = pd.read_csv(r_cn_TSV, sep="\t")
    print(cn_df)


    # define list of chromosomes
    nChrom    = 22
    cnChanges = [ [] for i in range(nChrom)]


    # store cn tuple info for each chromosome
    for i in range(len(cn_df)):
        # 0: chromID
        chromID = cn_df.iat[i,0]
        # 1: startPos; 2: endPos; 3: CN
        cnChanges[chromID-1].append(  (cn_df.iat[i,1], cn_df.iat[i,2], cn_df.iat[i,3])  ) # -1 for index


    # studies each chromosome's cn info
    test = []
    condition = 5
    for i in range(nChrom):
        # if 10 breakpoints on a given chromosome then simulation passes
        if len( cnChanges[i] ) >= condition:
            test.append(True)
        else: test.append(False)

    return



def main():

    # ensure output file exists
    src  = '../input/00'
    dest = '../output/00'
    copy(src,dest)


    # run simulation N times
    N = 100
    for i in range(N):

        # generate SVs
        sv_gen.main()

        # read copy number info
        readCN()

    return



if __name__ == '__main__':
    print(" Running ABC program.\n")
    main()
