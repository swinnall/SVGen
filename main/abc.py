" Approximate Bayesian Computation Analysis; Author: Samuel Winnall "

################################# Overview ###################################

# 1. run breakpoint_model

# 2. read and sort cn information

# 3. calculate distance from summary statistics/criteria

# 4. store distance and parameters in memory

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
    r_cn_TSV = '../output/00' + '/cn_data.tsv'
    cn_df = pd.read_csv(r_cn_TSV, sep="\t")


    # define list of chromosomes
    nChrom    = 22
    cnInfo = [ [] for i in range(nChrom)]


    # store cn tuple info for each chromosome
    for i in range(len(cn_df)):
        # 0: chromID
        chromID = cn_df.iat[i,0]
        # 1: startPos; 2: endPos; 3: CN
        cnInfo[chromID-1].append(  (cn_df.iat[i,1], cn_df.iat[i,2], cn_df.iat[i,3])  ) # -1 for index


    # sort cnInfo tuples by start position for each chromosome
    for i in range(nChrom):
        cnInfo[i].sort(key=lambda x:x[1])


    # read SVGen parameter info
    r_par_TSV = '../output/00' + '/parameters.tsv'
    par_df = pd.read_csv(r_par_TSV, sep="\t")
    mu    = par_df.iat[0,0]
    lmbda = par_df.iat[0,1]
    delta = par_df.iat[0,2]


    return cnInfo, mu, lmbda, delta



def genDist(cnInfo):

    ## Generate Summary Statistics ##
    nChrom  = 22
    sumStat = [ [] for i in range(nChrom)]


    # analyse each chromosome
    for i in range(nChrom):

        # determine number of breakpoints
        nbp = len(cnInfo[i])
        sumStat[i].append(nbp)


        # determine number of oscillating cn segments
        nOscSeg = 0
        if len(cnInfo[i]) > 0:
            for j in range(0,len(cnInfo[i])-1,1):

                if cnInfo[i][j][2] != cnInfo[i][j+1][2]:
                    nOscSeg += 1
        sumStat[i].append(nOscSeg)


    ## Generate Distance ##
    nOscCrit = 10
    nbpCrit  = 10

    d = 1


    return d



def main():

    # define memory matrix
    mem = []


    # ensure output file exists
    src  = '../input/00'
    dest = '../output/00'
    copy(src,dest)


    # run simulation N times
    N = 10000
    for i in range(N):

        # generate SVs
        sv_gen.main()

        # read copy number info
        cnInfo, mu, lmbda, delta = readCN()

        # generate distance to SS
        d = genDist(cnInfo)

        # append simulation info to memory
        mem.append( (d, mu, lmbda, delta) )

    print(mem)
    return



if __name__ == '__main__':
    print(" Running ABC program.\n")
    main()
