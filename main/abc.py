" Approximate Bayesian Computation Analysis; Author: SanDSBel Winnall "

################################# Overview ###################################

# 1. run SVGen

# 2. read and sort cn information

# 3. check chromothripsis or determine parameters

# 4. ?

##############################################################################

import shutil
import errno
import sv_gen
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys


def copy(src, dest):
    try:
        shutil.copytree(src, dest)
    except OSError as e:
        # If the error was caused because the source wasn't a directory
        if e.errno == errno.ENOTDIR:
            shutil.copy(src, dest)
        else:
            print('Directory not copied. Error: %s' % e)
    return


def readCN(nChroms):

    # read SVGen copy number info
    r_cn_TSV = '../output/00' + '/cn_data.tsv'
    cn_df = pd.read_csv(r_cn_TSV, sep="\t")


    # define list of chromosomes
    cnInfo = [ [] for i in range(nChroms)]


    # store cn tuple info for each chromosome
    for i in range(len(cn_df)):
        # 0: chromID
        chromID = cn_df.iat[i,0]
        # 1: startPos; 2: endPos; 3: CN
        cnInfo[chromID-1].append(  (cn_df.iat[i,1], cn_df.iat[i,2], cn_df.iat[i,3])  ) # -1 for index


    # sort cnInfo tuples by start position for each chromosome
    for i in range(nChroms):
        cnInfo[i].sort(key=lambda x:x[1])


    # read SVGen parameter info
    r_par_TSV = '../output/00' + '/parameters.tsv'
    par_df = pd.read_csv(r_par_TSV, sep="\t")
    nDSB  = par_df.iat[0,0]
    lmbda = par_df.iat[0,1]
    sigma = par_df.iat[0,2]

    # print("\n %s\n" %cnInfo)
    return cnInfo, nDSB, lmbda, sigma


def checkChromothripsis(cnInfo, nChroms):

    ## Generate Summary Statistics from Criteria ##
    nOscCriteria = 10
    nbpCriteria  = 10

    q = [[nOscCriteria, nbpCriteria] for i in range(nChroms)]

    p = calc_p(cnInfo, nChroms)

    d = calc_d(nChroms, p, q)

    return d


def calc_q(nChroms):

    ## SumStats obtained from raw simulation
    q = [ [] for i in range(nChroms)]

    # will read from input folder


    return q


def calc_p(cnInfo, nChroms):

    r_sumStats_TSV = '../output/00' + '/sumStats.tsv'
    sumStat_df = pd.read_csv(r_sumStats_TSV, sep="\t")

    ## Read Summary Statistics from Simulation ##
    p = [ [] for i in range(nChroms)]

    # analyse each chromosome
    for i in range(nChroms):

        # key: chr, haplo, nbp, nOsc, nDel, nIns, nInv

        # number of junctions
        nJ = sumStat_df.iat[i,2]
        p[i].append(nJ)

        # number of oscillating cn segments
        nOsc = sumStat_df.iat[i,3]
        p[i].append(nOsc)

    return p


def calc_d(nChroms, p, q):
    ## Generate Distance ##
    d = [ [] for i in range(nChroms)]

    # for each chromosome generate the distance between each summary statistic
    for i in range(nChroms):
        x = 0
        for j in range(len(p[i])):
            x += (q[i][j] - p[i][j])**2

            # if p is greater than threshold, zero distance therefore undo addition
            if p[i][j] > q[i][j]:
                x -= (q[i][j] - p[i][j])**2

        d[i].append( (x)**0.5 )

    return d


def acceptReject(d, nDSB, lmbda, sigma, nChroms):
    # scans d; if any chromosome is within acceptable range then sinDSBlation
    # parameters and d are saved in memory

    outcome    = False
    chromCount = 0
    for i in range(nChroms):

        if d[i][0] < 5:
            chromCount += 1
            outcome = True

    #print("distance: %s" %d)
    print("number of chromosomes affected: %s" %chromCount)
    return outcome


def analysis(mem):

    dsbData = []
    for i in range(len(mem)):
        print(mem[i][1])
        dsbData.append( mem[i][1] )

    dsb_df = pd.DataFrame(dsbData, columns=['nDSB'])

    sns_plot = sns.violinplot(y="nDSB", data=dsb_df)
    plt.savefig("../output/00/violinplot.png")

    return


def main():

    # outline purpose of program
    analysisType = 'check_chromothripsis' # 'determine_params'

    # define number of chromosomes
    nChroms = 22

    # define memory matrix
    mem = []

    # ensure output file exists
    src  = '../input/00'
    dest = '../output/00'
    copy(src,dest)

    # run sinDSBlation N times
    N = 1000
    for i in range(N):

        # generate SVs
        sv_gen.main()

        # read copy number info
        cnInfo, nDSB, lmbda, sigma = readCN(nChroms)


        if analysisType == 'check_chromothripsis':

            # generate distance to SS
            d = checkChromothripsis(cnInfo, nChroms)

            # determine validity of simulation
            outcome = acceptReject(d, nDSB, lmbda, sigma, nChroms)

            # append sinDSBlation info to memory
            if outcome == True:
                mem.append( (d, nDSB, lmbda, sigma) )
                print("Chromothripsis generated, d = %s." %min(d))
                sys.exit()


        elif analysisType == 'determine_params':

            # model/real data summary statistics
            q = calc_q(nChroms)

            # current simulation summary statistics
            p = calc_p(cnInfo, nChroms)

            # distance between statistics
            d = calc_d(p, q)

            # determine validity of simulation
            outcome = acceptReject(d, nDSB, lmbda, sigma, nChroms)

            # append sinDSBlation info to memory
            if outcome == True:
                mem.append( (d, nDSB, lmbda, sigma) )

    #print(mem)
    print("number of accepted sinDSBlations: %s" %(len(mem)/N))

    if len(mem) > 0:
        analysis(mem)

    print("End of simulation.")
    return



if __name__ == '__main__':
    print(" Running ABC program.\n")
    main()
