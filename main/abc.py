" Approximate Bayesian Computation Analysis; Author: Samuel Winnall "

################################### Notes #####################################


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


def readSumStatsTotal(nChroms):

    # read SVGen parameter info
    r_par_TSV = '../output/sumstats/sumStats_total.tsv'
    par_df    = pd.read_csv(r_par_TSV, sep="\t")
    nDSB      = par_df.iat[0,0]
    nDel      = par_df.iat[0,1]
    nInv      = par_df.iat[0,2]
    nIns      = par_df.iat[0,3]
    nDup      = par_df.iat[0,4]

    return nDSB, nDel, nInv, nIns, nDup


def checkChromothripsis(nChroms, analysisType):

    ## Generate Summary Statistics from Criteria ##
    nOscCriteria = 10
    nbpCriteria  = 10

    q = [[nOscCriteria, nbpCriteria] for i in range(nChroms)]

    # parameters from simulation
    p = calc_p(nChroms)

    # distance between the two
    d = calc_d(nChroms, p, q, analysisType)

    return d


def calc_p(nChroms):

    ## Read Summary Statistics from Simulation ##
    r_sumStats_TSV = '../output/sumstats/sumStats_chrom.tsv'
    sumStat_df = pd.read_csv(r_sumStats_TSV, sep="\t")

    p = [ [] for i in range(nChroms)]

    for i in range(nChroms):
        # key: chr, nDSB, nOsc, nDel, nIns, nInv, nDup

        # number of junctions
        nDSB = sumStat_df.iat[i,1]
        p[i].append( nDSB )

        # number of oscillating cn segments
        nOsc = sumStat_df.iat[i,2]
        p[i].append( nOsc )

    return p


def calc_q(nChroms, dataType):

    ## SumStats obtained from valid chromothripsis

    if dataType == 'model':

        q = [ [] for i in range(nChroms)]

        r_chromo_stats_TSV = '../input/model_chromo.tsv'
        model_df = pd.read_csv(r_chromo_stats_TSV, sep="\t")

        for i in range(nChroms):
            # key: chr, nDSB, nOsc, nDel, nIns, nInv, nDup

            # number of double strand breaks
            nDSB = model_df.iat[i,1]
            q[i].append( nDSB )

            # number of oscillating cn segments
            nOsc = model_df.iat[i,2]
            q[i].append( nOsc )


    elif dataType == 'real':

        q = [ [] for i in range(nChroms)]

        r_chromo_stats_TSV = '../input/real_chromo.tsv'
        real_df = pd.read_csv(r_chromo_stats_TSV, sep="\t")

        for i in range(nChroms):
            # key: id, Chr, Start, End, Intrchr. SVs, Total SVs, SVs in sample, Nb. DEL, Nb. DUP, Nb. CN segments, Nb. oscillating CN, chromo_label, ploidy, histo, type_chromothripsis, Nb. breakpoints in chromosome, Fraction SVs in chromothripsis

            # number of double strand breaks
            nDSB = real_df.iat[i,15]
            q[i].append(nDSB)

            # number of oscillating cn segments
            nOsc = real_df.iat[i,10]
            q[i].append(nOsc)

    return q


def calc_d(nChroms, p, q, analysisType):
    ## Generate Distance ##

    d = [ [] for i in range(nChroms)]

    # for each chromosome generate the distance between each summary statistic
    for i in range(nChroms):
        x = 0
        for j in range(len(p[i])):

            # choose non zero chromosomes for true comparison
            if p[i][j] != 0 or q[i][j] != 0:
                x += (q[i][j] - p[i][j])**2

                # if p is greater than threshold, zero distance therefore undo addition
                if analysisType == 'check_chromothripsis' and p[i][j] > q[i][j]:
                    x -= (q[i][j] - p[i][j])**2
                # parameter inference requires normal euclidian distance so no else statement

            else:
                x += 10

        d[i].append( (x)**0.5 )

    return d


def acceptReject(d, nChroms):
    # scans d; if any chromosome is within acceptable range then outcome == true
    # d and parameters are then saved in memory

    outcome    = False
    chromCount = 0
    for i in range(nChroms):

        if d[i][0] < 1:
            chromCount += 1
            outcome = True

    print("number of chromosomes affected: %s" %chromCount)
    return outcome


def plotAnalysis(analysisType, dataType, mem):

    if analysisType == 'countSV':

        sv_data = []
        for i in range(len(mem)):
            sv_data.append( (mem[i][0], mem[i][1], mem[i][2], mem[i][3], mem[i][4]) )

        sv_df = pd.DataFrame(sv_data, columns=['nDSB','nDel','nInv','nIns', 'nDup'])

        sns_plot = sns.violinplot(data=sv_df) #x="SV Type", y="Count",


        plt.savefig("../output/sumstats/sv_count.png")


    elif analysisType == 'determine_params':

        dsbData = []
        for i in range(len(mem)):
            dsbData.append( mem[i][1] )

        dsb_df = pd.DataFrame(dsbData, columns=['nDSB'])
        print('total number of dsbs per successful simulation: \n%s' %dsb_df)

        sns_plot = sns.violinplot(y="nDSB", data=dsb_df)
        plt.savefig("../output/sumstats/mutationRate_" + str(dataType) + "/.png")

    return


def main():

    ## outline purpose of program
    #analysisType = 'countSV'
    analysisType = 'check_chromothripsis'
    #analysisType = 'determine_params'

    ## state type of data being read
    dataType = 'model'
    #dataType = 'real'

    # define number of chromosomes
    nChroms = 22

    # define memory matrix for storing simulation statistics
    mem = []

    # run simulation N times
    N = 100
    for i in range(N):

        # generate SVs
        sv_gen.main()

        if analysisType == 'countSV':

            # import summary statistics of whole simulation
            nDSB, nDel, nInv, nIns, nDup = readSumStatsTotal(nChroms)

            # save to memory
            mem.append( (nDSB, nDel, nInv, nIns, nDup))


        elif analysisType == 'check_chromothripsis':

            # generate distance to SS
            d = checkChromothripsis(nChroms, analysisType)

            # determine validity of simulation
            outcome = acceptReject(d, nChroms)

            # append simulation info to memory
            if outcome == True:
                mem.append( (d, nDSB) )
                print("Chromothripsis generated, d = %s." %min(d))
                sys.exit()


        elif analysisType == 'determine_params':

            # model/real data summary statistics by chromosome
            q = calc_q(nChroms, dataType)

            # current simulation summary statistics by chromosome
            p = calc_p(nChroms)

            # distance between statistics
            d = calc_d(nChroms, p, q, analysisType)

            # determine validity of simulation
            outcome = acceptReject(d, nChroms)

            # append simulation info to memory
            if outcome == True:
                mem.append( (d, p[2]) )

            # clear variables
            p.clear()
            q.clear()
            d.clear()


    if len(mem) > 0 and (analysisType == 'countSV' or analysisType == 'determine_params'):
        print("number of accepted simulations: %s" %(len(mem)/N))
        plotAnalysis(analysisType, dataType, mem)


    mem.clear()
    print("End of simulation.")
    return



if __name__ == '__main__':
    print(" Running ABC program.\n")
    main()
