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


def readParams(nChroms):

    # read SVGen parameter info
    r_par_TSV = '../output/sw_out' + '/parameters.tsv'
    par_df  = pd.read_csv(r_par_TSV, sep="\t")
    nJ      = par_df.iat[0,0]

    return nJ


def checkChromothripsis(nChroms, analysisType):

    ## Generate Summary Statistics from Criteria ##
    nOscCriteria = 10
    nbpCriteria  = 10

    q = [[nOscCriteria, nbpCriteria] for i in range(nChroms)]

    p = calc_p(nChroms)

    d = calc_d(nChroms, p, q, analysisType)

    return d


def calc_q(nChroms, dataType):

    ## SumStats obtained from valid chromothripsis

    if dataType == 'model':

        q = [ [] for i in range(nChroms)]

        r_chromo_stats_TSV = '../input' + '/model_chromo.tsv'
        model_df = pd.read_csv(r_chromo_stats_TSV, sep="\t")

        ele = 0
        for i in range(0,2*nChroms,2):
            # key: chr, haplo, nbp, nOsc, nDel, nIns, nInv

            # number of junctions
            nJ1 = model_df.iat[i,2]
            nJ2 = model_df.iat[i+1,2]
            q[ele].append( nJ1+nJ2 )

            # number of oscillating cn segments
            nOsc1 = model_df.iat[i,3]
            nOsc2 = model_df.iat[i+1,3]
            q[ele].append( nOsc1+nOsc2 )

            ele += 1


    elif dataType == 'real':

        q = [ [] for i in range(nChroms)]

        r_chromo_stats_TSV = '../input' + '/real_chromo.tsv'
        real_df = pd.read_csv(r_chromo_stats_TSV, sep="\t")

        for i in range(nChroms):
            # key: id, Chr, Start, End, Intrchr. SVs, Total SVs, SVs in sample, Nb. DEL, Nb. DUP, Nb. CN segments, Nb. oscillating CN, chromo_label, ploidy, histo, type_chromothripsis, Nb. breakpoints in chromosome, Fraction SVs in chromothripsis

            # number of junctions
            nJ = 2*real_df.iat[i,15] # double as number of junctions not bp
            q[i].append(nJ)

            # number of oscillating cn segments
            nOsc = real_df.iat[i,10]
            q[i].append(nOsc)

    return q


def calc_p(nChroms):

    ## Read Summary Statistics from Simulation ##
    r_sumStats_TSV = '../output/sw_out' + '/sumStats.tsv'
    sumStat_df = pd.read_csv(r_sumStats_TSV, sep="\t")

    p = [ [] for i in range(nChroms)]

    # average number of junctions as real data isn't haplotype specific
    ele = 0
    for i in range(0,2*nChroms,2):
        # key: chr, haplo, nbp, nOsc, nDel, nIns, nInv

        # number of junctions
        nJ1 = sumStat_df.iat[i,2]
        nJ2 = sumStat_df.iat[i+1,2]
        p[ele].append( nJ1+nJ2 )

        # number of oscillating cn segments
        nOsc1 = sumStat_df.iat[i,3]
        nOsc2 = sumStat_df.iat[i+1,3]
        p[ele].append( nOsc1+nOsc2 )

        ele += 1

    return p


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
                # parameter inference requires normal euclidian distance so no else

            else:
                x += 10

        d[i].append( (x)**0.5 )

    return d


def acceptReject(d, nChroms):
    # scans d; if any chromosome is within acceptable range then simulation
    # parameters and d are saved in memory

    outcome    = False
    chromCount = 0
    for i in range(nChroms):

        if d[i][0] < 1:
            chromCount += 1
            outcome = True

    print("number of chromosomes affected: %s" %chromCount)
    return outcome


def analysis(mem):

    dsbData = []
    for i in range(len(mem)):
        dsbData.append( mem[i][1] )

    dsb_df = pd.DataFrame(dsbData, columns=['nJ'])
    print('total number of dsbs per successful simulation: \n%s' %dsb_df)

    sns_plot = sns.violinplot(y="nJ", data=dsb_df)
    plt.savefig("../output/sw_out/violinplot.png")

    return


def main():

    # outline purpose of program
    #analysisType = 'determine_params'
    analysisType = 'check_chromothripsis'

    # state type of data being read
    dataType = 'model'
    #dataType = 'real'

    # define number of chromosomes
    nChroms = 22

    # define memory matrix
    mem = []

    # ensure output file exists
    src  = '../input/00'
    dest = '../output/00'
    copy(src,dest)

    # import simulation parameters - total number of junctions
    nJ = readParams(nChroms)

    # run simulation N times
    N = 10000
    for i in range(N):

        # generate SVs
        sv_gen.main()

        if analysisType == 'check_chromothripsis':

            # generate distance to SS
            d = checkChromothripsis(nChroms, analysisType)

            # determine validity of simulation
            outcome = acceptReject(d, nChroms)

            # append simulation info to memory
            if outcome == True:
                mem.append( (d, nJ) )
                print("Chromothripsis generated, d = %s." %min(d))
                sys.exit()


        elif analysisType == 'determine_params':

            # model/real data summary statistics
            q = calc_q(nChroms, dataType)

            # current simulation summary statistics
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

    print("number of accepted simulations: %s" %(len(mem)/N))

    if len(mem) > 0 and analysisType == 'determine_params':
        analysis(mem)

    print("End of simulation.")
    return



if __name__ == '__main__':
    print(" Running ABC program.\n")
    main()
