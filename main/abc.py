" Approximate Bayesian Computation Analysis; Author: SanDSBel Winnall "

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


def checkChromothripsis(cnInfo, nChroms, analysisType, dataType):

    ## Generate Summary Statistics from Criteria ##
    nOscCriteria = 10
    nbpCriteria  = 10

    q = [[nOscCriteria, nbpCriteria] for i in range(2*nChroms)]

    p = calc_p(cnInfo, nChroms, dataType)

    d = calc_d(nChroms, p, q, analysisType, dataType)

    return d


def calc_q(nChroms, dataType):

    ## SumStats obtained from valid chromothripsis

    if dataType == 'model':

        q = [ [] for i in range(2*nChroms)]

        r_chromo_stats_TSV = '../input' + '/model_chromo.tsv'
        model_df = pd.read_csv(r_chromo_stats_TSV, sep="\t")

        for i in range(2*nChroms):
            # key: chr, haplo, nbp, nOsc, nDel, nIns, nInv

            # number of junctions
            nJ = model_df.iat[i,2]
            q[i].append(nJ)

            # number of oscillating cn segments
            nOsc = model_df.iat[i,3]
            q[i].append(nOsc)


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


    #print("\nq: \n%s\n" %q)
    return q


def calc_p(cnInfo, nChroms, dataType):

    ## Read Summary Statistics from Simulation ##
    r_sumStats_TSV = '../output/00' + '/sumStats.tsv'
    sumStat_df = pd.read_csv(r_sumStats_TSV, sep="\t")

    if dataType == 'model':

        # 2*nChroms for number of haplotypes
        p = [ [] for i in range(2*nChroms)]

        for i in range(0,2*nChroms,1):
            # key: chr, haplo, nbp, nOsc, nDel, nIns, nInv

            # number of junctions
            nJ = sumStat_df.iat[i,2]
            p[i].append(nJ)

            # number of oscillating cn segments
            nOsc = sumStat_df.iat[i,3]
            p[i].append(nOsc)

    if dataType == 'real':

        # average number of junctions as real data isn't haplotype specific
        for i in range(0,2*nChroms,2):
            # key: chr, haplo, nbp, nOsc, nDel, nIns, nInv

            # number of junctions
            nJ1 = sumStat_df.iat[i,2]
            nJ2 = sumStat_df.iat[i+1,2]
            p[i].append( (nJ1+nJ2)/2 )

            # number of oscillating cn segments
            nOsc1 = sumStat_df.iat[i,3]
            nOsc2 = sumStat_df.iat[i+1,3]
            p[i].append( (nOsc1+nOsc2)/2 )


    #print("\np: \n%s\n" %p)
    return p


def calc_d(nChroms, p, q, analysisType, dataType):
    ## Generate Distance ##

    if dataType == 'model':
        gma = 2*nChroms
    elif dataType == 'real':
        gma = nChroms

    d = [ [] for i in range(gma)]

    # for each chromosome generate the distance between each summary statistic
    for i in range(gma):
        x = 0
        for j in range(len(p[i])):

            # choose non zero chromosomes for true comparison
            if p[i][j] != 0 or q[i][j] != 0:
                x += (q[i][j] - p[i][j])**2

                # if p is greater than threshold, zero distance therefore undo addition
                if analysisType == 'check_chromothripsis' and p[i][j] > q[i][j]:
                    x -= (q[i][j] - p[i][j])**2
                # parameter inference requires normal euclidian distance

            else:
                x += 10

        d[i].append( (x)**0.5 )

    return d


def acceptReject(d, nDSB, lmbda, sigma, nChroms):
    # scans d; if any chromosome is within acceptable range then sinDSBlation
    # parameters and d are saved in memory

    outcome    = False
    chromCount = 0
    for i in range(nChroms):

        if d[i][0] < 1:
            chromCount += 1
            outcome = True

    #print("distance: %s" %d)
    print("number of chromosomes affected: %s" %chromCount)
    return outcome


def analysis(mem):

    dsbData = []
    for i in range(len(mem)):
        dsbData.append( mem[i][1] )

    dsb_df = pd.DataFrame(dsbData, columns=['nDSB'])
    print('total number of dsbs per successful simulation: \n%s' %dsb_df)

    sns_plot = sns.violinplot(y="nDSB", data=dsb_df)
    plt.savefig("../output/00/violinplot.png")

    return


def main():

    # outline purpose of program
    analysisType = 'determine_params'
    # analysisType = 'check_chromothripsis'

    # state type of data being read
    dataType = 'model'
    # dataType = 'real'

    # define number of chromosomes
    nChroms = 22

    # define memory matrix
    mem = []

    # ensure output file exists
    src  = '../input/00'
    dest = '../output/00'
    copy(src,dest)

    # run sinDSBlation N times
    N = 100

    for i in range(N):

        # generate SVs
        sv_gen.main()

        # read copy number info
        cnInfo, nDSB, lmbda, sigma = readCN(nChroms)


        if analysisType == 'check_chromothripsis':

            # generate distance to SS
            d = checkChromothripsis(cnInfo, nChroms, analysisType, dataType)

            # determine validity of simulation
            outcome = acceptReject(d, nDSB, lmbda, sigma, nChroms)

            # append sinDSBlation info to memory
            if outcome == True:
                mem.append( (d, nDSB, lmbda, sigma) )
                print("Chromothripsis generated, d = %s." %min(d))
                sys.exit()


        elif analysisType == 'determine_params':

            # model/real data summary statistics
            q = calc_q(nChroms, dataType)

            # current simulation summary statistics
            p = calc_p(cnInfo, nChroms, dataType)

            # distance between statistics
            d = calc_d(nChroms, p, q, analysisType, dataType)

            # determine validity of simulation
            outcome = acceptReject(d, nDSB, lmbda, sigma, nChroms)

            # append sinDSBlation info to memory
            if outcome == True:
                mem.append( (d, nDSB, lmbda, sigma) )


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
