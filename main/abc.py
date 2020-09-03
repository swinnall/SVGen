" Approximate Bayesian Computation Analysis; Author: Samuel Winnall "

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
import argparse
import os

def readSumStatsTotal(nChroms, sumStats_total, r_par_TSV=''):
    # read SVGen parameter info
    if r_par_TSV != '':
        r_par_TSV = '../output/sumstats/sumStats_total.tsv'
        par_df = pd.read_csv(r_par_TSV, sep="\t")
    else:
        par_df = sumStats_total

    nDSB      = par_df.iat[0,0]
    nDel      = par_df.iat[0,1]
    nInv      = par_df.iat[0,2]
    nIns      = par_df.iat[0,3]
    nDup      = par_df.iat[0,4]
    cycleID   = par_df.iat[0,5]
    nBiasedChroms = par_df.iat[0,6]
    nMitosisBreaks = par_df.iat[0,7]

    return nDSB, nDel, nInv, nIns, nDup, nBiasedChroms, nMitosisBreaks


def checkChromothripsis(nChroms, analysisType, sumStats_chrom, r_sumStats_TSV = ''):

    ## Generate Summary Statistics from Criteria ##
    nOscCriteria = 10
    nbpCriteria  = 10

    q = [[nOscCriteria, nbpCriteria] for i in range(nChroms)]

    # parameters from simulation
    p = calc_p(nChroms, sumStats_chrom, r_sumStats_TSV)

    # distance between the two
    d = calc_d(nChroms, p, q, analysisType)

    return d


def calc_p(nChroms, sumStats_chrom, r_sumStats_TSV = ''):

    ## Read Summary Statistics from Simulation ##
    # r_sumStats_TSV = '../output/sumstats/sumStats_chrom.tsv'
    if r_sumStats_TSV != '':
        sumStat_df = pd.read_csv(r_sumStats_TSV, sep="\t")
    else:
        sumStat_df = sumStats_chrom

    p = [ [] for i in range(nChroms)]

    for i in range(nChroms):
        # key: chr, nDSB, nOsc, nDel, nIns, nInv, nDup

        # number of DSBs
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

        r_chromo_stats_TSV = '../input/model/sumStats_chrom.tsv'
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

        r_chromo_stats_TSV = '../input/PCAWG/real_chromo.tsv'
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


def acceptReject(d, nChroms, verbose = 0):
    # scans d; if any chromosome is within acceptable range then outcome == true
    # d and parameters are then saved in memory

    outcome    = False
    chromCount = 0
    for i in range(nChroms):

        # set to 4 to consider low confidence chromothripsis cases
        if d[i][0] <= 4:
            chromCount += 1
            outcome = True

    if verbose > 0:
        print("number of chromosomes affected: %s" %chromCount)

    return outcome


def plotAnalysis(analysisType, dataType, mem):

    if analysisType == 'countSV':

        sv_data = []
        for i in range(len(mem)):
            # extract data from successful summary statistics stored in memory
            sv_data.append( (mem[i][0], mem[i][1], mem[i][2], mem[i][3], mem[i][4], mem[i][5]) )

        # create dataframe
        sv_df = pd.DataFrame(sv_data, columns=['nDSB','nDel','nInv','nIns', 'nDup', 'nMitosisBreaks'])

        sns_plot = sns.violinplot(data=sv_df)
        plt.savefig("../output/sumstats/sv_count.png")


    elif analysisType == 'determine_params':

        dsbData = []
        for i in range(len(mem)):
            # extract data from successful summary statistics stored in memory
            dsbData.append( mem[i][1] )

        # create dataframe
        dsb_df = pd.DataFrame(dsbData, columns=['nDSB'])
        print('total number of dsbs per successful simulation: \n%s' %dsb_df)

        sns_plot = sns.violinplot(y="nDSB", data=dsb_df)
        plt.savefig("../output/sumstats/mutationRate_" + str(dataType) + ".png")

    return


def main():
    parser = argparse.ArgumentParser(description='Running abc')

    default = 10000
    parser.add_argument('-N', '--num_sim', default=default, type=int, help='number of simulations [{}]'.format(default))

    default = 'countSV'
    parser.add_argument('-a', '--analysisType', default=default, choices=['countSV', 'check_chromothripsis', 'determine_params'], help='type of analysis [{}]'.format(default))

    default = 'biased'
    parser.add_argument('-b', '--chr_bias', default=default, choices=['', 'random', 'biased', 'fixed'], help='type of analysis [{}]'.format(default))

    default = 5
    parser.add_argument('-l', '--lmbda', default=default, type=int, help='lambda [{}]'.format(default))

    default = 0
    parser.add_argument('-m', '--minDSB', default=default, type=int, help='minimum number of DSBs [{}]'.format(default))

    default = 40
    parser.add_argument('-M', '--maxDSB', default=default, type=int, help='maximum number of DSBs [{}]'.format(default))

    default = 2
    parser.add_argument('-n', '--nCycles', default=default, type=int, help='number of cell cycles [{}]'.format(default))

    default = '../output/sumstats/prob_trends/nCycles/'
    parser.add_argument('-o', '--outdir', default=default, help='output directory [{}]'.format(default))

    default = '2cycle'
    parser.add_argument('-p', '--prefix', default=default, help='prefix of output file [{}]'.format(default))

    default = False
    parser.add_argument('-c', '--write_cicos', default=default, action = 'store_true', help='whether or not to write out files for circos plots [{}]'.format(default))

    default = False
    parser.add_argument('-k', '--write_shatterseek', default=default, action = 'store_true', help='whether or not to write out files for shatterseek plots[{}]'.format(default))

    default = False
    parser.add_argument('-s', '--write_sumStats', default=default, action = 'store_true', help='whether or not to write out summary statistics [{}]'.format(default))

    default = '../output/sumstats/sumStats_chrom.tsv'
    parser.add_argument('--r_sumStats_TSV', default=default, help='the file with summary statistics of each chromosome [{}]'.format(default))

    default = '../output/sumstats/sumStats_total.tsv'
    parser.add_argument('--r_par_TSV', default=default, help='the file with summary statistics of all chromosomes [{}]'.format(default))

    default = 0
    parser.add_argument('--verbose', default=default, type=int, help='detail level of output information [{}]'.format(default))


    args = parser.parse_args()

    analysisType = args.analysisType
    matType = args.chr_bias
    lmbda = args.lmbda
    minDSB = args.minDSB
    maxDSB = args.maxDSB
    nCycles = args.nCycles
    outdir = args.outdir
    prefix = args.prefix
    write_cicos = args.write_cicos
    write_shatterseek = args.write_shatterseek
    write_sumStats = args.write_sumStats
    r_sumStats_TSV = args.r_sumStats_TSV    # not used for now
    r_par_TSV = args.r_par_TSV   # not used for now
    verbose = args.verbose

    ## outline purpose of program
    # if analysisType == 'sv':
    #     analysisType = 'countSV'
    # elif analysisType == 'check':
    #     analysisType = 'check_chromothripsis'
    # else:
    #     analysisType = 'determine_params'

    print('type of analysis is {}\n'.format(analysisType))

    ## state type of data being read
    dataType = 'model'
    #dataType = 'real'

    # define number of chromosomes
    nChroms = 22

    # define memory matrix for storing simulation statistics
    mem = []

    # run simulation N times
    N = args.num_sim
    for i in range(N):
        print("simulation {}\n".format(i))
        # generate SVs
        sumStats_chrom, sumStats_total = sv_gen.main(matType, lmbda, minDSB, maxDSB, nCycles, write_cicos, write_shatterseek, write_sumStats, verbose)

        # print(sumStats_chrom)
        # print(sumStats_total)

        r_par_TSV = ''
        r_sumStats_TSV = ''

        # generating violin plot for counting SV accumulation
        if analysisType == 'countSV':

            # import summary statistics of whole simulation
            nDSB, nDel, nInv, nIns, nDup, nBiasedChroms, nMitosisBreaks = readSumStatsTotal(nChroms, sumStats_total, r_par_TSV)

            # save to memory
            mem.append( (nDSB, nDel, nInv, nIns, nDup, nMitosisBreaks))


        # runs SVGen until chromothripsis criteria is met,
        # either stops program for analysis of that particular case
        # or counts number of chromothripsis cases in N iterations
        elif analysisType == 'check_chromothripsis':

            # import summary statistics of SVGen
            nDSB, nDel, nInv, nIns, nDup, nBiasedChroms, nMitosisBreaks = readSumStatsTotal(nChroms, sumStats_total, r_par_TSV)

            # generate distance to summary statistics
            d = checkChromothripsis(nChroms, analysisType, sumStats_chrom, r_sumStats_TSV)

            # determine validity of simulation
            outcome = acceptReject(d, nChroms, verbose)

            # append simulation info to memory
            if outcome == True:
                mem.append( (d, nDSB, nBiasedChroms) )
                print("Chromothripsis generated, d = %s." %min(d))
                ## uncomment to stop immediately at chromothripsis
                ## leave commented for studying distribution of chromothripsis
                #sys.exit()


        # for determining parameters / number of DSBs in chromothripsis
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
                total_DSB = 0
                for i in range(nChroms):
                    total_DSB += p[i][0]
                mem.append( (d, total_DSB) )


    if len(mem) > 0 and (analysisType == 'countSV' or analysisType == 'determine_params'):
        plotAnalysis(analysisType, dataType, mem)
        print("number of accepted simulations: %d" %(len(mem)/N))

    if analysisType == 'check_chromothripsis':
        print("number of accepted simulations: %d" %(len(mem)/N))
        print("mem: %s" % mem)

        fout = os.path.join(outdir, prefix + '.tsv')
        print("Writing output to file {}".format(fout))
        with open(fout, 'w', newline='') as file:
            writer = csv.writer(file, delimiter = '\t')
            writer.writerow([len(mem)/N])


    mem.clear()
    print("\n~~ End of simulation ~~")
    return



if __name__ == '__main__':
    print(" Running ABC program.\n")
    main()
