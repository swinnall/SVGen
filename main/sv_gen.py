" Evolutionary Model of Structural Variation (EMoSV); Author: Samuel Winnall "

################################# comments ###################################

" analysis - insertions: currently not accycleIDing for BFB unconnected junctions "

" analysis - duplications: same problem as above "

" analysis - inversions: may be cycleIDed twice, not a big problem"

" do i need to return nodes everytime i change a key? "

################################# imports ####################################

import numpy as np
import pandas as pd
import random
import csv
import sys

################################# classes ####################################


class ConnectionsClass():
    chrom = 0
    side  = 0
    id_   = 0
    haplo = 0
    cnid  = 0
    cn    = 0

    def __init__(self, chrom, side, id_, haplo, cnid, cn):
        self.chrom
        self.side
        self.id_
        self.haplo
        self.cnid
        self.cn

connection = ConnectionsClass(0,0,0,0,0,0)


class cirosPlot():
    chr1     = 0
    coord1   = 0
    strand1  = '+'
    chr2     = 0
    coord2   = 0
    strand2  = '+'
    extra    = 0
    hap      = 0
    cycleNum  = 0

    def __init__(self, chr1, coord1, strand1, chr2, coord2, strand2, extra, hap, cycleNum):
        self.chr1
        self.coord1
        self.strand1
        self.chr2
        self.coord2
        self.strand2
        self.extra
        self.hap
        self.cycleNum

    def deletions(self, nodes, cycleID, dest):
        for i in range(len(nodes)):
            if nodes[i].get("cn") == 0:
                chr1     = nodes[i].get("chromID")
                coord1   = nodes[i].get("position")
                strand1  = '+'
                chr2     = chr1
                coord2   = coord1
                strand2  = '-'
                extra    = 'svtype=DEL'
                cycleNum  = cycleID

                # write to csv
                with open('../output/0' +  str(dest) + '/sv_data.tsv', 'a', newline='') as file:
                    writer = csv.writer(file, delimiter = '\t')
                    writer.writerow([chr1, coord1, strand1, chr2, coord2, strand2,	extra, cycleNum])
        return

    def insertions(self, nodes, cycleID, dest):
        for i in range(len(nodes)):

            if nodes[i].get("cn") > 0 and nodes[i].get("M") != 'none' and \
                nodes[i].get("position") != nodes[ int(nodes[i].get("M")) ].get("position"):

                chr1     = nodes[i].get("chromID")
                coord1   = nodes[i].get("position")
                strand1  = '+'
                chr2     = nodes[ int(nodes[i].get("M")) ].get("chromID")
                coord2   = nodes[ int(nodes[i].get("M")) ].get("position")
                strand2  = '+'
                extra    = 'svtype=INS'
                cycleNum  = cycleID

                # write to csv
                with open('../output/0' +  str(dest) + '/sv_data.tsv', 'a', newline='') as file:
                    writer = csv.writer(file, delimiter = '\t')
                    writer.writerow([chr1, coord1, strand1, chr2, coord2, strand2,	extra, cycleNum])

            else: pass

    def inversions(self, nodes, cycleID, dest, chromLengths):
        for i in range(len(nodes)):

            if nodes[i].get("cn") > 0 and nodes[i].get("inv") == True:

                if nodes[i].get("type") == 'nonTel':

                    chr1     = nodes[i].get("chromID")
                    coord1   = nodes[i].get("position")
                    strand1  = '+'
                    chr2     = chr1
                    coord2   = nodes[ int(nodes[i].get("WT")) ].get("position")
                    strand2  = '+'
                    extra    = 'svtype=INV'
                    cycleNum  = cycleID

                elif nodes[i].get("type") == 'pTel':

                    chr1     = nodes[i].get("chromID")
                    coord1   = nodes[i].get("position")
                    strand1  = '+'
                    chr2     = chr1
                    coord2   = 0
                    strand2  = '+'
                    extra    = 'svtype=INV'
                    cycleNum  = cycleID

                elif nodes[i].get("type") == 'qTel':

                    chr1     = nodes[i].get("chromID")
                    coord1   = nodes[i].get("position")
                    strand1  = '+'
                    chr2     = chr1
                    coord2   = chromLengths[1][chr1-1]
                    strand2  = '+'
                    extra    = 'svtype=INV'
                    cycleNum  = cycleID

                # write to csv
                with open('../output/0' +  str(dest) + '/sv_data.tsv', 'a', newline='') as file:
                    writer = csv.writer(file, delimiter = '\t')
                    writer.writerow([chr1, coord1, strand1, chr2, coord2, strand2,	extra, cycleNum])
        return

    def duplications(self, nodes, cycleID, dest, chromLengths):
        coveredNodes = []

        # debugging
        #print("\nAnalysis - Duplications Error:\n")

        for i in range(len(nodes)):
            nodeID = nodes[i].get("nodeID")
            coveredNodes.append(nodeID)

            if nodes[i].get("type") == 'nonTel' and nodes[i].get("M") != 'none':

                # debugging
                # print("Node causing error: %s" %nodes[i])


                AdjacentID = nodes[ int(nodes[i].get("WT")) ].get("nodeID")

                if AdjacentID not in coveredNodes:

                    chr1    = nodes[i].get("chromID")
                    start   = nodes[i].get("position")
                    end     = nodes[AdjacentID].get("position")
                    cn      = nodes[i].get("cn")
                    hap     = nodes[i].get("haplotype")
                    cycleNum = cycleID

                    # write to csv
                    with open('../output/0' +  str(dest) + '/cn_data.tsv', 'a', newline='') as file:
                        writer = csv.writer(file, delimiter = '\t')
                        writer.writerow([chr1, start, end, cn, hap, cycleNum])
                else: pass


            elif nodes[i].get("type") == 'pTel':

                chr1    = nodes[i].get("chromID")
                start   = 0
                end     = nodes[i].get("position")
                cn      = nodes[i].get("cn")
                hap     = nodes[i].get("haplotype")
                cycleNum = cycleID

                # write to csv
                with open('../output/0' +  str(dest) + '/cn_data.tsv', 'a', newline='') as file:
                    writer = csv.writer(file, delimiter = '\t')
                    writer.writerow([chr1, start, end, cn, hap, cycleNum])


            elif nodes[i].get("type") == 'qTel':

                chr1    = nodes[i].get("chromID")
                start   = nodes[i].get("position")
                end     = chromLengths[1][chr1-1]
                cn      = nodes[i].get("cn")
                hap     = nodes[i].get("haplotype")
                cycleNum = cycleID

                # write to csv
                with open('../output/0' +  str(dest) + '/cn_data.tsv', 'a', newline='') as file:
                    writer = csv.writer(file, delimiter = '\t')
                    writer.writerow([chr1, start, end, cn, hap, cycleNum])

            else: pass

        return

calcSVs = cirosPlot(0,0,'+',0,0,'+',0,0,0)


class CheckBool():
    def __init__(self):
        return

    def endGrowth(nodes, lmbda):
        tally = 0
        for i in range(len(nodes)):

            if nodes[i].get("M") == 'none' and nodes[i].get("cn") > 0:
                tally += 1
            else: pass

        print("%d remaining junctions.\n" %tally)

        if tally <= lmbda:
            endCondition = False
            print("End of growth function.\n")
        else:
            endCondition = True

        return endCondition

    def TelCheck(nodes,nodeID):
        telCondition = True

        if nodes[nodeID].get("type") == 'qTel' or nodes[nodeID].get("type") == 'pTel':
            telCondition = False

        else:
            telCondition = True

        return telCondition

    def checkCentromere(nodes,pathList,i):

        tally    = 0
        centList = []
        print('\nnext cent list info: ')
        for j in range(len(pathList.get(str(i)))):
            nodeID = pathList.get(str(i))[j]

            if nodes[nodeID].get("centromeric") == True and nodes[nodeID].get("type") != 'nonTel':
                tally += 1
                print(tally)
                centList.append(nodeID)

            # nonTel types have two junctions referring to the same centromere
            elif nodes[nodeID].get("centromeric") == True and nodes[nodeID].get("type") == 'nonTel':
                tally += 0.5
                print(tally)
                centList.append(nodeID)

        nCent = tally
        return nCent, centList



################################ functions ###################################


def generateDSBs(mu):

    # nDSB = int( np.random.poisson(mu, 1) )
    nDSB = int( random.randint(0,mu) )

    return nDSB


def generateNodes(nodes,nDSB,nChroms,chromLengths,firstEvent):
    print("\nEntering node generation\n")

    # reset prev path information
    for i in range(len(nodes)):
        if  nodes[i].get("cn") > 0:
            nodes[i]["type"]        = 'nonTel'
            nodes[i]["WT"]          = 'none'
            nodes[i]["centromeric"] = False
            nodes[i]["inv"]         = False

    # determine nbp on each chromosome prior to additional assignment; if 0 then cn = 1; else: call check func later
    # stores number of break points per chrom (row) for each haplotype (column)
    nbp_per_chrom = [ [0,0] for i in range(nChroms)]
    for i in range(len(nodes)):
        hapIndex = nodes[i].get("haplotype")
        chrIndex = nodes[i].get("chromID")-1

        nbp_per_chrom[chrIndex][hapIndex] = nbp_per_chrom[chrIndex][hapIndex] + 1

    # print bp information
    # print(len(nbp_per_chrom))
    # print(nbp_per_chrom)

    # local definitions
    uniqueID = len(nodes)
    hapChoice = [0,1]

    # nThreshold is the marker between nNodes in each cell cycle
    nThreshold = len(nodes)

    # twice the number of breakpoints added (LR junctions)
    for i in range(0, 2*nDSB, 2):

        # generate positional information
        chromosomeTarget = random.randint(1, nChroms)
        breakpointPos    = random.randint(0, chromLengths[1][chromosomeTarget-1])
        haplotype        = np.random.choice(hapChoice, 1)

        nodeData = {
            # identification:
            "nodeID":    uniqueID,
            "chromID":   chromosomeTarget,
            "haplotype": haplotype[0],
            "position":  breakpointPos,
            "side":      0,
            "cn":        'temp',
            "cnID":      'temp',
            # connections:
            "type":        'nonTel',
            "WT":          'none',
            "M":           'none',
            # properties:
            "centromeric": False,
            "inv":         False,
        }
        nodes.append(nodeData)
        uniqueID += 1

        nodeData = {
            # identification:
            "nodeID":    uniqueID,
            "chromID":   chromosomeTarget,
            "haplotype": haplotype[0],
            "position":  breakpointPos,
            "side":      1,
            "cn":        'temp',
            "cnID":      'temp',
            # connections:
            "type":        'nonTel',
            "WT":          'none',
            "M":           'none',
            # properties:
            "centromeric": False,
            "inv":         False,
        }
        nodes.append(nodeData)
        uniqueID += 1

    print("\nDetermining cn information\n")

    # initialise cn information
    if firstEvent == True:
        if len(nodes) == 0:
            print("nodes empty - why?!")
            sys.exit()

        for i in range(len(nodes)):
            nodes[i]["cnID"] = 1
            nodes[i]["cn"]   = 1

    # no change to system
    elif firstEvent == False and nDSB == 0:
        pass

    # non initial event
    elif firstEvent == False and nDSB > 0:
        # assumed no cn is greater than 5 in simulation
        cnidList   = [i for i in range(1,5)]

        # assign cn info for newly appended junctions
        for i in range(nThreshold, len(nodes), 2):

            # list of available cnIDs
            cnidChoice = []

            # if no breaks on a chromosome then it has cn = 1
            if nbp_per_chrom[nodes[i].get("chromID")-1][nodes[i].get("haplotype")] == 0:
                nodes[i]["cnID"]   = 1
                nodes[i]["cn"]     = 1
                nodes[i+1]["cnID"] = 1
                nodes[i+1]["cn"]   = 1

            # try different cnIDs, pick randomly between the ones that exist
            else:
                #print("number of prev junctions: %s" %nbp_per_chrom[nodes[i].get("chromID")-1][nodes[i].get("haplotype")])
                for j in cnidList:

                    # set index coords every cnID attempt
                    nodeID = i
                    nodes[nodeID]["cnID"]   = j
                    nodes[nodeID+1]["cnID"] = j

                    # those that exist will provide non np.pi return in >= 1 direction
                    for k in range(2):

                        adjID = findAdjacentJunction(nodes,nodeID)

                        if adjID == np.pi:
                            nodeID = i+1
                        else:
                            #print("node %s exists" %AdjID)
                            cnidChoice.append( (j, nodes[adjID].get("cn"), nodes[nodeID].get("side"), nodes[adjID].get("side")) )

                            # debugging
                            print("\nCurrent node:")
                            print(nodes[nodeID].get("nodeID"), nodes[nodeID].get("position"), nodes[nodeID].get("type"), nodes[nodeID].get("side"), nodes[nodeID].get("chromID"), nodes[nodeID].get("cnID"), nodes[nodeID].get("cn"), nodes[nodeID].get("haplotype"))
                            print("Adjacent node:")
                            print(nodes[adjID].get("nodeID"), nodes[adjID].get("position"), nodes[adjID].get("type"), nodes[adjID].get("side"), nodes[adjID].get("chromID"), nodes[adjID].get("cnID"), nodes[adjID].get("cn"), nodes[adjID].get("haplotype"))

                            break # prevents repeat assignment

                # choose cnID location based on which nodes define the segment (opposite directions/sides)
                for k in range(len(cnidChoice)):
                    idxChoice = cnidChoice[k]
                    if idxChoice[2] != idxChoice[3]:
                        cnRef = idxChoice
                        break

                # assign cnid
                nodes[i]["cnID"]   = cnRef[0]
                nodes[i+1]["cnID"] = cnRef[0]

                # assign cn
                nodes[i]["cn"]   = cnRef[1]
                nodes[i+1]["cn"] = cnRef[1]

                # debugging
                print("chosen cnID = %s" %cnRef[0])
                print("chosen cn   = %s" %cnRef[1])



    return nodes


def generateTelomeres(nodes):
    print("\nEntering telomere generation\n")

    for i in range(len(nodes)):
        nodeID = nodes[i].get("nodeID")
        print("scanning: %s" %nodeID)

        direction = nodes[nodeID].get("side")
        adjID = findAdjacentJunction(nodes,nodeID)

        if adjID == np.pi and direction == 1:
            nodes[nodeID]["type"] = 'qTel'
            print("qTel: %s" %nodeID)

        elif adjID == np.pi and direction == 0:
            nodes[nodeID]["type"] = 'pTel'
            print("pTel: %s" %nodeID)

        else:
            print("nonTel: %s" %nodeID)

    return nodes


def findAdjacentJunction(nodes,nodeID):
    rightAdjNodes = []
    leftAdjNodes  = []
    temp          = []

    direction = nodes[nodeID].get("side")
    location  = nodes[nodeID].get("position")

    connection.id_   = nodeID
    connection.chrom = nodes[nodeID].get("chromID")
    connection.cnid  = nodes[nodeID].get("cnID")
    #connection.cn    = nodes[nodeID].get("cn")
    connection.haplo = nodes[nodeID].get("haplotype")

    # debugging
    #print("\nConnection information: current junction information:")
    #print(connection.id_, location, nodes[nodeID].get("type"), direction, connection.chrom, connection.cnid, connection.haplo)
    #print(nodes[nodeID])


    # locates junctions along the chromosome in LR direction
    for i in range(len(nodes)):
        if nodes[i].get("chromID") == connection.chrom                    \
            and nodes[i].get("cnID") == connection.cnid                   \
            and nodes[i].get("haplotype") == connection.haplo             \
            and direction == 1 and nodes[i].get("position") > location:
                rightAdjNodes.append(nodes[i])

        elif nodes[i].get("chromID") == connection.chrom                  \
            and nodes[i].get("cnID") == connection.cnid                   \
            and nodes[i].get("haplotype") == connection.haplo             \
            and direction == 0 and nodes[i].get("position") < location:
                leftAdjNodes.append(nodes[i])


    # sorts list of potential junctions
    if len(rightAdjNodes) > 0:
        rightAdjNodes.sort(key = lambda x:x['position'])
        #print("RADJnodes: %s" %rightAdjNodes)

        # only telomere available
        if len(rightAdjNodes) == 1:
            nodeID = rightAdjNodes[0].get("nodeID")

        # closest position chosen
        elif len(rightAdjNodes) > 1:
            for i in range(len(rightAdjNodes)):

            # right direction: next junction is left
                if rightAdjNodes[i].get("side") == 0:
                    nodeID = rightAdjNodes[i].get("nodeID")
                    break # stops at first case


    # sorts list of potential junctions
    elif len(leftAdjNodes) > 0:
        leftAdjNodes.sort(key = lambda x:x['position'], reverse = True)
        #print("lADJnodes: %s" %leftAdjNodes)

        # only telomere available
        if len(leftAdjNodes) == 1:
            nodeID = leftAdjNodes[0].get("nodeID")

        # closest position chosen
        elif len(leftAdjNodes) > 1:
            for i in range(len(leftAdjNodes)):

            # right direction: next junction is right
                if leftAdjNodes[i].get("side") == 1:
                    nodeID = leftAdjNodes[i].get("nodeID")
                    break # stops at first case


    # no adjacent junctions means segment is telomeric, return pi
    elif len(leftAdjNodes) == 0 and len(rightAdjNodes) == 0:
        nodeID = np.pi

    return nodeID


def g1(nodes, lmbda):
    print('\nEntering G1\n')

    # list available wild type connections
    lWT = []
    for i in range(len(nodes)):
        if nodes[i].get("cn") > 0 and nodes[i].get("type") == 'nonTel' : #
            lWT.append(nodes[i].get("nodeID"))


    # check lWT has elements
    if len(lWT) > 0:
        WTcondition = True
    else: WTcondition = False

    print('WT Connections')
    print("Debug Key: nodeID, position, type, side, chromID, cnID, cn, haplo")

    # generate WT connections
    cc = 0
    while WTcondition:

        nodeID = int(np.random.choice(lWT, 1))
        adjID = findAdjacentJunction(nodes,nodeID)


        # debugging
        print("Current node:")
        print(nodes[nodeID].get("nodeID"), nodes[nodeID].get("position"), nodes[nodeID].get("type"), nodes[nodeID].get("side"), nodes[nodeID].get("chromID"), nodes[nodeID].get("cnID"), nodes[nodeID].get("cn"), nodes[nodeID].get("haplotype"))
        print("Adjacent node:")
        print(nodes[adjID].get("nodeID"), nodes[adjID].get("position"), nodes[adjID].get("type"), nodes[adjID].get("side"), nodes[adjID].get("chromID"), nodes[adjID].get("cnID"), nodes[adjID].get("cn"), nodes[adjID].get("haplotype"))


        # checks not connecting to telomere or repeated junctions
        if adjID in lWT:

            # assign connections
            nodes[nodeID]["WT"] = str(adjID)
            nodes[adjID]["WT"]  = str(nodeID)

            print("start J:  %s" %nodeID)
            print("WT J:     %s\n" %adjID)

            # prevent repeat assignment
            lWT.remove(nodeID)
            lWT.remove(adjID)


        # check end condition
        if len(lWT) > 0:
            WTcondition = True
        else: WTcondition = False

        cc += 1
        if cc == 20:
            sys.exit()

    print('M Connections\n')
    # list available mutant connections
    lM = []
    for i in range(len(nodes)):
        if nodes[i].get("cn") > 0 and nodes[i].get("M") == 'none':
            lM.append(nodes[i].get("nodeID"))


    # repairs breaks until endCondition is met
    endCondition = True
    while endCondition == True:

        ## start at random junction
        startNode = int(np.random.choice(lM, 1))

        # join random junction
        joinChoice = True
        while joinChoice:

            joinNode = int(np.random.choice(lM, 1))

            if joinNode != startNode:
                lM.remove(startNode)
                lM.remove(joinNode)
                joinChoice = False

        nodes[startNode]["M"] = str(joinNode)
        nodes[joinNode]["M"]  = str(startNode)

        print("start J:  %s" %startNode)
        print("joined J: %s" %joinNode)

        endCondition = CheckBool.endGrowth(nodes, lmbda)

    return nodes


def generateCentromeres(nodes, centromerePos):

    for i in range(len(nodes)):
        nodeID = nodes[i].get("nodeID")

        # check non telomeric segments
        if nodes[nodeID].get("type") == 'nonTel' and nodes[nodeID].get("cn") > 0:
            #print(nodeID)
            connID = int(nodes[nodeID].get("WT"))
            chromosome = nodes[nodeID].get("chromID")-1 # for index

            # left direction
            if nodes[nodeID].get("side") == 0 and centromerePos[chromosome] < nodes[nodeID].get("position") \
                and centromerePos[chromosome] > nodes[connID].get("position"):

                    nodes[nodeID]["centromeric"] = True
                    nodes[connID]["centromeric"] = True

            # right direction
            elif nodes[nodeID].get("side") == 1 and centromerePos[chromosome] > nodes[nodeID].get("position") \
                and centromerePos[chromosome] < nodes[connID].get("position"):

                    nodes[nodeID]["centromeric"] = True
                    nodes[connID]["centromeric"] = True
            else: pass

        # check telomeric segments
        elif nodes[nodeID].get("type") != 'nonTel' and nodes[nodeID].get("cn") > 0:
            chromosome = nodes[nodeID].get("chromID")-1 # for index

            if nodes[nodeID].get("type") == 'pTel' \
                and centromerePos[chromosome] < nodes[nodeID].get("position"):

                    nodes[nodeID]["centromeric"] = True

            elif nodes[nodeID].get("type") == 'qTel' \
                and centromerePos[chromosome] > nodes[nodeID].get("position"):

                    nodes[nodeID]["centromeric"] = True
            else: pass

    return nodes


def connectedPathConstruction(nodes,pathList):

    # locates all p telomeric junctions and sets as start of path
    pTel = []
    for i in range(len(nodes)):
        nodeID = nodes[i].get("nodeID")
        if nodes[nodeID].get("type") == 'pTel' and nodes[nodeID].get("cn") > 0:
            pTel.append(nodeID)


    # iterates along each path, checking for adjacent connections
    for i in pTel:
        temp = []; nodeID = i
        temp.append(nodeID)

        # condition stops at another telomere
        telCondition = True
        while telCondition:

            # check next connection (Mutant)
            if nodes[nodeID].get("M") != 'none':
                nodeID = int(nodes[nodeID].get("M"))
                temp.append(nodeID)
                #print("M Connection: %s" %nodeID)
                #print(nodes[nodeID])


                # debugging
                if nodes[nodeID].get("M") != 'none':
                    print(nodes[int(nodes[nodeID].get("M"))])
                if nodes[nodeID].get("WT") != 'none':
                    print(nodes[int(nodes[nodeID].get("WT"))])


                # check next connection (Wild Type)
                if nodes[nodeID].get("WT") != 'none':
                    print(nodes[int(nodes[nodeID].get("WT"))])

                    nodeID = int(nodes[nodeID].get("WT"))
                    temp.append(nodeID)
                    #print("WT Connection: %s" %nodeID)
                    #print(nodes[nodeID])


                    # debugging
                    if nodes[nodeID].get("M") != 'none':
                        print(nodes[int(nodes[nodeID].get("M"))])
                    if nodes[nodeID].get("WT") != 'none':
                        print(nodes[int(nodes[nodeID].get("WT"))])


                # checks if segment is telomeric
                else:
                    telCondition = CheckBool.TelCheck(nodes,nodeID)

            # no mutant connection; path unconnected & irrelevant
            else:
                temp.clear()
                break


         # check for paths defined by two identical pTelomeres
        if len(temp) > 0:
            for i in range(len(pathList)):
                # if end junction == start junction of a previous path: clear
                if temp[-1] in pathList.get(str(i)):
                    temp.clear()
                    break


        # append to end of pathList
        if len(temp) > 0:
            nPaths = len(pathList)
            pathList[str(nPaths)] = temp
            print("Connected path: %s\n" %temp)

    return pathList


def unconnectedPathConstruction(nodes):
    unconnected  = []
    telomeric    = []
    nonTelomeric = []

    # locates all unconnected junctions
    for i in range(len(nodes)):
        nodeID = nodes[i].get("nodeID")
        if nodes[nodeID].get("M") == 'none' and nodes[nodeID].get("cn") > 0:
            unconnected.append(nodeID)
    print("Unconnected junctions: %s" %unconnected)


    # identifies replication type for the unconnected segments
    for i in range(len(unconnected)):
        temp   = []
        nodeID = unconnected[i]
        temp.append(nodeID)

        # follows paths starting from unconnected junctions
        endPathCondition = True
        while endPathCondition:

            # checks next connection (WT)
            if nodes[nodeID].get("WT") != 'none':
                nodeID = int(nodes[nodeID].get("WT"))
                temp.append(nodeID)
                print("WT Connection: %s" %nodeID)
                print(nodes[nodeID])

                # checks next connection (M)
                if nodes[nodeID].get("M") != 'none':
                    nodeID = int(nodes[nodeID].get("M"))
                    temp.append(nodeID)
                    print("M Connection: %s" %nodeID)
                    print(nodes[nodeID])

                # if M == 'none': path ends on non-telomere
                else: break

            # if WT == 'none': path ends on telomere
            else:
                telCondition = CheckBool.TelCheck(nodes,nodeID)
                if telCondition == False:
                    break
                else:
                    print(nodes[nodeID])
                    sys.exit("Error in unconnected path construction")


        # unconnected path ends on telomere:
        if nodes[nodeID].get("type") == 'pTel' or nodes[nodeID].get("type") == 'qTel':
            telomeric.append(temp)

        # unconnected path ends on non-telomere:
        elif nodes[nodeID].get("type") == 'nonTel':

            # check for repeats
            for i in range(len(nonTelomeric)):
                if temp[-1] in nonTelomeric[i]:
                    temp = []
                    break

            if len(temp) > 0:
                nonTelomeric.append(temp)

    print("unconnected telomeric segments:     %s" %telomeric)
    print("unconnected non-telomeric segments: %s" %nonTelomeric)
    return telomeric, nonTelomeric


def syn_g2(nodes, pathList, telomeric, nonTelomeric):
    print('\nEntering S, G2\n')

    # lists for copied junctions
    telomericCopied    = [[] for i in range(len(telomeric))]
    nonTelomericCopied = [[] for i in range(len(nonTelomeric))]

    # unconnected paths
    uniqueID = len(nodes)
    for i in range(len(telomeric)):
        print("\nNew path, key i: %s" %i)

        # introduces new nodes to complete path
        for j in range(len(telomeric[i])):
            nodeID = telomeric[i][j]

            nodes[nodeID]["cn"] = nodes[nodeID].get("cn") + 1
            nodes.append({
                "nodeID":    uniqueID,
                "chromID":   nodes[nodeID].get("chromID"),
                "haplotype": nodes[nodeID].get("haplotype"),
                "position":  nodes[nodeID].get("position"),
                "side":      nodes[nodeID].get("side"),
                "cn":        nodes[nodeID].get("cn"),
                "cnID":      nodes[nodeID].get("cn"),

                "type":      nodes[nodeID].get("type"),
                "WT":        '',
                "M":         '',

                "centromeric": nodes[nodeID].get("centromeric"),
                "inv":         False,
            })
            telomericCopied[i].append(uniqueID)
            uniqueID += 1

        print("%d node(s) added (unconnected path length)" %len(telomeric[i]))

        # creates path
        temp = []
        for j in reversed(telomeric[i]):
            temp.append(j)

        for j in telomericCopied[i]:
            temp.append(j)

        # stores path information
        if len(temp) > 0:
            nPaths = len(pathList)
            pathList[str(nPaths)] = temp
        else: pass

        # path key:
        I = len(pathList)-1
        print("telomeric: %s" %telomeric[i])
        print("PathList:  %s" %pathList.get(str(I)))


        # connection information
        if len(pathList.get(str(I))) == 2:

            nodes[ pathList.get(str(I))[0] ]["M"]  = str( pathList.get(str(I))[1] )
            nodes[ pathList.get(str(I))[1] ]["M"]  = str( pathList.get(str(I))[0] )

            nodes[ pathList.get(str(I))[1] ]["WT"] = 'none'

        elif len(pathList.get(str(I))) > 2:

            for j in range(0, len(pathList.get(str(I)))-1, 2):

                if j == 0:

                    nodes[ pathList.get(str(I))[j]   ]["M"]  = str( pathList.get(str(I))[j+1] )
                    nodes[ pathList.get(str(I))[j+1] ]["M"]  = str( pathList.get(str(I))[j]   )

                else:
                    nodes[ pathList.get(str(I))[j]   ]["M"]  = str( pathList.get(str(I))[j+1] )
                    nodes[ pathList.get(str(I))[j+1] ]["M"]  = str( pathList.get(str(I))[j]   )

                    nodes[ pathList.get(str(I))[j]   ]["WT"]  = str( pathList.get(str(I))[j-1] )
                    nodes[ pathList.get(str(I))[j-1] ]["WT"]  = str( pathList.get(str(I))[j]   )

            nodes[ pathList.get(str(I))[ len(pathList.get(str(I)))-1 ] ]["WT"]  = 'none'


    # non telomeric segments:
    uniqueID = len(nodes)
    for i in range(len(nonTelomeric)):

        # introduces new nodes to complete path
        for j in range(len(nonTelomeric[i])):
            nodeID = nonTelomeric[i][j]

            nodes[nodeID]["cn"] = nodes[nodeID].get("cn") + 1
            nodes.append({
                "nodeID":    uniqueID,
                "chromID":   nodes[nodeID].get("chromID"),
                "haplotype": nodes[nodeID].get("haplotype"),
                "position":  nodes[nodeID].get("position"),
                "side":      nodes[nodeID].get("side"),
                "cn":        nodes[nodeID].get("cn"),
                "cnID":      nodes[nodeID].get("cn"),

                "type":        nodes[nodeID].get("type"),
                "WT":          '',
                "M":           '',
                "centromeric": nodes[nodeID].get("centromeric"),
                "inv":         False,
            })
            nonTelomericCopied[i].append(uniqueID)
            uniqueID += 1

        # creates path
        temp = []
        for j in nonTelomeric[i]:
            temp.append(j)

        for j in reversed(nonTelomericCopied[i]):
            temp.append(j)


        if len(temp) > 0:
            nPaths = len(pathList)
            pathList[str(nPaths)] = temp
        else: pass

        I = len(pathList)-1
        print("non-telomeric: %s" %nonTelomeric[i])
        print("PathList:      %s" %pathList.get(str(I)))

        # connection information
        for j in range(0, len(pathList.get(str(I)))-1, 2):

            if j == 0:

                nodes[ pathList.get(str(I))[j]   ]["M"]  = str( pathList.get(str(I))[-1] )
                nodes[ pathList.get(str(I))[-1]  ]["M"]  = str( pathList.get(str(I))[j]  )

                nodes[ pathList.get(str(I))[j]   ]["WT"]  = str( pathList.get(str(I))[j+1] )
                nodes[ pathList.get(str(I))[j+1] ]["WT"]  = str( pathList.get(str(I))[j]   )


            else:
                nodes[ pathList.get(str(I))[j]   ]["M"]  = str( pathList.get(str(I))[j-1] )
                nodes[ pathList.get(str(I))[j-1] ]["M"]  = str( pathList.get(str(I))[j]   )

                nodes[ pathList.get(str(I))[j]   ]["WT"]  = str( pathList.get(str(I))[j+1] )
                nodes[ pathList.get(str(I))[j+1] ]["WT"]  = str( pathList.get(str(I))[j]   )

    return nodes, pathList


def checkInv(nodes, pathList):

    for i in range(len(pathList)):
        for j in range(len(pathList.get(str(i)))):
            nodeID = pathList.get(str(i))[j]

            if nodes[ int(nodes[nodeID].get("M")) ].get("side") == 0 and nodes[nodeID].get("cn") > 0:
                nodes[nodeID]["inv"] = True
            else: pass

    return nodes


def cmplxSegregation(nodes, pathList, i, nCent, centList, centromerePos):

    # list of introduced junctions
    newJuncList = []

    # a threshold variable for iterating along PathList
    m = 0

    # a threshold variable for iterating along centList
    n = 0

    # nCent - 1 = number of breakpoints
    print('\n\n %s \n\n' %nCent)
    for j in range(int(nCent)-1):

        # segment options for breakpoints
        options = []

        for k in range(m, len(pathList.get(str(i))), 1):
            nodeID = pathList.get(str(i))[k]

            # appends first centromeric node
            if nodeID == centList[n] and len(options) == 0:
                options.append(nodeID)

                # debugging:
                print("\ncentList: %s\n " %centList)
                print("target 1st centromere %s" %centList[n])
                print("target 2nd centromere %s" %centList[n+1])

            # appends nodes between centromeres
            elif nodeID != centList[n+1] and len(options) != 0:
                options.append(nodeID)

            # appends second centromeric node, ending options
            elif nodeID == centList[n+1] and len(options) != 0:
                options.append(nodeID)
                # m becomes start point for next breakpoint
                m = k+1
                n = n+2
                break


        # choosing segment to insert breakpoint
        nodeChoice = int(np.random.choice(options,1))
        print("\noptions: %s" %options)
        print("decision: %s" %nodeChoice)


        # positional information
        jPos    = nodes[nodeChoice].get("position")
        jChrom  = nodes[nodeChoice].get("chromID")
        segType = nodes[nodeChoice].get("type")
        centPos = int(centromerePos[jChrom-1])

        print("key: nodeID, junction position, chromosome, type, centromere position")
        print("%s %s %s %s %s" %(nodeChoice, jPos, jChrom, segType, centPos))


        # choosing breakpoint position given the segment type
        if segType == 'pTel':
            pos = random.randint(centPos, jPos)


        elif segType == 'qTel':
            pos = random.randint(jPos, centPos)


        elif segType == 'nonTel' and nodes[nodeChoice].get("centromeric") == True:

            if jPos > centPos:
                pos = random.randint(centPos, jPos)

            elif jPos < centPos:
                pos = random.randint(jPos, centPos)


        elif segType == 'nonTel' and nodes[nodeChoice].get("centromeric") == False:
            connID = int(nodes[nodeChoice].get("WT"))

            if jPos > nodes[connID].get("position"):
                pos = random.randint(nodes[connID].get("position"), jPos)

            elif jPos < nodes[connID].get("position"):
                pos = random.randint(jPos, nodes[connID].get("position"))


        # appending new breakpoints to nodes
        uniqueID = len(nodes)
        nodes.append({
            "nodeID":    uniqueID,
            "chromID":   jChrom,
            "haplotype": nodes[nodeChoice].get("haplotype"),
            "position":  pos,
            "side":      0,
            "cn":        nodes[nodeChoice].get("cn"),
            "cnID":      nodes[nodeChoice].get("cnID"),

            "type":        'nonTel',
            "WT":          'none',
            "M":           'none',
            "centromeric": False,
            "inv":         nodes[nodeChoice].get("inv"),
        })
        nodes.append({
            "nodeID":    uniqueID+1,
            "chromID":   jChrom,
            "haplotype": nodes[nodeChoice].get("haplotype"),
            "position":  pos,
            "side":      1,
            "cn":        nodes[nodeChoice].get("cn"),
            "cnID":      nodes[nodeChoice].get("cnID"),

            "type":        'nonTel',
            "WT":          'none',
            "M":           'none',
            "centromeric": False,
            "inv":         nodes[nodeChoice].get("inv"),
        })
        newJuncList.append(uniqueID)
        newJuncList.append(uniqueID+1)

        # connect new junctions to path
        if segType == 'pTel':
            # segment becomes non telomeric as pos < jPos
            nodes[nodeChoice]["type"] = 'nonTel'
            nodes[nodeChoice]["WT"]   = str(uniqueID+1)

            # new WT connection
            nodes[uniqueID+1]["type"] = 'nonTel'
            nodes[uniqueID+1]["WT"]   = str(nodeChoice)

            # (new) telomeric segment unconnected
            nodes[uniqueID]["type"]        = 'pTel'
            nodes[uniqueID]["WT"]          = 'none'
            nodes[uniqueID]["M"]           = 'none'
            nodes[uniqueID]["centromeric"] = True

        elif segType == 'qTel':
            # segment becomes non telomeric as pos > jPos
            nodes[nodeChoice]["type"] = 'nonTel'
            nodes[nodeChoice]["WT"]   = str(uniqueID)

            # new WT connection
            nodes[uniqueID]["type"]   = 'nonTel'
            nodes[uniqueID]["WT"]     = str(nodeChoice)

            # (new) telomeric segment unconnected
            nodes[uniqueID+1]["type"]        = 'qTel'
            nodes[uniqueID+1]["WT"]          = 'none'
            nodes[uniqueID+1]["M"]           = 'none'
            nodes[uniqueID+1]["centromeric"] = True

        elif segType == 'nonTel':
            connID = int(nodes[nodeChoice].get("WT"))

            if jPos > nodes[connID].get("position"):
                nodes[connID]["WT"]     = str(uniqueID)
                nodes[uniqueID]["WT"]   = str(connID)

                nodes[uniqueID+1]["WT"] = str(nodeChoice)
                nodes[nodeChoice]["WT"] = str(uniqueID+1)

            elif jPos < nodes[connID].get("position"):
                nodes[connID]["WT"]     = str(uniqueID+1)
                nodes[uniqueID+1]["WT"] = str(connID)

                nodes[uniqueID]["WT"]   = str(nodeChoice)
                nodes[nodeChoice]["WT"] = str(uniqueID)

            # leave centromere assignment for next cell cycle

    print("\nEntering subpath construction\npath: %s\n" %pathList.get(str(i)))
    ## end of introducing breakpoints
    # define subpaths:
    subPaths = []
    # populating list of list from the newly introduced junctions
    # order does not matter as this is only for changing cn
    for ele in newJuncList:
        temp   = []
        temp.append(ele)
        nodeID = ele
        print("newJuncStart: %s" %ele)

        subPathCondition = True
        while subPathCondition:

            # check next connection (only WT)
            if nodes[nodeID].get("WT") != 'none':
                nodeID = int(nodes[nodeID].get("WT"))
                temp.append(nodeID)

                # debugging
                print("Next connection (WT): %s" %nodeID)

                # check next connection (M)
                if nodes[nodeID].get("M") != 'none':
                    nodeID = int(nodes[nodeID].get("M"))
                    temp.append(nodeID)

                    # debugging
                    print("Next connection (M): %s" %nodeID)

                # if no mutant connection, reached other junction
                else:
                    subPathCondition = False

            # no wild type; junction is telomeric
            else:
                # check
                telCondition = CheckBool.TelCheck(nodes,nodeID)
                if telCondition == True:
                    print("error: junction not telomeric but does not have WT")
                    sys.exit()

                subPathCondition = False
        print("temp: %s" %temp)
        subPaths.append(temp)

    print('Subpaths: ')
    print(subPaths)

    ## decide mitotic assignment and iterate along subpaths to delete
    if nCent == 2:

        # randomly chooses which subpath is being assigned to 0 (deleted)
        subChoice = int(random.randint(0,1))

        # deletes chosen subpath
        for q in range( len( subPaths[subChoice] ) ):
            nodeID = subPaths[subChoice][q]
            nodes[nodeID]["cn"] = 0

    elif nCent > 2:

        for p in range(len(subPaths)):
            daughterCell = int(random.randint(0,1))

            for q in range( len( subPaths[p] ) ):

                if daughterCell == 0:
                    nodeID = subPaths[p][q]
                    nodes[nodeID]["cn"] = 0

    return nodes


def mitosis(nodes, pathList, cycleID, delta, centromerePos):
    print('\nEntering M\n')

    for i in range(len(pathList)):
        daughterCell = int(random.randint(0,1))

        nCent, centList = CheckBool.checkCentromere(nodes, pathList, i)

        # stochastic assignment
        if nCent == 0:

            for j in range(len(pathList.get(str(i)))):
                nodeID  = pathList.get(str(i))[j]

                if daughterCell == 0:
                    nodes[nodeID]["cn"] = 0

        # balanced assignment
        elif nCent == 1:
            pass

        # breakage
        elif nCent > 1:
            nodes = cmplxSegregation(nodes, pathList, i, nCent, centList, centromerePos)
            
    return nodes


def analysis(nodes, cycleID, dest, mu, lmbda, delta, chromLengths):

    if cycleID == 0:

        with open('../output/0' +  str(dest) + '/parameters.tsv', 'w', newline='') as file:
            writer = csv.writer(file, delimiter = '\t')
            writer.writerow(["mu", "lmbda", "delta"])
            writer.writerow([mu, lmbda, delta])

        with open('../output/0' +  str(dest) + '/sv_data.tsv', 'w', newline='') as file:
            writer = csv.writer(file, delimiter = '\t')
            writer.writerow(["chr1", "coord1", "strand1", "chr2", "coord2", "strand2",	"extra", "cycleNum"])

        with open('../output/0' +  str(dest) + '/cn_data.tsv', 'w', newline='') as file:
            writer = csv.writer(file, delimiter = '\t')
            writer.writerow(["chr", "start", "end", "cn", "haplotype", "cycleNum"])

    # analyse data for SVs
    calcSVs.deletions(nodes,cycleID,dest)
    calcSVs.insertions(nodes,cycleID,dest)
    calcSVs.inversions(nodes,cycleID,dest,chromLengths)
    calcSVs.duplications(nodes,cycleID,dest,chromLengths)

    return


# generate SV:
def main():

    # variables
    nChroms  = 22
    dest     = 0
    nodes    = []
    pathList = {}

    # import chromosome lengths
    chromLengths = pd.read_csv("../input/hg38.size.tsv", header=None, sep='\t')

    # define centromere positions
    centromerePos = []
    for length in range(len(chromLengths)):
        centromerePos.append( chromLengths[1][length] / 2 )

    # parameters
    mu    = 10   # DSB
    lmbda = 5    # max number of unrepaired segments a cell can handle
    delta = 2    # nCycles

    # cell cycles
    cycleID = 0
    count   = 0

    for i in range(delta):
        print("\n############################ ")

        # generate breakpoints
        nDSB = generateDSBs(mu)
        print("Cycle: %d; Number of DSBs: %d" %(cycleID, nDSB))


        if nDSB > 0:
            count += 1

            # needed for junction cn initialisation
            if count == 1:
                firstEvent = True
            else: firstEvent = False


            # initialises nodes dictionary
            nodes = generateNodes(nodes, nDSB, nChroms, chromLengths, firstEvent)
            print("Total number of junctions in nucleus: %d\n" %len(nodes))


            # define path ends
            nodes = generateTelomeres(nodes)


            # growth phase
            nodes = g1(nodes, lmbda)


            # assign centromeres to segments
            nodes = generateCentromeres(nodes, centromerePos)


            # path construction
            pathList = connectedPathConstruction(nodes, pathList)
            telomeric, nonTelomeric = unconnectedPathConstruction(nodes)


            # synthesis & second growth phase
            nodes, pathList = syn_g2(nodes, pathList, telomeric, nonTelomeric)


            # path directionality
            nodes = checkInv(nodes, pathList)


            # mitosis phase
            nodes = mitosis(nodes, pathList, cycleID, delta, centromerePos)


            # open output file for writing SV data
            analysis(nodes, cycleID, dest, mu, lmbda, delta, chromLengths)


        else: print("Exception, no available junctions.")

        cycleID += 1
        pathList.clear()

    nodes.clear()
    return


if __name__ == '__main__':
    print(" Running SVGen.\n")
    main()
