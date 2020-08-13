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
        print("\nAnalysis - Duplications Error: ")

        for i in range(len(nodes)):
            nodeID = nodes[i].get("nodeID")
            coveredNodes.append(nodeID)

            if nodes[i].get("type") == 'nonTel' and nodes[i].get("M") != 'none':

                # debugging
                print("Node causing error: %s" %nodes[i])


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

        tally = 0
        centList = []

        for j in range(len(pathList.get(str(i)))):
            nodeID = pathList.get(str(i))[j]

            if nodes[nodeID].get("centromeric") == True:
                tally += 1
                centList.append(nodeID)
            else: pass

        nCent = tally
        return nCent, centList



################################ functions ###################################


def generateDSBs(mu):

    # nDSB = int( np.random.poisson(mu, 1) )
    nDSB = int( random.randint(0,mu) )

    return nDSB


def generateNodes(nodes,nDSB,nChroms,chromLengths,eventID):
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

        nodeData = {
            # identification:
            "nodeID":    uniqueID + 1,
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

        # update for next append pair
        uniqueID += 2


    # determine nbp on each chromosome; if !0 then check cn, else assign cn = 1
    nbp_per_chrom = [ 0 for i in range(nChroms)]
    for i in range(len(nodes)):
        nbp_per_chrom[nodes[i].get("chromID")-1] = nbp_per_chrom[nodes[i].get("chromID")-1] + 1


    # check copy number location
    if eventID == 0:
        for i in range(len(nodes)):
            nodes[i]["cnID"] = 1
            nodes[i]["cn"]   = 1

    else:
        # assumed no cn is greater than 5 in simulation
        cnidList   = [i for i in range(1,5)]

        # list of available cnIDs
        cnidChoice = []

        # assign cn info for newly appended junctions
        for i in range(nThreshold, len(nodes), 2):
            nodeID = i

            # if no breaks on a chromosome then it has cn = 1
            if nbp_per_chrom[nodes[i].get("chromID")-1] == 0:
                nodes[i]["cnID"]   = 1
                nodes[i]["cn"]     = 1
                nodes[i+1]["cnID"] = 1
                nodes[i+1]["cn"]   = 1

            # try different cnIDs, pick randomly between the ones that exist
            else:
                for j in cnidList:
                    nodes[i]["cnID"] = j

                    # those that exist will provide non np.pi return in >= 1 direction  - unless new chromosome...
                    for k in range(2):

                        AdjID = findAdjacentJunction(nodes,nodeID)
                        print(AdjID)

                        if AdjID == np.pi:
                            nodeID = i+1
                        else:
                            cnidChoice.append(j)
                            break # prevents repeat assignment of j

                # choose random cnid from available segments
                nodes[i]["cnID"]   = int(np.random.choice(cnidChoice,1))
                nodes[i+1]["cnID"] = int(np.random.choice(cnidChoice,1))

                # call associated cn for segment
                nodes[i]["cn"]   = max(cnidChoice)
                nodes[i+1]["cn"] = max(cnidChoice)


    # resets prev path information
    for i in range(len(nodes)):
        if  nodes[i].get("cn") > 0:
            nodes[i]["type"]        = 'nonTel'
            nodes[i]["WT"]          = 'none'
            nodes[i]["centromeric"] = False
            nodes[i]["inv"]         = False

    return nodes


def generateTelomeres(nodes):
    for i in range(0,len(nodes),1):
        nodeID = nodes[i].get("nodeID")

        direction = nodes[nodeID].get("side")
        AdjacentID = findAdjacentJunction(nodes,nodeID)

        if AdjacentID == np.pi and direction == 1:
            nodes[nodeID]["type"] = 'qTel'

        elif AdjacentID == np.pi and direction == 0:
            nodes[nodeID]["type"] = 'pTel'

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
    print("\nConnection information: current junction information:")
    print(connection.id_, location, nodes[nodeID].get("type"), direction, connection.chrom, connection.cnid, connection.haplo)
    print(nodes[nodeID])

# and nodes[i].get("cn") > 0                                    \
# and nodes[i].get("cn") == connection.cn                       \
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
        print("RADJnodes: %s" %rightAdjNodes)

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
        print("lADJnodes: %s" %leftAdjNodes)

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

    # list available wild type connections
    lWT = []
    for i in range(len(nodes)):
        if nodes[i].get("cn") > 0 and nodes[i].get("WT") == 'none' and nodes[i].get("type") == 'nonTel' : #
            lWT.append(nodes[i].get("nodeID"))


    # check lWT has elements
    if len(lWT) > 0:
        WTcondition = True
    else: WTcondition = False

    cc = 0
    # generate WT connections
    while WTcondition:

        nodeID = int(np.random.choice(lWT, 1))
        adjID = findAdjacentJunction(nodes,nodeID)


        # debugging
        print("Adjacent node information: ")
        print(nodes[adjID].get("nodeID"), nodes[adjID].get("position"), nodes[adjID].get("type"), nodes[adjID].get("side"), nodes[adjID].get("chromID"), nodes[adjID].get("cnID"), nodes[adjID].get("haplotype"))


        # checks not connecting to telomere or repeated junctions
        if adjID in lWT:

            # assign connections
            nodes[nodeID]["WT"] = str(adjID)
            nodes[adjID]["WT"]  = str(nodeID)

            print("start J:  %s" %nodeID)
            print("WT J:     %s" %adjID)

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
            print("\nConnected path: %s" %temp)

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

                "type":        nodes[nodeID].get("type"),
                "WT":          '',
                "M":           '',

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

            if nodes[ int(nodes[nodeID].get("M")) ].get("side") == 0 and nodes[nodeID].get("cn") > 0: # and j != 0:
                nodes[nodeID]["inv"] = True
            else: pass

    return nodes


def cmplxSegregation(nodes, pathList, i, nCent, centList, centromerePos):

            # nCent - 1 number of breakpoints
            for j in range(nCent-1):

                daughterCell = int(random.randint(0,1))
                if daughterCell == 1:
                    pass
                else:
                    marker = 0
                    for k in range(marker, len(pathList.get(str(i))), 1):

                        nodes[ pathList.get(str(i))[j] ]["cn"] = 0

                        if k == centList[j]:
                            marker = k
                            break
                        else: pass
                    break


                m = 0
                cycleID = 0
                options = []
                for k in range(m, len(pathList.get(str(i))), 1):
                    nodeID = pathList.get(str(i))[k]

                    if nodeID == centList[j] and len(options) == 0:
                        options.append(nodeID)

                    elif k != centList[j] and cycleID == 1:
                        options.append(nodeID)

                    elif nodeID == centList[j] and len(options) != 0:
                        options.append(nodeID)
                        m = k
                        break

                    else: pass

                print("options: %s" %options)
                nodeChoice = int(np.random.choice(options,1))
                jPos    = nodes[nodeChoice].get("position")
                jChrom  = nodes[nodeChoice].get("chromID")
                segType = nodes[nodeChoice].get("type")

                centPos = int(centromerePos[jChrom-1])

                print("\n %s %s %s %s %s" %(nodeChoice, jPos, jChrom, segType, centPos))

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


                nodes.append({
                    "nodeID":    len(nodes),
                    "chromID":   nodes[nodeChoice].get("chromID"),
                    "haplotype": nodes[nodeChoice].get("haplotype"),
                    "position":  pos,
                    "side":      0,
                    "cn":        nodes[nodeChoice].get("cn"),
                    "cnID":      nodes[nodeChoice].get("cnID"),

                    "type":        'nonTel',
                    "WT":          'none',
                    "M":           'none',
                    "centromeric": False,
                    "inv":         False,
                })
                nodes.append({
                    "nodeID":    len(nodes),
                    "chromID":   nodes[nodeChoice].get("chromID"),
                    "haplotype": nodes[nodeChoice].get("haplotype"),
                    "position":  pos,
                    "side":      1,
                    "cn":        nodes[nodeChoice].get("cn"),
                    "cnID":      nodes[nodeChoice].get("cnID"),

                    "type":        'nonTel',
                    "WT":          'none',
                    "M":           'none',
                    "centromeric": False,
                    "inv":         False,
                })


            return nodes


def mitosis(nodes, pathList, cycleID, delta, centromerePos):

    for i in range(len(pathList)):
        daughterCell = int(random.randint(0,1))

        nCent, centList = CheckBool.checkCentromere(nodes, pathList, i)

        # stochastic assignment
        if nCent == 0:

            for j in range(len(pathList.get(str(i)))):
                nodeID  = pathList.get(str(i))[j]

                if daughterCell == 1:
                    pass
                elif daughterCell == 0:
                    nodes[nodeID]["cn"] = 0

        # balanced assignment
        elif nCent == 1:
            pass

        # breakage
        elif nCent > 1:

            if cycleID < delta-1:
                nodes = cmplxSegregation(nodes, pathList, i, nCent, centList, centromerePos)
            elif cycleID == delta-1:

                for j in range(len(pathList.get(str(i)))):
                    nodeID  = pathList.get(str(i))[j]

                    if daughterCell == 1:
                        pass
                    elif daughterCell == 0:
                        nodes[nodeID]["cn"] = 0

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
    eventID = 0
    for i in range(delta):
        print("\n############################ ")

        # generate breakpoints
        nDSB = generateDSBs(mu)
        print("Cycle: %d; Number of DSBs: %d" %(cycleID, nDSB))


        if nDSB > 0:

            # initialises nodes dictionary
            nodes = generateNodes(nodes, nDSB, nChroms, chromLengths, eventID)
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

            eventID += 1

        else: print("Exception, no available junctions.")

        cycleID += 1
        pathList.clear()

    print(nodes)
    nodes.clear()
    return


if __name__ == '__main__':
    print(" Running SVGen.\n")
    main()
