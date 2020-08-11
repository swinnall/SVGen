" Evolutionary Model of Structural Variation (EMoSV); Author: Samuel Winnall "

################################# comments ###################################

" analysis - insertions: currently not accounting for BFB unconnected junctions "

" analysis - duplications: same problem as above "

" analysis - inversions: may be counted twice, not a big problem"

################################# imports ####################################

import numpy as np
import pandas as pd
import random
import csv

################################# classes ####################################


class ConnectionsClass():
    chrom = 0
    side  = 0
    id_   = 0
    haplo = 0
    cnid  = 0

    def __init__(self, chrom, side, id_, haplo, cnid):
        self.chrom
        self.side
        self.id_
        self.haplo
        self.cnid

connection = ConnectionsClass(0,0,0,0,0)


class cirosPlot():
    chr1     = 0
    coord1   = 0
    strand1  = '+'
    chr2     = 0
    coord2   = 0
    strand2  = '+'
    extra    = 0
    hap      = 0
    cycleID  = 0

    def __init__(self, chr1, coord1, strand1, chr2, coord2, strand2, extra, hap, cycleID):
        self.chr1
        self.coord1
        self.strand1
        self.chr2
        self.coord2
        self.strand2
        self.extra
        self.hap
        self.cycleID

    def deletions(self, nodes, count, dest):
        for i in range(len(nodes)):
            if nodes[i].get("cn") == 0:
                chr1     = nodes[i].get("chromID")
                coord1   = nodes[i].get("position")
                strand1  = '+'
                chr2     = chr1
                coord2   = coord1
                strand2  = '-'
                extra    = 'svtype=DEL'
                cycleID  = count

                # write to csv
                with open('../output/0' +  str(dest) + '/sv_data.tsv', 'a', newline='') as file:
                    writer = csv.writer(file, delimiter = '\t')
                    writer.writerow([chr1, coord1, strand1, chr2, coord2, strand2,	extra, cycleID])
        return

    def insertions(self, nodes, count, dest):
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
                cycleID  = count

                # write to csv
                with open('../output/0' +  str(dest) + '/sv_data.tsv', 'a', newline='') as file:
                    writer = csv.writer(file, delimiter = '\t')
                    writer.writerow([chr1, coord1, strand1, chr2, coord2, strand2,	extra, cycleID])

            else: pass

    def inversions(self, nodes, count, dest, chromLengths):
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
                    cycleID  = count

                elif nodes[i].get("type") == 'pTel':

                    chr1     = nodes[i].get("chromID")
                    coord1   = nodes[i].get("position")
                    strand1  = '+'
                    chr2     = chr1
                    coord2   = 0
                    strand2  = '+'
                    extra    = 'svtype=INV'
                    cycleID  = count

                elif nodes[i].get("type") == 'qTel':

                    chr1     = nodes[i].get("chromID")
                    coord1   = nodes[i].get("position")
                    strand1  = '+'
                    chr2     = chr1
                    coord2   = chromLengths[1][chr1-1]
                    strand2  = '+'
                    extra    = 'svtype=INV'
                    cycleID  = count

                # write to csv
                with open('../output/0' +  str(dest) + '/sv_data.tsv', 'a', newline='') as file:
                    writer = csv.writer(file, delimiter = '\t')
                    writer.writerow([chr1, coord1, strand1, chr2, coord2, strand2,	extra, cycleID])
        return

    def duplications(self, nodes, count, dest, chromLengths):
        coveredNodes = []
        for i in range(len(nodes)):
            nodeID = nodes[i].get("nodeID")
            coveredNodes.append(nodeID)

            if nodes[i].get("type") == 'nonTel' and nodes[i].get("M") != 'none':
                #print(nodes[i])
                AdjacentID = nodes[ int(nodes[i].get("WT")) ].get("nodeID")

                if AdjacentID not in coveredNodes:

                    chr1    = nodes[i].get("chromID")
                    start   = nodes[i].get("position")
                    end     = nodes[AdjacentID].get("position")
                    cn      = nodes[i].get("cn")
                    hap     = nodes[i].get("haplotype")
                    cycleID = count

                    # write to csv
                    with open('../output/0' +  str(dest) + '/cn_data.tsv', 'a', newline='') as file:
                        writer = csv.writer(file, delimiter = '\t')
                        writer.writerow([chr1, start, end, cn, hap, cycleID])
                else: pass


            elif nodes[i].get("type") == 'pTel':

                chr1    = nodes[i].get("chromID")
                start   = 0
                end     = nodes[i].get("position")
                cn      = nodes[i].get("cn")
                hap     = nodes[i].get("haplotype")
                cycleID = count

                # write to csv
                with open('../output/0' +  str(dest) + '/cn_data.tsv', 'a', newline='') as file:
                    writer = csv.writer(file, delimiter = '\t')
                    writer.writerow([chr1, start, end, cn, hap, cycleID])


            elif nodes[i].get("type") == 'qTel':

                chr1    = nodes[i].get("chromID")
                start   = nodes[i].get("position")
                end     = chromLengths[1][chr1-1]
                cn      = nodes[i].get("cn")
                hap     = nodes[i].get("haplotype")
                cycleID = count

                # write to csv
                with open('../output/0' +  str(dest) + '/cn_data.tsv', 'a', newline='') as file:
                    writer = csv.writer(file, delimiter = '\t')
                    writer.writerow([chr1, start, end, cn, hap, cycleID])

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


def initialiseNodes(nodes,nDSB):
    uniqueID = len(nodes)
    batchIDX = len(nodes)

    # twice the number of breakpoints added (LR junctions)
    for i in range(2*nDSB):
        nodeData = {
            "nodeID":   uniqueID,
            "chromID":   0,
            "haplotype": 0,
            "position":  0,
            "side":      0,
            "cn":        1,
            "cnID":      1,

            "type":        'nonTel',
            "WT":          'none',
            "M":           'none',
            "startPath":   False,
            "centromeric": False,
            "inv":         False,
        }
        nodes.append(nodeData)
        uniqueID += 1

    # resets prev path information
    for i in range(len(nodes)):
        if nodes[i].get("cn") > 0:
            nodes[i]["type"] = 'nonTel'
            nodes[i]["WT"]   = 'none'
            nodes[i]["startPath"] = False
            nodes[i]["centromeric"] = False
            nodes[i]["inv"] = False

    return nodes, batchIDX


def breakpointPos(nodes,batchIDX,nChroms,chromLengths):

    hapChoice = [0,1]
    for i in range(batchIDX,len(nodes),2):

        chromosomeTarget = random.randint(1, nChroms)
        breakpointPos    = random.randint(0, chromLengths[1][chromosomeTarget-1])
        haplotype        = np.random.choice(hapChoice, 1)
        cnID             = random.randint(1, nodes[i].get("cn"))

        nodes[i]["chromID"]   = chromosomeTarget
        nodes[i]["position"]  = breakpointPos
        nodes[i]["side"]      = 0
        nodes[i]["haplotype"] = haplotype[0]
        nodes[i]["cnID"]      = cnID

        nodes[i+1]["chromID"]   = chromosomeTarget
        nodes[i+1]["position"]  = breakpointPos
        nodes[i+1]["side"]      = 1
        nodes[i+1]["haplotype"] = haplotype[0]
        nodes[i+1]["cnID"]      = cnID

    return


def telomericJunctions(nodes):
    for i in range(0,len(nodes),1):
        nodeID = nodes[i].get("nodeID")

        direction = nodes[nodeID].get("side")
        AdjacentID = findAdjacentJunction(nodes,nodeID)

        if AdjacentID == np.pi and direction == 1:
            nodes[nodeID]["type"] = 'qTel'


        elif AdjacentID == np.pi and direction == 0:
            nodes[nodeID]["type"] = 'pTel'

        else: pass

    return


def findAdjacentJunction(nodes,nodeID):
    rightAdjNodes = []
    leftAdjNodes  = []
    temp = []

    connection.id_ = nodeID
    direction = nodes[nodeID].get("side")
    location  = nodes[nodeID].get("position")

    connection.chrom = nodes[nodeID].get("chromID")
    connection.cnid  = nodes[nodeID].get("cnID")
    connection.haplo = nodes[nodeID].get("haplotype")

    for i in range(len(nodes)):
        if nodes[i].get("chromID") == connection.chrom                    \
            and nodes[i].get("cn") > 0                                    \
            and nodes[i].get("cnID") == connection.cnid                   \
            and nodes[i].get("haplotype") == connection.haplo             \
            and direction == 1 and nodes[i].get("position") > location:
                rightAdjNodes.append(nodes[i].get("nodeID"))

        elif nodes[i].get("chromID") == connection.chrom                  \
            and nodes[i].get("cn") > 0                                    \
            and nodes[i].get("cnID") == connection.cnid                   \
            and nodes[i].get("haplotype") == connection.haplo             \
            and direction == 0 and nodes[i].get("position") < location:
                leftAdjNodes.append(nodes[i].get("nodeID"))

        else: pass


    if len(rightAdjNodes) > 0:
        rightAdjNodes.sort()

        if len(rightAdjNodes) == 1:
            nodeID = rightAdjNodes[0]

        elif len(rightAdjNodes) > 1:
            temp.append( (rightAdjNodes[0], rightAdjNodes[1]) )

            # right wild type connection: next junction is left
            for i in range(len(temp)):
                idx = temp[0][i]
                if nodes[idx].get("side") == 0:
                    nodeID = idx
                else: pass


    elif len(leftAdjNodes) > 0:
        leftAdjNodes.sort(reverse = True)

        if len(leftAdjNodes) == 1:
            nodeID = leftAdjNodes[0]

        elif len(leftAdjNodes) > 1:
            temp.append( (leftAdjNodes[0], leftAdjNodes[1]) )

            # left wild type connection: next junction is right
            for i in range(len(temp)):
                idx = temp[0][i]
                if nodes[idx].get("side") == 1:
                    nodeID = idx
                else: pass


    # telomeric junctions return pi
    elif len(leftAdjNodes) == 0 and len(rightAdjNodes) == 0:
        nodeID = np.pi

    rightAdjNodes.clear()
    leftAdjNodes.clear()
    temp.clear()
    return nodeID


def generateWT(nodes):

    for i in range(len(nodes)):
        nodeID = nodes[i].get("nodeID")

        if nodes[nodeID].get("type") == "nonTel" and nodes[nodeID].get("cn") > 0:
            adjID = findAdjacentJunction(nodes,nodeID)

            nodes[nodeID]["WT"] = str(adjID)
            nodes[adjID]["WT"]  = str(nodeID)

        elif nodes[nodeID].get("type") != 'nonTel' and nodes[nodeID].get("cn") > 0:
            nodes[nodeID]["WT"] = 'none'

    return


def g1(nodes, lmbda):

    l = []
    for i in range(len(nodes)):
        if nodes[i].get("cn") > 0 and nodes[i].get("M") == 'none':
            l.append(nodes[i].get("nodeID"))
        else:
            pass


    endCondition = True
    while endCondition == True:

        ## Start at Random Breakpoint
        startChoice = True
        while startChoice:
            startNode = int(np.random.choice(l, 1))

            if nodes[startNode].get("M") == 'none' and nodes[startNode].get("cn") > 0:
                startChoice = False
            else: pass


        joinChoice = True
        while joinChoice:
            joinNode = int(np.random.choice(l, 1))

            if nodes[joinNode].get("M") == 'none' and joinNode != startNode and nodes[joinNode].get("cn") > 0:
                joinChoice = False
            else: pass


        nodes[startNode]["M"] = str(joinNode)
        nodes[joinNode]["M"]  = str(startNode)

        print("start J: %s" %startNode)
        print("joined J: %s" %joinNode)


        endCondition = CheckBool.endGrowth(nodes, lmbda)

    return


def connectedPathConstruction(nodes,pathList):
    pTel = []

    # locates all left telomeric junctions
    for i in range(len(nodes)):
        nodeID = nodes[i].get("nodeID")

        if nodes[nodeID].get("type") == 'pTel' and nodes[nodeID].get("cn") > 0:
            nodes[nodeID]["startPath"] = True
            pTel.append(nodeID)
        else: pass


    for i in pTel:
        temp = []

        nodeID = i
        temp.append(nodeID)

        # condition stops at another telomere
        telCondition = True
        while telCondition:

            # if M = 'none': path unconnected & irrelevant
            if nodes[nodeID].get("M") != 'none':

                nodeID = int(nodes[nodeID].get("M"))
                temp.append(nodeID)
                print("M Connection: %s" %nodeID)

                # if WT == 'none': path ends as it is telomeric
                if nodes[nodeID].get("WT") == 'none':

                    telCondition = CheckBool.TelCheck(nodes,nodeID)
                    # print("telCondition: %s" %telCondition)

                elif nodes[nodeID].get("WT") != 'none':
                    nodeID = int(nodes[nodeID].get("WT"))
                    temp.append(nodeID)
                    print("WT Connection: %s" %nodeID)
                    # continue if path continues, else: break
                    if nodes[nodeID].get("M") != 'none':
                        pass

                        # check for cycles - might hide erroneous loops...
                        if temp[-1] == temp[0]:
                            break
                        else: pass

                    else:
                        temp.clear()
                        break

            else:
                temp.clear()
                break


         # check for repeats
        if len(temp) > 0:
            for i in range(len(pathList)):
                if temp[-1] in pathList.get(str(i)):
                    nodes[temp[0]]["startPath"] = False
                    temp.clear()
                    break
                else: pass
        else: pass


        # append to pathList
        if len(temp) > 0:
            nPaths = len(pathList)
            pathList[str(nPaths)] = temp
            print("\nConnected path: %s" %temp)
        else: pass


    pTel.clear()
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
        else: pass

    print("Unconnected junctions: %s" %unconnected)

    # identifies replication type for the unconnected segments
    for i in range(len(unconnected)):
        temp = []
        nodeID = unconnected[i]
        temp.append(nodeID)

        # goes through paths starting from unconnected junctions
        endPathCondition = True
        while endPathCondition:

            if nodes[nodeID].get("WT") != 'none':
                nodeID = int(nodes[nodeID].get("WT"))
                print("WT Connection: %s" %nodeID)
                temp.append(nodeID)

                if nodes[nodeID].get("M") != 'none':
                    nodeID = int(nodes[nodeID].get("M"))
                    temp.append(nodeID)
                else: break

            else:
                telCondition = CheckBool.TelCheck(nodes,nodeID)
                if telCondition == False:
                    break
                else:
                    print(nodes[nodeID])
                    print("\nWT is empty...\n")
                    pass
                    break


        if nodes[nodeID].get("type") == 'pTel' or nodes[nodeID].get("type") == 'qTel':
            nodes[temp[0]]["startPath"] = True
            telomeric.append(temp)

        elif nodes[nodeID].get("type") == 'nonTel':
            #print(temp)

            # checks for repeats
            for i in range(len(nonTelomeric)):
                if temp[-1] in nonTelomeric[i]:
                    temp = []
                    break
                else:
                    pass

            if len(temp) > 0:
                nodes[temp[0]]["startPath"] = True
                nonTelomeric.append(temp)
            else: pass

        # temp.clear()
    print("unconnected telomeric segments:     %s" %telomeric)
    print("unconnected non-telomeric segments: %s" %nonTelomeric)
    return telomeric, nonTelomeric


def syn_g2(nodes, pathList, telomeric, nonTelomeric):
    telomericCopied    = [[] for i in range(len(telomeric))]
    nonTelomericCopied = [[] for i in range(len(nonTelomeric))]

    uniqueID = len(nodes)
    for i in range(len(telomeric)):

        print("\nNew path, key i: %s" %i)

        # introduces new nodes
        for j in range(len(telomeric[i])):
            nodeID = telomeric[i][j]

            nodes[nodeID]["cn"] = nodes[nodeID].get("cn") + 1
            nodes.append( {
                "nodeID":    uniqueID,
                "chromID":   nodes[nodeID].get("chromID"),
                "haplotype": nodes[nodeID].get("haplotype"),
                "position":  nodes[nodeID].get("position"),
                "side":      nodes[nodeID].get("side"),
                "cn":        nodes[nodeID].get("cn"),
                "cnID":      nodes[nodeID].get("cn"),#nodes[nodeID].get("cnID") + 1,

                "type":        nodes[nodeID].get("type"),
                "WT":          '',
                "M":           '',
                "startPath":   False,
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


        if len(temp) > 0:
            nPaths = len(pathList)
            pathList[str(nPaths)] = temp
            nodes[ temp[0] ]["startPath"] = True
        else: pass

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


    uniqueID = len(nodes)
    for i in range(len(nonTelomeric)):

        # introduces new nodes
        for j in range(len(nonTelomeric[i])):
            nodeID = nonTelomeric[i][j]

            nodes[nodeID]["cn"] = nodes[nodeID].get("cn") + 1
            nodes.append( {
                "nodeID":    uniqueID,
                "chromID":   nodes[nodeID].get("chromID"),
                "haplotype": nodes[nodeID].get("haplotype"),
                "position":  nodes[nodeID].get("position"),
                "side":      nodes[nodeID].get("side"),
                "cn":        nodes[nodeID].get("cn"),
                "cnID":      nodes[nodeID].get("cn"),#nodes[nodeID].get("cnID") + 1, # what if segment 1/3 is being connected... needs to be 4 but currently would be 2..

                "type":        nodes[nodeID].get("type"),
                "WT":          '',
                "M":           '',
                "startPath":   False,
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
            nodes[ temp[0] ]["startPath"] = True
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

    return


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

    return


def cmplxSegregation(nodes, pathList, i, nCent, centList):

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
                count = 0
                options = []
                for k in range(m, len(pathList.get(str(i))), 1):
                    nodeID = pathList.get(str(i))[k]

                    if nodeID == centList[j] and len(options) == 0:
                        options.append(nodeID)

                    elif k != centList[j] and count == 1:
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
                    "startPath":   False,
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
                    "startPath":   False,
                    "centromeric": False,
                    "inv":         False,
                })


            return



def mitosis(nodes, pathList, count, delta):

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

            if count < delta-1:
                cmplxSegregation(i, nCent, centList)
            elif count == delta-1:

                for j in range(len(pathList.get(str(i)))):
                    nodeID  = pathList.get(str(i))[j]

                    if daughterCell == 1:
                        pass
                    elif daughterCell == 0:
                        nodes[nodeID]["cn"] = 0

    return


def analysis(nodes, count, dest, nDSB, lmbda, chromLengths):

    if count == 0:

        with open('../output/0' +  str(dest) + '/parameters.tsv', 'w', newline='') as file:
            writer = csv.writer(file, delimiter = '\t')
            writer.writerow(["nDSB", "max unrepaired", "cycleID"])
            writer.writerow([nDSB, lmbda, count])

        with open('../output/0' +  str(dest) + '/sv_data.tsv', 'w', newline='') as file:
            writer = csv.writer(file, delimiter = '\t')
            writer.writerow(["chr1", "coord1", "strand1", "chr2", "coord2", "strand2",	"extra", "cycleID"])

        with open('../output/0' +  str(dest) + '/cn_data.tsv', 'w', newline='') as file:
            writer = csv.writer(file, delimiter = '\t')
            writer.writerow(["chr", "start", "end", "cn", "haplotype", "cycleID"])

    elif count > 0:

        with open('../output/0' +  str(dest) + '/parameters.tsv', 'a', newline='') as file:
            writer = csv.writer(file, delimiter = '\t')
            writer.writerow([nDSB, lmbda, count])


    # analyse data for SVs
    calcSVs.deletions(nodes,count,dest)
    calcSVs.insertions(nodes,count,dest)
    calcSVs.inversions(nodes,count,dest,chromLengths)
    calcSVs.duplications(nodes,count,dest,chromLengths)

    return


# generate SV:
def main():

    # variables
    nChroms = 22
    dest    = 0
    nodes   = []
    pathList = { }


    # import chromosome lengths
    chromLengths = pd.read_csv("../input/hg38.size.tsv", header=None, sep='\t')


    # define centromere positions
    centromerePos = []
    for length in range(len(chromLengths)):
        centromerePos.append( chromLengths[1][length] / 2 )


    # parameters
    mu    = 10   # DSB
    lmbda = 5    # max number of unrepaired segments a cell can handle
    delta = 1    # nCycles


    # cell cycles
    count = 0
    for i in range(delta):
        print("\n############################ ")

        # generate breakpoints
        nDSB = generateDSBs(mu)
        print("Cycle: %d; Number of DSBs: %d" %(count, nDSB))


        # initialises nodes dictionary
        nodes, batchIDX = initialiseNodes(nodes,nDSB)
        print("Total number of junctions in nucleus: %d\n" %len(nodes))


        # assign breakpoint genetic positions
        breakpointPos(nodes, batchIDX, nChroms, chromLengths)


        # define path ends
        telomericJunctions(nodes)


        if nDSB > 0:

            # growth phase
            generateWT(nodes)
            g1(nodes, lmbda)

            # path construction
            pathList = connectedPathConstruction(nodes, pathList)
            telomeric, nonTelomeric = unconnectedPathConstruction(nodes)


            # synthesis & second growth phase
            nodes, pathList = syn_g2(nodes, pathList, telomeric, nonTelomeric)


            # path directionality
            checkInv(nodes, pathList)


            # assign centromeres to segments
            generateCentromeres(nodes, centromerePos)


            # mitosis phase
            mitosis(nodes, pathList, count, delta)


            # open output file for writing SV data
            analysis(nodes, count, dest, nDSB, lmbda, chromLengths)


        else: print("Exception, no available junctions.")
        count += 1

        pathList.clear()
    nodes.clear()
    return


if __name__ == '__main__':
    print(" Running as standalone programme.\n")
    main()
