" Approximate Bayesian Computation; Author: Samuel Winnall "

" canonical chromothripsis "

################################# Overview ###################################

# 1. run breakpoint_model N times
        # edit raw script to generate two output files: parameters, summary stats


# 2. check summary stats against chromothripsis criteria
        # specify confidence interval?


# 3. if pass: keep file; else: delete


# 4. plot successful parameters (poisson distribution?)


# 5. insert new parameters into original model script (via repeating main?)

##############################################################################

import shutil
import errno
import breakpoint_model_5                
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



def readCN(delete, keep, dest):
    
    r_cn_TSV = '../output/0' + str(dest) + '/cn_data.tsv'
    cn_df = pd.read_csv(r_cn_TSV, sep="\t")
    print(cn_df)
    
    
    
    for i in range(len(cn_df)):
        
        chromID = cn_df.iat[i,0]
        cnChanges[chromID-1].append(  (cn_df.iat[i,1], cn_df.iat[i,2], cn_df.iat[i,3])  )
            # -1 for index
    
    test = []
    for i in range(nChrom):
        # if 10 breakpoints on a given chromosome then simulation passes  
        if len( cnChanges[i] ) >= condition:
            test.append(True)
        else: test.append(False)
        
    if True in test: 
        keep.append(dest)
    else:
        delete.append(dest)
    
    #print(delete)
    #cnChanges.clear()
    
    return delete, keep



def deleteFails(delete):
    
    for i in range(len(delete)):
        
        fileID = delete[i]
        path   = '../output/0' + str(fileID)
        shutil.rmtree(path, ignore_errors=False, onerror=None)
    
    return


def plotParameters(keep):
    
    params = []
    for i in range(len(keep)):
        fileID = keep[i]
        
        r_filenameTSV = '../output/0' +  str(fileID) + '/parameters.tsv'
        p_df = pd.read_csv(r_filenameTSV, sep="\t")
  
        # fileID, nDSB, nRep
        params.append( (fileID, p_df.iat[0,0], p_df.iat[0,1]) )
        
    #print('\n %s' %params[:])

    x = []; y = []
    for i in range(len(params)):
        x.append( params[i][2] )
        y.append( params[i][1] )
    
    #print(x); print(y)
    plt.scatter( x, y)
    plt.show
    
    ax = plt.gca()
    ax.set_frame_on(True)  
    ax.set(xlim=(0,100), ylim=(0,100))
    plt.xlabel('nRepairAttempts')
    plt.ylabel('nDSBs')
    
    if len(x) and len(y) > 0:
        plt.savefig('../output/parameters.png', format='png')
    else: pass
    
    return 


def svGeneration(dest):
    
    nodes    = []
    pathList = {}
    
    mu    = 10
    lmbda = 5
    delta = 3
    
    breakpoint_model_5.main(mu, lmbda, delta, count, dest, nodes, pathList) #mu, lmbda, delta, nodes, pathList =
    
    print(nodes)
    print(pathList)
    
    return


################################# Main #######################################

N         = 100
count     = 0
nChrom    = 22
condition = 5
  
delete = []
keep   = []
cnChanges = [ [] for i in range(nChrom)]


src  = '../input/00'
for i in range(N):
        
    dest = '../output/0' + str(i)
    copy(src,dest)

   
for i in range(N):

    dest  = i
    svGeneration(dest)
    
    delete, keep = readCN(delete, keep, dest)
    cnChanges = [ [] for i in range(nChrom)]
    

    
deleteFails(delete)
    
plotParameters(keep)
    
