# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 12:42:25 2018

@author: Nicole Noelle Lama (verified by: )

DistCorrAnalaysis will calcuate distances between pairs and output
a file with pairs, correlations, and distances (uses C1*)

requirements: pdb file with polar coordinates of each residue and
correlation txt file output from RING. 
"""
#import time
#time_start = time.clock()

import numpy as np
import math
import os
import matplotlib.pyplot as plot
import scipy
import seaborn as sns
import sys
sys.path.append('D:\Documents2\Rotations\Weeks\PythonScripts') #add RNAtools dir
###############################################################################
#def findTertMatrix(ctF, contDist):
#    rnaObj = RNAtools.CT(ctF)  #create RNAtools object
#    ct = readCt(ctF)
#    terMat = np.zeros([len(ct),len(ct)]) #initialize matrix with 3* interactions
#    for i in range(len(ct)):
#        for j in range(i+1, len(ct)):
#            d = rnaObj.contactDistance(i,j)
#            if d > contDist:
#                #print(i,j)
#                terMat[i,j] = int(1)
#    return terMat
#    
#
#def readCt(ctF):
#    ct = open(ctF)
#    lines = []
#    for i in range(1):
#            ct.readline() #skip first line because it is a header
#    for line in ct:
#            lineList = line.split()
#            lines.append(lineList)
#    ct = pd.DataFrame(lines)
#    return ct
    
#def getTerts(ctF, contDist):
#    rnaObj = RNAtools.CT(ctF)  #create RNAtools object
#    terInts = [] #list containing tertiary interactions
#    ct = readCt(ctF)
#    for i in range(len(ct)):
#        for j in range(i+1, len(ct)):
#            d = rnaObj.contactDistance(i,j)
#            if d > contDist:
#                terInts.append([i,j])
#    return terInts
#            
#def filtTerts(cr,terInts):
#    if [int(cr[0]), int(cr[1])] in terInts or [cr[1], cr[0]] in terInts:
#        return cr
#    else:
#        cr = []
#        return cr

   
#####################FUNCTIONS#################################################

def corrPairFunc(pdbF, atom = 'C1*'):
    pyr = ['U', 'C']
    pur = ['A','G']
    c1List = [] #list stores full line with C1
    coorPair = {} #will store residue and its coordinates
    #atomPair = {}
    for line in pdbF:
        lineList = line.split()
        if atom == 'N1':
            if len(lineList)>5 and lineList[3] in pur:
                if lineList[2] == 'N1':
                    coorPair[lineList[5]]=lineList[6:9]
                    #('pur N3 appended to coorPair')
            elif len(lineList)>5 and lineList[3] in pyr: 
                if lineList[2] == 'N3':
                    coorPair[lineList[5]]=lineList[6:9]
                    #print('pyr N1 appended to coorPair')
        elif len(lineList)>5 and lineList[2] == atom: #only keep if C1 is 2nd eleme
            c1List.append(lineList) #if so, append entire line to c1List
            coorPair[lineList[5]]=lineList[6:9] #keep only residue and xyz
            #atomPair[lineList[5]]=lineList[3]
    return coorPair

#Look up residues in coorPair for computing distances
def resCorrDistFunc(coorPair, corrF):
    resCorrDist = [] #will store correlations and distances for residues i and j
    for i in range(2):
        corrF.readline() #skip first 2 lines because they are headers
    for line in corrF:
        coor = line.split() #split line into list
        #coor = filtTerts(coor, terInts)
        #print('returned coor')
        if len(coor) > 1: # and terMat[int(coor[0]), int(coor[1])] == 1:
            i = coor[0] #extract residue1 in pair
            j = coor[1] #extract residue2 in pair
            iCoords = coorPair[i] #find residue 1's coordinates
            jCoords = coorPair[j] #find residue 2's coordinates
            x = (float(iCoords[0]) - float(jCoords[0]))**2 #compute dist form x 
            y = (float(iCoords[1]) - float(jCoords[1]))**2 #compute dist form y
            z = (float(iCoords[2]) - float(jCoords[2]))**2 #compute dist form z
            distance = math.sqrt(x + y + z) # compute distance between residues
            resCorrDist.append([i , j , float(coor[2]), distance])
    if resCorrDist == []:
        raise ValueError("returned empty resCorrDist list")
    else:
        return resCorrDist


##Will create median!
##curveNum indicates which curve you want median to attach to. EX: if it's the
##first curve plotted on the graph, curveNum = 0   
def medianMaker(curveNum,distList,chooseColor='k', chooseLine = 'dashed'):
    x,y = distList.get_lines()[curveNum].get_data()
    cdf = scipy.integrate.cumtrapz(y, x, initial=0)    
    nearest_05 = np.abs(cdf-0.5).argmin()
    x_median = x[nearest_05]
    y_median = y[nearest_05] 
    #print(x_median, y_median)
    plot.vlines(x_median, 0, y_median, color = chooseColor, 
                linestyles= chooseLine)
      
def getDist(resCorrDist, cutoff=-10):
    distList = []
    for line in resCorrDist:
        if line[2] > cutoff:
            distList.append(line[3])
    return distList

def cutoffMaker(fileCorr):
    corrF = open(fileCorr)
    cutoffVec = []
    for i in range(2):
        corrF.readline() #skip first 2 lines because they are headers
    for line in corrF:
        coor = line.split() #split line into list
        if len(coor) > 1:
            cutoffVec.append(float(coor[2])) #extract correlation   
    q = [25,50,75,90,95,99,99.5]
    m = np.percentile(cutoffVec, q)
    return m
###################### MI + APC FILE ##########################################
#working directory
directory =  'D:\Documents2\Rotations\Weeks\Week1Data\CorrDistAnalysisRNaseP'
filePDB = '00952V_n2_sup.pdb' #pdb with xyz coordinates
#ctF = 'rnpb.ct'
f = 'rnasep.1.mi.apc.corr.txt' #RING output with residues and correlations
fileCorr = os.path.join(directory,f)
#outputFile = 'resCorrDistOutput.txt' #name of output file for distances
pdbF = open(filePDB) #open pdb file into memory
corrF = open(fileCorr) #same for RING output    
#contDist = 20
#terInts = getTerts(ctF, contDist)
#terMat = findTertMatrix(ctF, contDist)
coorPair = corrPairFunc(pdbF) #prepare corr file so it's NT and coordinate
resCorrDist = resCorrDistFunc(coorPair, corrF)#, ctF, terMat) #get coordinates and distances
m = cutoffMaker(fileCorr) #get quartile cutoffs for coors

distList25_apc = getDist(resCorrDist, cutoff = m[0]) #quartile 25
distList50_apc = getDist(resCorrDist, cutoff = m[1]) #quartile 50
distList75_apc = getDist(resCorrDist, cutoff = m[2])
distList90_apc = getDist(resCorrDist, cutoff = m[3])
distList95_apc = getDist(resCorrDist, cutoff = m[4])
distList99_apc = getDist(resCorrDist, cutoff = m[5])
distList995_apc = getDist(resCorrDist, cutoff = m[6])
binRange = range(0,160,3) #make bins all same for all plots

plot.figure(20)
plot.ylim(0, 0.105)
plot.xlim(-5,175)
sns.set(rc={"figure.figsize": (6, 6)})
sns.set_style('whitegrid')
flatui = ["#e74c3c",'#FF9933','#FACC2E', "#31B404", "#9b59b6", "#3498db",  "#34495e", "#2ecc71" ,  "#95a5a6"]
sns.set_palette(flatui)


#plot.hist(distList995_apc, alpha =1 , bins= binRange, histtype = 'step', linewidth = 2.3,
#          label = '{0:.2E} n = {1}'.format(m[6], str(len(distList995_apc))),
#          normed = True)
plot.hist(distList99_apc, alpha =1 , bins= binRange, histtype = 'step', linewidth = 2.3,
          label = '{0:.2E} n = {1}'.format(m[5], str(len(distList99_apc))),
          normed = True)
plot.hist(distList95_apc, alpha = 1,  bins= binRange, histtype = 'step', linewidth = 2.3,
          label = '{0:.2E} n = {1}'.format(m[4], str(len(distList95_apc))),
          normed = True)
plot.hist(distList75_apc, alpha = 1, bins= binRange, histtype = 'step',linewidth = 2.3,
          label = '{0:.2E} n = {1}'.format(m[2], str(len(distList75_apc))),
          normed = True)
plot.hist(distList25_apc, alpha = 1, bins= 20, histtype = 'step', linewidth = 2.3,
          label = '{0:.2E} n = {1}'.format(m[0], str(len(distList25_apc))),
          normed = True)


sns.set_palette(flatui)
plot.legend()
plot.xlabel('Distances')
plot.ylabel('Normalized Counts')
plot.title('{0}: contact distance'.format(fileCorr))
#plot.savefig('D:\Documents2\Rotations\Weeks\Week1Data\ContactDistanceAnalysis\{0}_histograms_tert{1}.png'.format(fileCorr, contDist)) #save plot to working dir

plot.figure(3)
plot.legend()
sns.set(rc={"figure.figsize": (6, 6)})
#flatui = ["#e74c3c", '#FA5858', '#FACC2E', "#31B404", "#9b59b6", "#3498db",  "#34495e", "#2ecc71" ,  "#95a5a6"]
#flatui = ["#e74c3c",'#FACC2E', "#31B404", "#9b59b6", "#3498db",  "#34495e", "#2ecc71" ,  "#95a5a6"]
sns.set_palette(flatui)
sns.set_style('whitegrid')
plot.ylim(0, 0.06)
plot.xlim(-5,160)

#a=sns.distplot(distList995_apc, bins=binRange, hist = False, kde_kws={"linewidth": 2.6, 'alpha' : 1}, 
#             label = '{0:.2E} n = {1}'.format(m[6], str(len(distList995_apc))))
a=sns.distplot(distList99_apc, bins=binRange, hist = False, kde_kws={"linewidth": 2.6, 'alpha' : 1}, 
             label = '{0:.2E} n = {1}'.format(m[5], str(len(distList99_apc))))
lo=sns.distplot(distList95_apc, bins=binRange,hist = False, kde_kws={"linewidth": 2.6, 'alpha' : 1}, 
             label = '{0:.2E} n = {1}'.format(m[4], str(len(distList95_apc))))
sns.distplot(distList75_apc, bins = binRange, hist = False, kde_kws={"linewidth": 2.6, 'alpha' : 1}, 
             label = '{0:.2E} n = {1}'.format(m[2], str(len(distList75_apc))))
sns.distplot(distList25_apc, bins=binRange, hist = False, kde_kws={"linewidth": 2.6, 'alpha' : 1}, 
             label = '{0:.4E} n = {1}'.format(m[0], str(len(distList25_apc))))

#medianMaker(0,a,chooseColor='r') #Add median line for plot
medianMaker(0,a,chooseColor='r') #Add median line for plot

plot.show()
plot.xlabel('Distances')
plot.ylabel('Normalized Counts')
#plot.title('{0} : contact distance {1}'.format(fileCorr, contDist))

#plot.savefig('D:\Documents2\Rotations\Weeks\Week1Data\ContactDistanceAnalysis\{0}_kde_tert{1}.png'.format(fileCorr, contDist))
#
#
###############################################################################
######################### PHI ONLY FILE #######################################

#directory =  'D:\Documents2\Rotations\Weeks\Week1Data\CorrDistAnalysisRNaseP'
filePDB = '00952V_n2_sup.pdb' #pdb with xyz coordinates
f = 'rnasep.1.corr.txt' #RING output with residues and correlations
fileCorr = os.path.join(directory,f)
#outputFile = 'resCorrDistOutput.txt' #name of output file for distances
pdbF = open(filePDB) #open pdb file into memory
corrF = open(fileCorr) #same for RING output    

coorPair = corrPairFunc(pdbF)
resCorrDist = resCorrDistFunc(coorPair, corrF)
m_phi = cutoffMaker(fileCorr)

distList25_phi = getDist(resCorrDist, cutoff = m_phi[0])
distList50_phi = getDist(resCorrDist, cutoff = m_phi[1])
distList75_phi = getDist(resCorrDist, cutoff = m_phi[2])
distList90_phi = getDist(resCorrDist, cutoff = m_phi[3])
distList95_phi = getDist(resCorrDist, cutoff = m_phi[4])
distList99_phi = getDist(resCorrDist, cutoff = m_phi[5])
distList995_phi = getDist(resCorrDist, cutoff = m_phi[6])

binRange = range(0,160,3)


######################### MI ONLY FILE ########################################

#directory =  'D:\Documents2\Rotations\Weeks\Week1Data\CorrDistAnalysisRNaseP'
filePDB = '00952V_n2_sup.pdb' #pdb with xyz coordinates
f = 'rnasep.1.mi.corr.txt' #RING output with residues and correlations
fileCorr = os.path.join(directory,f)
#outputFile = 'resCorrDistOutput.txt' #name of output file for distances
pdbF = open(filePDB) #open pdb file into memory
corrF = open(fileCorr) #same for RING output    

coorPair = corrPairFunc(pdbF)
resCorrDist = resCorrDistFunc(coorPair, corrF)
m_mi = cutoffMaker(fileCorr)

distList25_mi = getDist(resCorrDist, cutoff = m_mi[0])
distList50_mi = getDist(resCorrDist, cutoff = m_mi[1])
distList75_mi = getDist(resCorrDist, cutoff = m_mi[2])
distList90_mi = getDist(resCorrDist, cutoff = m_mi[3])
distList95_mi = getDist(resCorrDist, cutoff = m_mi[4])
distList99_mi = getDist(resCorrDist, cutoff = m_mi[5])
distList995_mi = getDist(resCorrDist, cutoff = m_mi[6])
binRange = range(0,160,3)


####### Figure comparing 99 or 95 percentile mi vs phi vs mi+phi #############

plot.figure(50)
plot.legend()
sns.set(rc={"figure.figsize": (8, 6)})
sns.set_style('whitegrid')
#sns.set_palette('husl')
plot.ylim(0, 0.05)
plot.xlim(-15,160)

a=sns.distplot(distList95_apc, bins=binRange, norm_hist = True,
               hist = False, kde_kws={"linewidth": 3, 'alpha' : 1}, 
             label = 'mi + apc')
b=sns.distplot(distList95_mi, bins=binRange, norm_hist = True,
               hist = False, kde_kws={"linewidth": 3, 'alpha' : 1}, 
             label = 'mi only')
c=sns.distplot(distList95_phi, bins=binRange, norm_hist = True,
               hist = False, kde_kws={"linewidth": 3, 'alpha' : 1}, 
             label = 'phi') 

medianMaker(0,a,'#FA5882')
medianMaker(1,b,'#FACC2E')
medianMaker(2,c,'g')

plot.show()
plot.xlabel('Distances')
plot.ylabel('Normalized Counts')
plot.title('KDE of PHI, MI and MI + APC (99)')
#plot.savefig('D:\Documents2\Rotations\Weeks\Week1Data\ContactDistanceAnalysis\{0}_kde_mi_apc_phi_99_tert{1}.png'.format(fileCorr, contDist))




#plot.figure(3)
#a=sns.distplot(distList99_phi, bins=binRange, hist = False, kde_kws={"linewidth": 2.6, 'alpha' : 1}, 
#             label = '{0:.2E} n = {1}'.format(m[5], str(len(distList99_apc))))
#lo=sns.distplot(distList95_phi, bins=binRange,hist = False, kde_kws={"linewidth": 2.6, 'alpha' : 1}, 
#             label = '{0:.2E} n = {1}'.format(m[4], str(len(distList95_apc))))
#sns.distplot(distList75_phi, bins = binRange, norm_hist = True,
#             hist = False, kde_kws={"linewidth": 2.6, 'alpha' : 1}, 
#             label = '{0:.2E} n = {1}'.format(m[2], str(len(distList75_apc))))
#sns.distplot(distList25_phi, bins=binRange, norm_hist = True,
#             hist = False, kde_kws={"linewidth": 2.6, 'alpha' : 1}, 
#             label = '{0:.4E} n = {1}'.format(m[0], str(len(distList25_apc))))
#
##medianMaker(0,a,chooseColor='r') #Add median line for plot
#medianMaker(0,a,chooseColor='r') #Add median line for plot
#
#plot.show()
#plot.xlabel('Distances')
#plot.ylabel('Normalized Counts')

