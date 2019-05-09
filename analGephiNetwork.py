# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 16:06:29 2018
Need CT file, csv from gephi, and node-edge list 
@author: nlama
"""
import pandas as pd
import os
import seaborn as sns
import random as rand
import numpy as np
import matplotlib.pyplot as plot
import RNAtools

condition = "incell"
rnaName = 'tmRNA'
gephiDat=pd.read_csv("D:/Weeks/Data/NetworkAnalysis/resCorrDistFiles/incell.tmRNA.3.50.csv")
nodeEdge=pd.read_csv("D:/Weeks/Data/NetworkAnalysis/resCorrDistFiles/incell4x.tmRNA.apc.win3.minCorr5.minCount50.50.nodist.resCorrDist.csv")
#ctHead = "D:/Documents2/Rotations/Weeks/RNAFiles/RNaseP/rnpb.ct"
ctHead = "D:/Weeks/Data/BasePairDetectionProject/SecondaryStructure/AmpliconData/tmRNA.ct"
ctHeaders = ['int', 'idx0', 'i_pair', 'j_pair', 'idx1']
#ct = open("D:/Documents2/Rotations/Weeks/RNAFiles/RNaseP/rnpb.noHeader.ct")#, names=ctHeaders)
ct = open("D:/Weeks/Data/BasePairDetectionProject/SecondaryStructure/AmpliconData/tmRNA_noHead.ct")

numClusters = 6
#Lists with stats about each module
numPairedNodes = np.zeros(numClusters)
sumCorrs = np.zeros(numClusters)
numNodes = np.zeros(numClusters)
numTriangles = np.zeros(numClusters)
numCorr = np.zeros(numClusters)
numDegrees = np.zeros(numClusters)
maxDegree = np.zeros(numClusters)
bp = {}
nodeCorrDict = {}
nodeEigenDict = {}
nodeDegDict = {}
oldS = nodeEdge.Source[0]
allC = 0.0
numC = 0.0
cListPerNode = []
medCorrDict = {}



for s,c in zip(nodeEdge.Source, nodeEdge.Correlation):
    newS = s
    if newS == oldS:
        allC += c
        cListPerNode.append(c)
        numC += 1
        continue
    if allC != 0 and numC != 0:
        avgC = allC/numC
        medC = np.median(cListPerNode)
        medCorrDict[s] = medC
        nodeCorrDict[s] = avgC
    oldS = newS
    allC = 0.0
    numC = 0.0

for line in ct:
    ctLine = line.split()
    i = ctLine[0]
    j = ctLine[4]
    bp[i] = j

rnaObj = RNAtools.CT(ctHead)  #create RNAtools object
totalKnown = 0 
totalKnownPredicted = 0

for position in rnaObj.ct:
    if position != 0:
        totalKnown += 1
        
rnaKnownPairs = rnaObj.pairList()
alreadyChecked = []
totalPredictedPairs = 0
nodeClustDict = {}

for x in range(numClusters):
    for y in range(len(gephiDat)):
        if gephiDat.modularity_class[y] != x:
            continue
        if gephiDat.modularity_class[y] == x:
            numDegrees[x] += int(gephiDat.Degree[y] )
            node = str(gephiDat.Id[y])
            nodeClustDict[node] = x
            nodeEigenDict[node] = gephiDat.eigencentrality[y]
            nodeDegDict[node] = gephiDat.Degree[y]
            if bp[node] != '0':
                numPairedNodes[x] += 1
            numTriangles[x] += int(gephiDat.triangles[y] )
            numCorr[x] += float(nodeCorrDict.get(int(node),0))
            numNodes[x] += 1
            

ppvSens = []           

for x in range(len(nodeEdge)):
    for y in range(3):#Accounting for way window 3 NTs are counted in RING-Map:  
            totalPredictedPairs += 1
            if nodeClustDict.get(str(nodeEdge['Source'][x]+y),50) == 1 and nodeClustDict.get(str(nodeEdge['Target'][x]+2-y),50) == 1:
                if (nodeEdge['Source'][x]+y, nodeEdge['Target'][x]+2-y) in rnaKnownPairs: #\
                    if (nodeEdge['Source'][x]+y, nodeEdge['Target'][x]+2-y) not in alreadyChecked: #\
                        #print(nodeEdge['Source'][x]+y, nodeEdge['Target'][x]+2-y)
                        totalKnownPredicted += 1
                        alreadyChecked.append((nodeEdge['Source'][x]+y, nodeEdge['Target'][x]+2-y))
                    
                        
totalKnown /= 2
#totalPredictedPairs /= 3
#        print(totalKnownPredicted, totalPredictedPairs, totalKnown)
sen = totalKnownPredicted/totalKnown
if totalPredictedPairs != 0:
    ppv = totalKnownPredicted/totalPredictedPairs
else:
    ppv = 0
print('sen', sen, 'ppv', ppv)
 

numParameters = 4 #change depending on num of interesting parameters
avgDegrees = [int(numDeg) / int(numNodes) for numDeg,numNodes in zip(numDegrees, numNodes)]
avgTriangles = [int(numTri) / int(numNodes) for numTri,numNodes in zip(numTriangles, numNodes)]
perPairedNodes = [int(numPair) / int(numNodes) for numPair,numNodes in zip(numPairedNodes, numNodes)]
avgCorr = [float(numCorr) / int(numNodes) for numCorr,numNodes in zip(numCorr, numNodes)]
clusterLab = [0,1,2,3]

#sumData = np.array(clusterLab+avgCorr+avgDegrees+avgTriangles+perPairedNodes)

#df = pd.DataFrame(sumData.reshape([numParameters,numClusters+1]).transpose(),
#                  columns=['Cluster','Avg Correlation','Avg Degree','Avg # Triangles','Percent of Paired Nodes in Cluster'])
##df = df.loc[df['Correlation'] >= self.cutoffs[3],]
#
#figNum = rand.randint(1,10000)
#sns.set(font_scale=4, style='whitegrid', 
#            rc={"figure.figsize": (20, 12),"lines.linewidth":4})
#
#plot.figure(figNum)
##plot.xlim([0,80])
#ax = sns.pointplot(x= 'Avg Degree', y= 'Avg Correlation', 
#                   data = df, hue='Cluster')
#
#figNum = rand.randint(1,10000)
#plot.figure(figNum)
##plot.xlim([0,80])
#ax = sns.pointplot(x= 'Avg Correlation', y= 'Percent of Paired Nodes in Cluster', 
#                   data = df, hue='Cluster')
#
#
#figNum = rand.randint(1,10000)
#plot.figure(figNum)
##plot.xlim([0,80])
#ax = sns.pointplot(x= 'Avg Degree', y= 'Percent of Paired Nodes in Cluster', 
#                   data = df, hue='Cluster')
if rnaName == "RNaseP":
    p18 = list(range(304,328))
    p13 = list(range(184,203))
    p14 = list(range(204,225))
    p12 = list(range(142,177))
    p9 = list(range(107,120))
    p8 = list(range(92,107))
    p6 = [82,83,84,85,279,278,277,276]
    p4 = [66, 67, 68, 69, 70, 71, 72, 73, 353, 
          354, 355, 356, 357, 358, 359, 360]
    p1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
          363, 364, 365, 366, 367, 368, 369, 370, 
          371, 372, 373, 374, 375, 376, 377]
        
    sDomain = list(range(82,187)) + list(range(222,246))
    sDomain.append(203)
    cDomain = list(range(1,76)) + list(range(247,263)) + list(range(287,306)) + list(range(324,380))


if rnaName == "tmRNA":
    tRNAlike_Domain = list(range(1,40)) + list(range(305,350))
    H2 = list(range(19,51)) + list(range(320,295))
    PK1 = list(range(52,82))
    PK4 = list(range(230,282))
    PK2 = list(range(130,182))
    PK3 = list(range(190,229))
outDir = 'D:/Weeks/Data/NetworkAnalysis/3DSuperPosition'
qVec = []
percentile = 90
clusterCut = np.zeros(numClusters)
for cluster in range(numClusters):
   for i,j,c in zip(nodeEdge['Source'],nodeEdge['Target'],nodeEdge['Correlation']):
       if nodeClustDict[str(i)] == cluster or nodeClustDict[str(j)] == cluster:
           qVec.append(c)
   clusterCut[cluster] = np.percentile(qVec,percentile)

for cluster in range(numClusters):
    with open(os.path.join(outDir,'{0}.rnaseP.{1}.{2}.txt'.format(condition,cluster,percentile)), "w+", newline = '') as outFile:
        for i,j,c in zip(nodeEdge['Source'], nodeEdge['Target'],nodeEdge['Correlation']):
            if (nodeClustDict[str(i)] == cluster and nodeClustDict[str(j)] == cluster) and c > clusterCut[cluster] and nodeDegDict[str(i)] < 15 and nodeDegDict[str(j)] < 15:
                outFile.write('{0}\t{1}\n'.format(i,j))    
#                if i in sDomain and j in cDomain:
#                    print("{0} and {1} have s to c domain interation".format(i,j))
#                if j in sDomain and i in cDomain:
#                    print("{0} and {1} have c to s domain interation".format(i,j))
#                if i not in sDomain and i not in cDomain and j in sDomain:
#                    print("{0} and {1} have s to none domain interation".format(i,j)) 
#                if i not in sDomain and i not in cDomain and j in cDomain:
#                    print("{0} and {1} have c to none domain interation".format(i,j)) 
#                if j not in sDomain and j not in cDomain and i in sDomain:
#                    print("{0} and {1} have s to none domain interation".format(i,j)) 
#                if j not in sDomain and j not in cDomain and i in cDomain:
#                    print("{0} and {1} have c to none domain interation".format(i,j)) 
    with open(os.path.join(outDir,'{0}.rnaseP.interclusterNodes.{1}.{2}.txt'.format(condition,cluster,percentile)), "w+", newline = '') as outFile:
        for i,j,c in zip(nodeEdge['Source'], nodeEdge['Target'],nodeEdge['Correlation']):
            if (nodeClustDict[str(i)] != cluster and nodeClustDict[str(j)] == cluster) and c > clusterCut[cluster]  and nodeDegDict[str(i)] < 15 and nodeDegDict[str(j)] < 15:
                outFile.write('{0}\t{1}\n'.format(i,j)) 
#                if i in sDomain and j in cDomain:
#                    print("{0} and {1} have s to c domain interation".format(i,j))
#                if j in sDomain and i in cDomain:
#                    print("{0} and {1} have c to s domain interation".format(i,j))
#                if i not in sDomain and i not in cDomain and j in sDomain:
#                    print("{0} and {1} have s to none domain interation".format(i,j)) 
#                if i not in sDomain and i not in cDomain and j in cDomain:
#                    print("{0} and {1} have c to none domain interation".format(i,j)) 
#                if j not in sDomain and j not in cDomain and i in sDomain:
#                    print("{0} and {1} have s to none domain interation".format(i,j)) 
#                if j not in sDomain and j not in cDomain and i in cDomain:
#                    print("{0} and {1} have c to none domain interation".format(i,j)) 
            if (nodeClustDict[str(i)] == cluster and nodeClustDict[str(j)] != cluster) and c > clusterCut[cluster]  and nodeDegDict[str(i)] < 15 and nodeDegDict[str(j)] < 15:
                outFile.write('{0}\t{1}\n'.format(i,j)) 
       

for cluster in range(numClusters):
    with open(os.path.join(outDir,'{0}.{1}.eig.{2}.txt'.format(condition,rnaName,cluster)), "w+", newline = '') as outFile:
        for i,j,c in zip(nodeEdge['Source'], nodeEdge['Target'],nodeEdge['Correlation']):
            if (nodeClustDict[str(i)] == cluster and nodeClustDict[str(j)] == cluster) and  c > clusterCut[cluster] and nodeEigenDict[str(i)]>0.2 and nodeDegDict[str(i)]<25 and nodeEigenDict[str(j)]>0.2 and nodeDegDict[str(j)]<25:
                outFile.write('{0}\t{1}\n'.format(i,j))    
            
    with open(os.path.join(outDir,'{0}.{1}.interclusterNodes.eig.{2}.txt'.format(condition,rnaName,cluster)), "w+", newline = '') as outFile:
        for i,j,c in zip(nodeEdge['Source'], nodeEdge['Target'],nodeEdge['Correlation']):
            if (nodeClustDict[str(i)] != cluster and nodeClustDict[str(j)] == cluster)  and  c > clusterCut[cluster] and nodeEigenDict[str(i)]>0.2 and nodeDegDict[str(i)]<25 and nodeEigenDict[str(j)]>0.2 and nodeDegDict[str(j)]<25:
                outFile.write('{0}\t{1}\n'.format(i,j)) 
            if (nodeClustDict[str(i)] == cluster and nodeClustDict[str(j)] != cluster) and  c > clusterCut[cluster] and nodeEigenDict[str(i)]>0.2 and nodeDegDict[str(i)]<25 and nodeEigenDict[str(j)]>0.2 and nodeDegDict[str(j)]<25:
                outFile.write('{0}\t{1}\n'.format(i,j)) 
                

                
      
'''      
         
for cluster in range(numClusters):
    with open(os.path.join(outDir,'{0}.{1}.{2}.{3}.txt'.format(condition,rnaName,cluster,percentile)), "w+", newline = '') as outFile:
        for i,j,c in zip(nodeEdge['Source'], nodeEdge['Target'],nodeEdge['Correlation']):
            if (nodeClustDict[str(i)] == cluster and nodeClustDict[str(j)] == cluster) and c > clusterCut[cluster]:
                outFile.write('{0}\t{1}\n'.format(i,j))    
#                if i in sDomain and j in cDomain:
#                    print("{0} and {1} have s to c domain interation".format(i,j))
#                if j in sDomain and i in cDomain:
#                    print("{0} and {1} have c to s domain interation".format(i,j))
#                if i not in sDomain and i not in cDomain and j in sDomain:
#                    print("{0} and {1} have s to none domain interation".format(i,j)) 
#                if i not in sDomain and i not in cDomain and j in cDomain:
#                    print("{0} and {1} have c to none domain interation".format(i,j)) 
#                if j not in sDomain and j not in cDomain and i in sDomain:
#                    print("{0} and {1} have s to none domain interation".format(i,j)) 
#                if j not in sDomain and j not in cDomain and i in cDomain:
#                    print("{0} and {1} have c to none domain interation".format(i,j)) 
    with open(os.path.join(outDir,'{0}.{1}.interclusterNodes.{2}.{3}.txt'.format(condition,rnaName,cluster,percentile)), "w+", newline = '') as outFile:
        for i,j,c in zip(nodeEdge['Source'], nodeEdge['Target'],nodeEdge['Correlation']):
            if (nodeClustDict[str(i)] != cluster and nodeClustDict[str(j)] == cluster) and c > clusterCut[cluster]:
                outFile.write('{0}\t{1}\n'.format(i,j)) 
#                if i in sDomain and j in cDomain:
#                    print("{0} and {1} have s to c domain interation".format(i,j))
#                if j in sDomain and i in cDomain:
#                    print("{0} and {1} have c to s domain interation".format(i,j))
#                if i not in sDomain and i not in cDomain and j in sDomain:
#                    print("{0} and {1} have s to none domain interation".format(i,j)) 
#                if i not in sDomain and i not in cDomain and j in cDomain:
#                    print("{0} and {1} have c to none domain interation".format(i,j)) 
#                if j not in sDomain and j not in cDomain and i in sDomain:
#                    print("{0} and {1} have s to none domain interation".format(i,j)) 
#                if j not in sDomain and j not in cDomain and i in cDomain:
#                    print("{0} and {1} have c to none domain interation".format(i,j)) 
            if (nodeClustDict[str(i)] == cluster and nodeClustDict[str(j)] != cluster) and c > clusterCut[cluster]:
                outFile.write('{0}\t{1}\n'.format(i,j)) 

                
'''             
                
                