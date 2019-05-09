# -*- coding: utf-8 -*-
"""
Created on Wed May 30 13:51:19 2018
Network Graph to test integrity of base pair detection. (inspired by
A.L. and N.L. meeting)
Note that a cutoff is imposed in the resCorr Function!!!
@author: nlama
"""

######################## IMPORT NECESSARY LIBRARIES ###########################
import os
import sys
import numpy as np
import csv
import pandas as pd

##### ADD PATHS ###############################################################
sys.path.append('D:\Documents2\Rotations\Weeks\PythonScripts') #add RNAtools dir
sys.path.append('D:\Documents2\Rotations\Weeks\RNAFiles\TPP')

############# CREATING NETWORK ANALYSIS CLASS #################################

class NetworkAnalysis(object):

    def __init__(self, directory = os.getcwd(), ct=None, 
                 fileCorr=None):
        
        #assign file names and appropriate directories
        self.directory = directory
        self.fileCorr = os.path.join(directory,fileCorr)
        self.rnaName = "champ.fill"
        if ct is not None:
            self.ct = self.readCt(ct) #returns ct, assigns rnaName and rnaLen
        self.fileName = os.path.splitext(os.path.basename(self.fileCorr))[0] #Get file name without extention
        self.corrF = open(self.fileCorr) #same for RING output   
        
        # prepare distance, coordinate, and nt lists
        self.cutoffs = self.cutoffMaker(self.fileCorr)
        self.resCorrDist = self.resCorrDistFunc(WC=False) #get coordinates and distances


################################ FUNCTIONS ####################################
        
    def readCt(self,ctF):
        ct = open(ctF)
        lines = []
        for i in range(1):
            l=ct.readline().split()#get header information
            self.rnaLen=int(l[0])
            self.rnaName=l[1]
        for line in ct:
                lineList = line.split()
                lines.append(lineList)
        ct.close()
        ct = pd.DataFrame(lines)      
        return ct
    
    def fillResCorr(self,res,ij):
        fillUpTo = int(res[0][0])
        x=1
        while x < fillUpTo:
            res.append([0,x,16])
            x+=1
        for num in range(x,self.rnaLen+1):
            ind = ij.get(num,0)
            if ind == 0:
                res.append([0,num,16])
        return res
        
#    Look up residues in coorPair for computing distances
    def resCorrDistFunc(self, WC = False):
        resCorrDist = [] #will store correlations and distances for residues i and j
        iJDict = {}
        for i in range(2):
            self.corrF.readline() #skip first 2 lines because they are headers

        for line in self.corrF:
            corrLine = line.split() #split line into list
            if len(corrLine) > 1:
                i = corrLine[0] #extract residue1 in pair
                j = corrLine[1] #extract residue2 in pair
                iJDict[int(i)]=int(j)
                iJDict[int(j)]=int(i)
                self.cutoff = self.cutoffs[1]
                if float(np.nanmean([float(corrLine[5]), float(corrLine[6])])) >= self.cutoff:
                    resCorrDist.append([i , j , 
                                   float(np.nanmean([float(corrLine[5]), float(corrLine[6])]))])  #i,j,mean zscore
        if resCorrDist == []:
            raise ValueError("returned empty resCorrDist list")
        else:
#            resCorrDist=self.fillResCorr(resCorrDist,iJDict)
            return resCorrDist
              
    
    def cutoffMaker(self,fileCorr):
        corrF = open(fileCorr)
        cutoffVec = []
        for i in range(2):
            corrF.readline() #skip first 2 lines because they are headers
        for line in corrF:
            coor = line.split() #split line into list
            if len(coor) > 1:
                cutoffVec.append(np.nanmean([float(coor[5]), float(coor[6])])) #extract correlation  
        self.q = [25,30,75,90]
        self.m = np.percentile(cutoffVec, self.q)
        return self.m

    

#######################  PLOTTING FUNCTIONS ###################################
    
###############################################################################

if __name__ == '__main__':


#    directory = 'D:/Weeks/Data/NetworkAnalysis/resCorrDistFiles'
    directory = 'D:/Weeks/Data/JE_TPP_TMO'
#    ct = 'D:/Weeks/Data/BasePairDetectionProject/SecondaryStructure/AmpliconData/ec16S.ct'
    ct = 'D:/Weeks/Data/TPP_pairMap_ringMap/TPP_+L_structure.ct'
#    directory = 'D:/Weeks/Data/ringRnaseP'
#    ct = 'D:/Weeks/Data/BasePairDetectionProject/SecondaryStructure/AmpliconData/rnpB.ct'

#    directory = 'D:/Weeks/Data/NetworkAnalysis/TPP'
#    ct = 'D:/Weeks/Data/NetworkAnalysis/TPP/TPP_+L_structure.ct'
##    

    fileCorr = 'tpp.l.win5.tmo.ring'
    net = NetworkAnalysis( 
                   directory = directory, ct=ct, 
                   fileCorr=fileCorr)
    
    header = ['Source','Target','Weight']
    resCorrFileName = net.rnaName + '.'+ fileCorr + '.' + str(net.q[1]) + '.csv'
    rowNum=1
    with open(os.path.join('D:/Weeks/Data/NetworkAnalysis/TPPResCorr',resCorrFileName), "w+", newline = '') as file1:
        fWriter = csv.writer(file1, delimiter=',')
        hWriter = csv.writer(file1, delimiter=',')
        hWriter.writerow(header)
        for line in net.resCorrDist:
            fWriter.writerow(line)
            rowNum+=1
#