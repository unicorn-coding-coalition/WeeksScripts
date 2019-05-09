# -*- coding: utf-8 -*-
"""
Score all decoy structures and native pdb to ring constraints,
check by constraints, calculate a Lama score. 
@author: nlama
"""
import numpy as np
import math
import os
import matplotlib.pyplot as plot
import scipy
import random as rand
import pandas as pd
import glob

#sys.path.append('D:\Documents2\Rotations\Weeks\PythonScripts') #add RNAtools dir

class DistCorr(object):

    def __init__(self, prefix = None, score = None,
                 filePDB=None, directory = None, decoyArr = None, fa = None, fileCorr=None, **kwargs):
        
        self.prefix = prefix
        self.pdbF = open(filePDB) #open pdb file into memory
        self.fileCorr = os.path.join(directory, fileCorr)
        self.corrF = open(self.fileCorr) #same for RING output 
        self.energies = score
        self.decoyArr = decoyArr
        
        self.cutoffs = self.cutoffMaker(self.fileCorr)
        nativeCoorPair, nativeNumbering = self.corrPairFunc(self.pdbF) #prepare corr file so it's NT and coordinate
        self.resCorrDist = self.resCorrDistFunc(nativeCoorPair,nativeNumbering) #get coordinates and distances
        self.ringConstraints = self.readCorr()
        self.kwargs = kwargs
        self.energyBonus = self.scoreDecoys() #lama bonus
        
        for x in range(len(self.decoyArr)):
            print(self.energyBonus[x], self.decoyArr[x])
        minIdx = min(range(len(self.energyBonus)), key=self.energyBonus.__getitem__)
        print('Max energy Bonus of' , min(self.energyBonus), 'awarded to' , 
        os.path.splitext(os.path.basename(self.decoyArr[minIdx]))[0])
        
      
##################### FUNCTIONS ###############################################
   
    def readFasta(self, fastapath):
        """Assign sequence from fastafile"""
        
        with open(fastapath) as inp:
            
            inp.readline() # pop off the header
            seq = inp.readline().strip()

        self.seq = list(seq.replace('T','U'))
        
    def corrPairFunc(self, pdbF, atom = "C1'", strand = 'X', decoy=False):
        coorPair = {} #will store residue and its coordinates
        self.atomPair = {}
        pdbNumbering = []
        pdbNumbering.append(0)
        if decoy == True: #Don't account for strand 
            for line in pdbF:
                if line[13:17].strip() == atom:
                    pdbNumbering.append(int(line[23:26]))
                    coordinateList=[line[31:38].strip(),line[39:46].strip(),line[47:54].strip()] #Get coords xyz from pdb
                    coorPair[int(line[23:26])]=coordinateList
                    self.atomPair[int(line[23:26])]=line[17:20].strip()
            return coorPair, pdbNumbering
        for line in pdbF:
            if line[13:17].strip() == atom and line[21].strip() == strand: #only keep if C1 is 2nd eleme
                pdbNumbering.append(int(line[23:26]))
                coordinateList=[line[31:38].strip(),line[39:46].strip(),line[47:54].strip()] #Get coords xyz from pdb
                coorPair[int(line[23:26])]=coordinateList
                self.atomPair[int(line[23:26])]=line[17:20].strip()
        return coorPair, pdbNumbering
    
    def watsonCrick(self, i,j, pdbNumbering):
        pdbI = pdbNumbering[int(i)]
        pdbJ = pdbNumbering[int(j)]
        nt1 = self.atomPair[pdbI] #Dict containing {'i': 'nucleotide (G,A,C,U)'}
        nt2 = self.atomPair[pdbJ]
        pair = nt1 + nt2
        
        if pair in ('AU', 'UA', 'GC', 'CG', 'GU', 'UG'):
            return True   
        else:
            return False
    
    def findDist(self, i, j, coorDict, pdbNumbering):       #Finds distances between residues
        pdbI = pdbNumbering[int(i)]
        pdbJ = pdbNumbering[int(j)]
        iCoords = coorDict[pdbI] #find residue 1's coordinates
        jCoords = coorDict[pdbJ] #find residue 2's coordinates
        x = (float(iCoords[0]) - float(jCoords[0]))**2 #compute dist form x 
        y = (float(iCoords[1]) - float(jCoords[1]))**2 #compute dist form y
        z = (float(iCoords[2]) - float(jCoords[2]))**2 #compute dist form z
        distance = math.sqrt(x + y + z) # compute distance between residues
        return distance
    
    
    
    #Look up residues in coorPair for computing distances
    def resCorrDistFunc(self, coorPair, pdbNumbering, WC = False):
        resCorrDist = [] #will store correlations and distances for residues i and j
        for i in range(2):
            self.corrF.readline() #skip first 2 lines because they are headers
        for line in self.corrF:
            corrLine = line.split() #split line into list
            if len(corrLine) > 1:
                i = corrLine[0] #extract residue1 in pair
                j = corrLine[1] #extract residue2 in pair
                distance = self.findDist(i,j, coorPair, pdbNumbering)
                if WC:
                    if self.watsonCrick(i,j, pdbNumbering):
                        resCorrDist.append([i , j , 
                                            np.mean([float(corrLine[5]), 
                                                     float(corrLine[6])]), 
                                                     distance])
                else: 
                    resCorrDist.append([i , j , np.mean([float(corrLine[5]), 
                                                         float(corrLine[6])]), 
                             distance])
#
        if resCorrDist == []:
            raise ValueError("returned empty resCorrDist list")
        else:
            return resCorrDist
    
    
    def cutoffMaker(self,fileCorr):
        corrF = open(fileCorr)
        cutoffVec = []
        for i in range(2):
            corrF.readline() #skip first 2 lines because they are headers
        for line in corrF:
            coor = line.split() #split line into list
            if len(coor) > 1:
                cutoffVec.append(np.mean([float(coor[5]), float(coor[6])])) #extract correlation  
        self.q = [25,50,75,90,99]
        m = np.percentile(cutoffVec, self.q)
        return m
           
    def readCorr(self): #returns ring constraint list based on correlation
        df = pd.DataFrame(self.resCorrDist,
                          columns=['i','j','Correlation','Distance'])
        df = df.loc[df['Correlation'] >= self.cutoffs[1],]
        i = list(df['i'])
        j = list(df['j'])
        ringCorrs=[]
        for nti, ntj in zip(i,j):
            ringCorrs.append([nti,ntj])
        return ringCorrs
    
    def scoreDecoys(self):
        scoreTracker = np.zeros(len(self.decoyArr))
        for x in range(len(self.decoyArr)):
            pdb = open(self.decoyArr[x])
            decoyCoorPair, decoyNumbering = self.corrPairFunc(pdb, decoy=True) #prepare corr file so it's NT and coordinate
            for constraint in self.ringConstraints:
                rI = constraint[0]
                rJ = constraint[1]
                if self.findDist(rI,rJ,decoyCoorPair,decoyNumbering) < 20:
                    scoreTracker[x] += -0.5
                elif self.findDist(rI,rJ,decoyCoorPair,decoyNumbering) < 30:
                    scoreTracker[x] += -0.001
            pdb.close()
        return scoreTracker
                    
                    
                
                        

    
###############################################################################

if __name__ == '__main__':
    ### CHANGE ALL
    nativePDB= 'D:/Weeks/Data/FARFAR/TPP/2GDIX.pdb'
    rosPDBArray= []
    for decoy in glob.glob('D:/Weeks/Data/FARFAR/TPP/650/*.pdb'):
        rosPDBArray.append(decoy)
#    fileFasta = 'D:/Weeks/Data/FARFAR/1HC8\\1HC8.omit1.fa'
    directory = 'D:/Documents2/Rotations/Weeks/RNAFiles/TPP/TPP_RINGMAP_OUTPUT'
    scoreArr = []
   ###############################################################################3
    # initialize the object and read in matrices

    Obj = DistCorr(prefix = 'TPP decoy', score = scoreArr,
                      filePDB=nativePDB, directory = directory, decoyArr= rosPDBArray, fa = None,
                      fileCorr='tpp.apc.3.minCount50.minCorr5.corr.txt', WC = 'WC')

    
