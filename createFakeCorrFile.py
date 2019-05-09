# -*- coding: utf-8 -*-
"""
 
@author: nlama
"""
import numpy as np
import math
import os
import scipy
import random as rand

class DistCorr(object):

    def __init__(self, filePDB=None, directory = None, fa = None, 
                 strand = None, **kwargs):
        
        self.pdbF = open(filePDB) #open pdb file into memory
        
        self.coorPair = self.corrPairFunc(self.pdbF, strand = strand) #prepare corr file so it's NT and coordinate
            
        self.kwargs = kwargs
        if 'ct' in kwargs: 
            self.ct = open(kwargs.get('ct'))
            self.backgroundDistList = self.backgroundDist()
            self.probRing()

#################################FUNCTIONS#####################################           
    
    def findDist(self, i, j):       #Finds distances between residues
#        print(self.coorPair)
#        print(i)
        pdbI = self.pdbNumbering[int(i)]
        pdbJ = self.pdbNumbering[int(j)]
        iCoords = self.coorPair[pdbI] #find residue 1's coordinates
        jCoords = self.coorPair[pdbJ] #find residue 2's coordinates
#        print(jCoords)
        x = (float(iCoords[0]) - float(jCoords[0]))**2 #compute dist form x 
        y = (float(iCoords[1]) - float(jCoords[1]))**2 #compute dist form y
        z = (float(iCoords[2]) - float(jCoords[2]))**2 #compute dist form z
        distance = math.sqrt(x + y + z) # compute distance between residues
        return distance
    
    def backgroundDist(self):
        residueDist = []
        self.backgroundDistances=[]
        #randPairs = []
        for x in range(1,len(self.coorPair)+1):
            for y in range(1,x):
                distance = self.findDist(str(x),str(y))
                residueDist.append([x , y , 0, distance])
                self.backgroundDistances.append(distance)
        return residueDist
    
                      
    def probRing(self):
        #Let's find random real data that matches shape distribution for highest correlations
        self.randBp = []
        ringProbs = [0.03646308, 0.2032817, 0.2315406, 0.2160438,
                 0.1599818, 0.07976299, 0.04056518, 
                 0.02096627, 0.009571559, 0.001823154, 0 ]
#        probs = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1] #uniform distribution
        cumProbs = ringProbs
        for i in range(len(ringProbs)):
            if i == 0:
                cumProbs[i] = ringProbs[i]
            else:
                cumProbs[i] += cumProbs[i-1]
        
        i = 0
        self.randRingPairs = []
        self.numCorrs = 20
        while i < self.numCorrs:
            randIdx1 = rand.randint(0,len(self.backgroundDistances)-1)
            indicator= self.randomPick(self.backgroundDistances[randIdx1], ringProbs, cumProbs)
            if indicator and self.backgroundDistList[randIdx1][0] - self.backgroundDistList[randIdx1][1] > 5:
                self.randBp.append(self.backgroundDistances[randIdx1])
                self.randRingPairs.append([self.backgroundDistList[randIdx1][0] , self.backgroundDistList[randIdx1][1]]) 
                i += 1

                         
    def readFasta(self, fastapath):
        """Assign sequence from fastafile"""
        
        with open(fastapath) as inp:
            
            inp.readline() # pop off the header
            seq = inp.readline().strip()

        self.seq = list(seq.replace('T','U'))
        
    def corrPairFunc(self,pdbF,strand, atom = "C1'"):
        coorPair = {} #will store residue and its coordinates
        self.atomPair = {}
        self.pdbNumbering = []
        self.pdbNumbering.append(0)
        for line in pdbF:
             if strand == None:
                 if line[13:17].strip() == atom: #only keep if C1 is 2nd eleme
    #                print('in if')
                    self.pdbNumbering.append(int(line[23:26]))
                    coordinateList=[line[31:38].strip(),line[39:46].strip(),line[47:54].strip()] #Get coords xyz from pdb
                    coorPair[int(line[23:26])]=coordinateList
                    self.atomPair[int(line[23:26])]=line[17:20].strip()
            
             else:
                 if line[13:17].strip() == atom and line[21].strip() == strand: #only keep if C1 is 2nd eleme
    #                print('in if')
                    self.pdbNumbering.append(int(line[23:26]))
                    coordinateList=[line[31:38].strip(),line[39:46].strip(),line[47:54].strip()] #Get coords xyz from pdb
                    coorPair[int(line[23:26])]=coordinateList
                    self.atomPair[int(line[23:26])]=line[17:20].strip()
#        print(coorPair)
#        print(self.atomPair)
        return coorPair
    
    def watsonCrick(self, i,j):
        pdbI = self.pdbNumbering[int(i)]
        pdbJ = self.pdbNumbering[int(j)]
        nt1 = self.atomPair[pdbI] #Dict containing {'i': 'nucleotide (G,A,C,U)'}
        nt2 = self.atomPair[pdbJ]
        pair = nt1 + nt2
        
        if pair in ('AU', 'UA', 'GC', 'CG', 'GU', 'UG'):
            return True   
        else:
            return False

    
    
    def rangeAssigner(self,dist):
        cumDist = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        for i in range(len(cumDist)):
#            print(i)
            if i+1 >= len(cumDist):
                break
            if dist > cumDist[10]:
                #print(dist)
                return 10
            elif cumDist[i] < dist <= cumDist[i+1]:
#                print("{0} is between {1} and {2}".format(dist, cumDist[i],cumDist[i+1]))
                return i
              
    def randomPick(self, dist, probabilities, cumulativeProbability):
       # cumulativeProbabilityTempy = probabilities
        x = rand.uniform(0,1)
        probIdx = self.rangeAssigner(dist)
        #print(probIdx)
        #cumulativeProbabilityTempy[probIdx] += probabilities[probIdx-1]
        if x >= cumulativeProbability[probIdx] and probabilities[probIdx] != 0: 
            #print(probabilities[probIdx])
            #print(cumulativeProbability[probIdx])
            return True#, cumulativeProbabilityTempy
        else:
            return False#, cumulativeProbability
        
        
    def createOutFile(self):
        fileNum = rand.randint(0,500)
        ringCorrFile = 'mockRingCorrellations.{0}.minCorr5.{1}.txt'.format(self.numCorrs, fileNum)
        with open(os.path.join('D:/Weeks/Data/FARFAR/2n3r',ringCorrFile), "w+", newline = '') as file1:
            file1.write("68\tWindow=3\tMetric=APC\n")
            file1.write("i\tj\tChi2\t+/-\tCorr\tZi\tZj\tDepth\tComuts\n")
            for ringPair in self.randRingPairs:
                i=ringPair[0]
                j=ringPair[1]
                file1.write("{0}\t{1}\t1\t1\t2\t2\t1\t1000\t50\n".format(i,j))


            
    
###############################################################################

if __name__ == '__main__':
    ### CHANGE ALL OF THESE TO FIT YOUR NEEDS
    filePDB= 'D:/Weeks/Data/FARFAR/2n3r\\2n3r.model1.pdb'
#    fileFasta = 'D:/Weeks/Data/FARFAR/1HC8\\1HC8.omit1.fa'
    fileCt = 'D:/Weeks/Data/FARFAR/2n3r\\2n3r.ct' #add if you want background distances of your molecule
    directory = 'D:/Weeks/Data/FARFAR/2n3r'
   ###############################################################################3
    # initialize the object and read in matrices

    mi_apc = DistCorr(prefix = '2n3r', 
                      filePDB=filePDB, directory = directory, fa = None,
                      ct = fileCt, WC = 'WC', strand = None)
    
    mi_apc.createOutFile()
#    
#    fName = "distListP546.txt"
#    with open(os.path.join('D:/Weeks/Data/NetworkAnalysis/resCorrDistFiles',fName), "w+", newline = '') as file1:
#        fWriter = csv.writer(file1, delimiter=',')
#        fWriter.writerow(mi_apc.distList[0])