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
import math
import numpy as np
import csv
import plotly.plotly as py
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import networkx as nx

##### ADD PATHS ###############################################################
sys.path.append('D:\Documents2\Rotations\Weeks\PythonScripts') #add RNAtools dir
sys.path.append('D:\Documents2\Rotations\Weeks\RNAFiles\TPP')

############# CREATING NETWORK ANALYSIS CLASS #################################

class NetworkAnalysis(object):

    def __init__(self, directory = os.getcwd(), 
                 filePDB=None, fileCorr=None):
        
        #assign file names and appropriate directories
        self.directory = directory
        self.fileCorr = os.path.join(directory,fileCorr)
        self.fileName = os.path.splitext(os.path.basename(self.fileCorr))[0] #Get file name without extention
        self.pdbF = os.path.join(directory,filePDB)
        self.pdbF = open(self.pdbF) #open pdb file into memory
        self.corrF = open(self.fileCorr) #same for RING output   
        
        # prepare distance, coordinate, and nt lists
        self.cutoffs = self.cutoffMaker(self.fileCorr)
        self.coorPair = self.corrPairFunc(self.pdbF) #prepare corr file so it's NT and coordinate
        self.resCorrDist = self.resCorrDistFunc(WC=False) #get coordinates and distances


################################ FUNCTIONS ####################################
        
    def corrPairFunc(self,pdbF):

        coorPair = {} #will store residue and its coordinates
        self.atomPair = {}
        for line in pdbF:
            lineList = line.split()
            if len(lineList)>5 and lineList[2] == "C1'": #only keep if C1' atom
                coorPair[lineList[5]]=lineList[6:9] #keep only residue and xyz
                self.atomPair[lineList[5]]=lineList[3]
        return coorPair
    
    def watsonCrick(self, i,j):
        nt1 = self.atomPair[i] #Dict containing {'i': 'nucleotide (G,A,C,U)'}
        nt2 = self.atomPair[j]
        pair = nt1 + nt2
        
        if pair in ('AU', 'UA', 'GC', 'CG', 'GU', 'UG'):
            return True   
        else:
            return False
        
    def findDist(self, i, j):       #Finds distances between residues
        iCoords = self.coorPair[i] #find residue 1's coordinates
        jCoords = self.coorPair[j] #find residue 2's coordinates
        self.sourceX=float(iCoords[0])
        self.sourceY=float(iCoords[1])
        self.sourceZ=float(iCoords[2])
        self.targetX=float(jCoords[0])
        self.targetY=float(jCoords[1])
        self.targetZ=float(jCoords[2])
        self.x1 = (float(iCoords[0]) - float(jCoords[0]))**2 #compute dist form x 
        self.y1 = (float(iCoords[1]) - float(jCoords[1]))**2 #compute dist form y
        self.z1 = (float(iCoords[2]) - float(jCoords[2]))**2 #compute dist form z
        distance = math.sqrt(self.x1 + self.y1 + self.z1) # compute distance between residues
        return distance
    
    #Look up residues in coorPair for computing distances
    def resCorrDistFunc(self, WC = False):
        resCorrDist = [] #will store correlations and distances for residues i and j
        for i in range(2):
            self.corrF.readline() #skip first 2 lines because they are headers
        for line in self.corrF:
            corrLine = line.split() #split line into list
            if len(corrLine) > 1:
                i = corrLine[0] #extract residue1 in pair
                j = corrLine[1] #extract residue2 in pair
#                print(i,j)
                distance = self.findDist(i,j)
#                print(distance)
#                if self.watsonCrick(i,j) and WC:
                self.cutoff = self.cutoffs[3]
                if np.mean([float(corrLine[5]), float(corrLine[6])]) >= self.cutoff:
                    resCorrDist.append([i , j , 
                                    np.mean([float(corrLine[5]), 
                                    float(corrLine[6])]), 
                                    distance])  #, self.sourceX, self.sourceY, self.sourceZ,
                                    #self.targetX, self.targetY, self.targetZ])
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
        self.q = [25,50,50,50]
        self.m = np.percentile(cutoffVec, self.q)
        return self.m

    def CreateNetwork(self,resCorrDist):
        G = nx.Graph()
        G.clear()
        correlations = []
        distances = []
        for x in range(len(resCorrDist)):
            G.add_edge(int(resCorrDist[x][0]), int(resCorrDist[x][1]))
            correlations.append(float(resCorrDist[x][2]))
            distances.append(int(resCorrDist[x][3]))
            
        for idx, (u,v,d) in enumerate(G.edges(data=True)):
            d['weight'] = correlations[idx]
            d['length'] = distances[idx]
#            print(u,v,d)
        
        edges, weights = zip(*nx.get_edge_attributes(G, 'weight').items())
        lengths = nx.get_edge_attributes(G, 'length')
        
        pos = nx.spring_layout(G)
        ## Drawing the graph
        options = {
                'node_color': '#2C3539', 
                'node_size': 1000,
                'width': 4,
                'alpha': 0.85,
                'linewidths': 4.0,
                'font_color': 'white',
                'font_weight': 'bold',
                'font_size': 20,
                'edge_color': weights,
                'edge_cmap': plt.cm.Blues,
                'edgelist': edges,
                }
        
        plt.clf()
        plt.figure(102)
        nx.draw(G, pos, with_labels=True , **options)
        nx.draw_networkx_edge_labels(G, pos, edge_labels= lengths, font_size=12)
        plt.savefig('metric_network.png')
        plt.show()
    

#######################  PLOTTING FUNCTIONS ###################################
    
###############################################################################

if __name__ == '__main__':

    #D:/Documents2/Rotations/Weeks/RNAFiles/TPP/2GDI.pdb
    filePDB= 'D:\Documents2\Rotations\Weeks\RNAFiles\RNaseP\\00742V_n2_sup.pdb'
#    filePDB= 'D:/Weeks/Data/BasePairDetectionProject/SecondaryStructure/AmpliconData/tmRNA.pdb'
    directory = 'D:/Weeks/Data/NetworkAnalysis/resCorrDistFiles'

    mi_apc = NetworkAnalysis( 
                   filePDB=filePDB,
                   directory = directory,
                   fileCorr='cellfree.RMRP.apc.win3.minCorr5.minCount50.txt')
    
    
    mi_apc.CreateNetwork(mi_apc.resCorrDist)
    header =["Source,Target,Correlation,Distance"]
    resCorrFileName = mi_apc.fileName + '.' + str(mi_apc.q[3]) + '.nodist.resCorrDist.csv'
    with open(os.path.join('D:/Weeks/Data/NetworkAnalysis/resCorrDistFiles',resCorrFileName), "w+", newline = '') as file1:
        fWriter = csv.writer(file1, delimiter=',')
        hWriter = csv.writer(file1, delimiter=' ')
        hWriter.writerow(header)
        for line in mi_apc.resCorrDist:
            fWriter.writerow(line)

    
    