# -*- coding: utf-8 -*-
"""
DistCorrAnalysis is a script used to create correlation - distance histograms
Created on Wed Mar 28 13:37:43 2018
Note that a cutoff is imposed in the resCorr Function!!!
line graph, histograms, kdes of all files. 
@author: nlama
"""
import numpy as np
from sklearn.metrics import r2_score
import math
import os
import matplotlib.pyplot as plot
import scipy
import seaborn as sns
import sys
import random as rand
import pandas as pd
sys.path.append('D:\Documents2\Rotations\Weeks\PythonScripts') #add RNAtools dir

class DistCorr(object):

    def __init__(self, prefix = None, score = 'Z',
                 filePDB=None, directory = None, fileCorr=None, **kwargs):
        
        self.prefix = prefix
        self.pdbF = open(filePDB) #open pdb file into memory
        self.fileCorr = os.path.join(directory, fileCorr)
        self.corrF = open(self.fileCorr) #same for RING output 
        self.score = score
        
        
        self.cutoffs = self.cutoffMaker(self.fileCorr)
        self.coorPair = self.corrPairFunc(self.pdbF) #prepare corr file so it's NT and coordinate
        self.resCorrDist = self.resCorrDistFunc() #get coordinates and distances
        if self.prefix == None:
            self.prefix = 'Correlation Histograms'
            
        self.kwargs = kwargs
        if 'ct' in kwargs: 
            self.ct = open(kwargs.get('ct'))
            self.backgroundDistList = self.backgroundDist()
            self.backgroundBasePairDistList = self.backgroundBasePairDist()
            
    
##################### Plotting Functions ######################################
            
    def plotter(self):
        distList = []
        self.m = self.cutoffMaker(self.fileCorr) #get quartile cutoffs for coors
        for i in range(len(self.m)):
            distList.append(self.getDist(self.resCorrDist, cutoff = self.m[i])) 
        self.distList=distList
        self.binRange = range(0,160,3) #make bins all same for all plots
        self.simpleFlatui = [ "#95a5a6", '#FF9933', "#31B404", '#FACC2E',  "#e74c3c", "#9b59b6", 
          "#3498db",  "#34495e", "#2ecc71" ]
        sns.set(font_scale=3.5, rc={"figure.figsize": (18, 13)})
        sns.set_style('whitegrid')
        sns.set_palette(self.simpleFlatui)
        figNum = rand.randint(1,10000) #Generate a random num for figure
        plot.xlim([0,175])
        plot.figure(figNum)
#        plot.ylim(0, 0.2)
        plot.xlim(-5,175) #Keep scales the same for easier comparisons
        self.plotHists() #Function to plot histograms
        if 'ct' in self.kwargs:
            self.plotBackgroundHist() #Function to plot background histograms
        plot.legend()
        plot.xlabel('Distances')
        plot.ylabel('Normalized Counts')
        plot.title('{0}'.format(self.prefix))
        plot.show()
#        plot.savefig('{0}_histograms.png'.format(self.fileCorr)) #save plot to working dir
    
    def plotHists(self):
        for i in range(len(self.distList)):
            plot.hist(self.distList[i], alpha =1 , bins= self.binRange, 
                      histtype = 'step', linewidth = 4,
                      label = '{2} percentile;   n = {1}'.format(self.m[i], 
                               str(len(self.distList[i])), self.q[i]),
                      normed = True)
                      
    def plotBackgroundHist(self):
        self.backgroundDistances = self.getDist(self.backgroundDistList, cutoff=-10)
        self.backgroundBasePairDistances = self.getDist(self.backgroundBasePairDistList, cutoff=-10)
        plot.hist(self.backgroundDistances, alpha =1, bins= self.binRange, 
                      histtype = 'step', linewidth = 4,
                      label = 'background', normed = True)
                      
    def genLineGraph(self, pf, title, save = 'no'):
        figNum = rand.randint(1,10000)
        sns.set(font_scale=4, style='whitegrid', 
                    rc={"figure.figsize": (20, 12),"lines.linewidth":4})
        plot.figure(figNum)
        plot.xlim([0,80])
        ax = sns.regplot(x= 'Distance', y= 'Correlation', 
                           data = pf, lowess=False)

        plot.title('Correlation v.s. Distance ' + '(' + title +')' )
        if save == 'yes':
            plot.savefig('D:\Documents2\Rotations\Weeks\Week7Structures/lineGraph_mi_apc_phi_allPercentiles_CD20.png')
            
    def comparativePlots(self,apc,mi,phi,percentile):
        self.flatui2 = ["#e74c3c", '#FF9933','#FACC2E', "#31B404", "#9b59b6", 
          "#3498db",  "#34495e", "#2ecc71" ,  "#95a5a6"]
        figNum = rand.randint(1,10000)
        sns.set(font_scale=2, rc={"figure.figsize": (11, 7)})
        sns.set_style('whitegrid')
        sns.set_palette(self.simpleFlatui)
        plot.figure(figNum)
        plot.legend()
    

#        plot.ylim(0, 0.04)
        plot.xlim(-5,160)
        
        a=sns.distplot(apc, bins= self.binRange, hist = False, kde_kws={"linewidth": 3, 'alpha' : 1}, 
                     label = 'mi + apc')
        b=sns.distplot(mi, bins= self.binRange, hist = False, kde_kws={"linewidth": 3, 'alpha' : 1}, 
                     label = 'mi only')
        c=sns.distplot(phi, bins= self.binRange, hist = False, kde_kws={"linewidth": 3, 'alpha' : 1}, 
                     label = 'phi') 
        
        self.medianMaker(0,a,"#95a5a6")
        self.medianMaker(1,b, '#FF9933')
        self.medianMaker(2,c, 'green')
        
        plot.show()
        plot.xlabel('Distances')
        plot.ylabel('Normalized Counts')
        plot.title('KDE of PHI, MI and MI + APC \n {} percentile'.format(percentile))
#        plot.savefig('{0}_kde_mi_apc_phi_{1}.png'.format(self.fileCorr, percentile))
         
##################### FUNCTIONS ###############################################
   
    def corrPairFunc(self,pdbF, atom = "C1'"):
        c1List = [] #list stores full line with C1
        coorPair = {} #will store residue and its coordinates
        self.atomPair = {}
#        residueNum = 1 #only switch on if pdb is numbered not 1-end
        for line in pdbF:
            lineList = line.split()
            if atom == 'N1':
                pyr = ['U', 'C']
                pur = ['A','G']
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
#                coorPair[str(residueNum)]=lineList[6:9] #Only do this for pdbs with weird numbering! (not 0-end)
                coorPair[lineList[5]]=lineList[6:9] #keep only residue and xyz
                self.atomPair[lineList[5]]=lineList[3]
#                residueNum += 1
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
        x = (float(iCoords[0]) - float(jCoords[0]))**2 #compute dist form x 
        y = (float(iCoords[1]) - float(jCoords[1]))**2 #compute dist form y
        z = (float(iCoords[2]) - float(jCoords[2]))**2 #compute dist form z
        distance = math.sqrt(x + y + z) # compute distance between residues
        return distance
    
    def backgroundDist(self):
        residueDist = []
        for x in range(1,len(self.coorPair)+1):
            for y in range(1,x):
                distance = self.findDist(str(x),str(y))
                residueDist.append([x , y , 0, distance])
        return residueDist
    
    def backgroundBasePairDist(self):
        residueBasePairDist = []
        for line in self.ct:
            ctLine = line.split()
            if len(ctLine) > 3:
                i = ctLine[0]
                j = ctLine[4]
                if j != '0':
                    distance = self.findDist(i,j)
                    residueBasePairDist.append([i , j , 0, distance])
        return residueBasePairDist
    
    #Look up residues in coorPair for computing distances
    def resCorrDistFunc(self, WC = True):
        resCorrDist = [] #will store correlations and distances for residues i and j
        for i in range(2):
            self.corrF.readline() #skip first 2 lines because they are headers
        for line in self.corrF:
            corrLine = line.split() #split line into list
            if len(corrLine) > 1:
                i = corrLine[0] #extract residue1 in pair
                j = corrLine[1] #extract residue2 in pair
                distance = self.findDist(i,j)
                if self.watsonCrick(i,j) and WC:
                    resCorrDist.append([i , j , 
                                        np.mean([float(corrLine[5]), 
                                                 float(corrLine[6])]), 
                                                 distance])
#                elif self.score == 'Chi2': 
#                    resCorrDist.append([i , j , float(corrLine[2]), distance])
#                elif self.score == 'Z':
#                    resCorrDist.append([i , j , np.mean([float(corrLine[5]), float(corrLine[6])]), distance])
        if resCorrDist == []:
            raise ValueError("returned empty resCorrDist list")
        else:
            return resCorrDist
    
    
    ##Will create median!
    ##curveNum indicates which curve you want median to attach to. EX: if it's the
    ##first curve plotted on the graph, curveNum = 0   
    def medianMaker(self,curveNum,distList,chooseColor='k', chooseLine = 'dashed'):
        x,y = distList.get_lines()[curveNum].get_data()
        cdf = scipy.integrate.cumtrapz(y, x, initial=0)    
        nearest_05 = np.abs(cdf-0.5).argmin()
        x_median = x[nearest_05]
        y_median = y[nearest_05] 
#        print(x_median, y_median)
        plot.vlines(x_median, 0, y_median, color = chooseColor, 
                    linestyles= chooseLine)
          
    def getDist(self,resCorrDist, cutoff=-10):
        distList = []
        for line in resCorrDist:
            if line[2] > cutoff:
                distList.append(line[3])
        return distList
    
    def cutoffMaker(self,fileCorr):
        corrF = open(fileCorr)
        cutoffVec = []
        for i in range(2):
            corrF.readline() #skip first 2 lines because they are headers
        for line in corrF:
            coor = line.split() #split line into list
            if len(coor) > 1 and self.score == 'Chi2':
                cutoffVec.append(float(coor[2])) #extract correlation  
            elif len(coor) > 1 and self.score == 'Z':
                cutoffVec.append(np.mean([float(coor[5]), float(coor[6])])) #extract correlation  
        self.q = [25,50,95,99]
        m = np.percentile(cutoffVec, self.q)
        return m
    

        
    def readCorr(self, resCorrDist):
        df = pd.DataFrame(np.float64(resCorrDist),
                          columns=['i','j','Correlation','Distance'])
        df = df.loc[df['Correlation'] >= self.cutoffs[3],]
        return df
        

    
###############################################################################

if __name__ == '__main__':
    ### CHANGE ALL OF THESE TO FIT YOUR NEEDS
#    directory = 'D:\Documents2\Rotations\Weeks\PythonScripts'
    filePDB= 'D:\Documents2\Rotations\Weeks\RNAFiles\RNaseP\\00742V_n2_sup.pdb'
    fileCt = 'D:\Documents2\Rotations\Weeks\PythonScripts\\rnpb.ct' #add if you want background distances of your molecule
#    directory = 'C:\Users\nicol\Desktop\Tempy
    directory = 'D:/Weeks/Data/BasePairDetectionProject/SecondaryStructure/FinalCorrFiles/MinCount50'
   ###############################################################################3
    # initialize the object and read in matrices
    score = 'Z'
#    phi = DistCorr(prefix = 'PHI (Window = 3), score =  ' + score + ')', score = score, 
#                      filePDB=filePDB,
#                      directory = directory,
#                      fileCorr='rnasep.phi.3.txt')
    mi_apc = DistCorr(prefix = 'APC  (Window = 3, score =  ' + score + ')', score = score,
                      filePDB=filePDB, directory = directory,
                      fileCorr='cellfree.rnpB.apc.win3.minCorr5.minCount50.txt', ct = fileCt, WC = 'WC')
#    mi = DistCorr(prefix = 'MI (Window = 3), score =  ' + score + ')', score = score,
#                      filePDB=filePDB,
#                      directory = directory,
#                      fileCorr='rnasep.mi.3.txt')

    
    
    
    mi_apc.plotter()
    mi.plotter()
    phi.plotter()

    distNum = 3
    mi.comparativePlots(mi_apc.distList[distNum], 
                        mi.distList[distNum], phi.distList[distNum],
                        percentile = mi.q[distNum])

    
    
pf = mi_apc.readCorr(mi_apc.resCorrDist) #This is using the last q in m. (99)
mi_apc.genLineGraph(pf, title='MI-APC', save = 'no')
mi_apcR2 = r2_score(pf['Correlation'], pf['Distance'])

pf = mi.readCorr(mi.resCorrDist)
mi.genLineGraph(pf,  title='MI', save = 'no')
miR2 = r2_score(pf['Correlation'], pf['Distance'])

pf = phi.readCorr(phi.resCorrDist)
phi.genLineGraph(pf, title='PHI', save = 'no')
phiR2 = r2_score(pf['Correlation'], pf['Distance'])
