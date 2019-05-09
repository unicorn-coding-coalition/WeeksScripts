# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 11:37:39 2018

@author: nlama
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 16:50:29 2018

@author: nlama
"""

import matplotlib.pyplot as plt
from math import isnan
import numpy as np
import plotTools 
import pandas as pd
from sklearn.metrics import  auc

import seaborn as sns
sns.set('talk', 'whitegrid', 'dark', font_scale=1.5, font='Arial',
        rc={"lines.linewidth": 2, 'grid.linestyle': '--'})

rnaNames = ['ec16S', 'ec23S']
incfs = ['1M7', 'cellfree', 'incell']
profiles = []
ctFs =[]
names = [] #used for plotting :D
for rna in rnaNames:
    for incf in incfs:
        names.append('{0} {1}'.format(incf,rna))
        profiles.append('D:/Weeks/Data/BasePairDetectionProject/SecondaryStructure/AmpliconData/defaultProfiles/{0}_{1}_profile.txt'.format(incf,rna))
        ctFs.append('D:/Weeks/Data/BasePairDetectionProject/SecondaryStructure/AmpliconData/{0}.ct'.format(rna))
        #fa = 'D:/Weeks/Data/BasePairDetectionProject/SecondaryStructure/AmpliconData/{}.fa'.format(rnaName)

################################################################################
class ROCobj(object):
    
    def __init__(self, profile=None, ctFile = None):
       
        self.profile = profile        
        self.rmrpReactivity = plotTools.ReactivityProfile(plusfile = self.profile)
        self.nmFactors = self.rmrpReactivity.normalize(DMS=True) #normalize reactivites for nts
        self.rmrpCt = self.readCt(ctFile)
        self.rmrpBp = self.wcMasker(self.rmrpCt) #Actually I end up not using this?
        self.rmrpMa = self.createNtMasks(self.rmrpReactivity) #Create array of booleans for nucleicAcids
        self.rmrpRc = self.rmrpReactivity.normprofile #store reactive nts in array
        self.rmrpNt = [] #initialize array that will contain T/F for nt
        self.ntIdx = [] #initialize array containing seq array for nucelic acids
        self.ntReact = []
        self.ntReactDict = []
        
        for x in range(len(self.rmrpMa)):
            self.rmrpNt.append(self.rmrpRc[self.rmrpMa[x]])
        
        self.tpSeq = []
        self.fpSeq = []
        
        
        for x in range(len(self.rmrpMa)):
            self.ntIdx.append(self.idxSaver(self.rmrpMa[x])) #get idx corresponding to reactivities. 
            
        
        for x in range(len(self.rmrpNt)):
            self.ntReact.append(self.zipper(self.ntIdx[x], self.rmrpNt[x])) #make idx reactivity dicts for each nucleic acid
            
        
        
        for idx, i in enumerate(self.rmrpCt[4]): #Get idx where base pairing occurs or not
            if i !="0":
                self.tpSeq.append(idx)
            else:
                self.fpSeq.append(idx)
        
        
        
        for nt in range(len(self.ntReact)):
            self.ntReactDict.append({k: self.ntReact[nt][k] for k in self.ntReact[nt] if not isnan(self.ntReact[nt][k])})
        
        
        threshReactivity=np.arange(-2, 5, 0.005)
        
        self.tprs = []
        self.fprs = []
        
        #get true and false positive rates for all nucleic acids
        for x in range(len(self.ntReactDict)):
            self.tprs.append(self.getTpr(threshReactivity, self.ntReactDict[x]))
            self.fprs.append(self.getFpr(threshReactivity,self.ntReactDict[x]))

        

    def readCt(self, ctF):
        ct = open(ctF)
        lines = []
        for i in range(1):
                ct.readline() #skip first line because it is a header
        for line in ct:
                lineList = line.split()
                lines.append(lineList)
        ct = pd.DataFrame(lines)
        return ct
    
    def wcMasker(self, ct):
        wcMask = ct[4] != '0' 
        return wcMask
    
    def createNtMasks(self,reactivity):
        aMask = reactivity.sequence == 'A'
        cMask = reactivity.sequence == 'C'
        gMask = reactivity.sequence == 'G'
        uMask = reactivity.sequence == 'U'
        return aMask, cMask, gMask, uMask
    
    def idxSaver(self,rmrpNt):
        idxSaver = []
        for idx,nt in enumerate(rmrpNt):
            if nt == True:
                idxSaver.append(idx)
        return idxSaver
    
    def zipper(self,ntIdx,reactivities): #makes dict with nt idxs and reactivities
        keys=ntIdx
        values=reactivities
        dictionary=dict(zip(keys,values))
        return dictionary
    
    
    def getTpr(self,threshReactivity,ntReact):
        tpr = []
        for threshold in threshReactivity:
            tpArray= 0
            allTpBpNts=0
            for key, value in ntReact.items():
                if key in self.tpSeq: #if the idx is base paired, add to tp list
                    allTpBpNts+=1
                if key in self.tpSeq and float(value) < threshold :
                    tpArray+=1 #array of low reactivites whose nts basepair
            tpr.append(tpArray/allTpBpNts)
        return tpr
        
        
            
            
    def getFpr(self,threshReactivity,ntReact):
        fpr = []
        for threshold in threshReactivity:
            allFpBpNts=0
            fpArray=0
            for key, value in ntReact.items():
                if key in self.fpSeq:#if the idx is single stranded, add to tp list
                    allFpBpNts+=1
                if key in self.fpSeq and float(value) < threshold:
                    fpArray+=1 #array of low reactivies that don't bp
            fpr.append(fpArray/allFpBpNts)
        return fpr

    
    
############################ END OF OBJECT ####################################  
def plotROC(fpr,tpr,nucleicAcid,names):
    
    roc_auc = []
    for x in range(len(fpr)):
        roc_auc.append(auc(fpr[x], tpr[x]))
       
    csfont = {'fontname':'Lucida Bright', 'fontsize':50}
    hfont = {'fontname':'Lucida Bright', 'fontsize':35}
    lw = 8
    colors = ['#9E70A5', '#704993', '#C5FFFD', '#076399', 
              '#FF96AB', '#AA1E4D', '#F9A852', '#BA5900',
              '#AFC2D5', '#374B4A', '#6C3F15', '#312509']
    plt.figure(figsize=(20,15))
    for x in range(len(fpr)):
        plt.hold(True)
        plt.plot(fpr[x], tpr[x], color=colors[x],
                 lw=lw, label='{0} (AUC = {1:0.2f})'.format(names[x],roc_auc[x]))
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', **hfont)
    plt.ylabel('True Positive Rate', **hfont)
    plt.title("ROC " + nucleicAcid + "\n(Norm Reactivites vs Base Pairing)", **csfont)
    plt.legend(loc="lower right")
    plt.show()
#    plt.savefig('ROC_{0}_{1}_{2}.png'.format(incf,rnaName,nucleicAcid))


###############################################################################
#Create plotTools reactivity object 
rocObjArray = []
nucleicAcids = ['Adenine','Cytosine','Guanine','Uracil']

for i,j in zip(profiles,ctFs):
        rocObjArray.append(ROCobj(i,j))

fprArray = []
tprArray = []

for x in range(len(rocObjArray)):
    fprArray.append(rocObjArray[x].fprs)
    tprArray.append(rocObjArray[x].tprs)


byNtfprs= []
byNttprs = []

for x in range(len(nucleicAcids)):
    byNtfprs_tempy = []
    byNttprs_tempy = []
    for y in range(len(fprArray)):
        byNtfprs_tempy.append(fprArray[y][x])
        byNttprs_tempy.append(tprArray[y][x])
    byNtfprs.append(byNtfprs_tempy)
    byNttprs.append(byNttprs_tempy)
    

for f,t,n in zip(byNtfprs,byNttprs,nucleicAcids):
        plotROC(f,t,n,names)

