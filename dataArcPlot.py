# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 20:34:47 2018

@author: nlama
"""

import numpy as np
import sys
import pandas as pd
sys.path.append('D:\Documents2\Rotations\Weeks\PythonScripts') #add RNAtools dir
import RNAtools #import RNAtools 

rnaObj = RNAtools.CT('rnpb.ct') #create RNAtools object

ct = open('rnpb.ct')
lines = []
for i in range(1):
        ct.readline() #skip first 2 lines because they are headers
for line in ct:
        lineList = line.split()
        lines.append(lineList)
ct = pd.DataFrame(lines)  

terInts = [] #list containing tertiary interactions
for i in range(len(ct)):
    for j in range(len(ct)):
        d = rnaObj.contactDistance(i,j)
        #print(d)
        if d > 20:
            terInts.append([i,j])

def filtTerts(cr,terInts):
    if [int(cr[0]), int(cr[1])] in terInts or [cr[1], cr[0]] in terInts:
        return cr
    else:
        return
    
def filtCorrsFunc(corrF,terInts):
    filtCorrs = [] #will store correlations and distances for residues i and j
    for i in range(2):
        corrF.readline() #skip first 2 lines because they are headers
    for line in corrF:
        coor = line.split() #split line into list
        filt = filtTerts(coor, terInts)
        if filt != None and filt[2] != '0.00e+00':
            filtCorrs.append(filt)
    df = pd.DataFrame(filtCorrs, columns = ['i' , 'j', 'mi', 'depth'])
    df[['mi']] = df[['mi']].apply(pd.to_numeric)
    return df

def topCorrs(table, n): #Returns list of top n correlations
    table = table.sort_values('mi', ascending=False)
    table = table.head(n)
    return table

directory =  'D:\Documents2\Rotations\Weeks\Week1Data\CorrDistAnalysisRNaseP'
filePDB = '00952V_n2_sup.pdb' #pdb with xyz coordinates
ctF = 'rnpb.ct'
fileCorr = 'rnasep.1.corr.txt' #RING output with residues and correlations
pdbF = open(filePDB) #open pdb file into memory
corrF = open(fileCorr) #same for RING output
filtCorrs = filtCorrsFunc(corrF, terInts)    
topList = topCorrs(filtCorrs, 50)
np.savetxt(r'D:\Documents2\Rotations\Weeks\Week3MIAPC\filtPHICorr_Top50_CD20.txt', topList, fmt='%s')

#terInts = getTerts(ctF)