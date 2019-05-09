import pandas as pd
import os
import random as rand
import numpy as np
import csv
import matplotlib.pyplot as plot

os.chdir("D:/Weeks/Data/NetworkAnalysis/reWeighting")

degDict={}
edgeList=pd.read_csv("cf.rnaseP.apc.win3.edeList.csv") #contains weights
#nodeList=pd.read_csv("cf.rnaseP.apc.win3.nodelist.csv") #contains degrees
#self.fileName = os.path.splitext(os.path.basename(fileCorr))[0] #Get file name without extention


with open('cf.rnaseP.apc.win3.nodeList.csv', mode='r') as infile:
    reader = csv.reader(infile)
    next(reader, None)  # skip the headers
    degDict = {int(rows[0]):float(rows[3]) for rows in reader}


meanDeg=np.mean(list(degDict.values()))

newWeight=[]


for i,j,w in zip(edgeList["Source"],edgeList["Target"],edgeList["Weight"]):
    if degDict[i] > 26 or degDict[j] >26:
        newWeight.append(0)
    elif  np.mean((degDict[i],degDict[j])) >= meanDeg:
        newWeight.append(w/4)
    else:
        newWeight.append(w)
        
        
    
    
se = pd.Series(newWeight)
edgeList['New_Weight'] = se.values

edgeList.to_csv("output.strict.deg26.csv", index=False)


