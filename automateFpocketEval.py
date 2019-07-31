# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 13:33:18 2019

@author: nlama
"""
from ast import literal_eval
import pandas as pd
from glob import iglob
import math
import os

#Change directory to where all the data is
os.chdir("D://Weeks//PocketFinding//dataSetA")

#Functions
def calcDistFromList(x1,y1,z1,L):
    inL=0
    for nt,d in L.items():
        x = (x1 - d[0])**2 #compute dist form x 
        y = (y1 - d[1])**2 #compute dist form y
        z = (z1 - d[2])**2 #compute dist form z
        distance = math.sqrt(x + y + z) # compute distance between dues
        if distance <= 6:
            inL+=1
            if inL >= 2:
                return True


def getPDBCoords(f):  
    pdbCoords = []
    f=open(f)
    nt = 1
    for row in f:
        if "STP" in row[17:21]:
            nt=int(row[22:26])
            x2=float(row[30:38])
            y2=float(row[38:46])
            z2=float(row[46:54])
            pdbCoords.append([nt,x2,y2,z2])           
    f.close()  
    return pdbCoords

def getListPDBCoords(f,L):  
    pdbCoords = {}
    f=open(f)
    nt = 1
    for row in f:
        if int(row[22:26]) in L:
            nt=int(row[22:26])
            x2=float(row[30:38])
            y2=float(row[38:46])
            z2=float(row[46:54])
            pdbCoords[nt]=[x2,y2,z2]       
    f.close()  
    return pdbCoords

#Global Variables
wL2NZ4 = [2,3,31,32,33,34,35,42,43,44,45,57,58,59,60,61,83,84,85,86,125,126,127,128,4,28,52,54,55,5,27,51]
bL2NZ4=[17,72,73,136,137,138,139,140,141,89,90,91,92,93,129,130,131,132,133,71,74,75,76]
numPockets_2NZ4 = 15
wL3MXH=[14,15,16,17,18,19,20,21,46,47,48,49,50,91,92,93,95]
bL3MXH=[75,660,661,662,663,664,665,666,667,668,669,9,10,11,12,96,97,98]
numPockets_3MXH = 9

wL = wL2NZ4
bL = bL2NZ4
numPockets = numPockets_2NZ4

pdbs = [f for f in iglob("fpocketOutput.Riboswitches*/*/2NZ4_Clean_out.pdb", recursive=True) if os.path.isfile(f)]
#pdbs = ["fpocketOutput.Riboswitches.M.3.5.i.20.m.1.D.0.5.s.0.5.r.2.5.n.10/3MXH_Clean_out/3MXH_Clean_out.pdb","fpocketOutput.Riboswitches.M.3.5.i.20.m.1.D.0.5.s.0.5.r.2.5.n.3/3MXH_Clean_out/3MXH_Clean_out.pdb"]
stage2 = []
reject = set()

whiteBlackListDict=[]

for L in [wL,bL]:
    whiteBlackListDict.append(getListPDBCoords(pdbs[1],L))

with open("fPocketEval_2NZ4_da.txt",'w') as file:
    for pdb in pdbs:
        inB = False
        inW = False
        p=os.path.basename(pdb)
        opt=os.path.split((os.path.split(pdb))[0])[0]
        pdbCoords=getPDBCoords(pdb)
        for x in pdbCoords:
    #        inW=calcDistFromList(x[1],x[2],x[3], whiteBlackListDict[0])
            inB=calcDistFromList(x[1],x[2],x[3], whiteBlackListDict[1])
            if inB == True:
                reject.add(opt)
                break
        idx = len(pdbCoords)
        if idx < 1 or inB == True:
            reject.add(opt)
        elif int(pdbCoords[idx-1][0])-2 > numPockets:
            reject.add(opt)
        else:
            file.write(opt+"\n")
            stage2.append(opt)
            
#

#l3 = [x for x in stage2 if x not in reject]
#
    









