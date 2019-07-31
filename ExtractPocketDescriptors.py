# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 14:34:41 2019
Go through fpocket output and extract all pocket descriptors
@author: nlama
"""
from ast import literal_eval
import pandas as pd
from glob import iglob
import os

os.chdir("D:/Weeks/PocketFinding/dataSetA")

desc=["PDB","Pocket","Score ","Number of Alpha Spheres ","Polarity score",
      "Mean local hydrophobic density ","Hydrophobicity score",
      "Volume score","Apolar alpha sphere proportion ","Alpha sphere density ",
      "Volume ", "Proportion of polar atoms", "Druggability Score "]

pockDescFiles = [f for f in iglob("Analyze/Final/fpocketOutput.Riboswitches.M.4.i.20.m.3.D.1.s.4.r.3.n.1/*/*Clean_info.*", 
                                  recursive=True) if os.path.isfile(f)]
df = pd.read_csv("dataSetARiboswitches.csv")
all = {}
pDesDict = {}
for file in pockDescFiles:
    pdbName = os.path.basename(file)[0:4]
    ma = df['PDB'] == pdbName
    kp=df[ma]['pocketRank']
    for d in desc:
        pDesDict.setdefault(d, []) 
    
    with open(file) as f:
        for line in f:
            l=line.split(":")
            if "Pocket" in l[0]:
                pnum=l[0].split()[1]
                pDesDict[l[0].split()[0]].append(float(pnum))
                pDesDict['PDB'].append(pdbName)
            else:
                try:
                    pDesDict[l[0].strip('\t')].append(float(l[1].strip('\t').strip("\n").lstrip()))
                except:
                    continue
                
with open('defaultFpocketDescriptors_FinalFilter.csv','w') as outF:
    outF.write("PDB,Pocket,Volume,Proportion of Polar Atoms,Druggability Score,Hydrophobicity Score,Mean Local Hydrophobic Density,Number of Alpha Spheres,Polarity Score,Proportion of Apolar Alpha Sphere,Density of Alpha Sphere,Score\n")
    for x in range(len(pDesDict['PDB'])):
        outF.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n".format(pDesDict['PDB'][x],
                   pDesDict['Pocket'][x],pDesDict['Volume '][x],pDesDict['Proportion of polar atoms'][x],
                   pDesDict['Druggability Score '][x],pDesDict['Hydrophobicity score'][x],
                   pDesDict['Mean local hydrophobic density '][x],pDesDict['Number of Alpha Spheres '][x],
                   pDesDict['Polarity score'][x],pDesDict['Apolar alpha sphere proportion '][x],
                   pDesDict['Alpha sphere density '][x],pDesDict['Score '][x]))
    

                
        