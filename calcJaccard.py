'''
READ!!
Only works if the cluster ID is important. If cluster ID
is randomized for each graph, then do not use this script.
Only use if cluster ID is maintained for each partition. 
'''


import pandas as pd
import numpy as np


fileCorr1=pd.read_csv('D:/Weeks/Data/NetworkAnalysis/resCorrDistFiles/cf.RMRP.Id.Mod.csv')
fileCorr2=pd.read_csv('D:/Weeks/Data/NetworkAnalysis/resCorrDistFiles/in.RMRP.Id.Mod.csv')

comNum = max(fileCorr1['modularity_class']) +1

if comNum != max(fileCorr2['modularity_class']) +1:
    raise ValueError('two networks must have same number of communities')

nodeClustDict1={}
nodeClustDict2={}

totalClust = np.zeros(comNum)

for nt,clus in zip(fileCorr1['Id'], fileCorr1['modularity_class']):
    nodeClustDict1[nt]=clus
    totalClust[clus] += 1
    
for nt,clus in zip(fileCorr2['Id'], fileCorr2['modularity_class']):
    nodeClustDict2[nt]=clus
    totalClust[clus] += 1

union=0.0
total=len(nodeClustDict1)+len(nodeClustDict2)

clustUnion = np.zeros(comNum)

for nt in fileCorr1['Id']:
    try:
        if nodeClustDict1[nt] == nodeClustDict2[nt]:
            union+=1
            clustUnion[nodeClustDict1[nt]] += 1
    except:
        continue
    
jaccardEntire=union/total 
jaccardClust=np.zeros(comNum)

for x in range(comNum):
    jaccardClust[x]=clustUnion[x]/totalClust[x]
    
      





        
