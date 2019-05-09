import numpy as np
import h5py

"""Generates jaccard index for all partitions in hdf5 graph files"""

#############################FUNCTIONS####################################
    
def jaccard(list1, list2):

    intersection = len(set(list1).intersection(list2))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection / union)


def getComDict(array):
    comDict={}
    for pos, comId in enumerate(array):
        try:
            comDict[comId].append(pos)
        except:
            comDict[comId]=[pos]
    return comDict
    
                    
def getJaccard(d1,d2):
   
    jIdx=[]
    jIdxTempy=[]
    for i in d1.keys():
        for j in d2.keys():
            jIdxTempy.append(jaccard(d1[i],d2[j]))
        jIdx.append(max(jIdxTempy))
    return np.mean(jIdx)

def findSim(part1,part2):
    progress=0
    jaccardIndPerPartition=[]
    for p1 in part1:
        #print(len(p1))
        progress+=1
        p1Dict=getComDict(p1)
        jac=[]
        for p2 in part2:
            #print(len(p2))
            p2Dict=getComDict(p2)
            jac.append(getJaccard(p1Dict,p2Dict))
        jaccardIndPerPartition.append(max(jac))
        print("{0} out of {1} partitions left".format(progress,len(part1)))
#        if progress == 10:
#            break
    return jaccardIndPerPartition
 

############################################DATA##########################
#
#ensemble1=h5py.File("D:/Weeks/Data/NetworkAnalysis/cf_rnpB2_ensemble.fin2.numruns200.hdf5","r")
#ensemble2=h5py.File("D:/Weeks/Data/NetworkAnalysis/cf_rnpB_ensemble.fin2.numruns200.hdf5","r")

ensemble1=h5py.File("D:/Weeks/tempyNicole/networkAnalysis/in4x.tmRNA.40.graphml.mincom50.hdf5","r")
ensemble2=h5py.File("D:/Weeks/tempyNicole/networkAnalysis/in1x.tmRNA.40.graphml.mincom50.hdf5","r")


partitions1=ensemble1['_partitions']
partitions2=ensemble2['_partitions']


a=findSim(partitions1,partitions2)
print("Median Jaccard Index: ", np.median(a))
print("Mean Jaccard Index: ", np.mean(a))



