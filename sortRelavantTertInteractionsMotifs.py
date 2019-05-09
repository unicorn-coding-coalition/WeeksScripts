# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 14:05:50 2018

@author: nlama
"""
import numpy as np
import os

p18 = list(range(304,328))
p13 = list(range(184,203))
p14 = list(range(204,225))
p12 = list(range(142,177))
p9 = list(range(107,120))
p8 = list(range(92,107))
p6 = [82,83,84,85,279,278,277,276]
p4 = [66, 67, 68, 69, 70, 71, 72, 73, 353, 
      354, 355, 356, 357, 358, 359, 360]
p1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
      363, 364, 365, 366, 367, 368, 369, 370, 
      371, 372, 373, 374, 375, 376, 377]
knownDoms = p18 + p13 + p14 + p12 + p9 + p8 + p4 + p1

singleStrandedNTs = [62, 63,22, 59, 343, 344, 345, 346, 
                     347, 348, 349, 350, 351, 352, 328, 
                     329, 330, 331, 332, 333, 334, 335,
                     225, 226, 227, 228, 229, 230, 231, 
                     232, 233, 234, 235, 236, 237, 238,
                     275, 247, 248, 249, 300, 301, 302, 
                     303, 291, 292, 293, 294, 295, 266,
                     267, 268, 269, 254, 255, 256, 257, 
                     258, 259, 260, 177, 178, 179, 180, 
                     181, 182, 183, 184, 185, 127, 128, 
                     129, 130, 131, 132, 133, 134, 135, 
                     136, 137, 138, 139, 140, 141, 122,
                     123, 124]

p3 = [39, 40, 41, 42, 33, 52]

##Single Strands on Known Tert Helices
pS12 = [147, 148, 157, 158, 159, 160, 166, 173, 174]
pS13 = [191, 192, 193, 194, 195, 196]
pS14 = [212, 213, 214, 215, 219]
pS9 = [111, 113, 114, 118, 119, 200]
pS8 = [97, 98, 99, 100, 101]
pS18 = [314, 315, 316, 317]

################# BULGES ######################################################
B1 = [21,22,23]
B2 = [58,59,60]
B3 = [121,122,123,124,125]
B4 = [146,147,148,149]
B6 = list(range(265,271))
################ Single Strands Only ##########################################
S1 = list(range(61,65))
S2 = [127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 
      137, 138, 139, 140, 141]
S3 = [177, 178, 179, 180, 181, 182, 183, 184, 185]
S4 = [225, 226, 227, 228, 229, 230, 231,232, 233, 234, 
      235, 236, 237, 238]
S5 = list(range(253,262))
S6 = list(range(290,297))
S7 = list(range(300,305))
S8 = list(range(327,334))
S9 = [343,344,345,346,347,348,349,350,351]


SS = S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8 + S9
allSingleStrands = singleStrandedNTs + pS12 + pS13 + pS14 + pS9 + pS8 + pS18
pSAll = pS12 + pS13 + pS14 + pS9 + pS8 + pS18
bulges = B1 + B2 + B3 + B4 + B6


def getCorrF(directory, f):
    corrF = []
    for i in f:
        fileCorr = os.path.join(directory,i)
        corrF.append(open(fileCorr)) #same for RING output 
    return corrF

def filtTertCorr(line, foundNts):
    if [int(line[0]), int(line[1])] in foundNts or [int(line[1]), int(line[0])] in foundNts:
        return line
    else:
        line = []
        return line

def readCorrFilts(corrF, foundNts):
    tables = []
    for count, line in enumerate(corrF):
        line = line.split()
        print(count)
        print([line[0],line[1]],[int(line[0]),int(line[1])] in foundNts)
        filtLine = filtTertCorr(line, foundNts) 
        if filtLine != []:
            tables.append(filtLine)
        else:
            print('unrecognized')
    return tables
                


def tertRNAMotifsList(prepTable):
    foundNts = []
    tertMotifs = [0,0,0,0,0]
#    print(prepTable[0])
    for i, j in zip(prepTable['i'], prepTable['j']):
#        if [i,j] in already_checked:
#            continue
        if i in pSAll and j in pSAll:
            tertMotifs[0] += 1
            foundNts.append([i,j])
        if j in pSAll and i in pSAll:
            tertMotifs[0] += 1
            foundNts.append([i,j])
        if i in pSAll and j in bulges:
            tertMotifs[1] += 1
            foundNts.append([i,j])
        if j in pSAll and i in bulges:
            tertMotifs[1] += 1
            foundNts.append([i,j])
        if i in p4 and j in SS:
            tertMotifs[2] += 1
            foundNts.append([i,j])
        if j in p4 and i in SS:
            tertMotifs[2] += 1
            foundNts.append([i,j])
        if i in SS and j in pSAll:
            tertMotifs[3] += 1
            foundNts.append([i,j])
        if j in SS and i in pSAll:
            tertMotifs[3] += 1
            foundNts.append([i,j])
        if i in SS and j in bulges:
            tertMotifs[4] += 1
            foundNts.append([i,j])
#            print(i,j)
        if j in SS and i in bulges:
            tertMotifs[4] += 1
            foundNts.append([i,j])
#            print(i,j)       
    return foundNts

def getTertsTables(corrTables, fmt):
    tertDictList = []
    if fmt == 'list':
        for i in corrTables:
            tertDictList.append(tertRNAMotifsList(i))
    return tertDictList


fileNames = ['rnasep.mi.apc.top100.CD20_test.txt', 
             'rnasep.mi.top100.CD20_test.txt', 
             'rnasep.phi.top100.CD20_test.txt'] 
corrs = getCorrF('D:\Documents2\Rotations\Weeks\PythonScripts', fileNames)
#corrTables = readCorr(corrs)
#tertDictsList = getTertsTables(corrTables, 'dict')
foundNts = getTertsTables(corrTables, 'list')


def countTertsList(prepTable):
    tertIntsMat = np.array([pS14+pS13, pS8+p4,
                            pS9+p1, pS12+pS13, pS8+pS18, pS8+pS14, p6+pS18
                            ])
    tertMotifs = [0,0,0,0,0,0,0]
    foundNts = []
    for i, j in zip(prepTable['i'], prepTable['j']):
        if i in tertIntsMat[0] and j in tertIntsMat[0]:
            tertMotifs[0] += 1
            foundNts.append([i,j])
        if i in tertIntsMat[1] and j in tertIntsMat[1]:
            tertMotifs[1] += 1
            foundNts.append([i,j])
        if i in tertIntsMat[2] and j in tertIntsMat[2]:
            tertMotifs[2] += 1
            foundNts.append([i,j])
        if i in tertIntsMat[3] and j in tertIntsMat[3]:
            tertMotifs[3] += 1
            foundNts.append([i,j])
        if i in tertIntsMat[4] and j in tertIntsMat[4]:
            tertMotifs[4] += 1
            foundNts.append([i,j])
        if i in tertIntsMat[5] and j in tertIntsMat[5]:
            tertMotifs[5] += 1
            foundNts.append([i,j])
        if i in tertIntsMat[6] and j in tertIntsMat[6]:
            tertMotifs[6] += 1
            foundNts.append([i,j])
#        if i in tertIntsMat[7] and j in tertIntsMat[7]:
#            tertMotifs[7] += 1
#        if i in tertIntsMat[8] and j in tertIntsMat[8]:
#            tertMotifs[8] += 1
#        if i in tertIntsMat[9] and j in tertIntsMat[9]:
#            tertMotifs[9] += 1
#        if i in tertIntsMat[10] and j in tertIntsMat[10]:
#            tertMotifs['B0-P1'] += 1        
#    print(tertMotifs)
    return foundNts 
#
#
fileNames = ['rnasep.mi.apc.top225.CD20_test.txt', 
             'rnasep.mi.top225.CD20_test.txt', 
             'rnasep.phi.top225.CD20_test.txt'] 
corrF = getCorrF('D:\Documents2\Rotations\Weeks\PythonScripts', fileNames)
foundNts2 = []
for i in corrTables:
    foundNts2.append(countTertsList(i))
      
apcNts = foundNts[0] + foundNts2[0]
miNts = foundNts[1] + foundNts2[1]
phiNts = foundNts[2] + foundNts2[2]

#apcNts = foundNts2[0]
#miNts = foundNts2[1]
#phiNts = foundNts2[2]
#
#apcNts = foundNts[0]
#miNts = foundNts[1]
#phiNts = foundNts[2]

corrTableapc = readCorrFilts(corrF[0], apcNts)
np.savetxt(r'rnasep.mi.apc.tertints.cd20.txt', corrTableapc, fmt = '%s')

corrTablemi = readCorrFilts(corrF[1], miNts)
np.savetxt(r'rnasep.mi.tertints.cd20.txt', corrTablemi, fmt = '%s')

corrTablephi = readCorrFilts(corrF[2], phiNts)
np.savetxt(r'rnasep.phi.tertints.cd20.txt', corrTablephi, fmt = '%s')

for i,j in zip(foundNts[0], foundNts2[0]):
    if i == j:
        print(i)
        
        
