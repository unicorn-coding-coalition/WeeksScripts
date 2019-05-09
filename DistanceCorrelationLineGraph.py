# -*- coding: utf-8 -*-
"""
Created on Thu May 31 13:15:43 2018

@author: nlama
"""

def getCorrF(directory, f):
    corrF = []
    for i in f:
        fileCorr = os.path.join(directory,i)
        corrF.append(open(fileCorr)) #same for RING output 
    return corrF

def readCorr(corrF):
    tables = []
    for j in range(len(corrF)): 
        for i in range(2):
            corrF[j].readline() #skip first 2 lines because they are headers
    for t in corrF:
        df = pd.read_table(t, delim_whitespace=True, names=('i', 'j', 'mi', 'depth'))
        map(float, df['mi'])
        tables.append(df)
    return tables

def corrPairFunc(pdbF):

        coorPair = {} #will store residue and its coordinates
#        residueNum = 1 #only switch on if pdb is numbered not 1-end
        for line in pdbF:
            lineList = line.split()
            if len(lineList)>5 and lineList[2] == "C1'": #only keep if C1' atom
#                coorPair[str(residueNum)]=lineList[6:9] #Only do this for pdbs with weird numbering! (not 0-end)
                coorPair[lineList[5]]=lineList[6:9] #keep only residue and xyz
#                residueNum += 1
        return coorPair
    
#Look up residues in coorPair for computing distances
def resCorrDistFunc(coorPair, corrF):
    resCorrDist = [] #will store correlations and distances for residues i and j
    for i in range(2):
        corrF.readline() #skip first 2 lines because they are headers
    for line in corrF:
        coor = line.split() #split line into list
        if len(coor) > 1 and coor[2] != '0.00e+00':
            i = coor[0] #extract residue1 in pair
            j = coor[1] #extract residue2 in pair
            iCoords = coorPair[i] #find residue 1's coordinates
            jCoords = coorPair[j] #find residue 2's coordinates
            x = (float(iCoords[0]) - float(jCoords[0]))**2 #compute dist form x 
            y = (float(iCoords[1]) - float(jCoords[1]))**2 #compute dist form y
            z = (float(iCoords[2]) - float(jCoords[2]))**2 #compute dist form z
            distance = math.sqrt(x + y + z) # compute distance between residues
            if float(coor[2]) > cutoff:
                resCorrDist.append([i , j , float(coor[2]), distance])
    if resCorrDist == []:
        raise ValueError("returned empty resCorrDist list")
    else:
        return resCorrDist


        
        
        
        
fileNames = ['rnasep.1.mi.apc.corr.txt', 'rnasep.1.mi.corr.txt', 'rnasep.1.corr.txt'] 
corrs = getCorrF('D:\Documents2\Rotations\Weeks\Week1Data\CorrDistAnalysisRNaseP', fileNames)
corrTables = readCorr(corrs)
#q =[25, 40, 50, 70, 85, 90, 95, 99.9]
#q = [98.5, 98.5]
#cutoffs = getCutoffs(corrTables, q)
#pKnowns = percentCalculator(corrTables, cutoffs)
correlations = []

correlations.append() 

#### Create DataFrame for Plotting #####################################    
m = ['mi + apc',] * len(q), ['mi',] * len(q), ['phi',] * len(q)
m = m[0] + m[1] + m[2] #concatanate all strings for data frame
pTable = pd.DataFrame({'percentile cutoffs': q*3, 'percentage of relevant tertiary interactions': pKnowns, 'metric': m })
genLineGraph(pTable, save = 'yes')