import numpy as np
import yaml

fileName = 'factoredLT.yaml'

inFile = open(fileName,'r')
matDat = yaml.safe_load(inFile)
inFile.close()

maxRow = 0
maxCol = 0
for e in matDat['matrix']:
    if(e[0] > maxRow):
        maxRow = e[0]
    if(e[1] > maxCol):
        maxCol = e[1]

lMat = np.zeros((maxRow+1,maxCol+1))
dMat = np.zeros((maxRow+1,maxCol+1))

for i in range(0,maxRow):
    lMat[i,i] = 1.0
    
for e in matDat['matrix']:
    r = e[0]
    c = e[1]
    if(r == c):
        dMat[r,c] = e[2]
    else:
        lMat[r,c] = e[2]

LD = np.matmul(lMat,dMat)
LDL = np.matmul(LD,np.transpose(lMat))

prodList = list();

for r in range(0,maxRow):
    for c in range(0,maxRow):
        ent = str([r,c,float(LDL[r,c])])
        prodList.append(ent)

outDat = dict()
outDat['matrix'] = prodList

outFile = open('ldlProd.yaml','w')
yaml.dump(outDat,stream=outFile,sort_keys=False)
outFile.close()
