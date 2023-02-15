import numpy as np
from pythonUtilities.Segment2DClass import *

class Boundary2D():

    def __init__(self,segList=[]):
        self.segList = segList
        
    def addSegment(self,segType,keyPts,numEls):
        self.segList.append(Segment2D(segType,keyPts,numEls))
        
    def getBoundaryMesh(self):
        allNds = list()
        allEds = list()
        totNds = 0
        for seg in self.segList:
            segMesh = seg.getNodesEdges()
            allNds.extend(segMesh['nodes'])
            allEds.extend(segMesh['edges'] + totNds)
            totNds = len(allNds)
        allNds = np.array(allNds)
        allEds = np.array(allEds)
        ndElim = -np.ones(totNds,dtype=int)
        ndNewInd = -np.ones(totNds,dtype=int)
        for n1i in range(0,totNds):
            if(ndElim[n1i] == -1):
                for n2i in range(n1i+1,totNds):
                    if(ndElim[n2i] == -1):
                        proj = allNds[n2i] - allNds[n1i]
                        dist = np.linalg.norm(proj)
                        if(dist < 1e-12):
                            ndElim[n2i] = n1i
        ndi = 0
        nodesFinal = list()
        for n1i in range(0,totNds):
            if(ndElim[n1i] == -1):
                nodesFinal.append(allNds[n1i])
                ndNewInd[n1i] = ndi
                ndi = ndi + 1
        nodesFinal = np.array(nodesFinal)
        for edi in range(0,len(allEds)):
            for j in range(0,2):
                nd = allEds[edi,j]
                if(ndElim[nd] == -1):
                    allEds[edi,j] = ndNewInd[nd]
                else:
                    allEds[edi,j] = ndNewInd[ndElim[nd]]
        output = dict()
        output['nodes'] = nodesFinal
        output['edges'] = allEds
        return output