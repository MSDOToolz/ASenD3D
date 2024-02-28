import numpy as np
import asendUtils.meshing.MeshTools as mt
from asendUtils.meshing.Segment2D import *

class Boundary2D():

    def __init__(self,segList=None,segNames=None,meshList=None,meshNames=None):
        self.segList = list()
        if(segList != None):
            self.segList.extend(segList)
        self.segNames = list()
        if(segNames != None):
            self.segNames.extend(segNames)
        self.meshList = list()
        if(meshList != None):
            self.meshList.extend(meshList)
        self.meshNames = list()
        if(meshNames != None):
            self.meshNames.extend(meshNames)
        
    def addSegment(self,segType,keyPts,numEls,name=None):
        self.segList.append(Segment2D(segType,keyPts,numEls))
        if(name == None):
            numSeg = len(self.segList)
            segName = 'segment_' + str(numSeg)
            self.segNames.append(segName)
        else:
            self.segNames.append(name)
        
    def addMesh(self,meshData,name=None):
        self.meshList.append(meshData)
        if(name == None):
            numMsh = len(self.meshList)
            meshName = 'Sub-Mesh_' + str(numMsh)
            self.meshNames.append(meshName)
        else:
            self.meshNames.append(name)
        
    def getBoundaryMesh(self):
        allNds = list()
        allEds = list()
        elSets = list()
        totNds = 0
        totEls = 0
        for i, seg in enumerate(self.segList):
            segMesh = seg.getNodesEdges()
            allNds.extend(segMesh['nodes'])
            allEds.extend(segMesh['edges'] + totNds)
            newSet = dict()
            newSet['name'] = self.segNames[i]
            totElsNext = totEls + len(segMesh['edges'])
            newSet['labels'] = list(range(totEls,totElsNext))
            elSets.append(newSet)
            totNds = len(allNds)
            totEls = totElsNext
        for i, mesh in enumerate(self.meshList):
            allNds.extend(mesh['nodes'])
            allEds.extend(mesh['edges'] + totNds)
            newSet = dict()
            newSet['name'] = self.meshNames[i]
            totElsNext = totEls + len(mesh['edges'])
            newSet['labels'] = list(range(totEls,totElsNext))
            elSets.append(newSet)
            totNds = len(allNds)
            totEls = totElsNext   
        allNds = np.array(allNds)
        allEds = np.array(allEds)
        
        meshData = dict()
        meshData['nodes'] = allNds
        meshData['elements'] = allEds
        meshData['sets'] = dict()
        meshData['sets']['element'] = elSets
        
        output = mt.mergeDuplicateNodes(meshData)
        
        return output