import numpy as np
import MeshTools as mt
from ShellRegionClass import *

class Surface():

    def __init__(self,regionList=[],regionNames=[],meshList=[],meshNames=[]):
        self.shellRegions = regionList
        self.regionNames = regionNames
        self.meshes = meshList
        self.meshNames = meshNames
        
    def addShellRegion(self,regType,keyPts,numEls,name='',natSpaceCrds=[],elType='quad',meshMethod='free'):
        self.shellRegions.append(ShellRegion(regType,keyPts,numEls,natSpaceCrds,elType,meshMethod))
        if(name == ''):
            numReg = len(self.shellRegions)
            regName = 'Sub-Region_' + str(numReg)
            self.regionNames.append(regName)
        else:
            self.regionNames.append(name)
        
    def addMesh(self,meshData,name=''):
        self.meshes.append(meshData)
        if(name == ''):
            numMsh = len(self.meshes)
            meshName = 'Sub-Mesh_' + str(numMsh)
            self.meshNames.append(meshName)
        else:
            self.meshNames.append(name)
        
    def getSurfaceMesh(self):
        allNds = list()
        allEls = list()
        elSetList = list()
        numNds = 0
        numEls = 0
        regi = 0
        for reg in self.shellRegions:
            regMesh = reg.createShellMesh()
            allNds.extend(regMesh['nodes'])
            setList = list()
            eli = 0
            for el in regMesh['elements']:
                for i in range(0,4):
                    if(el[i] != -1):
                        el[i] = el[i] + numNds
                allEls.append(el)
                setList.append((eli + numEls))
                eli = eli + 1
            thisSet = dict()
            thisSet['name'] = self.regionNames[regi]
            thisSet['labels'] = setList
            elSetList.append(thisSet)
            numNds = len(allNds)
            numEls = len(allEls)
            regi = regi + 1
        mshi = 0
        for msh in self.meshes:
            allNds.extend(msh['nodes'])
            setList = list()
            eli = 0
            for el in msh['elements']:
                for i in range(0,4):
                    if(el[i] != -1):
                        el[i] = el[i] + numNds
                allEls.append(el)
                setList.append((eli + numEls))
                eli = eli + 1
            thisSet = dict()
            thisSet['name'] = self.meshNames[mshi]
            thisSet['labels'] = setList
            elSetList.append(thisSet)
            numNds = len(allNds)
            numEls = len(allEls)
            mshi = mshi + 1
        mData = dict()
        mData['nodes'] = np.array(allNds)
        mData['elements'] = np.array(allEls)
        mData = mt.mergeDuplicateNodes(mData)
        mData['sets'] = dict()
        mData['sets']['element'] = elSetList
        return mData