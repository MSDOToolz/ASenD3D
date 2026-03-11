import numpy as np
import asendUtils.meshing.MeshTools as mt
from asendUtils.meshing.ShellRegion import *

class Surface():

    def __init__(self,regionList=None,regionNames=None,meshList=None,meshNames=None):
        self.shellRegions = list()
        if(regionList != None):
            self.shellRegions.extend(regionList)
        self.regionNames = list()
        if(regionNames != None):
            self.regionNames.extend(regionNames)
        self.meshes = list()
        if(meshList != None):
            self.meshes.extend(meshList)
        self.meshNames = list()
        if(meshNames != None):
            self.meshNames.extend(meshNames)
        
    def addShellRegion(self,regType,keyPts,numEls,name=None,natSpaceCrds=None,elType='quad',meshMethod='free'):
        self.shellRegions.append(ShellRegion(regType,keyPts,numEls,natSpaceCrds,elType,meshMethod))
        if(name == None):
            numReg = len(self.shellRegions)
            regName = 'Sub-Region_' + str(numReg)
            self.regionNames.append(regName)
        else:
            self.regionNames.append(name)
        
    def addMesh(self,meshData,name=None):
        self.meshes.append(meshData)
        self.meshNames.append(name)
        # if(name == None):
        #     numMsh = len(self.meshes)
        #     meshName = 'Sub-Mesh_' + str(numMsh)
        #     self.meshNames.append(meshName)
        # else:
        #     self.meshNames.append(name)
        
    def getSurfaceMesh(self):
        allNds = list()
        allEls = list()
        elSetList = dict()
        numNds = 0
        numEls = 0
        for regi, reg in enumerate(self.shellRegions):
            regMesh = reg.createShellMesh()
            setList = list()
            eli = 0
            for el in regMesh['elements']:
                for i in range(0,4):
                    if(el[i] != -1):
                        el[i] = el[i] + numNds
                allEls.append(el)
                setList.append(int(eli + numEls))
                eli = eli + 1
            elSetList[self.regionNames[regi]] = setList
            allNds.extend(regMesh['nodes'])
            numNds = len(allNds)
            numEls = len(allEls)
        for mshi, msh in enumerate(self.meshes):
            setList = list()
            for eli, el in enumerate(msh['elements']):
                newEl = -1*np.ones(4,dtype=int)
                for i in range(0,4):
                    if(el[i] != -1):
                        newEl[i] = el[i] + numNds
                allEls.append(newEl)
                setList.append(int(eli + numEls))
            nm = self.meshNames[mshi]
            if nm != None:
                if nm in elSetList:
                    elSetList[nm].extend(setList)
                else:
                    elSetList[nm] = setList
            try:
                for es in msh['sets']['element']:
                    if es in elSetList:
                        for el in msh['sets']['element'][es]:
                            elSetList[es].append(el + numEls)
                    else:
                        setList = list()
                        for el in msh['sets']['element'][es]:
                            setList.append(el + numEls)
                        elSetList[es] = setList
            except:
                pass
            allNds.extend(msh['nodes'])
            numNds = len(allNds)
            numEls = len(allEls)
        mData = dict()
        mData['nodes'] = np.array(allNds)
        mData['elements'] = np.array(allEls)
        mData = mt.mergeDuplicateNodes(mData)
        mData['sets'] = dict()
        mData['sets']['element'] = elSetList
        mData = mt.getAllMatchingNodeSets(mData)
        return mData