import numpy as np
import os
import yaml
from asendUtils.model.Section import Section
from asendUtils.model.Material import Material
from asendUtils.model.Constraint import Constraint
from asendUtils.syst.pathTools import *

class Model():

    def __init__(self):
        self.fileName = ""
        self.modelData = dict()
        self.modelData['nodes'] = list()
        self.modelData['elements'] = list()
        self.modelData['sets'] = dict()
        self.modelData['sets']['node'] = list()
        self.modelData['sets']['element'] = list()
        self.modelData['sections'] = list()
        self.modelData['materials'] = list()
        self.totNds = 0
        self.totEls = 0
        
    def addMeshData(self,meshData,meshType='solid'):
        nodes = meshData['nodes']
        allEls = meshData['elements']
        nNds = self.totNds
        nEls = self.totEls
        nInpNds = len(nodes)
        nInpEls = len(allEls)
        newNds = list()
        for i in range(0,nInpNds):
            nd = list()
            nd.append(nNds+i)
            nd.extend(nodes[i])
            newNds.append(str(nd))
        self.modelData['nodes'].extend(newNds)
        if(meshType == 'solid'):
            hexList = list()
            wedList = list()
            tetList = list()
            for i in range(0,nInpEls):
                if(allEls[i,4] == -1):
                    el = list()
                    el.append(nEls+i)
                    el.extend(allEls[i,0:4]+nNds)
                    tetList.append(str(el))
                elif(allEls[i,6] == -1):
                    el = list()
                    el.append(nEls+i)
                    el.extend(allEls[i,0:6]+nNds)
                    wedList.append(str(el))
                else:
                    el = list()
                    el.append(nEls+i)
                    el.extend(allEls[i,0:8]+nNds)
                    hexList.append(str(el))
            if(len(tetList) > 0):
                tetDic = dict()
                tetDic['type'] = 'tet4'
                tetDic['connectivity'] = tetList
                self.modelData['elements'].append(tetDic)
            if(len(wedList) > 0):
                wedDic = dict()
                wedDic['type'] = 'wedge6'
                wedDic['connectivity'] = wedList
                self.modelData['elements'].append(wedDic)
            if(len(hexList) > 0):
                hexDic = dict()
                hexDic['type'] = 'brickIM'
                hexDic['connectivity'] = hexList
                self.modelData['elements'].append(hexDic)
        elif(meshType == 'shell'):
            quadList = list()
            triList = list()
            for i in range(0,nInpEls):
                if(allEls[i,3] == -1):
                    el = list()
                    el.append(nEls+i)
                    el.extend(allEls[i,0:3]+nNds)
                    triList.append(str(el))
                else:
                    el = list()
                    el.append(nEls+i)
                    el.extend(allEls[i,0:4]+nNds)
                    quadList.append(str(el))
            if(len(triList) > 0):
                triDic = dict()
                triDic['type'] = 'shell3'
                triDic['connectivity'] = triList
                self.modelData['elements'].append(triDic)
            if(len(quadList) > 0):
                quadDic = dict()
                quadDic['type'] = 'shell4'
                quadDic['connectivity'] = quadList
                self.modelData['elements'].append(quadDic)
        elif(meshType == 'beam'):
            beamList = list()
            for i in range(0,nInpEls):
                el = list()
                el.append(nEls+i)
                el.extend(allEls[i,0:2]+nNds)
                beamList.append(str(el))
            beamDic = dict()
            beamDic['type'] = 'beam2'
            beamDic['connectivity'] = beamList
            self.modelData['elements'].append(beamDic)
        try:
            inpSets = meshData['sets']
            try:
                nodeSets = inpSets['node']
                newNdSets = list()
                for ns in nodeSets:
                    newLabs = list()
                    for nd in ns['labels']:
                        newLabs.append(nd+nNds)
                    ns['labels'] = newLabs
                    newNdSets.append(ns)
                self.modelData['sets']['node'].extend(newNdSets)
            except:
                pass
            try:
                elSets = inpSets['element']
                newElSets = list()
                for es in elSets:
                    newLabs = list()
                    for el in es['labels']:
                        newLabs.append(el+nEls)
                    es['labels'] = newLabs
                    newElSets.append(es)
                self.modelData['sets']['element'].extend(newElSets)
            except:
                pass
        except:
            pass
        self.totNds = nNds + nInpNds
        self.totEls = nEls + nInpEls
        
    def addNodeSet(self,name,labelList):
        newSet = dict()
        newSet['name'] = name
        newSet['labels'] = labelList
        self.modelData['sets']['node'].append(newSet)
    
    def addElementSet(self,name,labelList):
        newSet = dict()
        newSet['name'] = name
        newSet['labels'] = labelList
        self.modelData['sets']['element'].append(newSet)
            
    def addSection(self,newSection):
        self.modelData['sections'].append(newSection.secData)
            
    def addMaterial(self,newMaterial):
        self.modelData['materials'].append(newMaterial.matData)
        
    def addConstraint(self,newConstraint):
        try:
            self.modelData['constraints'].append(newConstraint.constData)
        except:
            constraints = list()
            constraints.append(newConstraint.constData)
            self.modelData['constraints'] = constraints
            
    def addAnyLoad(self,newLd):
       try:
           self.modelData['loads'].append(newLd)
       except:
           loads = list()
           loads.append(newLd)
           self.modelData['loads'] = loads       
            
    def addNodalForce(self,nodeSet,F,M,stTime=0.0,endTime=1e+100):
        newLd = dict()
        newLd['type'] = 'nodalForce'
        newLd['activeTime'] = str([stTime,endTime])
        ld = list()
        ld.extend(F)
        ld.extend(M)
        newLd['load'] = str(ld)
        newLd['nodeSet'] = nodeSet
        self.addAnyLoad(newLd) 
            
    def addBodyForce(self,elementSet,F,M,stTime=0.0,endTime=1e+100):
        newLd = dict()
        newLd['type'] = 'bodyForce'
        newLd['activeTime'] = str([stTime,endTime])
        ld = list()
        ld.extend(F)
        ld.extend(M)
        newLd['load'] = str(ld)
        newLd['elementSet'] = elementSet
        self.addAnyLoad(newLd)    
            
    def addGravityForce(self,elementSet,G,stTime=0.0,endTime=1e+100):
        newLd = dict()
        newLd['type'] = 'gravitational'
        newLd['activeTime'] = str([stTime,endTime])
        newLd['load'] = str(G)
        newLd['elementSet'] = elementSet
        self.addAnyLoad(newLd)
        
    def addCentifugalForce(self,elementSet,center,axis,angularVelocity,stTime=0.0,endTime=1e+100):
        newLd = dict()
        newLd['type'] = 'centrifugal'
        newLd['activeTime'] = str([stTime,endTime])
        newLd['center'] = str(center)
        newLd['axis'] = str(axis)
        newLd['angularVelocity'] = angularVelocity
        newLd['elementSet'] = elementSet
        self.addAnyLoad(newLd)
        
    def addSurfaceTraction(self,elementSet,T,normDir,normTol=5.0,stTime=0.0,endTime=1e+100):
        newLd = dict()
        newLd['type'] = 'surfaceTraction'
        newLd['activeTime'] = str([stTime,endTime])
        newLd['normDir'] = str(normDir)
        newLd['normTolerance'] = normTol
        newLd['load'] = str(T)
        newLd['elementSet'] = elementSet
        self.addAnyLoad(newLd)
        
    def addSurfacePressure(self,elementSet,P,normDir,normTol=5.0,stTime=0.0,endTime=1e+100):
        newLd = dict()
        newLd['type'] = 'surfacePressure'
        newLd['activeTime'] = str([stTime,endTime])
        newLd['normDir'] = str(normDir)
        newLd['normTolerance'] = normTol
        newLd['load'] = P
        newLd['elementSet'] = elementSet
        self.addAnyLoad(newLd)
        
    def addNodalHeatGen(self,nodeSet,stTime=0.0,endTime=1e+100,Q=0.0):
        newLd = dict()
        newLd['type'] = 'nodalHeatGen'
        newLd['activeTime'] = str([stTime,endTime])
        newLd['load'] = Q
        newLd['nodeSet'] = nodeSet
        self.addAnyLoad(newLd)
        
    def addBodyHeatGen(self,elementSet,stTime=0.0,endTime=1e+100,specQ=0.0):
        newLd = dict()
        newLd['type'] = 'bodyHeatGen'
        newLd['activeTime'] = str([stTime,endTime])
        newLd['load'] = specQ
        newLd['elementSet'] = elementSet
        self.addAnyLoad(newLd)
        
    def addSurfaceFlux(self,elementSet,flux,normDir,normTol=5.0,stTime=0.0,endTime=1e+100):
        newLd = dict()
        newLd['type'] = 'surfaceFlux'
        newLd['activeTime'] = str([stTime,endTime])
        newLd['normDir'] = str(normDir)
        newLd['normTolerance'] = normTol
        newLd['load'] = flux
        newLd['elementSet'] = elementSet
        self.addAnyLoad(newLd)
        
    def addInitialState(self,field,state):
        allFields = 'displacement velocity acceleration temperature tdot'
        if(field not in allFields):
            errstr = 'Error: ' + field + ' is not a currently supported solution variable for initial state. Valid fields: ' + allFields
            raise TypeError(errstr)
        strState = list()
        for s in state:
            strState.append(str(s))
        try:
            self.modelData['initialState'][field] = strState
        except:
            initialState = dict()
            initialState[field] = strState
            self.modelData['initialState'] = initialState
        
    def writeModelInput(self,fileName):
        self.fileName = makeAbsolute(fileName)

        outFile = open('temp.yaml','w')
        yaml.dump(self.modelData,stream=outFile,sort_keys=False)
        outFile.close()
        
        inFile = open('temp.yaml','r')
        outFile = open(fileName,'w')

        fLine = inFile.readline()
        while(fLine != ''):
            newSt = fLine.replace("'","")
            newSt = newSt.replace('"','')
            outFile.write(newSt)
            fLine = inFile.readline()

        inFile.close()
        outFile.close()
        
        os.remove('temp.yaml')