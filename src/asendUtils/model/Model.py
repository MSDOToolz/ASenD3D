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
        self.massElements = list()
        self.forceElements = list()
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
            if(len(allEls[0]) == 8):
                hexList = list()
                wedList = list()
                tetList = list()
                for i, eli in enumerate(allEls):
                    if(eli[4] == -1):
                        el = list()
                        el.append(nEls+i)
                        el.extend(eli[0:4]+nNds)
                        tetList.append(str(el))
                    elif(eli[6] == -1):
                        el = list()
                        el.append(nEls+i)
                        el.extend(eli[0:6]+nNds)
                        wedList.append(str(el))
                    else:
                        el = list()
                        el.append(nEls+i)
                        el.extend(eli[0:8]+nNds)
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
            else:
                tetList = list()
                for i, eli in enumerate(allEls):
                    el = [nEls+i]
                    el.extend(eli[0:10]+nNds)
                    tetList.append(str(el))
                if(len(tetList) > 0):
                    tetDic = dict()
                    tetDic['type'] = 'tet10'
                    tetDic['connectivity'] = tetList
                    self.modelData['elements'].append(tetDic)
                
        elif(meshType == 'shell'):
            quadList = list()
            triList = list()
            for i, eli in enumerate(allEls):
                if(eli[3] == -1):
                    el = list()
                    el.append(nEls+i)
                    el.extend(eli[0:3]+nNds)
                    triList.append(str(el))
                else:
                    el = list()
                    el.append(nEls+i)
                    el.extend(eli[0:4]+nNds)
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
            for i, eli in enumerate(allEls):
                el = list()
                el.append(nEls+i)
                el.extend(eli[0:2]+nNds)
                beamList.append(str(el))
            beamDic = dict()
            beamDic['type'] = 'beam2'
            beamDic['connectivity'] = beamList
            self.modelData['elements'].append(beamDic)
        elif(meshType == 'frcFld'):
            elList = list()
            for i, eRow in enumerate(allEls,nEls):
                el = [i]
                el.extend(eRow)
                elList.append(str(el))
            elDic = dict()
            elDic['type'] = 'frcFld'
            elDic['connectivity'] = elList
            self.modelData['elements'].append(elDic)
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
        try:
            self.massElements.extend(meshData['massElements'])
        except:
            pass
        try:
            self.forceElements.extend(meshData['forceElements'])
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
            
    def fixDisplacement(self,nodeSet,ux=None,uy=None,uz=None,rx=None,ry=None,rz=None):
        if(ux != None):
            newConst = Constraint('displacement')
            newConst.addTerm(nodeSet,1,1.0)
            newConst.setRHS(ux)
            self.addConstraint(newConst)
        if(uy != None):
            newConst = Constraint('displacement')
            newConst.addTerm(nodeSet,2,1.0)
            newConst.setRHS(uy)
            self.addConstraint(newConst)
        if(uz != None):
            newConst = Constraint('displacement')
            newConst.addTerm(nodeSet,3,1.0)
            newConst.setRHS(uz)
            self.addConstraint(newConst)
        if(rx != None):
            newConst = Constraint('displacement')
            newConst.addTerm(nodeSet,4,1.0)
            newConst.setRHS(rx)
            self.addConstraint(newConst)
        if(ry != None):
            newConst = Constraint('displacement')
            newConst.addTerm(nodeSet,5,1.0)
            newConst.setRHS(ry)
            self.addConstraint(newConst)
        if(rz != None):
            newConst = Constraint('displacement')
            newConst.addTerm(nodeSet,6,1.0)
            newConst.setRHS(rz)
            self.addConstraint(newConst)
            
    def fixTemperature(self,nodeSet,T):
        newConst = Constraint('temperature')
        newConst.addTerm(nodeSet,1,1.0)
        newConst.setRHS(T)
        self.addConstraint(newConst)
        
    def periodicDisplacement(self,ndSets=None):
        if(ndSets is None):
            sN = ['periodicXMin','periodicXMax',
                  'periodicYMin','periodicYMax',
                  'periodicZMin','periodicZMax',
                  'xMinRef','xMaxRef',
                  'yMinRef','yMaxRef',
                  'zMinRef','zMaxRef']
        else:
            sN = ndSets
        for i in range(0,3):
            dof = i + 1
            # x periodic constraint
            newConst = Constraint('displacement')
            newConst.addTerm(sN[1],dof,1.0)
            newConst.addTerm(sN[6],dof,1.0)
            newConst.addTerm(sN[0],dof,-1.0)
            newConst.addTerm(sN[7],dof,-1.0)
            newConst.setRHS(0.0)
            self.addConstraint(newConst)
            # y periodic constraint
            newConst = Constraint('displacement')
            newConst.addTerm(sN[3],dof,1.0)
            newConst.addTerm(sN[8],dof,1.0)
            newConst.addTerm(sN[2],dof,-1.0)
            newConst.addTerm(sN[9],dof,-1.0)
            newConst.setRHS(0.0)
            self.addConstraint(newConst)
            # z periodic constraint
            newConst = Constraint('displacement')
            newConst.addTerm(sN[5],dof,1.0)
            newConst.addTerm(sN[10],dof,1.0)
            newConst.addTerm(sN[4],dof,-1.0)
            newConst.addTerm(sN[11],dof,-1.0)
            newConst.setRHS(0.0)
            self.addConstraint(newConst)
            
            
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
        
    def addCentrifugalForce(self,elementSet,center,axis,angularVelocity,stTime=0.0,endTime=1e+100):
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
    
    def integrateMassElements(self):
        for me in self.massElements:
            setLabs = list()
            for ns in self.modelData['sets']['node']:
                if(ns['name'] == me[0]):
                    setLabs = ns['labels']
            ei = self.totEls
            newESLabs = list()
            newEls = list()
            for nd in setLabs:
                newEl = [ei,nd]
                newEls.append(str(newEl))
                newESLabs.append(ei)
                ei = ei + 1
            eList = dict()
            eList['type'] = 'mass'
            eList['connectivity'] = newEls
            self.modelData['elements'].append(eList)
            eSet = dict()
            eSet['name'] = me[1]
            eSet['labels'] = newESLabs
            self.modelData['sets']['element'].append(eSet)
            self.totEls = ei
        self.massElements = list()

    def integrateForceElements(self):
        for fe in self.forceElements:
            set1Labs = list()
            set2Labs = list()
            for ns in self.modelData['sets']['node']:
                if(ns['name'] == fe[0]):
                    set1Labs = ns['labels']
                if(ns['name'] == fe[1]):
                    set2Labs = ns['labels']
            ei = self.totEls
            newEls = list()
            newES = list()
            for s1 in set1Labs:
                for s2 in set2Labs:
                    if(s2 != s1):
                        newEl = [ei,s1,s2]
                        newEls.append(str(newEl))
                        newES.append(ei)
                        ei = ei + 1
            eList = dict()
            eList['type'] = 'frcFld'
            eList['connectivity'] = newEls
            self.modelData['elements'].append(eList)
            eSet = dict()
            eSet['name'] = fe[2]
            eSet['labels'] = newES
            self.modelData['sets']['element'].append(eSet)
            self.totEls = ei
        self.forceElements = list()
    
    def writeModelInput(self,fileName):
        self.fileName = makeAbsolute(fileName)

        self.integrateForceElements()
        self.integrateMassElements()
        
        fileStr = yaml.dump(self.modelData,width=200,sort_keys=False)
        
        fileStr = fileStr.replace("'","")
        fileStr = fileStr.replace('"','')
        
        outFile = open(self.fileName,'w')
        outFile.write(fileStr)
        outFile.close()
