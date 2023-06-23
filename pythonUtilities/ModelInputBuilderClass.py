import numpy as np
import os
##from ruamel.yaml import YAML
import yaml

class ModelInputBuilder():

    def __init__(self):
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
        
    def addSet(self,setType,name,labelList):
        newSet = dict()
        newSet['name'] = name
        newSet['labels'] = labelList
        if(setType == 'node'):
            self.modelData['sets']['node'].append(newSet)
        else:
            self.modelData['sets']['element'].append(newSet)
            
    def addSection(self,elementSet,secType='solid',material='',orientation=[],zOffset=0.0,layup=[],area=0.0,areaMoment=[],polarMoment=0.0,stiffnessMat=[],massMat=[],expLoadCoef=[],conductivity=[],specHeat=0.0):
        newSec = dict()
        newSec['type'] = secType
        newSec['elementSet'] = elementSet
        if(material != ''):
            newSec['material'] = material
        if(len(orientation) > 0):
            newSec['orientation'] = str(orientation)
        if(len(layup) > 0):
            newSec['layup'] = dict()
            newSec['layup']['zOffset'] = zOffset
            lst = list()
            for ly in layup:
                lst.append(str(ly))
            newSec['layup']['layers'] = lst
        if(len(areaMoment) > 0):
            newSec['beamProps'] = dict()
            newSec['beamProps']['area'] = area
            newSec['beamProps']['I'] = areaMoment
            newSec['beamProps']['J'] = polarMoment
        if(len(stiffnessMat) > 0):
            newSec['beamProps'] = dict()
            lst = list()
            for i in stiffnessMat:
                lst.append(str(i))
            newSec['beamProps']['stiffness'] = lst
            if(len(massMat) > 0):
                lst = list()
                for i in stiffnessMat:
                    lst.append(str(i))
                newSec['beamProps']['mass'] = lst
            if(len(expLoadCoef) > 0):
                newSec['beamProps']['expLoadCoef'] = str(expLoadCoeff)
            if(len(conductivity) > 0):
                newSec['beamProps']['conductivity'] = str(conductivity)
            if(specHeat > 0.0):
                newSec['beamProps']['specHeat'] = specHeat
        self.modelData['sections'].append(newSec)
            
    def addMaterial(self,name,density=0.0,modulus=[],poissonRatio=[],shearModulus=[],stiffnessMat=[],conductivity=[],expansion=[],specHeat=0.0,maxTensileStress=[],maxCompressiveStress=[],maxShearStress=[],maxTensileStrain=[],maxCompressiveStrain=[],maxShearStrain=[],maxStrainEnergy=0.0,maxMises=0.0):
        newMat = dict()
        newMat['name'] = name
        if(density > 0.0):
            newMat['density'] = density
        try:
            lm = len(modulus)
        except:
            lm = 1
        if(lm > 0):
            newMat['elastic'] = dict()
            newMat['elastic']['E'] = str(modulus)
            newMat['elastic']['nu'] = str(poissonRatio)
            newMat['elastic']['G'] = str(shearModulus)
        if(len(stiffnessMat) > 0):
            newMat['elastic'] = dict()
            lst = list()
            for i in stiffnessMat:
                lst.append(str(i))
            newMat['elastic']['stiffness'] = lst
        try:
            lc = len(conductivity)
        except:
            lc = 1
        if(lc > 0):
            newMat['thermal'] = dict()
            newMat['thermal']['conductivity'] = str(conductivity)
        try:
            le = len(expansion)
        except:
            le = 1
        if(le > 0):
            newMat['thermal']['expansion'] = str(expansion)
        if(specHeat > 0.0):
            newMat['thermal']['specHeat'] = specHeat
        
        stressDic = dict()
        strainDic = dict()
        failDic = dict()
        
        try:
            le = len(maxTensileStress)
        except:
            le = 1
        if(le > 0):
            stressDic['tensile'] = str(maxTensileStress)
            
        try:
            le = len(maxCompressiveStress)
        except:
            le = 1
        if(le > 0):
            stressDic['compressive'] = str(maxCompressiveStress)
            
        try:
            le = len(maxShearStress)
        except:
            le = 1
        if(le > 0):
            stressDic['shear'] = str(maxShearStress)
        
        if(stressDic):
            failDic['maxStress'] = stressDic
            
        try:
            le = len(maxTensileStrain)
        except:
            le = 1
        if(le > 0):
            strainDic['tensile'] = str(maxTensileStrain)
            
        try:
            le = len(maxCompressiveStrain)
        except:
            le = 1
        if(le > 0):
            strainDic['compressive'] = str(maxCompressiveStrain)
            
        try:
            le = len(maxShearStrain)
        except:
            le = 1
        if(le > 0):
            strainDic['shear'] = str(maxShearStrain)
            
        if(strainDic):
            failDic['maxStrain'] = strainDic
            
        if(maxStrainEnergy > 0.0):
            failDic['maxStrainEnergy'] = maxStrainEnergy
        
        if(maxMises > 0.0):
            failDic['maxMises'] = maxMises
        
        if(failDic):
            newMat['failure'] = failDic
        
        self.modelData['materials'].append(newMat)
        
    def writeModelInput(self,fileName):
        #yamlReader = YAML() ## or YAML(typ='safe'), default is 'rt' for round trip

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