#from ruamel.yaml import YAML
import numpy as np
import yaml
from yaml import CLoader as Loader
from asendUtils.visualization.plotlyUtils import *

class ResultsProcessor:

    def __init__(self,modelFile,nodeResFile='none',elementResFile='none',modalResFile='none',objResFile='none'):

        inFile = open(modelFile,'r')
        self.modelData = yaml.load(inFile,Loader=Loader)
        inFile.close()

        if(not nodeResFile == 'none'):
            self.loadNodeResults(nodeResFile)
        else:
            self.nodeData = dict()

        if(not elementResFile == 'none'):
            self.loadElementResults(elementResFile)
        else:
            self.elementData = dict()
            
        if(not modalResFile == 'none'):
            self.loadModalResults(modalResFile)
        else:
            self.modalData = dict()
            
        if(not objResFile == 'none'):
            self.loadObjectiveResults(objResFile)
        else:
            self.objectiveData = dict()

    def loadNodeResults(self,nodeResFile):
        inFile = open(nodeResFile,'r')
        self.nodeData = yaml.load(inFile,Loader=Loader)
        inFile.close()

    def loadElementResults(self,elementResFile):
        inFile = open(elementResFile,'r')
        self.elementData = yaml.load(inFile,Loader=Loader)
        inFile.close()
        
    def loadModalResults(self,modalResFile):
        inFile = open(modalResFile,'r')
        self.modalData = yaml.load(inFile,Loader=Loader)
        inFile.close()
        
    def loadObjectiveResults(self,objResFile):
        inFile = open(objResFile,'r')
        self.objectiveData = yaml.load(inFile,Loader=Loader)
        inFile.close()
        
    def getPlotNdElSet(self,elementSet):
        numNds = len(self.modelData['nodes'])
        if(elementSet == 'all'):
            ndSet = set(range(0,numNds))
            elSet = set()
            for et in self.modelData['elements']:
                for el in et['connectivity']:
                    elSet.add(el[0])
        else:
            for es in self.modelData['sets']['element']:
                if(es['name'] == elementSet):
                    elSet = set(es['labels'])
                    ndSet = set()
                    for et in self.modelData['elements']:
                        for el in et['connectivity']:
                            lab = el[0]
                            if(lab in elSet):
                                for ni, nd in enumerate(el):
                                    if(ni > 0):
                                        ndSet.add(nd)
        return ndSet, elSet
        
    def buildNodalPlotCrd(self,ndSet,deformed=False,defScaleFact=1.0):
        ## ndSet a python set() of the desired node labels
        allNds = self.modelData['nodes']
        numNds = len(allNds)
        crdAr = np.zeros((numNds,3),dtype=float)
        for nd in allNds:
            lab = nd[0]
            if(lab in ndSet):
                crdAr[lab] = np.array(nd[1:4])
        if(deformed):
            for ndU in self.nodeData['nodeResults']['displacement']:
                lab = ndU[0]
                if(lab in ndSet):
                    crdAr[lab] = crdAr[lab] + defScaleFact*np.array(ndU[1:4])
        xLst = list()
        yLst = list()
        zLst = list()
        for nd in allNds:
            lab = nd[0]
            if(lab in ndSet):
                xLst.append(crdAr[lab,0])
                yLst.append(crdAr[lab,1])
                zLst.append(crdAr[lab,2])
            
        return {'xLst': xLst, 'yLst': yLst, 'zLst': zLst}

    def buildElementVertexList(self,elSet):
        # elSet = python set() with the desired element labels
        v1 = []
        v2 = []
        v3 = []
        for et in self.modelData['elements']:
            if('brick' in et['type']):
                for el in et['connectivity']:
                    eli = el[0]
                    if(eli in elSet):
                        n1 = el[1]
                        n2 = el[2]
                        n3 = el[3]
                        n4 = el[4]
                        n5 = el[5]
                        n6 = el[6]
                        n7 = el[7]
                        n8 = el[8]
                        
                        v1.append(n1)
                        v2.append(n4)
                        v3.append(n5)
    
                        v1.append(n4)
                        v2.append(n5)
                        v3.append(n8)
    
                        v1.append(n2)
                        v2.append(n3)
                        v3.append(n6)
    
                        v1.append(n3)
                        v2.append(n6)
                        v3.append(n7)
    
                        v1.append(n1)
                        v2.append(n2)
                        v3.append(n5)
    
                        v1.append(n2)
                        v2.append(n5)
                        v3.append(n6)
    
                        v1.append(n3)
                        v2.append(n4)
                        v3.append(n7)
    
                        v1.append(n4)
                        v2.append(n7)
                        v3.append(n8)
    
                        v1.append(n1)
                        v2.append(n2)
                        v3.append(n3)
    
                        v1.append(n1)
                        v2.append(n3)
                        v3.append(n4)
    
                        v1.append(n5)
                        v2.append(n6)
                        v3.append(n7)
    
                        v1.append(n5)
                        v2.append(n7)
                        v3.append(n8)
            elif('tet' in et['type']):
                for el in et['connectivity']:
                    eli = el[0]
                    if(eli in elSet):
                        n1 = el[1]
                        n2 = el[2]
                        n3 = el[3]
                        n4 = el[4]
                        
                        v1.append(n1)
                        v2.append(n2)
                        v3.append(n3)
    
                        v1.append(n1)
                        v2.append(n2)
                        v3.append(n4)
    
                        v1.append(n1)
                        v2.append(n3)
                        v3.append(n4)
    
                        v1.append(n2)
                        v2.append(n3)
                        v3.append(n4)
            elif('shell4' in et['type']):
                for el in et['connectivity']:
                    eli = el[0]
                    if(eli in elSet):
                        n1 = el[1]
                        n2 = el[2]
                        n3 = el[3]
                        n4 = el[4]
                        
                        v1.append(n1)
                        v2.append(n2)
                        v3.append(n3)
    
                        v1.append(n1)
                        v2.append(n3)
                        v3.append(n4)
            elif('shell3' in et['type']):
                for el in et['connectivity']:
                    eli = el[0]
                    if(eli in elSet):
                        n1 = el[1]
                        n2 = el[2]
                        n3 = el[3]
                        
                        v1.append(n1)
                        v2.append(n2)
                        v3.append(n3)
                    
        return {'v1': v1, 'v2': v2, 'v3': v3}
    
    def getFaceValues(self,elSet,elVal):
        fcVal = list()
        for et in self.modelData['elements']:
            eTp = et['type']
            for el in et['connectivity']:
                eli = el[0]
                if(eli in elSet):
                    si = elVal[eli]
                    if('brick' in eTp):
                        fcVal.extend([si,si,si,si,si,si,si,si,si,si,si,si])
                    elif('wedge' in eTp):
                        fcVal.extend([si,si,si,si,si,si,si,si])
                    elif('tet' in eTp):
                        fcVal.extend([si,si,si,si])
                    elif('shell4' in eTp):
                        fcVal.extend([si,si])
                    elif('shell3' in eTp):
                        fcVal.extend([si])
                        
        return fcVal
    
    def plotElementProperty(self,prop='section',elementSet='all'):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
        ndCrd = self.buildNodalPlotCrd(ndSet)
        verts = self.buildElementVertexList(elSet)
        
        if(prop == 'section'):
            numEls = 0
            for et in self.modelData['elements']:
                numEls = numEls + len(et['connectivity'])
            elVal = np.zeros(numEls,dtype=int)
            
            elSets = self.modelData['sets']['element']
            setDic = dict()
            for si, es in enumerate(elSets):
                setDic[es['name']] = si
            
            for si, sec in enumerate(self.modelData['sections']):
                setNm = sec['elementSet']
                seti = setDic[setNm]
                for eli in elSets[seti]['labels']:
                    elVal[eli] = si
        
        fcVals = self.getFaceValues(elSet,elVal)
        cbTitle = prop
        plotMeshSolution(ndCrd,fcVals,verts,valMode='cell',title=cbTitle)

    def plotNodeResults(self,field,component=1,elementSet='all',deformed=False,defScaleFact=1.0):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
        ndCrd = self.buildNodalPlotCrd(ndSet,deformed,defScaleFact)
        numNds = len(self.modelData['nodes'])
        
        valAr = np.zeros(numNds,dtype=float)
        for nd in self.nodeData['nodeResults'][field]:
            lab = nd[0]
            if(component == 'mag'):
                vec = np.array(nd[1:4])
                val = np.linalg.norm(vec)
            else:
                val = nd[component]
            valAr[lab] = val
        values = list()
        for nd in self.modelData['nodes']:
            lab = nd[0]
            if(lab in ndSet):
                values.append(valAr[lab])
        
        verts = self.buildElementVertexList(elSet)
        cbTitle = field + str(component)
        plotMeshSolution(ndCrd,values,verts,valMode='vertex',title=cbTitle)
        
    def plotElementResults(self,field,component=1,elementSet='all',deformed=False,defScaleFact=1.0):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
                                        
        ndCrd = self.buildNodalPlotCrd(ndSet,deformed,defScaleFact)
        numNds = len(self.modelData['nodes'])
        
        ndValues = np.zeros(numNds,dtype=float)
        numHits = np.zeros(numNds,dtype=int)
        numEls = 0
        for et in self.modelData['elements']:
            numEls = numEls + len(et['connectivity'])
        elValues = np.zeros(numEls,dtype=float)
        for el in self.elementData['elementResults'][field]:
            lab = el[0]
            i = component + 2
            val = el[i]
            elValues[lab] = val
        for et in self.modelData['elements']:
            elst = et['connectivity']
            for el in elst:
                lab = el[0]
                v = elValues[lab]
                for i, nd in enumerate(el):
                    if(i > 0):
                        ndValues[nd] = ndValues[nd] + v
                        numHits[nd] = numHits[nd] + 1
        nVLst = list()
        for nd in self.modelData['nodes']:
            lab = nd[0]
            if(lab in ndSet):
                nVLst.append(ndValues[lab]/numHits[lab])
        
        verts = self.buildElementVertexList(elSet)
        cbTitle = field + str(component)
        plotMeshSolution(ndCrd,nVLst,verts,valMode='vertex',title=cbTitle)
        
    def plotModalResults(self,mode,elementSet='all',defScaleFact=1.0):
        nodeCopy = self.nodeData.copy()
        for md in self.modalData['modalResults']['modes']:
            if(md['mode'] == mode):
                self.nodeData['nodeResults'] = dict()
                self.nodeData['nodeResults']['displacement'] = md['displacement']
        self.plotNodeResults('displacement',component='mag',elementSet=elementSet,deformed=True,defScaleFact=defScaleFact)
        self.nodeData = nodeCopy
        
    def animateNodeResults(self,fileName,field,timeSteps,component=1,elementSet='all',deformed=False,defScaleFact=1.0,frameDuration=1000):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
        
        allNdCrd = list()
        allNdValues = list()
        numNds = len(self.modelData['nodes'])
        fnLst = fileName.split('.')
        verts = self.buildElementVertexList(elSet)
        valAr = np.zeros(numNds,dtype=float)
        for ts in timeSteps:
            fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
            self.loadNodeResults(fn)
            ndCrd = self.buildNodalPlotCrd(ndSet,deformed,defScaleFact)
            for nd in self.nodeData['nodeResults'][field]:
                lab = nd[0]
                valAr[lab] = nd[component]
            ndValues = list()
            for nd in self.modelData['nodes']:
                lab = nd[0]
                if(lab in ndSet):
                    ndValues.append(valAr[lab])
            allNdCrd.append(ndCrd)
            allNdValues.append(ndValues)
        cbTitle = field + str(component)
        animateMeshSolution(allNdCrd,allNdValues,verts,frameDuration,cbTitle)
        
    def animateElementResults(self,fileName,field,timeSteps,component=1,elementSet='all',deformed=False,defScaleFact=1.0,nodeResFile=None,frameDuration=1000):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
        
        allNdCrd = list()
        allNdValues = list()
        numNds = len(self.modelData['nodes'])
        fnLst = fileName.split('.')
        if(nodeResFile != None):
            ndFnLst = nodeResFile.split('.')
        verts = self.buildElementVertexList(elSet)
        
        numEls = 0
        for et in self.modelData['elements']:
            numEls = numEls + len(et['connectivity'])
        elValues = np.zeros(numEls,dtype=float)
        
        for ts in timeSteps:
            if(deformed and nodeResFile != None):
                fn = ndFnLst[0] + '_timestep' + str(ts) + '.' + ndFnLst[1]
                self.loadNodeResults(fn)
            ndCrd = self.buildNodalPlotCrd(ndSet,deformed,defScaleFact)
            fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
            self.loadElementResults(fn)
            for el in self.elementData['elementResults'][field]:
                lab = el[0]
                i = component + 2
                val = el[i]
                elValues[lab] = val
            valAr = np.zeros(numNds,dtype=float)
            numHits = np.zeros(numNds,dtype=int)
            for et in self.modelData['elements']:
                elst = et['connectivity']
                for el in elst:
                    lab = el[0]
                    v = elValues[lab]
                    for i, nd in enumerate(el):
                        if(i > 0):
                            valAr[nd] = valAr[nd] + v
                            numHits[nd] = numHits[nd] + 1
            ndValues = list()
            for nd in self.modelData['nodes']:
                lab = nd[0]
                if(lab in ndSet):
                    ndValues.append(valAr[lab]/numHits[lab])
            allNdCrd.append(ndCrd)
            allNdValues.append(ndValues)
        cbTitle = field + str(component)
        animateMeshSolution(allNdCrd,allNdValues,verts,frameDuration,cbTitle)
        
    def animateModalSolution(self,mode,elementSet='all',defScaleFact=1.0):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
        nodeCopy = self.nodeData.copy()
        for md in self.modalData['modalResults']['modes']:
            if(md['mode'] == mode):
                self.nodeData['nodeResults'] = dict()
                self.nodeData['nodeResults']['displacement'] = md['displacement']
        
        allNdCrd = list()
        allNdValues = list()
        numNds = len(self.modelData['nodes'])
        valAr = np.zeros(numNds,dtype=float)
        verts = self.buildElementVertexList(elSet)
        
        for theta in range(0,360,30):
            tRad = 0.0174533*theta
            sinTh = np.math.sin(tRad)
            sf = sinTh*defScaleFact
            ndCrd = self.buildNodalPlotCrd(ndSet,deformed=True,defScaleFact=sf)
            for nd in self.nodeData['nodeResults']['displacement']:
                lab = nd[0]
                if(lab in ndSet):
                    vec = np.array(nd[1:4])
                    valAr[lab] = np.linalg.norm(vec)
            ndValues = list()
            for nd in self.modelData['nodes']:
                lab = nd[0]
                if(lab in ndSet):
                    ndValues.append(valAr[lab])
            allNdCrd.append(ndCrd)
            allNdValues.append(ndValues)
        cbTitle = 'displacement'
        animateMeshSolution(allNdCrd,allNdValues,verts,title=cbTitle)
        self.nodeData = nodeCopy
        
    def plotNodeHistory(self,fileName,field,timeSteps,nodeSet,component=1):
        fnLst = fileName.split('.')
        try:
            ndI = int(nodeSet)
            vals = list()
            timePts = list()
            for ts in timeSteps:
                fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
                self.loadNodeResults(fn)
                for nd in self.nodeData['nodeResults'][field]:
                    if(nd[0] == ndI):
                        vals.append(nd[component])
                        timePts.append(self.nodeData['nodeResults']['time'])
            series = dict()
            lab = 'node_' + str(ndI)
            series[lab] = vals
            ytitle = field + str(component)
            plotTimeHistory(series,timePts,ytitle)
        except:
            series = dict()
            for ns in self.modelData['sets']['node']:
                if(ns['name'] == nodeSet):
                    for nd in ns['labels']:
                        lab = 'node_' + str(nd)
                        series[lab] = list()
            timePts = list()
            for ts in timeSteps:
                fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
                self.loadNodeResults(fn)
                timePts.append(self.nodeData['nodeResults']['time'])
                for nd in self.nodeData['nodeResults'][field]:
                    try:
                        lab = 'node_' + str(nd[0])
                        series[lab].append(nd[component])
                    except:
                        pass
            ytitle = field + str(component)
            plotTimeHistory(series,timePts,ytitle)
        
    def plotElementHistory(self,fileName,field,timeSteps,elementSet,component=1):
        fnLst = fileName.split('.')
        try:
            elI = int(elementSet)
            vals = list()
            timePts = list()
            for ts in timeSteps:
                fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
                self.loadElementResults(fn)
                for el in self.elementData['elementResults'][field]:
                    if(el[0] == elI):
                        vals.append(el[component+2])
                        timePts.append(self.nodeData['elementResults']['time'])
            series = dict()
            lab = 'element_' + str(ndI)
            series[lab] = vals
            ytitle = field + str(component)
            plotTimeHistory(series,timePts,ytitle)
        except:
            series = dict()
            for es in self.modelData['sets']['element']:
                if(es['name'] == elementSet):
                    for el in es['labels']:
                        lab = 'element_' + str(el)
                        series[lab] = list()
            timePts = list()
            for ts in timeSteps:
                fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
                self.loadElementResults(fn)
                timePts.append(self.elementData['elementResults']['time'])
                for el in self.elementData['elementResults'][field]:
                    try:
                        lab = 'element_' + str(el[0])
                        series[lab].append(el[component+2])
                    except:
                        pass
            ytitle = field + str(component)
            plotTimeHistory(series,timePts,ytitle)