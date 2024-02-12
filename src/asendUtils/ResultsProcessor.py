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
        
    def buildNodalPlotCrd(self,component=1,deformed=False):
        xLst = []
        yLst = []
        zLst = []      
        for nd in self.modelData['nodes']:
            xLst.append(nd[1])
            yLst.append(nd[2])
            zLst.append(nd[3])
        if(deformed):
            for ndU in self.nodeData['nodeResults']['displacement']:
                i = ndU[0]
                xLst[i] = xLst[i] + ndU[1]
                yLst[i] = yLst[i] + ndU[2]
                zLst[i] = zLst[i] + ndU[3]
            
        return {'xLst': xLst, 'yLst': yLst, 'zLst': zLst}

    def buildElementVertexList(self):
        v1 = []
        v2 = []
        v3 = []
        for et in self.modelData['elements']:
            if('brick' in et['type']):
                for el in et['connectivity']:
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
                    n1 = el[1]
                    n2 = el[2]
                    n3 = el[3]
                    
                    v1.append(n1)
                    v2.append(n2)
                    v3.append(n3)
                    
        return {'v1': v1, 'v2': v2, 'v3': v3}

    def plotNodeResults(self,field,component=1,deformed=False):
        ndCrd = self.buildNodalPlotCrd(component,deformed)
        
        values = []
        for nd in self.nodeData['nodeResults'][field]:
            val = nd[component]
            values.append(val)
        
        verts = self.buildElementVertexList()
        cbTitle = field + str(component)
        plotMeshSolution(ndCrd,values,verts,cbTitle)
        
    def plotElementResults(self,field,component=1,deformed=False):
        ndCrd = self.buildNodalPlotCrd(component,deformed)
        
        numNds = len(self.modelData['nodes'])
        ndValues = np.zeros(numNds,dtype=float)
        numHits = np.zeros(numNds,dtype=int)
        numEls = len(self.elementData['elementResults'][field])
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
        for i, nh in enumerate(numHits):
            if(nh > 0):
                ndValues[i] = ndValues[i]/nh
        
        verts = self.buildElementVertexList()
        cbTitle = field + str(component)
        plotMeshSolution(ndCrd,ndValues,verts,cbTitle)
        
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