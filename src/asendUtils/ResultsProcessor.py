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

    def plotNodeResults(self,field,component=1,deformed=False):
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
        xMax = np.max(xLst)
        xMin = np.min(xLst)
        xLen = xMax - xMin
        xMid = 0.5*(xMax+xMin)
        yMax = np.max(yLst)
        yMin = np.min(yLst)
        yLen = yMax - yMin
        yMid = 0.5*(yMax+yMin)
        zMax = np.max(zLst)
        zMin = np.min(zLst)
        zLen = zMax - zMin
        zMid = 0.5*(zMax+zMin)
        
        maxLen = np.max([xLen,yLen,zLen])
        midpt = np.array([xMid,yMid,zMid])
        for i in range(0,3):
            proj = np.zeros(3)
            proj[i] = 0.5*maxLen
            newPt = midpt + proj
            xLst.append(newPt[0])
            yLst.append(newPt[1])
            zLst.append(newPt[2])
            newPt = midpt - proj
            xLst.append(newPt[0])
            yLst.append(newPt[1])
            zLst.append(newPt[2])
        values = []
        for nd in self.nodeData['nodeResults'][field]:
            val = nd[component]
            values.append(val)
        values.extend([0.,0.,0.,0.,0.,0.])
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
        cbTitle = field + str(component)
        plotMeshSolution(xLst,yLst,zLst,values,v1,v2,v3,cbTitle)
