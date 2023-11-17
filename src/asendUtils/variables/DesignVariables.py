import os
import yaml
from asendUtils.syst.pathTools import *

class DesignVariables():
    
    def __init__(self):
        self.fileName = ''
        self.desVarData = dict()
        self.desVarData['designVariables'] = list()
        self.categories = 'density massMat modulus shearModulus poissonRatio stiffnessMat thermalCond thermalExp specHeat initialStrain orientation nodeCoord elasticLoad thermalLoad thickness angle zOffset area areaMoment polarMoment dampingMat potFldCoef dampFldCoef massPerEl'
        
    def addDesignVariable(self,category,component=1,layer=-1,nodeSet='',elementSet='',activeTime=0.0,coefficients=1.0):
        newVar = dict()
        if(category in self.categories):
            newVar['category'] = category
        else:
            erStr = 'Error: ' + category + ' is not a recognized category for a design variable.  Applicable categories are: \n' + self.categories
            raise Exception(erStr)
            return
        newVar['component'] = component
        if(layer > -1):
            newVar['layer'] = layer
        if(nodeSet != ''):
            newVar['nodeSet'] = nodeSet
        if(elementSet != ''):
            newVar['elementSet'] = elementSet
        newVar['activeTime'] = str(activeTime)
        newVar['coefficients'] = coefficients
        self.desVarData['designVariables'].append(newVar)
        
    def writeInput(self,fileName):
        self.fileName = makeAbsolute(fileName)
        
        outFile = open('temp.yaml','w')
        yaml.dump(self.desVarData,stream=outFile,sort_keys=False)
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