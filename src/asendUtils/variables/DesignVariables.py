import os
import yaml

class DesignVariables():
    
    def __init__(self):
        self.desVarData = dict()
        self.desVarData['designVariables'] = list()
        self.categories = 'density massMat modulus shearModulus poissonRatio stiffnessMat thermalCond thermalExp specHeat orientation nodeCoord elasticLoad thermalLoad thickness angle zOffset area areaMoment polarMoment'
        
    def addDesignVariable(self,category,component=1,layer=-1,nodeSet='',elementSet='',activeTime=0.0,coefficients=[]):
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
        try:
            lc = len(coefficients)
            if(lc > 0):
                newVar['coefficients'] = coefficients
        except:
            newVar['coefficients'] = coefficients
        self.desVarData['designVariables'].append(newVar)
        
    def writeInput(self,fileName):
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