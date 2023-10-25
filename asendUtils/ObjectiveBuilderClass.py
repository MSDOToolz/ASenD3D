import os
import yaml

class ObjectiveBuilder():

    def __init__(self):
        self.objData = dict()
        self.objData['objectiveTerms'] = list()
        self.categories = 'displacement velocity acceleration temperature tDot stress strain strainEnergy shellDef shellFrcMom beamDef beamFrcMom flux tempGradient mass volume massDisp'
        self.operators = 'powerNorm volumeIntegral volumeAverage'
        
    def addObjectiveTerm(self,category,operator='powerNorm',component=1,layer=-1,elementSet='',nodeSet='',activeTime=0.0,coefficient=1.0,exponent=2.0,targetValue=[]):
        newTerm = dict()
        if(category in self.categories):
            newTerm['category'] = category
        else:
            erStr = 'Error: ' + str(category) + ' is not a recognized category for an objective function term.  Applicable categories are:\n' + self.categories
            raise Exception(erStr)
            return
        if(operator in self.operators):
            newTerm['operator'] = operator
        else:
            erStr = 'Error: ' + str(operator) + ' is not a recognized category for an objective function term.  Applicable categories are:\n' + self.operators
            raise Exception(erStr)
            return
        newTerm['component'] = component
        if(layer > -1):
            newTerm['layer'] = layer
        if(elementSet != ''):
            newTerm['elementSet'] = elementSet
        if(nodeSet != ''):
            newTerm['nodeSet'] = nodeSet
        newTerm['activeTime'] = str(activeTime)
        newTerm['coefficient'] = coefficient
        newTerm['exponent'] = exponent
        try:
            tvl = len(targetValue)
            if(tvl == 0):
                newTerm['targetValue'] = 0.0
            else:
                newTerm['targetValue'] = targetValue
        except:
            newTerm['targetValue'] = targetValue
        self.objData['objectiveTerms'].append(newTerm)
        
    def writeInput(self,fileName)
        outFile = open('temp.yaml','w')
        yaml.dump(self.objData,stream=outFile,sort_keys=False)
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