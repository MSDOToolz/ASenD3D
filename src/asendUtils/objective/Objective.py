import os
import yaml
from asendUtils.syst.pathTools import *

class Objective():

    def __init__(self):
        self.fileName = ''
        self.objData = dict()
        self.objData['objectiveTerms'] = list()
        self.categories = 'displacement velocity acceleration temperature tDot stress strain strainEnergy shellDef shellFrcMom beamDef beamFrcMom flux tempGradient mass volume massDisp'
        self.operators = 'powerNorm volumeIntegral volumeAverage'
        
    def addObjectiveTerm(self,category,operator='powerNorm',component=1,layer=-1,elementSet='',nodeSet='',stTime=0.0,endTime=1e+100,coefficient=1.0,exponent=2.0,targetValue=0.0):
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
        newTerm['activeTime'] = str([stTime,endTime])
        newTerm['coefficient'] = coefficient
        newTerm['exponent'] = exponent
        newTerm['targetValue'] = targetValue
        self.objData['objectiveTerms'].append(newTerm)
        
    def writeInput(self,fileName):
        self.fileName = makeAbsolute(fileName)
        
        fileStr = yaml.dump(self.objData,sort_keys=False)
        
        fileStr = fileStr.replace("'","")
        fileStr = fileStr.replace('"','')
        
        outFile = open(self.fileName,'w')
        outFile.write(fileStr)
        outFile.close()
        