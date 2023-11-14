# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 09:18:03 2023

@author: evans
"""

import os
import sys
import subprocess
import platform
from datetime import datetime
import yaml
from asendUtils.syst.pathTools import *

class ASenDJob:
    
    def __init__(self):
        self.fileName = ''
        self.jobData = dict()
        self.jobData['jobCommands'] = list()
        
    def readModelInput(self,fileName):
        newCmd = dict()
        newCmd['command'] = 'readModelInput'
        newCmd['fileName'] = makeAbsolute(fileName)
        self.jobData['jobCommands'].append(newCmd)
        
    def readLoads(self,fileName):
        newCmd = dict()
        newCmd['command'] = 'readLoads'
        newCmd['fileName'] = makeAbsolute(fileName)
        self.jobData['jobCommands'].append(newCmd)
        
    def readDesignVarInput(self,fileName):
        newCmd = dict()
        newCmd['command'] = 'readDesignVarInput'
        newCmd['fileName'] = makeAbsolute(fileName)
        self.jobData['jobCommands'].append(newCmd)
        
    def readObjectiveInput(self,fileName):
        newCmd = dict()
        newCmd['command'] = 'readObectiveInput'
        newCmd['fileName'] = makeAbsolute(fileName)
        self.jobData['jobCommands'].append(newCmd)
        
    def solve(self,elastic=True,thermal=False,nonlinearGeom=False,staticLoadTime=0.0,
              loadRampSteps=1,dynamic=False,timeStep=1.0,newmarkBeta=0.25,newmarkGamma=0.5,
              simPeriod=1.0,saveSolnHist=False,solnHistDir='',solverMethod='direct',
              solverBandwidth=1e+10,solverBlockDim=1e+10):
        newCmd = dict()
        newCmd['command'] = 'solve'
        if(not elastic):
            newCmd['elastic'] = 'no'
        if(thermal):
            newCmd['thermal'] = 'yes'
        if(nonlinearGeom):
            newCmd['nonlinearGeom'] = 'yes'
        newCmd['staticLoadTime'] = staticLoadTime
        newCmd['loadRampSteps'] = loadRampSteps
        if(dynamic):
            newCmd['dynamic'] = 'yes'
        newCmd['timeStep'] = timeStep
        newCmd['newmarkBeta'] = newmarkBeta
        newCmd['newmarkGamma'] = newmarkGamma
        newCmd['simPeriod'] = simPeriod
        if(saveSolnHist):
            newCmd['saveSolnHist'] = 'yes'
        newCmd['solnHistDir'] = makeAbsolute(solnHistDir)
        newCmd['solverMethod'] = solverMethod
        newCmd['solverBandwidth'] = solverBandwidth
        newCmd['solverBlockDim'] = solverBlockDim
        self.jobData['jobCommands'].append(newCmd)
        
    def modalAnalysis(self,analysisType='buckling',numModes=10,solverMethod='direct'):
        newCmd = dict()
        newCmd['command'] = 'modalAnalysis'
        newCmd['type'] = analysisType
        newCmd['numModes'] = numModes
        newCmd['solverMethod'] = solverMethod
        self.jobData['jobCommands'].append(newCmd)
        
    def setSolnToMode(self,solnField='displacement',mode=0,maxAmplitude = 1.0):
        newCmd = dict()
        newCmd['command'] = 'setSolnToMode'
        newCmd['solnField'] = solnField
        newCmd['mode'] = mode
        newCmd['maxAmplitude'] = maxAmplitude
        self.jobData['jobCommands'].append(newCmd)
        
    def calcObjective(self):
        newCmd = dict()
        newCmd['command'] = 'calcObjective'
        
    def calcObjGradient(self):
        newCmd = dict()
        newCmd['command'] = 'calcObjGradient'
        
    def writeNodeResults(self,fileName,fields,timeSteps=None,nodeSet='all'):
        newCmd = dict()
        newCmd['command'] = 'writeNodeResults'
        newCmd['fileName'] = makeAbsolute(fileName)
        newCmd['nodeSet'] = nodeSet
        newCmd['fields'] = fields
        if(timeSteps != None):
            newCmd['timeSteps'] = timeSteps
        self.jobData['jobCommands'].append(newCmd)
        
    def writeElementResults(self,fileName,fields,timeSteps=None,elementSet='all'):
        newCmd = dict()
        newCmd['command'] = 'writeElementResults'
        newCmd['fileName'] = makeAbsolute(fileName)
        newCmd['elementSet'] = elementSet
        newCmd['fields'] = fields
        if(timeSteps != None):
            newCmd['timeSteps'] = timeSteps
        self.jobData['jobCommands'].append(newCmd)
    
    def writeModalResults(self,fileName,writeModes=True):
        newCmd = dict()
        newCmd['command'] = 'writeModalResults'
        newCmd['fileName'] = makeAbsolute(fileName)
        if(not writeModes):
            newCmd['writeModes'] = 'no'
        self.jobData['jobCommands'].append(newCmd)
    
    def writeObjective(self,fileName,include=None,writeGradient=True):
        newCmd = dict()
        newCmd['command'] = 'writeObjective'
        newCmd['fileName'] = makeAbsolute(fileName)
        if(include != None):
            newCmd['include'] = include
        if(not writeGradient):
            newCmd['writeGradient'] = 'no'
    
    def writeJobInput(self,fileName):
        self.fileName = makeAbsolute(fileName)
        
        outFile = open(self.fileName,'w')
        yaml.dump(self.jobData,stream=outFile,sort_keys=False)
        outFile.close()
        
    def executeJob(self,solverPath='default'):
        if(self.fileName == ''):
            dt = str(datetime.now())
            dt = dt.split('.')[0]
            dt = dt.replace(':','_')
            dt = dt.replace(' ','_')
            fn = 'ASenDJob_' + dt
            self.writeJobInput(fn)
        slvPth = getSolverPath()
        if(not os.path.exists(slvPth)):
            errSt = 'Error: solver path ' + slvPth + 'does not exist. Double check and reset with asendUtils.syst.pathTools.setSolverPath()'
            raise Exception(errStr)
        try:
            subprocess.run([slvPth,self.fileName])
        except:
            errSt = 'Error: problem executing job: ' + self.fileName + '.  Could not complete analysis.\n'
            raise Exception(errSt)