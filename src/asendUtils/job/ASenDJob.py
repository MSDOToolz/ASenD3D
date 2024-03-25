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
        
    def readConstraints(self,fileName):
        newCmd = dict()
        newCmd['command'] = 'readConstraints'
        newCmd['fileName'] = makeAbsolute(fileName)
        self.jobData['jobCommands'].append(newCmd)
        
    def readLoads(self,fileName):
        newCmd = dict()
        newCmd['command'] = 'readLoads'
        newCmd['fileName'] = makeAbsolute(fileName)
        self.jobData['jobCommands'].append(newCmd)
        
    def readInitialState(self,fileName):
        newCmd = dict()
        newCmd['command'] = 'readInitialState'
        newCmd['fileName'] = makeAbsolute(fileName)
        self.jobData['jobCommands'].append(newCmd)
        
    def readNodeResults(self,fileName):
        newCmd = dict()
        newCmd['command'] = 'readNodeResults'
        newCmd['fileName'] = makeAbsolute(fileName)
        self.jobData['jobCommands'].append(newCmd)
        
    def readDesignVarInput(self,fileName):
        newCmd = dict()
        newCmd['command'] = 'readDesignVarInput'
        newCmd['fileName'] = makeAbsolute(fileName)
        self.jobData['jobCommands'].append(newCmd)
        
    def readDesignVarValues(self,fileName):
        newCmd = dict()
        newCmd['command'] = 'readDesignVarValues'
        newCmd['fileName'] = makeAbsolute(fileName)
        self.jobData['jobCommands'].append(newCmd)
        
    def readObjectiveInput(self,fileName):
        newCmd = dict()
        newCmd['command'] = 'readObjectiveInput'
        newCmd['fileName'] = makeAbsolute(fileName)
        self.jobData['jobCommands'].append(newCmd)
        
    def solvePrep(self,elastic=True,thermal=False,nonlinearGeom=False,staticLoadTime=0.0,
              loadRampSteps=1,dynamic=False,timeStep=1.0,newmarkBeta=0.25,newmarkGamma=0.5,
              simPeriod=1.0,saveSolnHist=False,solnHistDir='',solverMethod='direct',
              solverBlockDim=2000000000,maxIt=0,convTol=1.0e-12):
        newCmd = dict()
        newCmd['command'] = 'setSolveOptions'
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
        newCmd['solverBlockDim'] = solverBlockDim
        if(maxIt != 0):
            newCmd['maxIterations'] = maxIt
        if(convTol != 1.0e-12):
            newCmd['convergenceTol'] = convTol
        self.jobData['jobCommands'].append(newCmd)
        
    def solve(self,elastic=True,thermal=False,nonlinearGeom=False,staticLoadTime=0.0,
              loadRampSteps=1,dynamic=False,timeStep=1.0,newmarkBeta=0.25,newmarkGamma=0.5,
              simPeriod=1.0,saveSolnHist=False,solnHistDir='',solverMethod='direct',
              solverBlockDim=2000000000,maxIt=0,convTol=1.0e-12):
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
            newCmd['dynamic'] = "yes"
        newCmd['timeStep'] = timeStep
        newCmd['newmarkBeta'] = newmarkBeta
        newCmd['newmarkGamma'] = newmarkGamma
        newCmd['simPeriod'] = simPeriod
        if(saveSolnHist):
            newCmd['saveSolnHist'] = 'yes'
        newCmd['solnHistDir'] = makeAbsolute(solnHistDir)
        newCmd['solverMethod'] = solverMethod
        newCmd['solverBlockDim'] = solverBlockDim
        if(maxIt != 0):
            newCmd['maxIterations'] = maxIt
        if(convTol != 1.0e-12):
            newCmd['convergenceTol'] = convTol
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
        self.jobData['jobCommands'].append(newCmd)
        
    def calcObjGradient(self):
        newCmd = dict()
        newCmd['command'] = 'calcObjGradient'
        self.jobData['jobCommands'].append(newCmd)
        
    def writeNodeResults(self,fileName,fields,timeSteps=None,nodeSet='all'):
        newCmd = dict()
        newCmd['command'] = 'writeNodeResults'
        newCmd['fileName'] = makeAbsolute(fileName)
        newCmd['nodeSet'] = nodeSet
        newCmd['fields'] = fields
        if(timeSteps != None):
            newCmd['timeSteps'] = timeSteps
        self.jobData['jobCommands'].append(newCmd)
        
    def writeElementResults(self,fileName,fields,position='centroid',timeSteps=None,elementSet='all'):
        newCmd = dict()
        newCmd['command'] = 'writeElementResults'
        newCmd['fileName'] = makeAbsolute(fileName)
        newCmd['elementSet'] = elementSet
        newCmd['fields'] = fields
        if(position != 'centroid'):
            newCmd['position'] = position
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
        self.jobData['jobCommands'].append(newCmd)
    
    def writeJobInput(self,fileName):
        self.fileName = makeAbsolute(fileName)
        
        fileStr = yaml.dump(self.jobData,sort_keys=False)
        
        fileStr = fileStr.replace("'","")
        fileStr = fileStr.replace('"','')
        
        outFile = open(self.fileName,'w')
        outFile.write(fileStr)
        outFile.close()
        
    def executeJob(self,solverPath='default'):
        if(self.fileName == ''):
            dt = str(datetime.now())
            dt = dt.split('.')[0]
            dt = dt.replace(':','_')
            dt = dt.replace(' ','_')
            fn = 'ASenDJob_' + dt
            self.fileName = fn
            self.writeJobInput(fn)
        if(solverPath != 'default'):
            slvPth = makeAbsolute(solverPath)
            setEnvPath('solverpath',slvPth)
        else:
            slvPth = getEnvPath('solverpath')
        if(not os.path.exists(slvPth)):
            errSt = 'Error: solver path ' + slvPth + 'does not exist. Double check and reset with asendUtils.syst.pathTools.setEnvPath()'
            raise Exception(errSt)
        try:
            ptStr = 'Executing job file ' + self.fileName + ' ...'
            print(ptStr)
            procRes = subprocess.run([slvPth,self.fileName],capture_output=True,text=True)
            ptStr = 'Job ' + self.fileName
            print(ptStr)
            ptStr = 'return code: ' + str(procRes.returncode)
            print(ptStr)
            lfName = self.fileName.split('.')[0] + '.log'
            ptStr = 'see log file, ' + lfName + ' for details.'
            print(ptStr)
            outFile = open(lfName,'w')
            outFile.write(procRes.stdout)
            outFile.write(procRes.stderr)
            outFile.close()
        except:
            errSt = 'Error: problem executing job: ' + self.fileName + '.  Could not complete analysis.\n'
            raise Exception(errSt)