# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 16:57:06 2024

@author: evans
"""
import os
import sys
from asendUtils.model.Model import *
from asendUtils.syst.pathTools import *
from asendUtils.job.ASenDJob import *
from asendUtils.ResultsProcessor import *

if(not os.path.exists('shellDiskImpact')):
    os.mkdir('shellDiskImpact')
    
if(not os.path.exists('shellDiskImpact/results')):
    os.mkdir('shellDiskImpact/results')
    
rtDir = getEnvPath('rootpath')
modFile = rtDir + '/examples/common/shellDiskProjectile.yaml'
modScrDir = rtDir + '/examples/modelGeneration'
if(not os.path.exists(modFile)):
    sys.path.append(modScrDir)
    import shellDiskProjectile
    
constMod = Model()
constMod.fixDisplacement('boundaryNodes',ux=0.0,uy=0.0,uz=0.0)
constFile = 'shellDiskImpact/constraints.yaml'
constMod.writeModelInput(constFile)

initMod = Model()
state = [['projectileNode',0.0,0.0,-10.0,0.0,0.0,0.0]]
initMod.addInitialState('velocity',state)
initFile = 'shellDiskImpact/initialState.yaml'
initMod.writeModelInput(initFile)

job = ASenDJob()
job.readModelInput(modFile)
job.readConstraints(constFile)
job.readInitialState(initFile)
job.solve(nonlinearGeom=True,dynamic=True,timeStep=0.000625,simPeriod=0.2,solnHistDir='shellDiskImpact/results/',lumpMass=True,solverMethod='iterative',solverBlockDim=2)
#job.solve(nonlinearGeom=True,dynamic=True,timeStep=0.005,simPeriod=0.5,saveSolnHist=True,solnHistDir='shellDiskImpact/results/')
resFile = 'shellDiskImpact/results/nodeResults.yaml'
ts = list(range(0,40))
job.writeNodeResults(resFile,['displacement'],timeSteps=ts)
jobFile = 'shellDiskImpact/job.yaml'
job.writeJobInput(jobFile)
job.executeJob()

rp = ResultsProcessor(modFile)
ndResFile = 'shellDiskImpact/results/nodeResults.yaml'
#rp.animateNodeResults(ndResFile,'displacement',ts,component=3,elementSet='allDiskEls',deformed=True,defScaleFact=50.0)
rp.plotNodeHistory(ndResFile,'displacement',ts,'projectileNode',component=3)

inFile = open('shellDiskImpact/job.log')
fileLine = inFile.readline()
totIt = 0
while(fileLine != ''):
    if('Total CG iterations:' in fileLine):
        lst = fileLine.split(':')
        totIt += int(lst[1])
    fileLine = inFile.readline()
print('Total CG iterations: ' + str(totIt))