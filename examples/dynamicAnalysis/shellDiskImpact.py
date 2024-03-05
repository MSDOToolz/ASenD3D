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

if(not os.path.exists('shellDiskImpact')):
    os.mkdir('shellDiskImpact')
    
if(not os.path.exists('shellDiskImpact/results')):
    os.mkdir('shellDiskImpact/results')
    
rtDir = getEnvPath('rootpath')
modFile = rtDir + '/examples/common/shellDisk.yaml'
modScrDir = rtDir + '/examples/modelGeneration'
if(not os.path.exists(modFile)):
    sys.path.append(modScrDir)
    import shellDisk
    
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
job.solve(nonlinearGeom=True,dynamic=True,timeStep=0.005,simPeriod=0.2,saveSolnHist=True)
resFile = 'shellDiskImpact/results/nodeResults.yaml'
job.writeNodeResults(resFile,['displacement'])
jobFile = 'shellDiskImpact/job.yaml'
job.writeJobInput(jobFile)
#job.executeJob()