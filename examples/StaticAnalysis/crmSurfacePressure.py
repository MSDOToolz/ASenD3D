# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 20:50:57 2024

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.syst.pathTools import *
from asendUtils.job.ASenDJob import *
from asendUtils.ResultsProcessor import *

if(not os.path.exists('crmSurfacePressure')):
    os.mkdir('crmSurfacePressure')
    
if(not os.path.exists('crmSurfacePressure/results')):
    os.mkdir('crmSurfacePressure/results')
    
rtDir = getEnvPath('rootpath')
modFile = rtDir + '/examples/common/crmWing.yaml'
modScrDir = rtDir + '/examples/modelGeneration'
if(not os.path.exists(modFile)):
    sys.path.append(modScrDir)
    import crmWing

constMod = Model()
constMod.fixDisplacement('rootNodes',ux=0.0,uy=0.0,uz=0.0)
constFile = 'crmSurfacePressure/constraints.yaml'
constMod.writeModelInput(constFile)

loadMod = Model()
loadMod.addSurfacePressure('pressureSideEls',2000.0,[-0.5,-1.0,0.0],normTol=20.0)
loadFile = 'crmSurfacePressure/loads.yaml'
loadMod.writeModelInput(loadFile)

job = ASenDJob()
job.readModelInput(modFile)
job.readConstraints(constFile)
job.readLoads(loadFile)
job.solve()
resFile = 'crmSurfacePressure/results/nodeResults.yaml'
job.writeNodeResults(resFile, ['displacement'])
job.writeJobInput('crmSurfacePressure/job.yaml')
job.executeJob()

rp = ResultsProcessor(modFile,nodeResFile=resFile)
rp.plotNodeResults('displacement',component='mag',deformed=True,defScaleFact=50.0)