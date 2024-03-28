# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 11:17:38 2024

@author: evans
"""

import os
import sys
from asendUtils.syst.pathTools import *
from asendUtils.model.Model import *
from asendUtils.job.ASenDJob import *
from asendUtils.ResultsProcessor import *

if(not os.path.exists('biaxialComposite')):
    os.mkdir('biaxialComposite')
    
if(not os.path.exists('biaxialComposite/results')):
    os.mkdir('biaxialComposite/results')

rtDir = getEnvPath('rootpath')
modFile = rtDir + '/examples/common/biaxialComposite.yaml'
modScr = rtDir + '/examples/modelGeneration'
if(not os.path.exists(modFile)):
    sys.path.append(modScr)
    import biaxialComposite
    
constMod = Model()
constMod.fixDisplacement('xMinRef',ux=0.,uy=0.,uz=0.)
constMod.fixDisplacement('xMaxRef',ux=0.01,uy=0.,uz=0.)
constMod.fixDisplacement('yMinRef',uz=0.)
constMod.periodicDisplacement()
constFile = 'biaxialComposite/constraints.yaml'
constMod.writeModelInput(constFile)

# loadMod = Model()
# loadMod.addNodalForce('periodicXMax',F=[1000.,0.,0.],M=[0.,0.,0.])
# loadFile = 'biaxialComposite/loads.yaml'
# loadMod.writeModelInput(loadFile)

job = ASenDJob()
job.readModelInput(modFile)
job.readConstraints(constFile)
#job.readLoads(loadFile)
job.solve(solverMethod='iterative',solverBlockDim=30)
resFile = 'biaxialComposite/results/nodeResults.yaml'
job.writeNodeResults(resFile,['displacement'])
jobFile = 'biaxialComposite/job.yaml'
job.writeJobInput(jobFile)
job.executeJob()

rp = ResultsProcessor(modFile,nodeResFile=resFile)
rp.plotNodeResults('displacement',component=3 , elementSet='all', deformed=True,defScaleFact=10.0)