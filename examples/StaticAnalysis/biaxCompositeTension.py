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

constFile = 'biaxialComposite/constraints.yaml'
loadFile = 'biaxialComposite/loads.yaml'
nodeResFile = 'biaxialComposite/results/nodeResults.yaml'
elResFile = 'biaxialComposite/results/elementResults.yaml'
jobFile = 'biaxialComposite/job.yaml'

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
    
# constMod = Model()
# constMod.fixDisplacement('xMinRef',ux=0.,uy=0.,uz=0.)
# constMod.fixDisplacement('xMaxRef',uy=0.,uz=0.)
# constMod.fixDisplacement('yMinRef',uz=0.)
# constMod.periodicDisplacement()
# constMod.writeModelInput(constFile)

# loadMod = Model()
# #loadMod.addNodalForce('xMaxRef',F=[1000000.,0.,0.],M=[0.,0.,0.])
# loadMod.addNodalForce('yMaxRef',F=[1000000.,0.,0.],M=[0.,0.,0.])
# loadMod.addNodalForce('yMinRef',F=[-1000000.,0.,0.],M=[0.,0.,0.])
# loadMod.writeModelInput(loadFile)

# job = ASenDJob()
# job.readModelInput(modFile)
# job.readConstraints(constFile)
# job.readLoads(loadFile)
# job.solve()
# job.writeNodeResults(nodeResFile,['displacement'])
# job.writeElementResults(elResFile,['stress'])
# job.writeJobInput(jobFile)
# job.executeJob()

rp = ResultsProcessor(modFile,nodeResFile=nodeResFile,elementResFile=elResFile)
rp.plotNodeResults('displacement',component=1,elementSet='all', deformed=True,defScaleFact=1000.0)
rp.plotElementResults('stress',component=4,deformed=True,defScaleFact=1000.0)