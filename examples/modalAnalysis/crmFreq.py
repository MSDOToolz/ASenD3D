# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 09:43:52 2024

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.syst.pathTools import *
from asendUtils.job.ASenDJob import *
from asendUtils.ResultsProcessor import *

if(not os.path.exists('crmFreq')):
    os.mkdir('crmFreq')
    
if(not os.path.exists('crmFreq/results')):
    os.mkdir('crmFreq/results')
    
rtDir = getEnvPath('rootpath')
modFile = rtDir + '/examples/common/crmWing.yaml'
modScrDir = rtDir + '/examples/modelGeneration'
if(not os.path.exists(modFile)):
    sys.path.append(modScrDir)
    import crmWing

constMod = Model()
constMod.fixDisplacement('rootNodes',ux=0.0,uy=0.0,uz=0.0)
constFile = 'crmFreq/constraints.yaml'
constMod.writeModelInput(constFile)

job = ASenDJob()
job.readModelInput(modFile)
job.readConstraints(constFile)
job.solve()
job.modalAnalysis('freq')
resFile = 'crmFreq/results/modalResults.yaml'
job.writeModalResults(resFile)
job.writeJobInput('crmFreq/job.yaml')
job.executeJob()

rp = ResultsProcessor(modFile,modalResFile=resFile)
for md in range(0,5):
    rp.animateModalSolution(md,defScaleFact=50.0)