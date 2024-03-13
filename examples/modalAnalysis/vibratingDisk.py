# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 18:50:28 2024

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.syst.pathTools import *
from asendUtils.job.ASenDJob import *
from asendUtils.ResultsProcessor import *

if(not os.path.exists('vibratingDisk')):
    os.mkdir('vibratingDisk')
    
if(not os.path.exists('vibratingDisk/results')):
    os.mkdir('vibratingDisk/results')
    
rtDir = getEnvPath('rootpath')
modFile = rtDir + '/examples/common/shellDisk.yaml'
modScrDir = rtDir + '/examples/modelGeneration'
if(not os.path.exists(modFile)):
    sys.path.append(modScrDir)
    import shellDisk

constMod = Model()
constMod.fixDisplacement('boundaryNodes',ux=0.0,uy=0.0,uz=0.0)
constFile = 'vibratingDisk/constraints.yaml'
constMod.writeModelInput(constFile)

job = ASenDJob()
job.readModelInput(modFile)
job.readConstraints(constFile)
job.solve()
job.modalAnalysis('freq')
resFile = 'vibratingDisk/results/modalResults.yaml'
job.writeModalResults(resFile)
job.writeJobInput('vibratingDisk/job.yaml')
job.executeJob()

rp = ResultsProcessor(modFile,modalResFile=resFile)
for md in range(0,10):
    rp.animateModalSolution(md,defScaleFact=5.0)
