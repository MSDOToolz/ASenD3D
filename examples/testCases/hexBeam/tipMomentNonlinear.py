# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 17:05:14 2024

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *
import numpy as np

if(not os.path.exists('tipMomentNonlinear')):
    os.mkdir('tipMomentNonlinear')
    
if(not os.path.exists('tipMomentNonlinear/results')):
    os.mkdir('tipMomentNonlinear/results')

## Define loads

myMod = Model()
ld = (65.4498469/(11*np.sqrt(2)))*np.array([1.0,0.0,1.0])
myMod.addNodalForce('xMaxzMin',F=ld,M=[0.0,0.0,0.0])
ld = -ld
myMod.addNodalForce('xzMax', F=ld, M=[0.0,0.0,0.0])

## Write Load file
myMod.writeModelInput('tipMomentNonlinear/loads.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('hexBeam.yaml')
myJob.readLoads('tipMomentNonlinear/loads.yaml')

myJob.solve(nonlinearGeom=True,loadRampSteps=5)

myJob.writeNodeResults('tipMomentNonlinear/results/nodeResults.yaml',['displacement'])
myJob.writeElementResults('tipMomentNonlinear/results/elementResults.yaml',['strain','stress'])

myJob.writeJobInput('tipMomentNonlinear/job.yaml')

myJob.executeJob()