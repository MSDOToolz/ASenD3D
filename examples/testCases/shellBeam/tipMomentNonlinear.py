# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 16:40:55 2023

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *

if(not os.path.exists('tipMomentLoading')):
    os.mkdir('tipMomentLoading')
    
if(not os.path.exists('tipMomentLoading/results')):
    os.mkdir('tipMomentLoading/results')

## Define loads

myMod = Model()
myMod.addNodalForce('xMax', [0.0,0.0,0], [0.0,-3.272492345,0.0])

## Write Load file
myMod.writeModelInput('tipMomentLoading/loads.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('shellBeam.yaml')
myJob.readLoads('tipMomentLoading/loads.yaml')

myJob.solve(nonlinearGeom=True,loadRampSteps=5)

myJob.writeNodeResults('tipMomentLoading/results/nodeResults.yaml',['displacement'])

myJob.writeJobInput('tipMomentLoading/job.yaml')

myJob.executeJob()