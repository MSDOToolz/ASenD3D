# -*- coding: utf-8 -*-
"""
Created on Sat Dec  2 19:38:45 2023

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *

if(not os.path.exists('buckling')):
    os.mkdir('buckling')
    
if(not os.path.exists('buckling/results')):
    os.mkdir('buckling/results')

## Define loads

myMod = Model()
myMod.addNodalForce('xMax', [-1.028088,0.0,0.0], [0.0,0.0,0.0])

## Write Load file
myMod.writeModelInput('buckling/loads.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('shellBeam.yaml')
myJob.readLoads('buckling/loads.yaml')

myJob.solve()
myJob.modalAnalysis()

myJob.writeModalResults('buckling/results/bucklingResults.yaml')

myJob.writeJobInput('buckling/job.yaml')

#myJob.executeJob()