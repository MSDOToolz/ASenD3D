# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 19:40:02 2024

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *
import numpy as np

if(not os.path.exists('buckling')):
    os.mkdir('buckling')
    
if(not os.path.exists('buckling/results')):
    os.mkdir('buckling/results')

## Define loads

myMod = Model()
ld = 20.56176/2
myMod.addSurfacePressure('all',ld,[1.0,0.0,0.0])

## Write Load file
myMod.writeModelInput('buckling/loads.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('hexBeam.yaml')
myJob.readLoads('buckling/loads.yaml')

myJob.solve()
myJob.modalAnalysis()

myJob.writeModalResults('buckling/results/bucklingResults.yaml')

myJob.writeJobInput('buckling/job.yaml')

myJob.executeJob()