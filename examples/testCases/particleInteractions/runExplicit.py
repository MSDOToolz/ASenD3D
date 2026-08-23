# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 09:15:42 2026

@author: evaande
"""

import os
from asendUtils.job.ASenDJob import *

if not os.path.exists('runExplicit'):
    os.mkdir('runExplicit')
    
if not os.path.exists('runExplicit/results'):
    os.mkdir('runExplicit/results')

## Define job

myJob = ASenDJob()

myJob.readModelInput('collidingParticles.yaml')
myJob.solve(nonlinearGeom=True, dynamic=True, explicit=True, timeStep=0.025, simPeriod=3, solnHistFreq=4, solnHistDir='runExplicit/results')
numts = int(3/0.025)
ts = list(range(0, numts, 4))
myJob.writeNodeResults('runExplicit/results/node_results.csv', ['displacement'], timeSteps=ts)

myJob.writeJobInput('runExplicit/job.yaml')

myJob.executeJob()