# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 16:08:41 2026

@author: evaande
"""

import os
from asendUtils.job.ASenDJob import *
from asendUtils.ResultsProcessor import *

if not os.path.exists('runImplicit'):
    os.mkdir('runImplicit')
    
if not os.path.exists('runImplicit/results'):
    os.mkdir('runImplicit/results')

## Define job

myJob = ASenDJob()

myJob.readModelInput('collidingParticles.yaml')
myJob.solve(nonlinearGeom=True, dynamic=True, timeStep=0.025, simPeriod=3, solnHistFreq=4, solnHistDir='runImplicit/results')
numts = int(3/0.025)
ts = list(range(0, numts, 4))
myJob.writeNodeResults('runImplicit/results/node_results.csv', ['displacement'], timeSteps=ts)

myJob.writeJobInput('runImplicit/job.yaml')

myJob.executeJob()

rp = ResultsProcessor('collidingParticles.yaml')
rp.animateNodeResults('runImplicit/results/node_results.csv', 'displacement', ts, deformed=True, massElOptns={'showAsDots': True, 'refSize': 0.1, 'refMass': 1.0})