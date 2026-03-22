# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 18:35:24 2026

@author: evans
"""

import os
import sys

from asendUtils.syst.pathTools import *
from asendUtils.job.ASenDJob import *
from asendUtils.ResultsProcessor import *

rt = getEnvPath('rootpath')

modelFile = rt + '/examples/common/bouncingBall.yaml'

if not os.path.exists(modelFile):
    scrptDir = rt + '/examples/modelGeneration'
    sys.path.append(scrptDir)
    import bouncingBall
    
if not os.path.exists('bouncingBallExplicit'):
    os.mkdir('bouncingBallExplicit')

if not os.path.exists('bouncingBallExplicit/results'):
    os.mkdir('bouncingBallExplicit/results')
    
myJob = ASenDJob()
myJob.readModelInput(modelFile)

myJob.solve(nonlinearGeom=True, dynamic=True, explicit=True, timeStep=5.0e-6, simPeriod=0.65, solnHistFreq=500, solnHistDir='bouncingBallExplicit/results')
ts = list(range(0,130000,500))
#myJob.solve(nonlinearGeom=True, dynamic=True, explicit=True, timeStep=5.0e-6, simPeriod=5.0e-4, solnHistFreq=1, solnHistDir='bouncingBallExplicit/results')
#ts = list(range(0,100))

nodeFile = 'bouncingBallExplicit/results/node_results.csv'
myJob.writeNodeResults(nodeFile, ['displacement'], timeSteps=ts)

myJob.writeJobInput('bouncingBallExplicit/job.yaml')

myJob.executeJob()

rp = ResultsProcessor(modelFile)
rp.animateNodeResults(nodeFile, 'displacement', ts, component=3, deformed=True)