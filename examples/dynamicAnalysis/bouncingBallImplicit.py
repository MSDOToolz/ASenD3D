# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 07:08:39 2026

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
    
if not os.path.exists('bouncingBallImplicit'):
    os.mkdir('bouncingBallImplicit')

if not os.path.exists('bouncingBallImplicit/results'):
    os.mkdir('bouncingBallImplicit/results')
    
myJob = ASenDJob()
myJob.readModelInput(modelFile)
myJob.solve(nonlinearGeom=True, dynamic=True, timeStep=0.005, simPeriod=0.65, solnHistFreq=1, solnHistDir='bouncingBallImplicit/results', 
            lumpMass=True, solverMethod='iterative', solverBlockDim=10)
ts = list(range(0,130))
nodeFile = 'bouncingBallImplicit/results/node_results.csv'
myJob.writeNodeResults(nodeFile, ['displacement'], timeSteps=ts)

myJob.writeJobInput('bouncingBallImplicit/job.yaml')

myJob.executeJob()

rp = ResultsProcessor(modelFile)
rp.animateNodeResults(nodeFile, 'displacement', ts, component=3, deformed=True)