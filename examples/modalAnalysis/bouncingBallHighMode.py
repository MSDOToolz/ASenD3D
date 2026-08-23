# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 10:35:45 2026

@author: evaande
"""

import os
import sys

from asendUtils.syst.pathTools import *
from asendUtils.job.ASenDJob import *

rt = getEnvPath('rootpath')

modelFile = rt + '/examples/common/bouncingBall.yaml'

if not os.path.exists(modelFile):
    scrptDir = rt + '/examples/modelGeneration'
    sys.path.append(scrptDir)
    import bouncingBall
    
if not os.path.exists('bouncingBall'):
    os.mkdir('bouncingBall')

if not os.path.exists('bouncingBall/results'):
    os.mkdir('bouncingBall/results')
    
myJob = ASenDJob()
myJob.readModelInput(modelFile)
myJob.solvePrep(nonlinearGeom=True, dynamic=True, explicit=True, solnHistDir='bouncingBall/results')
myJob.modalAnalysis(analysisType='highestFreq')
myJob.writeModalResults('bouncingBall/results/modal_results.csv')

myJob.writeJobInput('bouncingBall/job.yaml')

myJob.executeJob()