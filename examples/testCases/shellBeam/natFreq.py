# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 16:41:54 2023

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *
# from asendUtils.ResultsProcessor import *

if(not os.path.exists('natFreq')):
    os.mkdir('natFreq')
    
if(not os.path.exists('natFreq/results')):
    os.mkdir('natFreq/results')

## Define job
myJob = ASenDJob()
myJob.readModelInput('shellBeam.yaml')

##myJob.solve()
myJob.solvePrep()
myJob.modalAnalysis(analysisType='frequency')

myJob.writeModalResults('natFreq/results/natFreqResults.csv')

myJob.writeJobInput('natFreq/job.yaml')

myJob.executeJob()

# rp = ResultsProcessor('shellBeam.yaml')
# rp.loadModalVals('natFreq/results/natFreqResults.csv')
# rp.loadModalVec('natFreq/results/natFreqResults.csv', mode=0)
# rp.plotModalResults()
# rp.animateModalSolution()