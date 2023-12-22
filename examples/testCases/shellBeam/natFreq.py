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

if(not os.path.exists('natFreq')):
    os.mkdir('natFreq')
    
if(not os.path.exists('natFreq/results')):
    os.mkdir('natFreq/results')

## Define job
myJob = ASenDJob()
myJob.readModelInput('shellBeam.yaml')

myJob.solve()
myJob.modalAnalysis(analysisType='frequency')

myJob.writeModalResults('natFreq/results/natFreqResults.yaml')

myJob.writeJobInput('natFreq/job.yaml')

myJob.executeJob()