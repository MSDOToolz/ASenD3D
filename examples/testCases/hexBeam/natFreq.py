# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 08:35:28 2024

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *
import numpy as np

if(not os.path.exists('natFreq')):
    os.mkdir('natFreq')
    
if(not os.path.exists('natFreq/results')):
    os.mkdir('natFreq/results')

## Define job
myJob = ASenDJob()
myJob.readModelInput('hexBeam.yaml')

myJob.solve()
myJob.modalAnalysis(analysisType='freq')

myJob.writeModalResults('natFreq/results/natFreqResults.yaml')

myJob.writeJobInput('natFreq/job.yaml')

myJob.executeJob()