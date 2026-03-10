# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 11:29:44 2026

@author: evaande
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.model.Constraint import *
from asendUtils.job.ASenDJob import *
# from asendUtils.ResultsProcessor import *

if(not os.path.exists('highestMode')):
    os.mkdir('highestMode')
    
if(not os.path.exists('highestMode/results')):
    os.mkdir('highestMode/results')

## Define constraints
myMod = Model()

myMod.fixDisplacement('xMin', ux=0.0, uy=0.0, uz=0.0)

## Write constraint file
myMod.writeModelInput('highestMode/elasticConstraints.yaml')

## Define job

myJob = ASenDJob()
myJob.readModelInput('singleHex.yaml')
myJob.readConstraints('highestMode/elasticConstraints.yaml')
myJob.solvePrep(dynamic=True, explicit=True)
myJob.modalAnalysis(analysisType='highestFreq')
#myJob.modalAnalysis(analysisType='frequency', numModes=2)
myJob.writeModalResults('highestMode/results/modalRes.csv')

myJob.writeJobInput('highestMode/job.yaml')

myJob.executeJob()