# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 20:11:26 2023

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.model.Constraint import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *
# from asendUtils.ResultsProcessor import *

if(not os.path.exists('dynamicElastic')):
    os.mkdir('dynamicElastic')
    
if(not os.path.exists('dynamicElastic/results')):
    os.mkdir('dynamicElastic/results')

## Define constraints
myMod = Model()

myMod.fixDisplacement('xMin', ux=0.0, uy=0.0, uz=0.0)

## Write constraint file
myMod.writeModelInput('dynamicElastic/elasticConstraints.yaml')


## Define initial state

initialVel = [[1, 5.47723, 0.0, 0.0, 0.0, 0.0, 0.0],
               [3, 5.47723, 0.0, 0.0, 0.0, 0.0, 0.0],
               [5, 5.47723, 0.0, 0.0, 0.0, 0.0, 0.0],
               [7, 5.47723, 0.0, 0.0, 0.0, 0.0, 0.0]]

myMod = Model()
myMod.addInitialState('velocity',initialVel)
myMod.writeModelInput('dynamicElastic/initialState.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('singleHex.yaml')
myJob.readConstraints('dynamicElastic/elasticConstraints.yaml')
myJob.readInitialState('dynamicElastic/initialState.yaml')
myJob.readDesignVarInput('singleHexDVars.yaml')

myJob.solve(dynamic=True,timeStep=0.00573573,simPeriod=0.115,saveSolnHist=True,
            solnHistDir='dynamicElastic/results/')


tSteps = list(range(0,20))
myJob.writeNodeResults('dynamicElastic/results/nodeResults.csv',['displacement'],timeSteps=tSteps)

myJob.writeJobInput('dynamicElastic/dynamicElasticJob.yaml')

myJob.executeJob()

# rp = ResultsProcessor('singleHex.yaml')
# rp.animateNodeResults('dynamicElastic/results/nodeResults.csv','displacement',tSteps,deformed=True)
