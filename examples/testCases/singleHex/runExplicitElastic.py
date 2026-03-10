# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 13:03:00 2026

@author: evaande
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.model.Constraint import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *
# from asendUtils.ResultsProcessor import *

if(not os.path.exists('explicitElastic')):
    os.mkdir('explicitElastic')
    
if(not os.path.exists('explicitElastic/results')):
    os.mkdir('explicitElastic/results')

## Define constraints
myMod = Model()

myMod.fixDisplacement('xMin', ux=0.0, uy=0.0, uz=0.0)

## Write constraint file
myMod.writeModelInput('explicitElastic/elasticConstraints.yaml')


## Define initial state

initialVel = [[1, 5.47723, 0.0, 0.0, 0.0, 0.0, 0.0],
               [3, 5.47723, 0.0, 0.0, 0.0, 0.0, 0.0],
               [5, 5.47723, 0.0, 0.0, 0.0, 0.0, 0.0],
               [7, 5.47723, 0.0, 0.0, 0.0, 0.0, 0.0]]

myMod = Model()
myMod.addInitialState('velocity',initialVel)
myMod.writeModelInput('explicitElastic/initialState.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('singleHex.yaml')
myJob.readConstraints('explicitElastic/elasticConstraints.yaml')
myJob.readInitialState('explicitElastic/initialState.yaml')
myJob.readDesignVarInput('singleHexDVars.yaml')

myJob.solve(dynamic=True,explicit=True,timeStep=0.00573573,simPeriod=0.115,saveSolnHist=True,
            solnHistDir='explicitElastic/results/')


tSteps = list(range(0,20))
myJob.writeNodeResults('explicitElastic/results/nodeResults.csv',['displacement'],timeSteps=tSteps)

myJob.writeJobInput('explicitElastic/dynamicElasticJob.yaml')

myJob.executeJob()
