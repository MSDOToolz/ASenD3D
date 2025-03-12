# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 16:36:16 2023

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *
# from asendUtils.ResultsProcessor import *

if(not os.path.exists('transverseTipLoading')):
    os.mkdir('transverseTipLoading')
    
if(not os.path.exists('transverseTipLoading/results')):
    os.mkdir('transverseTipLoading/results')

## Define loads

myMod = Model()
myMod.addNodalForce('xMax',F=[0.0,0.0,0.0909091],M=[0.0,0.0,0.0])

## Write Load file
myMod.writeModelInput('transverseTipLoading/loads.yaml')

## Define objective
myObj = Objective()
myObj.addObjectiveTerm('displacement',operator='powerNorm',
                       component=3,nodeSet='xMax',coefficient=0.0454545,exponent=1.0)

myObj.writeInput('transverseTipLoading/objective.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('hexBeam.yaml')
myJob.readLoads('transverseTipLoading/loads.yaml')
myJob.readDesignVarInput('hexBeamDVars.yaml')
myJob.readObjectiveInput('transverseTipLoading/objective.yaml')

myJob.solve(solnHistDir='transverseTipLoading/results')
myJob.calcObjGradient()

myJob.writeNodeResults('transverseTipLoading/results/nodeResults.csv',['displacement'])
myJob.writeElementResults('transverseTipLoading/results/elementResults.csv',['strain','stress'])
myJob.writeObjective('transverseTipLoading/results/objectiveResults.csv')

myJob.writeJobInput('transverseTipLoading/job.yaml')

myJob.executeJob()

# rp = ResultsProcessor('hexBeam.yaml')
# rp.plotNodeResults('displacement',3,deformed=True)