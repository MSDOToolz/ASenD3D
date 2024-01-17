# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 16:23:17 2024

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *

if(not os.path.exists('gravitationalLoading')):
    os.mkdir('gravitationalLoading')
    
if(not os.path.exists('gravitationalLoading/results')):
    os.mkdir('gravitationalLoading/results')

## Define loads

myMod = Model()
myMod.addGravityForce('all',[0.0,0.0,1.0])

## Write Load file
myMod.writeModelInput('gravitationalLoading/loads.yaml')

## Define objective
myObj = Objective()
myObj.addObjectiveTerm('displacement',operator='powerNorm',
                       component=3,nodeSet='xMax',coefficient=0.0454545,exponent=1.0)

myObj.writeInput('gravitationalLoading/objective.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('hexBeam.yaml')
myJob.readLoads('gravitationalLoading/loads.yaml')
myJob.readDesignVarInput('hexBeamDVars.yaml')
myJob.readObjectiveInput('gravitationalLoading/objective.yaml')

myJob.solve()
myJob.calcObjGradient()

myJob.writeNodeResults('gravitationalLoading/results/nodeResults.yaml',['displacement'])
myJob.writeElementResults('gravitationalLoading/results/elementResults.yaml',['strain','stress'])
myJob.writeObjective('gravitationalLoading/results/objectiveResults.yaml')

myJob.writeJobInput('gravitationalLoading/job.yaml')

myJob.executeJob()