# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 10:05:44 2023

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *

if(not os.path.exists('distributedLoading')):
    os.mkdir('distributedLoading')
    
if(not os.path.exists('distributedLoading/results')):
    os.mkdir('distributedLoading/results')

## Define loads

myMod = Model()
myMod.addSurfacePressure('all',0.1,[0.0,0.0,-1.0])

## Write Load file
myMod.writeModelInput('distributedLoading/loads.yaml')

## Define objective
myObj = Objective()
myObj.addObjectiveTerm('displacement',operator='powerNorm',
                       component=3,nodeSet='xMax',coefficient=0.5,exponent=1.0)

myObj.writeInput('distributedLoading/objective.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('shellBeam.yaml')
myJob.readLoads('distributedLoading/loads.yaml')
myJob.readDesignVarInput('shellBeamDVars.yaml')
myJob.readObjectiveInput('distributedLoading/objective.yaml')

myJob.solve()
myJob.calcObjGradient()

myJob.writeNodeResults('distributedLoading/results/nodeResults.yaml',['displacement'])
myJob.writeElementResults('distributedLoading/results/elementResults.yaml',['strain','stress'])
myJob.writeObjective('distributedLoading/results/objectiveResults.yaml')

myJob.writeJobInput('distributedLoading/job.yaml')

myJob.executeJob()