# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 15:40:49 2023

@author: evans
"""

import os
import sys
from asendUtils.model.Model import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *

if(not os.path.exists('centrifugalLoading')):
    os.mkdir('centrifugalLoading')
    
if(not os.path.exists('centrifugalLoading/results')):
    os.mkdir('centrifugalLoading/results')

## Define loads

myMod = Model()
myMod.addCentrifugalForce('all', [0.0,0.5,0.0], [0.0,0.0,1.0], 1.0)

## Write Load file
myMod.writeModelInput('centrifugalLoading/loads.yaml')

## Define objective
myObj = Objective()
myObj.addObjectiveTerm('displacement',operator='powerNorm',
                       component=1,nodeSet='xMax',coefficient=0.5,exponent=1.0)

myObj.writeInput('centrifugalLoading/objective.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('shellBeam.yaml')
myJob.readLoads('centrifugalLoading/loads.yaml')
myJob.readDesignVarInput('shellBeamDVars.yaml')
myJob.readObjectiveInput('centrifugalLoading/objective.yaml')

myJob.solve()
myJob.calcObjGradient()

myJob.writeNodeResults('centrifugalLoading/results/nodeResults.yaml',['displacement'])
myJob.writeElementResults('centrifugalLoading/results/elementResults.yaml',['strain','stress'])
myJob.writeObjective('centrifugalLoading/results/objectiveResults.yaml')

myJob.writeJobInput('centrifugalLoading/job.yaml')

myJob.executeJob()