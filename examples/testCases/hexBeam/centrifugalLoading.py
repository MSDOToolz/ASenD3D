# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 16:48:04 2024

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
                       component=1,nodeSet='xMax',coefficient=0.0454545,exponent=1.0)

myObj.writeInput('centrifugalLoading/objective.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('hexBeam.yaml')
myJob.readLoads('centrifugalLoading/loads.yaml')
myJob.readDesignVarInput('hexBeamDVars.yaml')
myJob.readObjectiveInput('centrifugalLoading/objective.yaml')

myJob.solve()
myJob.calcObjGradient()

myJob.writeNodeResults('centrifugalLoading/results/nodeResults.yaml',['displacement'])
myJob.writeElementResults('centrifugalLoading/results/elementResults.yaml',['strain','stress'])
myJob.writeObjective('centrifugalLoading/results/objectiveResults.yaml')

myJob.writeJobInput('centrifugalLoading/job.yaml')

myJob.executeJob()