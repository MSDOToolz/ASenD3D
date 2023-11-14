# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 20:08:59 2023

@author: evans
"""
import os
import sys
from asendUtils.model.Model import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *

## Define loads

myMod = Model()
myMod.addNodalForce('xMax',F=[0.25,0.0,0.0],M=[0.0,0.0,0.0])

## Write Load file
myMod.writeModelInput('staticNodalLoads.yaml')

## Define objective
myObj = Objective()
myObj.addObjectiveTerm('displacement',operator='powerNorm',
                       component=1,nodeSet='xMax',coefficient=0.25,exponent=1.0)
myObj.addObjectiveTerm('strain',operator='powerNorm',component=1,
                       elementSet='all',coefficient=1.0,exponent=1.0)
myObj.addObjectiveTerm('stress',operator='powerNorm',component=1,
                       elementSet='all',coefficient=1.0,exponent=1.0)

myObj.writeInput('staticObjective.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('singleHex.yaml')
myJob.readLoads('staticNodalLoads.yaml')
myJob.readDesignVarInput('singleHexDVars.yaml')
myJob.readObjectiveInput('staticObjective.yaml')

myJob.solve()

if(not os.path.exists('Results')):
    os.mkdir('Results')

if(not os.path.extist('staticElastic')):
    os.mkdir('staticElastic')
    
myJob.writeNodeResults('Results/staticElastic/nodeResults.yaml',['displacement'])
myJob.writeElementResults('Results/staticElastic/elementResults.yaml',['strain','stress'])

myJob.writeJobInput('staticElastic.yaml')
