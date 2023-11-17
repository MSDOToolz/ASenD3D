# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 20:08:59 2023

@author: evans
"""
import os
import sys
from asendUtils.model.Model import *
from asendUtils.model.Constraint import *
from asendUtils.objective.Objective import *
from asendUtils.job.ASenDJob import *

if(not os.path.exists('staticElastic')):
    os.mkdir('staticElastic')
    
if(not os.path.exists('staticElastic/results')):
    os.mkdir('staticElastic/results')

## Define constraints
myMod = Model()
blkConst = Constraint('displacement')
blkConst.addTerm('xMin', 1, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

blkConst = Constraint('displacement')
blkConst.addTerm('xMin', 2, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

blkConst = Constraint('displacement')
blkConst.addTerm('xMin', 3, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

blkConst = Constraint('displacement')
blkConst.addTerm('xMin', 4, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

blkConst = Constraint('displacement')
blkConst.addTerm('xMin', 5, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

blkConst = Constraint('displacement')
blkConst.addTerm('xMin', 6, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

## Write constraint file
myMod.writeModelInput('staticElastic/elasticConstraints.yaml')

## Define loads

myMod = Model()
myMod.addNodalForce('xMax',F=[0.5,0.0,0.0],M=[0.0,0.0,0.0])

## Write Load file
myMod.writeModelInput('staticElastic/staticNodalLoads.yaml')

## Define objective
myObj = Objective()
myObj.addObjectiveTerm('displacement',operator='powerNorm',
                       component=1,nodeSet='xMax',coefficient=0.5,exponent=1.0)
myObj.addObjectiveTerm('strain',operator='powerNorm',component=1,layer=0,
                       elementSet='all',coefficient=1.0,exponent=1.0)
myObj.addObjectiveTerm('stress',operator='powerNorm',component=1,layer=0,
                       elementSet='all',coefficient=1.0,exponent=1.0)

myObj.writeInput('staticElastic/staticObjective.yaml')

## Define job
myJob = ASenDJob()
myJob.readModelInput('singleShell.yaml')
myJob.readConstraints('staticElastic/elasticConstraints.yaml')
myJob.readLoads('staticElastic/staticNodalLoads.yaml')
myJob.readDesignVarInput('singleShellDVars.yaml')
myJob.readObjectiveInput('staticElastic/staticObjective.yaml')

myJob.solve()
myJob.calcObjGradient()

myJob.writeNodeResults('staticElastic/results/nodeResults.yaml',['displacement'])
myJob.writeElementResults('staticElastic/results/elementResults.yaml',['strain','stress'])
myJob.writeObjective('staticElastic/results/objectiveResults.yaml')

myJob.writeJobInput('staticElastic/staticElasticJob.yaml')

myJob.executeJob()
