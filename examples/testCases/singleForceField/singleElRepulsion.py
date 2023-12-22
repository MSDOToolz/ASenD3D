# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 10:31:53 2023

@author: evans
"""

import os
from asendUtils.model.Model import *
from asendUtils.model.Section import *
from asendUtils.model.Constraint import *
from asendUtils.job.ASenDJob import *
from asendUtils.meshing.MeshTools import *

if(not os.path.exists('singleElRepulsion')):
    os.mkdir('singleElRepulsion')
    
if(not os.path.exists('singleElRepulsion/results')):
    os.mkdir('singleElRepulsion/results')
    

frcNodes = [[0.0,0.0,0.0],
            [1.0,0.0,0.0]]

frcEls = [[0,1]]

meshData = meshFromScratch(frcNodes, frcEls)
meshData = getNodeSetInRadius(meshData, [0.5,0.0,0.0], 0.6, 'allNodes')
meshData = getElementSetInRadius(meshData, [0.5,0.0,0.0], 0.6, 'forceElement')
meshData = addMassElements(meshData, 'allNodes', 'massEls')

myMod = Model()
myMod.addMeshData(meshData,meshType='frcFld')

newSec = Section('frcFld')
newSec.setElementSet('forceElement')
newSec.setPotentialField(-1.0,2)
newSec.setDampingField(1.0,2)
myMod.addSection(newSec)

newSec = Section('mass')
newSec.setElementSet('massEls')
newSec.setMassPerElement(1.0)
myMod.addSection(newSec)

const = Constraint('displacement')
const.addTerm(0,1,1.0)
const.setRHS(0.0)
myMod.addConstraint(const)

const = Constraint('displacement')
const.addTerm(0,2,1.0)
const.setRHS(0.0)
myMod.addConstraint(const)

const = Constraint('displacement')
const.addTerm(0,3,1.0)
const.setRHS(0.0)
myMod.addConstraint(const)

myMod.writeModelInput('singleElRepulsion/model.yaml')

myJob = ASenDJob()
myJob.readModelInput('singleElRepulsion/model.yaml')
myJob.solve(nonlinearGeom=True,dynamic=True,timeStep=0.121,simPeriod=12.1,saveSolnHist=True,
            solnHistDir='singleElRepulsion/results')
myJob.writeNodeResults('singleElRepulsion/results/nodeResults.yaml',['displacement'],timeSteps=[0,25,50,75,100])

myJob.writeJobInput('singleElRepulsion/job.yaml')