# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 21:15:29 2024

@author: evans
"""

import numpy as np
from asendUtils.meshing.Boundary2D import *
from asendUtils.meshing.Mesh2D import *
from asendUtils.meshing.Mesh3D import *
from asendUtils.meshing.MeshTools import *
from asendUtils.meshing.Surface import *
from asendUtils.visualization.plotlyUtils import *
from asendUtils.model.Model import *
from asendUtils.model.Section import *
from asendUtils.model.Material import *
from asendUtils.syst.pathTools import *
from asendUtils.ResultsProcessor import *

volumeFraction = 0.6

faceBd = Boundary2D()
kp = [[0.,-0.5],
      [0.,0.5],
      [0.,-0.5]]
nEls = 24
faceBd.addSegment('arc',kp,nEls,'fiberBoundary')
boundMesh = faceBd.getBoundaryMesh()
fiberMesher = Mesh2D(boundMesh['nodes'],boundMesh['elements'])
fibMesh = fiberMesher.createUnstructuredMesh('quad')
fibMesh = make3D(fibMesh)

kp = [[-1.,-1.],
      [1.,-1.]]
nEls = 12
faceBd.addSegment('line',kp,nEls,'bottom')
kp = [[1.,-1.],
      [1.,1.]]
nEls = 12
faceBd.addSegment('line',kp,nEls,'right')
kp = [[1.,1.],
      [-1.,1.]]
nEls = 12
faceBd.addSegment('line',kp,nEls,'top')
kp = [[-1.,1.],
      [-1.,-1.]]
nEls = 12
faceBd.addSegment('line',kp,nEls,'left')

boundMesh = faceBd.getBoundaryMesh()

matrixMesher = Mesh2D(boundMesh['nodes'],boundMesh['elements'])
matMesh = matrixMesher.createUnstructuredMesh('quad')
matMesh = make3D(matMesh)

faceSurf = Surface()
faceSurf.addMesh(fibMesh,name='fiberElements')
faceSurf.addMesh(matMesh,name='matrixElements')

faceMesh = faceSurf.getSurfaceMesh()

#plotShellMesh(faceMesh)

cellMesher = Mesh3D(faceMesh['nodes'],faceMesh['elements'])
cellMesh = cellMesher.createSweptMesh('inDirection',2,sweepDistance=0.1,axis=[0.,0.,1.])

#plotSolidMesh(cellMesh)

extSets = getExtrudedSets(faceMesh,2)
cellMesh['sets'] = extSets

cellMesh = getPeriodicSets(cellMesh,2.,2.,0.1)

myMod = Model()
myMod.addMeshData(cellMesh,meshType='solid')

newSec = Section('solid')
newSec.setElementSet('fiberElements')
newSec.setMaterial('carbon')
myMod.addSection(newSec)

newSec = Section('matrix')
newSec.setElementSet('matrixElements')
newSec.setMaterial('epoxy')
myMod.addSection(newSec)

newMat = Material('carbon')
newMat.setIsotropic(2.0e+11,0.2)
myMod.addMaterial(newMat)

newMat = Material('epoxy')
newMat.setIsotropic(5.0e+9,0.3)
myMod.addMaterial(newMat)

modFile = getEnvPath('rootpath') + '/examples/common/fiberUnitCell.yaml'

myMod.writeModelInput(modFile)

rp = ResultsProcessor(modFile)
rp.plotElementProperty()
