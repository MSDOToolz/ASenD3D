# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 21:15:29 2024

@author: evans
"""

import numpy as np
import copy as cp
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

## Input parameters

volumeFraction = 0.6
xDim = 0.05
yDim = 1.0
zDim = yDim*np.sqrt(3)
elSize = 0.05
elLayers = 1

## Calculated parameters

totArea = yDim*zDim
fibArea = totArea*volumeFraction
fibRad = np.sqrt(0.5*fibArea/np.pi)
hW = 0.5*yDim
hHt = 0.5*zDim

faceSurf = Surface()

## Mesh center fiber

fiberBd = Boundary2D()
kp = [[0.,-fibRad],
      [0.,fibRad],
      [0.,-fibRad]]
nEl = int(np.ceil(2.0*np.pi*fibRad/elSize))
fiberBd.addSegment('arc',kp,nEl)
bdData = fiberBd.getBoundaryMesh()

mesher = Mesh2D(bdData['nodes'],bdData['elements'])
mesh = mesher.createUnstructuredMesh('quad')
mesh = make3D(mesh)

faceSurf.addMesh(mesh,name='centerFiber')

## Mesh quarter fiber corners

qtrBd = Boundary2D()

kp = [[0.,0.],
      [fibRad,0.]]
nEl = int(np.ceil(fibRad/elSize))
qtrBd.addSegment('line',kp,nEl)

sr2 = fibRad/np.sqrt(2.0)
kp = [[fibRad,0.], 
      [sr2,sr2],
      [0.,fibRad]]
nEl = int(np.ceil(fibRad*0.5*np.pi/elSize))
qtrBd.addSegment('arc',kp,nEl)

kp = [[0.,fibRad],
      [0.,0.]]
nEl = int(np.ceil(fibRad/elSize))
qtrBd.addSegment('line',kp,nEl)

bdData = qtrBd.getBoundaryMesh()
mesher = Mesh2D(bdData['nodes'],bdData['elements'])
mesh = mesher.createUnstructuredMesh('quad')
mesh = make3D(mesh)

c1 = cp.deepcopy(mesh)
c1 = translateMesh(c1,[-hW,-hHt,0.])
faceSurf.addMesh(c1,'corner1Fiber')

c2 = cp.deepcopy(mesh)
c2 = translateMesh(c2,[hW,-hHt,0.])
c2 = rotateMesh(c2,[hW,-hHt,0.],[0.,0.,1.],90)
faceSurf.addMesh(c2,'corner2Fiber')

c3 = cp.deepcopy(mesh)
c3 = translateMesh(c3,[hW,hHt,0.])
c3 = rotateMesh(c3,[hW,hHt,0.],[0.,0.,1.],180)
faceSurf.addMesh(c3,'corner3Fiber')

c4 = cp.deepcopy(mesh)
c4 = translateMesh(c4,[-hW,hHt,0.])
c4 = rotateMesh(c4,[-hW,hHt,0.],[0.,0.,1.],270)
faceSurf.addMesh(c4,'corner4Fiber')

## Mesh matrix

matBd = Boundary2D()

kp = [[0.,-fibRad],
      [0.,fibRad],
      [0.,-fibRad]]
nEl = int(np.ceil(2.0*np.pi*fibRad/elSize))
matBd.addSegment('arc',kp,nEl)

kp = [[-hW,-hHt+fibRad],
      [-hW+sr2,-hHt+sr2],
      [-hW+fibRad,-hHt]]
nEl = int(np.ceil(0.5*np.pi*fibRad/elSize))
matBd.addSegment('arc',kp,nEl)

kp = [[-hW+fibRad,-hHt],
      [hW-fibRad,-hHt]]
nEl = int(np.ceil((yDim-2.0*fibRad)/elSize))
matBd.addSegment('line',kp,nEl)

kp = [[hW-fibRad,-hHt],
      [hW-sr2,-hHt+sr2],
      [hW,-hHt+fibRad]]
nEl = int(np.ceil(0.5*np.pi*fibRad/elSize))
matBd.addSegment('arc',kp,nEl)

kp = [[hW,-hHt+fibRad],
      [hW,hHt-fibRad]]
nEl = int(np.ceil((zDim-2.0*fibRad)/elSize))
matBd.addSegment('line',kp,nEl)

kp = [[hW,hHt-fibRad],
      [hW-sr2,hHt-sr2],
      [hW-fibRad,hHt]]
nEl = int(np.ceil(0.5*np.pi*fibRad/elSize))
matBd.addSegment('arc',kp,nEl)

kp = [[hW-fibRad,hHt],
      [-hW+fibRad,hHt]]
nEl = int(np.ceil((yDim-2.0*fibRad)/elSize))
matBd.addSegment('line',kp,nEl)

kp = [[-hW+fibRad,hHt],
      [-hW+sr2,hHt-sr2],
      [-hW,hHt-fibRad]]
nEl = int(np.ceil(0.5*np.pi*fibRad/elSize))
matBd.addSegment('arc',kp,nEl)

kp = [[-hW,hHt-fibRad],
      [-hW,-hHt+fibRad]]
nEl = int(np.ceil((zDim-2.0*fibRad)/elSize))
matBd.addSegment('line',kp,nEl)

bdData = matBd.getBoundaryMesh()
mesher = Mesh2D(bdData['nodes'],bdData['elements'])
mesh = mesher.createUnstructuredMesh('quad')
mesh = make3D(mesh)

faceSurf.addMesh(mesh,name='matrix')

## Generate face mesh

faceMesh = faceSurf.getSurfaceMesh()

faceMesh = rotateMesh(faceMesh,[0.,0.,0.],[0.,1.,0.],90.0)
faceMesh = rotateMesh(faceMesh,[0.,0.,0.],[1.,0.,0.],90.0)

## Generate 3D cell mesh

cellMesher = Mesh3D(faceMesh['nodes'],faceMesh['elements'])
cellMesh = cellMesher.createSweptMesh('inDirection',elLayers,sweepDistance=xDim,axis=[1.,0.,0.])

## Form additional sets

extSets = getExtrudedSets(faceMesh,elLayers)
cellMesh['sets'] = extSets

cellMesh = getPeriodicSets(cellMesh,xDim,yDim,zDim)

unSets = ['centerFiber',
          'corner1Fiber',
          'corner2Fiber',
          'corner3Fiber',
          'corner4Fiber']
cellMesh = getElementSetUnion(cellMesh,unSets,'allFiber')

## Build model

myMod = Model()
myMod.addMeshData(cellMesh,meshType='solid')

newSec = Section('solid')
newSec.setElementSet('allFiber')
newSec.setMaterial('carbon')
myMod.addSection(newSec)

newSec = Section('matrix')
newSec.setElementSet('matrix')
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
