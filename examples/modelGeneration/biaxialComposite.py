# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 15:13:44 2024

@author: evans
"""

import numpy as np
import copy
from asendUtils.meshing.Mesh2D import *
from asendUtils.meshing.Mesh3D import *
from asendUtils.meshing.Boundary2D import *
from asendUtils.meshing.Surface import *
from asendUtils.meshing.MeshTools import *
from asendUtils.visualization.plotlyUtils import *
from asendUtils.model.Model import *
from asendUtils.model.Section import *
from asendUtils.model.Material import *
from asendUtils.syst.pathTools import *
from asendUtils.ResultsProcessor import *

## Set laminate specifications

xDim = [-0.5,0.5]
yDim = [-0.5,0.5]
zDim = [-0.15,0.15]
numXtows = 1
numYtows = 1
thkRatio = 2
numStacks = 1
volFract = 0.6
elSize = 0.035

## Calculate additional specs

xLen = xDim[1] - xDim[0]
yLen = yDim[1] - yDim[0]
zLen = zDim[1] - zDim[0]
layThk = (zDim[1] - zDim[0])/numStacks
sqvf = np.math.sqrt(volFract)
xWid = sqvf*(yDim[1] - yDim[0])/numXtows
yWid = sqvf*(xDim[1] - xDim[0])/numYtows
yThk = sqvf*layThk/(thkRatio + 1.0)
xThk = yThk*thkRatio

## Draw x-tow cell cross section

xSurf = Surface()

radCrv = 0.25*xThk
rt2 = radCrv*(1.0 - 1.0/np.math.sqrt(2.0))
hW = 0.5*xWid
hT = 0.5*xThk

towBound = Boundary2D()

kp = [[-hW,-hT+radCrv],
      [-hW+rt2,-hT+rt2],
      [-hW+radCrv,-hT]]
nEls = 2*int(np.math.ceil(radCrv/elSize))
towBound.addSegment('arc',kp,nEls,'c1')

kp = [[-hW+radCrv,-hT],
      [hW-radCrv,-hT]]
nEls = int(np.math.ceil(xWid/elSize))
towBound.addSegment('line',kp,nEls,'e1')

kp = [[hW-radCrv,-hT],
      [hW-rt2,-hT+rt2],
      [hW,-hT+radCrv]]
nEls = 2*int(np.math.ceil(radCrv/elSize))
towBound.addSegment('arc',kp,nEls,'c2')

kp = [[hW,-hT+radCrv],
      [hW,hT-radCrv]]
nEls = int(np.math.ceil(xThk/elSize))
towBound.addSegment('line',kp,nEls,'e2')

kp = [[hW,hT-radCrv],
      [hW-rt2,hT-rt2],
      [hW-radCrv,hT]]
nEls = 2*int(np.math.ceil(radCrv/elSize))
towBound.addSegment('arc',kp,nEls,'c3')

kp = [[hW-radCrv,hT],
      [-hW+radCrv,hT]]
nEls = int(np.math.ceil(xWid/elSize))
towBound.addSegment('line',kp,nEls,'e3')

kp = [[-hW+radCrv,hT],
      [-hW+rt2,hT-rt2],
      [-hW,hT-radCrv]]
nEls = 2*int(np.math.ceil(radCrv/elSize))
towBound.addSegment('arc',kp,nEls,'c4')

kp = [[-hW,hT-radCrv],
      [-hW,-hT+radCrv]]
nEls = int(np.math.ceil(xThk/elSize))
towBound.addSegment('line',kp,nEls,'e4')

tbData = towBound.getBoundaryMesh()

towMesher = Mesh2D(tbData['nodes'],tbData['elements'])
towMesh = towMesher.createUnstructuredMesh('quad')

towMesh = make3D(towMesh)

hW = 0.5*(yDim[1] - yDim[0])/numXtows
hT = 0.5*xThk/sqvf

kp = [[-hW,-hT],
      [hW,-hT]]
nEls = int(np.ceil(2*hW/elSize))
towBound.addSegment('line',kp,nEls,'outerE1')
yExtrudeEls = numXtows*nEls

kp = [[hW,-hT],
      [hW,hT]]
nEls = int(np.ceil(2*hT/elSize))
towBound.addSegment('line',kp,nEls,'outerE2')

kp = [[hW,hT],
      [-hW,hT]]
nEls = int(np.ceil(2*hW/elSize))
towBound.addSegment('line',kp,nEls,'outerE3')

kp = [[-hW,hT],
      [-hW,-hT]]
nEls = int(np.ceil(2*hT/elSize))
towBound.addSegment('line',kp,nEls,'outerE4')

tbData = towBound.getBoundaryMesh()

matrixMesher = Mesh2D(tbData['nodes'],tbData['elements'])
matrixMesh = matrixMesher.createUnstructuredMesh('quad')
matrixMesh = make3D(matrixMesh)

xShft = 2.0*hW
yShft = layThk

ct = 0
for i in range(0,numXtows):
    for j in range(0,numStacks):
        sct = str(ct)
        transVec = [i*xShft,j*yShft,0.0]
        newTow = copy.deepcopy(towMesh)
        newTow = translateMesh(newTow,transVec)
        nm = 'xTow' + sct
        xSurf.addMesh(newTow,nm)
        newMat = copy.deepcopy(matrixMesh)
        newMat = translateMesh(newMat,transVec)
        nm = 'xMatrix' + sct
        xSurf.addMesh(newMat,nm)
        ct = ct + 1
        
xTFaceMesh = xSurf.getSurfaceMesh()

xTFaceMesh = rotateMesh(xTFaceMesh,[0.,0.,0.],[1.,0.,0.],90.0)
xTFaceMesh = rotateMesh(xTFaceMesh,[0.,0.,0.],[0.,0.,1.],90.0)
nodes = xTFaceMesh['nodes']
maxX = np.max(nodes[:,0])
maxY = np.max(nodes[:,1])
maxZ = np.max(nodes[:,2])
transVec = [(xDim[1]-maxX),(yDim[1]-maxY),(zDim[1]-maxZ)]
xTFaceMesh = translateMesh(xTFaceMesh, transVec)

## Draw y-tow cell cross section

ySurf = Surface()

radCrv = 0.25*yThk
rt2 = radCrv*(1.0 - 1.0/np.math.sqrt(2.0))
hW = 0.5*yWid
hT = 0.5*yThk

towBound = Boundary2D()

kp = [[-hW,-hT+radCrv],
      [-hW+rt2,-hT+rt2],
      [-hW+radCrv,-hT]]
nEls = 2*int(np.math.ceil(radCrv/elSize))
towBound.addSegment('arc',kp,nEls,'c1')

kp = [[-hW+radCrv,-hT],
      [hW-radCrv,-hT]]
nEls = int(np.math.ceil(yWid/elSize))
towBound.addSegment('line',kp,nEls,'e1')

kp = [[hW-radCrv,-hT],
      [hW-rt2,-hT+rt2],
      [hW,-hT+radCrv]]
nEls = 2*int(np.math.ceil(radCrv/elSize))
towBound.addSegment('arc',kp,nEls,'c2')

kp = [[hW,-hT+radCrv],
      [hW,hT-radCrv]]
nEls = int(np.math.ceil(yThk/elSize))
towBound.addSegment('line',kp,nEls,'e2')

kp = [[hW,hT-radCrv],
      [hW-rt2,hT-rt2],
      [hW-radCrv,hT]]
nEls = 2*int(np.math.ceil(radCrv/elSize))
towBound.addSegment('arc',kp,nEls,'c3')

kp = [[hW-radCrv,hT],
      [-hW+radCrv,hT]]
nEls = int(np.math.ceil(yWid/elSize))
towBound.addSegment('line',kp,nEls,'e3')

kp = [[-hW+radCrv,hT],
      [-hW+rt2,hT-rt2],
      [-hW,hT-radCrv]]
nEls = 2*int(np.math.ceil(radCrv/elSize))
towBound.addSegment('arc',kp,nEls,'c4')

kp = [[-hW,hT-radCrv],
      [-hW,-hT+radCrv]]
nEls = int(np.math.ceil(yThk/elSize))
towBound.addSegment('line',kp,nEls,'e4')

tbData = towBound.getBoundaryMesh()

towMesher = Mesh2D(tbData['nodes'],tbData['elements'])
towMesh = towMesher.createUnstructuredMesh('quad')

towMesh = make3D(towMesh)

hW = 0.5*(xDim[1] - xDim[0])/numYtows
hT = 0.5*yThk/sqvf

kp = [[-hW,-hT],
      [hW,-hT]]
nEls = int(np.ceil(2*hW/elSize))
towBound.addSegment('line',kp,nEls,'outerE1')
xExtrudeEls = numYtows*nEls

kp = [[hW,-hT],
      [hW,hT]]
nEls = int(np.ceil(2*hT/elSize))
towBound.addSegment('line',kp,nEls,'outerE2')

kp = [[hW,hT],
      [-hW,hT]]
nEls = int(np.ceil(2*hW/elSize))
towBound.addSegment('line',kp,nEls,'outerE3')

kp = [[-hW,hT],
      [-hW,-hT]]
nEls = int(np.ceil(2*hT/elSize))
towBound.addSegment('line',kp,nEls,'outerE4')

tbData = towBound.getBoundaryMesh()

matrixMesher = Mesh2D(tbData['nodes'],tbData['elements'])
matrixMesh = matrixMesher.createUnstructuredMesh('quad')
matrixMesh = make3D(matrixMesh)

xShft = 2.0*hW
yShft = layThk

ct = 0
for i in range(0,numYtows):
    for j in range(0,numStacks):
        sct = str(ct)
        transVec = [i*xShft,j*yShft,0.0]
        newTow = copy.deepcopy(towMesh)
        newTow = translateMesh(newTow,transVec)
        nm = 'yTow' + sct
        ySurf.addMesh(newTow,nm)
        newMat = copy.deepcopy(matrixMesh)
        newMat = translateMesh(newMat,transVec)
        nm = 'yMatrix' + sct
        ySurf.addMesh(newMat,nm)
        ct = ct + 1
        
yTFaceMesh = ySurf.getSurfaceMesh()

yTFaceMesh = rotateMesh(yTFaceMesh,[0.,0.,0.],[1.,0.,0.],90.0)
nodes = yTFaceMesh['nodes']
maxX = np.max(nodes[:,0])
maxY = np.max(nodes[:,1])
minZ = np.min(nodes[:,2])
transVec = [(xDim[1]-maxX),(yDim[1]-maxY),(zDim[0]-minZ)]
yTFaceMesh = translateMesh(yTFaceMesh, transVec)

# Extrude the x and y tow faces to make a 3D mesh

xTMesher = Mesh3D(xTFaceMesh['nodes'],xTFaceMesh['elements'])
xLen = xDim[1] - xDim[0]
xTMesh = xTMesher.createSweptMesh('inDirection',xExtrudeEls,sweepDistance=xLen,axis=[-1.0,0.0,0.0])
exSets = getExtrudedSets(xTFaceMesh,xExtrudeEls)
xTMesh['sets'] = exSets

yTMesher = Mesh3D(yTFaceMesh['nodes'],yTFaceMesh['elements'])
yLen = yDim[1] - yDim[0]
yTMesh = yTMesher.createSweptMesh('inDirection',yExtrudeEls,yLen,axis=[0.,-1.,0.])
exSets = getExtrudedSets(yTFaceMesh,yExtrudeEls)
yTMesh['sets'] = exSets

wholeMesh = mergeMeshes(xTMesh,yTMesh)

## Define additional sets

xTowSets = list()
yTowSets = list()
matrixSets = list()

for es in wholeMesh['sets']['element']:
    nm = es['name']
    if('xTow' in nm):
        xTowSets.append(nm)
    elif('yTow' in nm):
        yTowSets.append(nm)
    elif('Matrix' in nm):
        matrixSets.append(nm)

wholeMesh = getElementSetUnion(wholeMesh,xTowSets,'allXTows')
wholeMesh = getElementSetUnion(wholeMesh,yTowSets,'allYTows')
wholeMesh = getElementSetUnion(wholeMesh,matrixSets,'allMatrix')

wholeMesh = getPeriodicSets(wholeMesh,xLen,yLen,zLen)

## Create model and add data

mod = Model()

mod.addMeshData(wholeMesh)

newSec = Section('solid')
newSec.setElementSet('allXTows')
newSec.setMaterial('carbonFiber')
newSec.setOrientation([1.,0.,0.],[0.,1.,0.])
mod.addSection(newSec)

newSec = Section('solid')
newSec.setElementSet('allYTows')
newSec.setMaterial('glassFiber')
newSec.setOrientation([0.,1.,0.],[-1.,0.,0.])
mod.addSection(newSec)

newSec = Section('solid')
newSec.setElementSet('allMatrix')
newSec.setMaterial('epoxy')
mod.addSection(newSec)

newMat = Material('carbonFiber')
newMat.setIsotropic(2.0e+11,0.2)
mod.addMaterial(newMat)

newMat = Material('glassFiber')
newMat.setIsotropic(1.0e+11,0.2)
mod.addMaterial(newMat)

newMat = Material('epoxy')
newMat.setIsotropic(5.0e+9,0.3)
mod.addMaterial(newMat)

rtDir = getEnvPath('rootpath')
modFile = rtDir + '/examples/common/biaxialComposite.yaml'
mod.writeModelInput(modFile)

rp = ResultsProcessor(modFile)
rp.plotElementProperty()