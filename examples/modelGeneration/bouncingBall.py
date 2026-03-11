# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 13:44:46 2026

@author: evaande
"""

import numpy as np

from asendUtils.meshing.Surface import *
from asendUtils.meshing.MeshTools import *
from asendUtils.model.Model import *
from asendUtils.model.Section import *
from asendUtils.model.Material import *
from asendUtils.model.Interaction import *
from asendUtils.visualization.plotlyUtils import *
from asendUtils.syst.pathTools import *

## Input parameters

plateDim = 1.0
ballDiam = 0.5
centerHt = 0.75
elSize = 0.02

## end input

# Floor plate

plateSurf = Surface()

hd = 0.5*plateDim
kp = [[-hd, -hd, 0.0],
      [hd, -hd, 0.0],
      [hd, hd, 0.0],
      [-hd, hd, 0.0]]
ne = int(plateDim/elSize + 1)
nels = [ne, ne, ne, ne]

plateSurf.addShellRegion('quad1', kp, nels, name='plate', meshMethod='structured')

plateMesh = plateSurf.getSurfaceMesh()

# plotShellMesh(plateMesh)

## Ball

ballSurf = Surface()

rad = 0.5*ballDiam
kp = np.array([[0.0, 0.0, centerHt],
      [rad, 0.0, centerHt],
      [0., rad, centerHt]])
nels = [int(np.pi*ballDiam/elSize)]
ballSurf.addShellRegion('sphere', kp, nels, name='upperBall')

kp[2,1] = -rad

ballSurf.addShellRegion('sphere', kp, nels, name='lowerBall')

ballMesh = ballSurf.getSurfaceMesh()

#plotShellMesh(ballMesh)
#wholeDeal = mergeMeshes(plateMesh, ballMesh)
#plotShellMesh(wholeDeal)

ballMesh = getElementSetUnion(ballMesh, ['upperBall','lowerBall'], 'wholeBall')
ballMesh = getMatchingNodeSet(ballMesh, 'wholeBall', 'wholeBall')

## Build model

myMod = Model()
myMod.addMeshData(plateMesh, meshType='shell')
myMod.addMeshData(ballMesh, meshType='shell')

## Sections
newSec = Section('shell')
newSec.setElementSet('plate')
newSec.addLayer('nylon', elSize)
newSec.setOrientation([1., 0., 0.], [0., 1., 0.])
myMod.addSection(newSec)

newSec = Section('shell')
newSec.setElementSet('wholeBall')
newSec.addLayer('nylon', 0.5*elSize)
newSec.setOrientation([1., 0., 0.], [0., 1., 0.])
myMod.addSection(newSec)

## Material
newMat = Material('nylon')
newMat.setIsotropic(2.0e+9, 0.3)
newMat.setDensity(1200.0)
myMod.addMaterial(newMat)

## Constraint
myMod.fixDisplacement('plate', ux=0., uy=0., uz=0.)

## Load

myMod.addGravityForce('wholeBall', G3=-9.8)

## Interaction

intctn = contactInteraction(1.0e+6, 0.5, elSize, name='plate_ball_cont', nodeSet1='plate', nodeSet2='wholeBall', maxDistance=2*elSize)
myMod.addInteraction(intctn)

## Write input

rt = getEnvPath('rootpath')
fn = rt + '/examples/common/bouncingBall.yaml'
myMod.writeModelInput(fn)

