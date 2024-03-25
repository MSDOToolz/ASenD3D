# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:28:25 2024

@author: evans
"""

from asendUtils.meshing.Boundary2D import *
from asendUtils.meshing.Mesh2D import *
from asendUtils.meshing.MeshTools import *
from asendUtils.visualization.plotlyUtils import *
from asendUtils.model.Section import *
from asendUtils.model.Model import *
from asendUtils.syst.pathTools import *

diskBound = Boundary2D()

keypts = [[-1.0,0.0],
          [1.0,0.0],
          [-1.0,0.0]]
numEls = 80
diskBound.addSegment('arc',keypts,numEls,'diskBoundary')

boundMesh = diskBound.getBoundaryMesh()

diskMesher = Mesh2D(boundMesh['nodes'],boundMesh['elements'])
meshData = diskMesher.createUnstructuredMesh('quad')

meshData = make3D(meshData)

plotShellMesh(meshData)

meshData = getElementSetInRadius(meshData,[0.0,0.0,0.0],1.01,'allDiskEls')

meshData = getNearestNodes(meshData,[0.0,0.0,0.0],1,'centerNode')
meshData = getNodeSetInRadius(meshData,[0.0,0.0,0.0],0.2,'impactNodes')
meshData = getNodeSetInRadius(meshData,[0.0,0.0,0.0],1.01,'allDiskNodes')
meshData = getNodeSetInRadius(meshData,[0.0,0.0,0.0],0.99,'innerNodes')
meshData = subtractNodeSet(meshData,'allDiskNodes','innerNodes','boundaryNodes')

meshData = addFreeNodes(meshData,[[0.0,0.0,1.0]],'projectileNode')
meshData = addMassElements(meshData,'projectileNode','projectileElement')
meshData = addForceElements(meshData,'projectileNode','impactNodes','forceElements')

model = Model()
model.addMeshData(meshData,'shell')

diskSec = Section('shell')
diskSec.setElementSet('allDiskEls')
diskSec.addLayer('nylon',0.001)
model.addSection(diskSec)

forceSec = Section('forceField')
forceSec.setElementSet('forceElements')
potCoef = -1.0*getPotFrcCoef(0.01,10.0,6,0.1)
forceSec.setPotentialField(potCoef,6)
model.addSection(forceSec)

massSec = Section('mass')
massSec.setElementSet('projectileElement')
massSec.setMassPerElement(0.01)
model.addSection(massSec)

mat = Material('nylon')
mat.setDensity(1140.0)
E = 2800000000.0
nu = 0.3
mat.setIsotropic(E,nu)
model.addMaterial(mat)

rtDir = getEnvPath('rootpath')
modFile = rtDir + '/examples/common/shellDiskProjectile.yaml'
model.writeModelInput(modFile)