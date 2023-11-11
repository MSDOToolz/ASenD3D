# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 09:46:34 2023

@author: evans
"""

from asendUtils.meshing.Surface import Surface
from asendUtils.meshing.Mesh3D import Mesh3D
from asendUtils.meshing.MeshTools import *
from asendUtils.model.Model import Model
from asendUtils.model.Section import Section
from asendUtils.model.Material import Material
from asendUtils.model.Constraint import Constraint

## Temporary Imports
# import os
# import sys

# currDir = os.getcwd()

# if('ASenD3D' in currDir):
#     dirLst = currDir.split('ASenD3D')
#     mainDir = dirLst[0] + 'ASenD3D'
#     meshDir = dirLst[0] + 'ASenD3D/src/asendUtils/meshing'
#     utilsDir = dirLst[0] + 'ASenD3D/src/asendUtils'
    
# sys.path.insert(0,utilsDir)
# sys.path.insert(0,meshDir)

# from Surface import Surface
# from Mesh3D import Mesh3D
# from MeshTools import *

## End temp imports

testCases = ['staticElastic','staticThermal','dynamicElastic','dynamicThermal']

## Build model geometry

surf = Surface()
kp = [[0.0,0.0,0.0],
      [1.0,0.0,0.0],
      [1.0,1.0,0.0],
      [0.0,1.0,0.0]]
ne = [1,1,1,1]
surf.addShellRegion('quad1',kp,ne,elType='quad',meshMethod='structured')
surfMesh = surf.getSurfaceMesh()
print(surfMesh)

block = Mesh3D(surfMesh['nodes'],surfMesh['elements'])
blkMesh = block.createSweptMesh('inDirection',1,sweepDistance=1.0,axis=[0.0,0.0,1.0])
print(blkMesh)

blkMesh = getNodeSetInXYZRange(blkMesh,'xMin',xRange=[-0.1,0.1])
blkMesh = getNodeSetInXYZRange(blkMesh,'xMax',xRange=[0.9,1.1])

## Create model and add in mesh
myMod = Model()
myMod.addMeshData(blkMesh)

## Define section
blkSec = Section('solid')
blkSec.setMaterial('myMat')
blkSec.setOrientation()
blkSec.setElementSet('all')

myMod.addSection(blkSec)

## Define material
blkMat = Material('myMat')
blkMat.setOrthotropic(E1=1000000.0,E2=1000000.0,E3=1000000.0,nu12=0.3,nu13=0.3,nu23=0.3,G12=500000.0,G13=5000000.0,G23=500000.0)

myMod.addMaterial(blkMat)

## Define constraints
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

## Define loads

myMod.addNodalForce('xMax',F=[0.25,0.0,0.0])

myMod.writeModelInput('singleHex.yaml')
