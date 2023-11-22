# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 16:18:08 2023

@author: evans
"""

from asendUtils.meshing.Surface import Surface
from asendUtils.meshing.Mesh3D import Mesh3D
from asendUtils.meshing.MeshTools import *
from asendUtils.model.Model import Model
from asendUtils.model.Constraint import Constraint
from asendUtils.model.Section import Section
from asendUtils.model.Material import Material


## Build model geometry

surf = Surface()
kp = [[0.0,0.0,0.0],
      [10.0,0.0,0.0],
      [10.0,1.0,0.0],
      [0.0,1.0,0.0]]
ne = [10,1,10,1]
surf.addShellRegion('quad1',kp,ne,elType='quad',meshMethod='structured')
surfMesh = surf.getSurfaceMesh()

surfMesh = getNodeSetInXYZRange(surfMesh,'xMin',xRange=[-0.1,0.1])
surfMesh = getNodeSetInXYZRange(surfMesh,'xMax',xRange=[9.9,10.1])

## Create model and add in mesh
myMod = Model()
myMod.addMeshData(surfMesh,meshType='shell')

## Define section
shSec = Section('shell')
shSec.setElementSet('all')
shSec.setOrientation(xDir=[1.0,0.0,0.0],xyVec=[0.0,1.0,0.0])
shSec.setZOffset(0.0)
shSec.addLayer(material='myMat',thickness=0.1,angle=0.0)

myMod.addSection(shSec)

## Define material
blkMat = Material('myMat')
blkMat.setOrthotropic(E1=1000000.0,E2=1000000.0,E3=1000000.0,nu12=0.0,nu13=0.0,nu23=0.0,G12=500000.0,G13=500000.0,G23=500000.0)
blkMat.setDensity(1.0)

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

## Write Input file

myMod.writeModelInput('shellBeam.yaml')
