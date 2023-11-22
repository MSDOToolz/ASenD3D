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
ne = [100,10,100,10]
surf.addShellRegion('quad1',kp,ne,elType='quad',meshMethod='structured')
surfMesh = surf.getSurfaceMesh()

beam = Mesh3D(surfMesh['nodes'],surfMesh['elements'])
beamMesh = beam.createSweptMesh('inDirection',1,sweepDistance=0.1,axis=[0.0,0.0,1.0])

beamMesh = getNodeSetInXYZRange(beamMesh,'xMin',xRange=[-0.1,0.1])
beamMesh = getNodeSetInXYZRange(beamMesh,'xMax',xRange=[9.9,10.1])

## Create model and add in mesh
myMod = Model()
myMod.addMeshData(beamMesh,meshType='solid')

## Define section
beamSec = Section('solid')
beamSec.setElementSet('all')
beamSec.setOrientation(xDir=[1.0,0.0,0.0],xyVec=[0.0,1.0,0.0])
beamSec.setMaterial('myMat')

myMod.addSection(beamSec)

## Define material
myMat = Material('myMat')
myMat.setOrthotropic(E1=1000000.0,E2=1000000.0,E3=1000000.0,nu12=0.0,nu13=0.0,nu23=0.0,G12=500000.0,G13=500000.0,G23=500000.0)
myMat.setDensity(1.0)

myMod.addMaterial(myMat)

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

## Write Input file

myMod.writeModelInput('hexBeam.yaml')
