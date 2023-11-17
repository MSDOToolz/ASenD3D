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


## Build model geometry

surf = Surface()
kp = [[0.0,0.0,0.0],
      [1.0,0.0,0.0],
      [1.0,1.0,0.0],
      [0.0,1.0,0.0]]
ne = [1,1,1,1]
surf.addShellRegion('quad1',kp,ne,elType='quad',meshMethod='structured')
surfMesh = surf.getSurfaceMesh()

block = Mesh3D(surfMesh['nodes'],surfMesh['elements'])
blkMesh = block.createSweptMesh('inDirection',1,sweepDistance=1.0,axis=[0.0,0.0,1.0])

blkMesh = getNodeSetInXYZRange(blkMesh,'xMin',xRange=[-0.1,0.1])
blkMesh = getNodeSetInXYZRange(blkMesh,'xMax',xRange=[0.9,1.1])

## Create model and add in mesh
myMod = Model()
myMod.addMeshData(blkMesh)

## Define section
blkSec = Section('solid')
blkSec.setMaterial('myMat')
blkSec.setElementSet('all')

myMod.addSection(blkSec)

## Define material
blkMat = Material('myMat')
blkMat.setOrthotropic(E1=1000000.0,E2=1000000.0,E3=1000000.0,nu12=0.0,nu13=0.0,nu23=0.0,G12=500000.0,G13=500000.0,G23=500000.0)
blkMat.setDensity(1000.0)

myMod.addMaterial(blkMat)

## Write Input file

myMod.writeModelInput('singleHex.yaml')
