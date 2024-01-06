# -*- coding: utf-8 -*-
"""
Created on Sun Dec 31 21:29:30 2023

@author: evans
"""

from asendUtils.meshing.Surface import Surface
from asendUtils.meshing.Mesh3D import Mesh3D
from asendUtils.meshing.MeshTools import *
from asendUtils.model.Model import Model
from asendUtils.model.Section import Section
from asendUtils.model.Material import Material
from asendUtils.model.Constraint import Constraint
from asendUtils.job.ASenDJob import *

if(not os.path.exists('twoPlyPlusMinus45')):
    os.mkdir('twoPlyPlusMinus45')
    
if(not os.path.exists('twoPlyPlusMinus45/results')):
    os.mkdir('twoPlyPlusMinus45/results')

## Build model geometry

surf = Surface()
kp = [[-1.0,-1.0,0.0],
      [1.0,-1.0,0.0],
      [1.0,1.0,0.0],
      [-1.0,1.0,0.0]]
ne = [1,1,1,1]
surf.addShellRegion('quad1',kp,ne,elType='quad',meshMethod='structured')
surfMesh = surf.getSurfaceMesh()

surfMesh = getNodeSetInXYZRange(surfMesh,'xyMin',xRange=[-1.1,-0.9],yRange=[-1.1,-0.9])
surfMesh = getNodeSetInXYZRange(surfMesh, 'xMinyMax',xRange=[-1.1,-0.9],yRange=[0.9,1.1])
surfMesh = getNodeSetInXYZRange(surfMesh,'xMax',xRange=[0.9,1.1])

## Create model and add in mesh
myMod = Model()
myMod.addMeshData(surfMesh,meshType='shell')

## Define section
shSec = Section('shell')
shSec.setElementSet('all')
shSec.setOrientation(xDir=[1.0,0.0,0.0],xyVec=[0.0,1.0,0.0])
shSec.addLayer(material='myMat',thickness=0.05,angle=45.0)
shSec.addLayer(material='myMat',thickness=0.05,angle=-45.0)

myMod.addSection(shSec)

## Define material
blkMat = Material('myMat')
blkMat.setOrthotropic(E1=10000000.0,E2=1000000.0,E3=1000000.0,nu12=0.0,nu13=0.0,nu23=0.0,G12=5000000.0,G13=5000000.0,G23=5000000.0)
blkMat.setDensity(1000.0)

myMod.addMaterial(blkMat)

## Define constraints
blkConst = Constraint('displacement')
blkConst.addTerm('xyMin', 1, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

blkConst = Constraint('displacement')
blkConst.addTerm('xyMin', 2, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

blkConst = Constraint('displacement')
blkConst.addTerm('xyMin', 3, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

blkConst = Constraint('displacement')
blkConst.addTerm('xyMin', 4, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

blkConst = Constraint('displacement')
blkConst.addTerm('xyMin', 5, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

blkConst = Constraint('displacement')
blkConst.addTerm('xyMin', 6, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

blkConst = Constraint('displacement')
blkConst.addTerm('xMinyMax', 1, 1.0)
blkConst.setRHS(0.0)
myMod.addConstraint(blkConst)

## Define Loads

myMod.addNodalForce('xMax',F=[1.0,0.0,0.0],M=[0.0,0.0,0.0])

## Write Input file

myMod.writeModelInput('twoPlyPlusMinus45/model.yaml')

myJob = ASenDJob()
myJob.readModelInput('twoPlyPlusMinus45/model.yaml')
myJob.solve()
myJob.writeNodeResults('twoPlyPlusMinus45/results/nodeResults.yaml',['displacement'])

myJob.writeJobInput('twoPlyPlusMinus45/job.yaml')