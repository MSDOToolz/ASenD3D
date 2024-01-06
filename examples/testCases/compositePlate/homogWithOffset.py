# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 15:59:53 2023

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

if(not os.path.exists('homogWithOffset')):
    os.mkdir('homogWithOffset')
    
if(not os.path.exists('homogWithOffset/results')):
    os.mkdir('homogWithOffset/results')

## Build model geometry

surf = Surface()
kp = [[-1.0,-1.0,0.0],
      [1.0,-1.0,0.0],
      [1.0,1.0,0.0],
      [-1.0,1.0,0.0]]
ne = [1,1,1,1]
surf.addShellRegion('quad1',kp,ne,elType='quad',meshMethod='structured')
surfMesh = surf.getSurfaceMesh()

surfMesh = getNodeSetInXYZRange(surfMesh,'xMin',xRange=[-1.1,-0.9])
surfMesh = getNodeSetInXYZRange(surfMesh,'xMax',xRange=[0.9,1.1])

## Create model and add in mesh
myMod = Model()
myMod.addMeshData(surfMesh,meshType='shell')

## Define section
shSec = Section('shell')
shSec.setElementSet('all')
shSec.setOrientation(xDir=[1.0,0.0,0.0],xyVec=[0.0,1.0,0.0])
shSec.setZOffset(-1.0)
shSec.addLayer(material='myMat',thickness=0.1,angle=0.0)

myMod.addSection(shSec)

## Define material
blkMat = Material('myMat')
blkMat.setOrthotropic(E1=1000000.0,E2=1000000.0,E3=1000000.0,nu12=0.0,nu13=0.0,nu23=0.0,G12=500000.0,G13=500000.0,G23=500000.0)
blkMat.setDensity(1000.0)

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

## Define Loads

myMod.addNodalForce('xMax',F=[1.0,0.0,0.0],M=[0.0,0.0,0.0])

## Write Input file

myMod.writeModelInput('homogWithOffset/model.yaml')

myJob = ASenDJob()
myJob.readModelInput('homogWithOffset/model.yaml')
myJob.solve()
myJob.writeNodeResults('homogWithOffset/results/nodeResults.yaml',['displacement'])

myJob.writeJobInput('homogWithOffset/job.yaml')