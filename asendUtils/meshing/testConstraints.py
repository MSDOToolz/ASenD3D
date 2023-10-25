# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:10:42 2023

@author: evans
"""

from SurfaceClass import *
from Mesh3DClass import *
from MeshTools import *

## With 2 separate meshes
newSurf = Surface()
keyPts = [[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0],[0.0,1.0,0.0]]
numEls = [1,1,1,1]
newSurf.addShellRegion('quad1',keyPts,numEls,name='block1',elType='quad',meshMethod='structured')
surfMesh = newSurf.getSurfaceMesh()

newMesh = Mesh3D(surfMesh['nodes'],surfMesh['elements'])
block1Mesh = newMesh.createSweptMesh('inDirection',1,sweepDistance=1.0,axis=[0,0,1])

newSurf = Surface()
keyPts = [[0.0,0.0,1.0],[0.5,0.0,1.0],[0.5,0.5,1.0],[0.0,0.5,1.0]]
numEls = [1,1,1,1]
newSurf.addShellRegion('quad1',keyPts,numEls,name='block2',elType='quad',meshMethod='structured')
surfMesh = newSurf.getSurfaceMesh()

newMesh = Mesh3D(surfMesh['nodes'],surfMesh['elements'])
block2Mesh = newMesh.createSweptMesh('inDirection',1,sweepDistance=1.0,axis=[0,0,1])

constraints = tie2MeshesConstraints(block2Mesh,block1Mesh,0.01)

## With 2 sets of a single mesh
newSurf = Surface()
keyPts = [[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0],[0.0,1.0,0.0]]
numEls = [1,1,1,1]
newSurf.addShellRegion('quad1',keyPts,numEls,name='block1',elType='quad',meshMethod='structured')

keyPts = [[0.0,0.0,1.0],[0.5,0.0,1.0],[0.5,0.5,1.0],[0.0,0.5,1.0]]
numEls = [1,1,1,1]
newSurf.addShellRegion('quad1',keyPts,numEls,name='block2',elType='quad',meshMethod='structured')

surfMesh = newSurf.getSurfaceMesh()

surfMesh = getMatchingNodeSets(surfMesh)
extSets = getExtrudedSets(surfMesh,1)

newMesh = Mesh3D(surfMesh['nodes'],surfMesh['elements'])
bothBlocksMesh = newMesh.createSweptMesh('inDirection',1,sweepDistance=1.0,axis=[0,0,1])
bothBlocksMesh['sets'] = extSets

constraints2 = tie2SetsConstraints(bothBlocksMesh,'block2','block1',0.01)