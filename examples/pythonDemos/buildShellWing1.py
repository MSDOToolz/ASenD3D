import os
import sys
import numpy as np
import plotly.express as px

currDir = os.getcwd()

if('ASenD3D' in currDir):
    dirLst = currDir.split('ASenD3D')
    mainDir = dirLst[0] + 'ASenD3D'
    meshDir = dirLst[0] + 'ASenD3D/pythonUtilities/meshing'
    utilsDir = dirLst[0] + 'ASenD3D/pythonUtilities'
    
sys.path.insert(0,utilsDir)
sys.path.insert(0,meshDir)

from Mesh2DClass import *
from Mesh3DClass import *
from Boundary2DClass import *
import MeshTools as mt
from ModelInputBuilderClass import *


## Draw out cross section boundary of wing

boundary = Boundary2D()
boundary.addSegment('arc',[[0.1,-0.1],[0.0,0.0],[0.1,0.1]],16)
boundary.addSegment('curve',[[0.1,0.1],[0.55,0.085],[1.0,0.0]],24)
boundary.addSegment('curve',[[1.0,0.0],[0.55,-0.085],[0.1,-0.1]],24)
bData = boundary.getBoundaryMesh()

bData = mt.make3D(bData)
mesh2 = Mesh2D(bData['nodes'],bData['elements'])
sData = mesh2.createSweptMesh('inDirection',3,sweepDistance=0.1,axis=[0.0,0.0,-1.0])
mt.plotShellMesh(sData)

modInp = ModelInputBuilder()

modInp.addMeshData(sData,meshType='shell')
nEls = len(sData['elements'])
modInp.addSet(setType='element',name='allEls',labelList=list(range(0,nEls)))
modInp.addSection(elementSet='allEls',secType='shell',orientation=[0.0,0.0,-1.0,-1.0,0.0,0.0],layup=[['myMat',0.01,0.0]])
modInp.addMaterial(name='myMat',modulus=1000000.0,poissonRatio=0.3,shearModulus=500000.0)
modInp.writeModelInput('shellWing1.yaml')
