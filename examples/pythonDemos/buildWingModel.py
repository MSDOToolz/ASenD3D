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

nds = np.transpose(bData['nodes'])

#fig = px.scatter(x=nds[0], y=nds[1])
#fig.show()

cnt = input('continue?\n')
if(cnt == 'y'):
    mesh2 = Mesh2D(bData['nodes'],bData['elements'])
    sData = mesh2.createUnstructuredMesh('quad')
    #mt.plotShellMesh(sData)
    cnt = input('continue?\n')

if(cnt == 'y'):
    sData = mt.make3D(sData)
    mesh3 = Mesh3D(sData['nodes'],sData['elements'])
    mData = mesh3.createSweptMesh('inDirection',3,sweepDistance=0.1,axis=[0,0,-1.0])
    #mt.plotSolidMesh(mData)
    cnt = input('continue?\n')
    
if(cnt == 'y'):
    modInp = ModelInputBuilder()
    modInp.addMeshData(mData,meshType='solid')
    nEls = len(mData['elements'])
    modInp.addSet(setType='element',name='allEls',labelList=list(range(0,nEls)))
    modInp.addSection(elementSet='allEls',material='myMat')
    modInp.addMaterial(name='myMat',modulus=1000000.0,poissonRatio=0.3,shearModulus=500000.0)
    modInp.writeModelInput('solidWing.yaml')

