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
from SurfaceClass import *
import MeshTools as mt
from ModelInputBuilderClass import *

surf = Surface()

kp = [[0.1,0.1,0.0],[1.0,0.0,0.0],[1.0,0.0,-0.1],[0.1,0.1,-0.1],[0.5,0.085,0.0],[1.0,0.0,-0.05],[0.5,0.085,-0.1],[0.1,0.1,-0.05],[0.5,0.085,-0.05]]
surf.addShellRegion(regType='quad2',keyPts=kp,numEls=[24,3,24,3],name='suctionSide',elType='quad',meshMethod='structured')

kp = [[0.0,0.0,0.0],[0.1,0.1,0.0],[0.1,0.1,-0.1],[0.0,0.0,-0.1],[0.02,0.08,0.0],[0.1,0.1,-0.05],[0.02,0.08,-0.1],[0.0,0.0,-0.05],[0.02,0.08,-0.05]]
surf.addShellRegion(regType='quad2',keyPts=kp,numEls=[6,3,6,3],name='leadingEdge1',elType='quad',meshMethod='structured')

kp = [[0.0,0.0,0.0],[0.1,-0.1,0.0],[0.1,-0.1,-0.1],[0.0,0.0,-0.1],[0.02,-0.08,0.0],[0.1,-0.1,-0.05],[0.02,-0.08,-0.1],[0.0,0.0,-0.05],[0.02,-0.08,-0.05]]
surf.addShellRegion(regType='quad2',keyPts=kp,numEls=[6,3,6,3],name='leadingEdge2',elType='quad',meshMethod='structured')

kp = [[0.1,-0.1,0.0],[1.0,0.0,0.0],[1.0,0.0,-0.1],[0.1,-0.1,-0.1],[0.5,-0.085,0.0],[1.0,0.0,-0.05],[0.5,-0.085,-0.1],[0.1,-0.1,-0.05],[0.5,-0.085,-0.05]]
surf.addShellRegion(regType='quad2',keyPts=kp,numEls=[24,3,24,3],name='pressureSide',elType='quad',meshMethod='structured')

mData = surf.getSurfaceMesh()

#mt.plotShellMesh(mData)

modInp = ModelInputBuilder()

modInp.addMeshData(mData,meshType='shell')
nEls = len(mData['elements'])
modInp.addSet(setType='element',name='allEls',labelList=list(range(0,nEls)))
modInp.addSection(elementSet='allEls',secType='shell',orientation=[0.0,0.0,-1.0,-1.0,0.0,0.0],layup=[['myMat',0.01,0.0]])
modInp.addMaterial(name='myMat',modulus=1000000.0,poissonRatio=0.3,shearModulus=500000.0)
modInp.writeModelInput('shellWing2.yaml')
