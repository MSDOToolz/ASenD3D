# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 15:41:41 2024

@author: evans
"""

import numpy as np
import plotly.express as px
import yaml
from scipy import interpolate
from asendUtils.meshing.Surface import *
from asendUtils.meshing.MeshTools import *
from asendUtils.model.Model import *
from asendUtils.model.Section import *
from asendUtils.model.Material import *
from asendUtils.syst.pathTools import *
from asendUtils.visualization.plotlyUtils import *
from asendUtils.ResultsProcessor import *

## Set basic parameters
span = 26.49
rootChrd = 11.84
breakChrd = 7.243
tipChrd = 2.745
sweep = 31.0
sweepSt = 0.05
elSize = 0.2
nonDimSpan = [0.0, 0.05, 0.35, 1.0]

## End parameters

## Local functions

def transformKp(kp,chord,yrot,x,z):
    tKp = chord*kp
    for pi, pt in enumerate(tKp):
        rPt = rotateVector(pt,[0.0,1.0,0.0],yrot)
        rPt = rPt + np.array([x,0.0,z])
        tKp[pi] = rPt
    return tKp

def getRegNumEls(kp,eSz):
    nEls = list()
    n = int(np.ceil(np.linalg.norm(kp[0] - kp[1])/elSize))
    nEls.append(n)
    n = int(np.ceil(np.linalg.norm(kp[1] - kp[2])/elSize))
    nEls.append(n)
    n = int(np.ceil(np.linalg.norm(kp[2] - kp[3])/elSize))
    nEls.append(n)
    n = int(np.ceil(np.linalg.norm(kp[3] - kp[0])/elSize))
    nEls.append(n)
    return nEls

## Main procedure

chord = [rootChrd,rootChrd,breakChrd,tipChrd]
tanSw = np.math.tan(sweep*0.0174533)
xShft = list()
zShft = list()
for sp in nonDimSpan:
    if(sp < sweepSt):
        x = 0.0
    else:
        x = (sp - sweepSt)*span*tanSw
    xShft.append(x)
    zShft.append(sp*span)
yRot = [0.0,0.0,31.0,31.0]

iSpan = np.linspace(0.0,1.0,41)
chrdFun = interpolate.interp1d(nonDimSpan,chord,kind='linear')
iChrd = chrdFun(iSpan)
xFun = interpolate.interp1d(nonDimSpan,xShft,kind='linear')
iX = xFun(iSpan)
zFun = interpolate.interp1d(nonDimSpan,zShft,kind='linear')
iZ = zFun(iSpan)
rotFun = interpolate.interp1d(nonDimSpan,yRot,kind='linear')
iRot = rotFun(iSpan)

rtDir = getEnvPath('rootpath')

afPath = rtDir + '/examples/common/crmAfKeyPts.yaml'
inFile = open(afPath,'r')
afData = yaml.safe_load(inFile)
inFile.close()
lenKp = len(afData['keypts'])
refKp = np.zeros((lenKp,3),dtype=float)
refKp[:,0:2] = np.array(afData['keypts'])

wingSurf = Surface()
lenSpan = len(iSpan)
for si in range(0, (lenSpan-1), 2):
    thisKp = transformKp(refKp,iChrd[si],iRot[si],iX[si],iZ[si])
    nextKp = transformKp(refKp,iChrd[si+1],iRot[si+1],iX[si+1],iZ[si+1])
    next2Kp = transformKp(refKp,iChrd[si+2],iRot[si+2],iX[si+2],iZ[si+2])
    for ci in range(0, (lenKp-1), 2):
        kp = [thisKp[ci],
              thisKp[ci+2],
              next2Kp[ci+2],
              next2Kp[ci],
              thisKp[ci+1],
              nextKp[ci+2],
              next2Kp[ci+1],
              nextKp[ci],
              nextKp[ci+1]]
        nEl = getRegNumEls(kp, elSize)
        regNm = 'oml_sp' + str(si) + '_ch' + str(ci)
        wingSurf.addShellRegion('quad2',kp,nEl,name=regNm,meshMethod='structured')
    
    kp = np.zeros((9,3),dtype=float)
    kp[0] = thisKp[8]
    kp[1] = next2Kp[8]
    kp[2] = next2Kp[24]
    kp[3] = thisKp[24]
    kp[4] = nextKp[8]
    kp[6] = nextKp[24]
    kp[5] = 0.5*(kp[1] + kp[2])
    kp[7] = 0.5*(kp[0] + kp[3])
    kp[8] = 0.5*(kp[4] + kp[6])
    nEl = getRegNumEls(kp, elSize)
    regNm = 'tespar_' + str(si)
    wingSurf.addShellRegion('quad2',kp,nEl,name=regNm,meshMethod='structured')
    
    kp = np.zeros((9,3),dtype=float)
    kp[0] = thisKp[12]
    kp[1] = next2Kp[12]
    kp[2] = next2Kp[20]
    kp[3] = thisKp[20]
    kp[4] = nextKp[12]
    kp[6] = nextKp[20]
    kp[5] = 0.5*(kp[1] + kp[2])
    kp[7] = 0.5*(kp[0] + kp[3])
    kp[8] = 0.5*(kp[4] + kp[6])
    nEl = getRegNumEls(kp, elSize)
    regNm = 'lespar_' + str(si)
    wingSurf.addShellRegion('quad2',kp,nEl,name=regNm,meshMethod='structured')
    
    kp = np.zeros((9,3),dtype=float)
    kp[0] = thisKp[8]
    kp[1] = thisKp[10]
    kp[2] = thisKp[22]
    kp[3] = thisKp[24]
    kp[4] = thisKp[9]
    kp[6] = thisKp[23]
    kp[5] = 0.5*(kp[1] + kp[2])
    kp[7] = 0.5*(kp[0] + kp[3])
    kp[8] = 0.5*(kp[4] + kp[6])
    nEl = getRegNumEls(kp, elSize)
    regNm = 'terib_' + str(si)
    wingSurf.addShellRegion('quad2',kp,nEl,name=regNm,meshMethod='structured')
    
    kp = np.zeros((9,3),dtype=float)
    kp[0] = thisKp[10]
    kp[1] = thisKp[12]
    kp[2] = thisKp[20]
    kp[3] = thisKp[22]
    kp[4] = thisKp[11]
    kp[6] = thisKp[21]
    kp[5] = 0.5*(kp[1] + kp[2])
    kp[7] = 0.5*(kp[0] + kp[3])
    kp[8] = 0.5*(kp[4] + kp[6])
    nEl = getRegNumEls(kp, elSize)
    regNm = 'lerib_' + str(si)
    wingSurf.addShellRegion('quad2',kp,nEl,name=regNm,meshMethod='structured')
    
    kp = np.zeros((9,3),dtype=float)
    kp[0] = thisKp[0]
    kp[1] = next2Kp[0]
    kp[2] = next2Kp[32]
    kp[3] = thisKp[32]
    kp[4] = nextKp[0]
    kp[6] = nextKp[32]
    kp[5] = 0.5*(kp[1] + kp[2])
    kp[7] = 0.5*(kp[0] + kp[3])
    kp[8] = 0.5*(kp[4] + kp[6])
    nEl = getRegNumEls(kp, elSize)
    regNm = 'tecon_' + str(si)
    wingSurf.addShellRegion('quad2',kp,nEl,name=regNm,meshMethod='structured')
    
    # fullMesh = wingSurf.getSurfaceMesh()
    # plotShellMesh(fullMesh)
    # inp = input('continue?')
    # if(inp == 'no'):
    #     break

fullMesh = wingSurf.getSurfaceMesh()
#plotShellMesh(fullMesh)
secSets = fullMesh['sets']['element'].copy()

pressureSets = list()
for es in fullMesh['sets']['element']:
    if('oml_' in es['name']):
        stLst = es['name'].split('_ch')
        if(int(stLst[1]) < 0.5*lenKp):
            pressureSets.append(es['name'])
            
fullMesh = getElementSetUnion(fullMesh,pressureSets,'pressureSideEls')

zMrg = 0.5*elSize
fullMesh = getNodeSetInXYZRange(fullMesh,'rootNodes',zRange=[-zMrg,zMrg])

mod = Model()
mod.addMeshData(fullMesh,'shell')

for es in secSets:
    newSec = Section('shell')
    newSec.setElementSet(es['name'])
    newSec.setOrientation([0.0,0.0,1.0],[1.0,0.0,0.0])
    newSec.addLayer('aluminum',0.01)
    mod.addSection(newSec)
    
mat = Material('aluminum')
mat.setDensity(2800.0)
mat.setIsotropic(70.0e+9,0.3)

mod.addMaterial(mat)

rtDir = getEnvPath('rootpath')

fn = rtDir + '/examples/common/crmWing.yaml'
mod.writeModelInput(fn)

rp = ResultsProcessor(fn)
rp.plotElementProperty()