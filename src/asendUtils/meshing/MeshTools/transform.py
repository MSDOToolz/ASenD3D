# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 06:37:24 2026

@author: evans
"""

import numpy as np
from asendUtils.meshing.ElementUtils import *

def rotateVector(vec,axis,angle):
    if(angle < 0.0000000001):
        return vec.copy()
    else:
        axAr = np.array(axis)
        mag = np.linalg.norm(axis)
        unitAxis = (1.0/mag)*axAr
        alp1 = np.zeros((3,3),dtype=float)
        alp1[0] = unitAxis
        i1 = 0
        if(abs(unitAxis[1]) < abs(unitAxis[0])):
            i1 = 1
        if(abs(unitAxis[2]) < abs(unitAxis[i1])):
            i1 = 2
        alp1[1,i1] = np.sqrt(1.0 - alp1[0,i1]*alp1[0,i1])
        for i2 in range(0,3):
            if(i2 != i1):
                alp1[1,i2] = -alp1[0,i1]*alp1[0,i2]/alp1[1,i1]
        alp1[2] = crossProd(alp1[0], alp1[1])
        theta = angle*np.pi/180.0
        cs = np.cos(theta)
        sn = np.sin(theta)
        alp2 = np.array([[1.0,0.0,0.0],
                         [0.0,cs,-sn],
                         [0.0,sn,cs]])
        rV = np.matmul(alp1,vec)
        rV = np.matmul(alp2,rV)
        rV = np.matmul(rV,alp1)
        return rV

def translateMesh(meshData,tVec):
    tAr = np.array(tVec)
    nLen = len(meshData['nodes'])
    newNds = np.zeros((nLen,3),dtype=float)
    for i, nd in enumerate(meshData['nodes']):
        newNds[i] = nd + tAr
    meshData['nodes'] = newNds
    return meshData

def rotateMesh(meshData,pt,axis,angle):
    ptAr = np.array(pt)
    nLen = len(meshData['nodes'])
    newNds = np.zeros((nLen,3),dtype=float)
    for i, nd in enumerate(meshData['nodes']):
        tCrd = nd - ptAr
        rCrd = rotateVector(tCrd,axis,angle)
        newNds[i] = ptAr + rCrd
    meshData['nodes'] = newNds
    return meshData