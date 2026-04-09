# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 16:41:32 2026

@author: evans
"""

import numpy as np

from asendUtils.meshing.MeshTools import *

def splitToTri(meshData):
    newEList = list()
    newLabels = list()
    newEi = 0
    for el in meshData['elements']:
        if(el[3] == -1):
            newEList.append(el)
            newLabels.append([newEi])
            newEi += 1
        else:
            newE = np.array([el[0],el[1],el[2],-1])
            newEList.append(newE)
            newE = np.array([el[0],el[2],el[3],-1])
            newEList.append(newE)
            newLabels.append([newEi, newEi + 1])
            newEi += 2
    meshData['elements'] = np.array(newEList)
    try:
        eSets = meshData['sets']['element']
        newESets = dict()
        for es in eSets:
            newLst = list()
            for el in eSets[es]:
                newLst.extend(newLabels[el])
            newESets[es] = newLst
        meshData['sets']['element'] = newESets
    except:
        pass
    return meshData

def convertToQuadratic(meshData):
    nodes = meshData['nodes']
    elements = meshData['elements']
    numEls = len(elements)
    elNds = len(elements[0])
    numNds = len(nodes)
    newNds = list(nodes)
    if(elNds == 8):
        newEls = -1*np.ones((numEls,10),dtype=int)
        for i, el in enumerate(elements):
            if(el[4] == -1):
                newEls[i,0:4] = el[0:4]
                newNds.append(0.5*(nodes[el[0]] + nodes[el[1]]))
                newNds.append(0.5*(nodes[el[1]] + nodes[el[2]]))
                newNds.append(0.5*(nodes[el[2]] + nodes[el[0]]))
                newNds.append(0.5*(nodes[el[0]] + nodes[el[3]]))
                newNds.append(0.5*(nodes[el[1]] + nodes[el[3]]))
                newNds.append(0.5*(nodes[el[2]] + nodes[el[3]]))
                newEls[i,4:10] = np.array(range(numNds,numNds+6))
                numNds = numNds + 6
    elif(elNds == 4):
        newEls = -1*np.ones((numEls,8),dtype=int)
        for i, el in enumerate(elements):
            if(el[3] == -1):
                newEls[i,0:3] = el[0:3]
                newNds.append(0.5*(nodes[el[0]] + nodes[el[1]]))
                newNds.append(0.5*(nodes[el[1]] + nodes[el[2]]))
                newNds.append(0.5*(nodes[el[2]] + nodes[el[0]]))
                newEls[i,3:6] = np.array(range(numNds,numNds+3))
                numNds = numNds + 3
            else:
                newEls[i,0:4] = el[0:4]
                newNds.append(0.5*(nodes[el[0]] + nodes[el[1]]))
                newNds.append(0.5*(nodes[el[1]] + nodes[el[2]]))
                newNds.append(0.5*(nodes[el[2]] + nodes[el[3]]))
                newNds.append(0.5*(nodes[el[3]] + nodes[el[0]]))
                newEls[i,4:8] = np.array(range(numNds,numNds+4))
                numNds = numNds + 4
    newData = dict()
    newData['nodes'] = np.array(newNds)
    newData['elements'] = newEls
    return mergeDuplicateNodes(newData)

def addFreeNodes(meshData,ndList,setName):
    stLen = len(meshData['nodes'])
    newLen = len(ndList)
    totLen = stLen + newLen
    newNds = np.zeros((totLen,3),dtype=float)
    newNds[0:stLen] = meshData['nodes'].copy()
    newNds[stLen:totLen] = ndList.copy()
    newSet = {setName: list(range(stLen,totLen))}
    meshData['nodes'] = newNds
    meshData = addNodeSet(meshData,newSet)
    return meshData

def addMassElements(meshData,nodeSet,elSetName):
    newSet = [nodeSet,elSetName]
    try:
        meshData['massElements'].append(newSet)
    except:
        meshData['massElements'] = [newSet]
    return meshData
        
def addForceElements(meshData,nodeSet1,nodeSet2,elSetName):
    newSet = [nodeSet1,nodeSet2,elSetName]
    try:
        meshData['forceElements'].append(newSet)
    except:
        meshData['forceElements'] = [newSet]
    return meshData

def make3D(meshData):
    numNodes = len(meshData['nodes'])
    nodes3D = np.zeros((numNodes,3))
    nodes3D[:,0:2] = meshData['nodes']
    dataOut = dict()
    dataOut['nodes'] = nodes3D
    dataOut['elements'] = meshData['elements']
    try:
        inSets = meshData['sets']
        dataOut['sets'] = inSets
    except:
        pass
    return dataOut