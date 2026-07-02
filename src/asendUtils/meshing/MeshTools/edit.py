# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 16:41:32 2026

@author: evans
"""

import numpy as np

from asendUtils.meshing.MeshTools.extract import *
from asendUtils.meshing.MeshTools.merge import *
from asendUtils.meshing.MeshTools.find import *
from asendUtils.meshing.ElementUtils import *

from asendUtils.visualization.plotlyUtils import *

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

def removeNodeSet(meshData, nodeSet):
    nodes = meshData['nodes']
    elements = meshData['elements']
    nset = set(meshData['sets']['node'][nodeSet])
    
    ndNewLab = -1*np.ones(len(nodes), dtype=int)
    i = 0
    for j in range(0, len(nodes)):
        if j not in nset:
            ndNewLab[j] = i
            i += 1
    
    elElim = -1*np.ones(len(elements), dtype=int)
    for i, el in enumerate(elements):
        for nd in el:
            if nd != -1:
                if ndNewLab[nd] == -1:
                    elElim[i] = 1
                    
    elNewLab = -1*np.ones(len(elements), dtype=int)
    i = 0
    for j in range(0, len(elements)):
        if elElim[j] == -1:
            elNewLab[j] = i
            i += 1
            
    newNds = list()
    for i, nd in enumerate(nodes):
        if ndNewLab[i] != -1:
            newNds.append(nd)
            
    newEls = list()
    for i, el in enumerate(elements):
        if elElim[i] == -1:
            newEl = list()
            for nd in el:
                if nd == -1:
                    newEl.append(-1)
                else:
                    newEl.append(ndNewLab[nd])
            newEls.append(newEl)
    
    sets = meshData['sets']
    newSets = {'node': dict(), 'element': dict()}
    
    for ns in sets['node']:
        if ns != nodeSet:
            newLabs = list()
            for nd in sets['node'][ns]:
                if ndNewLab[nd] != -1:
                    newLabs.append(ndNewLab[nd])
            if len(newLabs) > 0:
                newSets['node'][ns] = newLabs
    
    try:
        for es in sets['element']:
            newLabs = list()
            for el in sets['element'][es]:
                if elNewLab[el] != -1:
                    newLabs.append(elNewLab[el])
            if len(newLabs) > 0:
                newSets['element'][es] = newLabs
    except:
        pass
            
    meshData['nodes'] = np.array(newNds)
    meshData['elements'] = np.array(newEls)
    meshData['sets'] = newSets
    
    return meshData
    
def removeElementSet(meshData, elSet, removeNds=True):
    nodes = meshData['nodes']
    elements = meshData['elements']
    eset = set(meshData['sets']['element'][elSet])
    
    elNewLab = -1*np.ones(len(elements), dtype=int)
    i = 0
    for j in range(0, len(elements)):
        if j not in eset:
            elNewLab[j] = i 
            i += 1
    
    if removeNds:
        ndElim = np.ones(len(nodes), dtype=int)
        for i, el in enumerate(elements):
            if elNewLab[i] != -1:
                for nd in el:
                    if nd != -1:
                        ndElim[nd] = -1
        
        ndNewLab = -1*np.ones(len(nodes), dtype=int)
        i = 0
        for j in range(0, len(nodes)):
            if ndElim[j] == -1:
                ndNewLab[j] = i
                i += 1
                
        newNds = list()
        for i, nd in enumerate(nodes):
            if ndElim[i] == -1:
                newNds.append(nd)
        newNds = np.array(newNds)
    else:
        ndElim = -1*np.ones(len(nodes), dtype=int)
        ndNewLab = np.array(range(0, len(nodes)))
        newNds = nodes
        
    newEls = list()
    for i, el in enumerate(elements):
        if elNewLab[i] != -1:
            newEl = list()
            for nd in el:
                if nd == -1:
                    newEl.append(-1)
                else:
                    newEl.append(ndNewLab[nd])
            newEls.append(newEl)
            
    sets = meshData['sets']
    newSets = {'node': dict(), 'element': dict()}
    
    try:
        for ns in sets['node']:
            newLabs = list()
            for nd in sets['node'][ns]:
                if ndNewLab[nd] != -1:
                    newLabs.append(ndNewLab[nd])
            if len(newLabs) > 0:
                newSets['node'][ns] = newLabs
    except:
        pass
    
    for es in sets['element']:
        if es != elSet:
            newLabs = list()
            for el in sets['element'][es]:
                if elNewLab[el] != -1:
                    newLabs.append(elNewLab[el])
            if len(newLabs) > 0:
                newSets['element'][es] = newLabs
            
    meshData['nodes'] = newNds
    meshData['elements'] = np.array(newEls)
    meshData['sets'] = newSets
    
    return meshData

def cutMesh(meshData, surfaceData, criteria='elOverlap', distance=None):
    ## criteria : 'elOverlap', 'nodeDistance'
    meshNds = meshData['nodes']
    meshEls = meshData['elements']
    surfNds = surfaceData['nodes']
    
    if len(meshEls[0]) == 8:
        meshSolid = True
    else:
        meshSolid = False
    
    if criteria == 'elOverlap':
        meshGL = getMeshSpatialList(meshNds, meshEls)
        for i, el in enumerate(meshEls):
            crd = getElCoord(el, meshNds)
            cent = getElCentroid(crd)
            meshGL.addEntry(i, cent)
        
        if distance == None:
            dist = 0.25*meshGL.xGSz
        else:
            dist = distance
        srcRad = dist + 0.5*meshGL.xGSz
        
        elElim = -1*np.ones(len(meshEls), dtype=int)
        for sn in surfNds:
            nearEls = meshGL.findInRadius(sn, srcRad)
            for ei in nearEls:
                crd = getElCoord(meshEls[ei], meshNds)
                nn = len(crd[0])
                if meshSolid:
                    if nn == 8:
                        eltp = 'brick8'
                    elif nn == 6:
                        eltp = 'wedge6'
                    else:
                        eltp = 'tet4'
                else:
                    if nn == 4:
                        eltp = 'shell4'
                    else:
                        eltp = 'shell3'
                pO = getProjDist(crd, eltp, sn)
                if pO['distance'] < dist:
                    elElim[ei] = 1
        
        ndElim = np.ones(len(meshNds), dtype=int)
        for i, el in enumerate(meshEls):
            if elElim[i] == -1:
                for nd in el:
                    if nd != -1:
                        ndElim[nd] = -1
                        
    elif criteria == 'nodeDistance':
        meshGL = getMeshSpatialList(meshNds, meshEls)
        for i, nd in enumerate(meshNds):
            meshGL.addEntry(i, nd)
        
        if distance == None:
            dist = 0.25*meshGL.xGSz
        else:
            dist = distance
        srcRad = dist + 0.5*meshGL.xGSz
        
        ## Eliminate nodes matching distance criteria
        ndElim = -1*np.ones(len(meshNds), dtype=int)
        for sn in surfNds:
            nearNds = meshGL.findInRadius(sn, srcRad)
            for ni in nearNds:
                dvec = meshNds[ni] - sn
                di = np.linalg.norm(dvec)
                if di < dist:
                    ndElim[ni] = 1
        
        ## Eliminate connected elements
        elElim = -1*np.ones(len(meshEls), dtype=int)
        for i, el in enumerate(meshEls):
            for nd in el:
                if nd != -1:
                    if ndElim[nd] == 1:
                        elElim[i] = 1
                        
        ## Eliminate hanging nodes
        ndElim = np.ones(len(meshNds), dtype=int)
        for i, el in enumerate(meshEls):
            if elElim[i] == -1:
                for nd in el:
                    if nd != -1:
                        ndElim[nd] = -1
                        
    
    ## New node list
    ndNewLab = -1*np.ones(len(meshNds), dtype=int)
    newNds = list()
    i = 0
    for j, nd in enumerate(meshNds):
        if ndElim[j] == -1:
            ndNewLab[j] = i
            newNds.append(nd)
            i += 1
            
    ## New element list
    elNewLab = -1*np.ones(len(meshEls), dtype=int)
    newEls = list()
    i = 0
    for j, el in enumerate(meshEls):
        if elElim[j] == -1:
            elNewLab[j] = i
            newEl = list()
            for nd in el:
                if nd != -1:
                    newEl.append(ndNewLab[nd])
                else:
                    newEl.append(-1)
            newEls.append(newEl)
            i += 1
    
    ## New sets list
    newSets = {'node': dict(), 'element': dict()}
    
    try:
        for sn in meshData['sets']['node']:
            newLst = list()
            for nd in meshData['sets']['node'][sn]:
                if ndNewLab[nd] != -1:
                    newLst.append(ndNewLab[nd])
            if len(newLst) > 0:
                newSets['node'][sn] = newLst
    except:
        pass
    
    try:
        for sn in meshData['sets']['element']:
            newLst = list()
            for el in meshData['sets']['element'][sn]:
                if elNewLab[el] != -1:
                    newLst.append(elNewLab[el])
            if len(newLst) > 0:
                newSets['element'][sn] = newLst
    except:
        pass
    
    ## Update and return
    meshData['nodes'] = np.array(newNds)
    meshData['elements'] = np.array(newEls)
    meshData['sets'] = newSets
    
    return meshData

def eliminateHangingNodes(meshData):
    nds = meshData['nodes']
    newEls = meshData['elements']
    
    ndElim = np.ones(len(nds), dtype=int)
    for el in newEls:
        for nd in el:
            if nd != -1:
                ndElim[nd] = -1
    
    ndNewLab = -1*np.ones(len(nds), dtype=int)
    newNds = list()
    j = 0
    for i, nd in enumerate(nds):
        if ndElim[i] == -1:
            ndNewLab[i] = j
            newNds.append(nd)
            j += 1
    newNds = np.array(newNds)
    
    for i, el in enumerate(newEls):
        for j, nd in enumerate(el):
            if nd != -1:
                newEls[i,j] = ndNewLab[nd]
    
    meshData['nodes'] = newNds
    meshData['elements'] = newEls
    
    try:
        newSets = dict()
        for ns in meshData['sets']['node']:
            newSet = list()
            for nd in meshData['sets']['node'][ns]:
                newLab = ndNewLab[nd]
                if newLab != -1:
                    newSet.append(newLab)
            newSets[ns] = newSet
        meshData['sets']['node'] = newSets
    except:
        pass
    
    return meshData   

def smoothQuadCorners(meshData, elSet='all'):
    els = meshData['elements']
    nds = meshData['nodes']
    if elSet == 'all':
        eSet = set(range(0, len(els)))
    else:
        eSet = set(meshData['sets']['element'][elSet])
        
    surfEdges = getSurfaceEdges(meshData)
    surfNodes = set()
    for fk in surfEdges:
        for nd in surfEdges[fk]:
            surfNodes.add(nd)
            
    newEls = list()
    for i, el in enumerate(els):
        if i in eSet and el[3] != -1:
            chstr = ''
            ct = 0
            for nd in el:
                if nd in surfNodes:
                    chstr += '1'
                    ct += 1
                else:
                    chstr += '0'
            if ct == 3:
                if chstr == '1101': # 0
                    newEl = np.array([el[1], el[2], el[3], -1])
                elif chstr == '1110': # 1
                    newEl = np.array([el[0], el[2], el[3], -1])
                elif chstr == '0111': # 2
                    newEl = np.array([el[0], el[1], el[3], -1])
                elif chstr == '1011': # 3
                    newEl = np.array([el[0], el[1], el[2], -1])
                else:
                    newEl = np.copy(el)
            else:
                newEl = np.copy(el)
        else:
            newEl = np.copy(el)
            
        newEls.append(newEl)
        
    meshData['elements'] = np.array(newEls)
    
    return eliminateHangingNodes(meshData)


def smoothHexEdges(meshData, elSet='all'):
    els = meshData['elements']
    nds = meshData['nodes']
    if elSet == 'all':
        eSet = set(range(0, len(els)))
    else:
        eSet = set(meshData['sets']['element'][elSet])
    
    surfFaces = getSurfaceFaces(meshData)
    surfNodes = set()
    for fk in surfFaces:
        for nd in surfFaces[fk]:
            surfNodes.add(nd)
    
    newEls = list()
    for i, el in enumerate(els):
        if i in eSet and el[6] != -1:
            chstr = ''
            ct = 0
            for nd in el:
                if nd in surfNodes:
                    chstr += '1'
                    ct += 1
                else:
                    chstr += '0'
            ## Edges
            if ct == 6:
                if chstr == '11111100': # 0 - 1
                    newEl = np.array([el[2], el[5], el[6], el[3], el[4], el[7], -1, -1])
                elif chstr == '11110110': # 1 - 2
                    newEl = np.array([el[0], el[4], el[5], el[3], el[7], el[6], -1, -1])
                elif chstr == '11110011': # 2 - 3
                    newEl = np.array([el[0], el[7], el[4], el[1], el[6], el[5], -1, -1])
                elif chstr == '11111001': # 3 - 0
                    newEl = np.array([el[1], el[4], el[5], el[2], el[7], el[6], -1, -1])
                elif chstr == '11001111': # 4 - 5
                    newEl = np.array([el[0], el[3], el[7], el[1], el[2], el[6], -1, -1])
                elif chstr == '01101111': # 5 - 6
                    newEl = np.array([el[0], el[4], el[1], el[3], el[7], el[2], -1, -1])
                elif chstr == '00111111': # 6 - 7
                    newEl = np.array([el[0], el[3], el[4], el[1], el[2], el[5], -1, -1])
                elif chstr == '10011111': # 7 - 4
                    newEl = np.array([el[0], el[5], el[1], el[3], el[6], el[2], -1, -1])
                elif chstr == '11011101': # 0 - 4
                    newEl = np.array([el[1], el[2], el[3], el[5], el[6], el[7], -1, -1])
                elif chstr == '11101110': # 1 - 5
                    newEl = np.array([el[0], el[2], el[3], el[4], el[6], el[7], -1, -1])
                elif chstr == '01110111': # 2 - 6
                    newEl = np.array([el[0], el[1], el[3], el[4], el[5], el[7], -1, -1])
                elif chstr == '10111011': # 3 - 7
                    newEl = np.array([el[0], el[1], el[2], el[4], el[5], el[6], -1, -1])
                else:
                    newEl = np.copy(el)
            elif ct == 7:
                if chstr == '01111111':
                    newEl = np.array([el[0], el[1], el[3], el[4], -1, -1, -1, -1])
                elif chstr == '10111111':
                    newEl = np.array([el[1], el[2], el[0], el[5], -1, -1, -1, -1])
                elif chstr == '11011111':
                    newEl = np.array([el[2], el[3], el[1], el[6], -1, -1, -1, -1])
                elif chstr == '11101111':
                    newEl = np.array([el[3], el[0], el[2], el[7], -1, -1, -1, -1])
                elif chstr == '11110111':
                    newEl = np.array([el[4], el[0], el[7], el[5], -1, -1, -1, -1])
                elif chstr == '11111011':
                    newEl = np.array([el[5], el[1], el[4], el[6], -1, -1, -1, -1])
                elif chstr == '11111101':
                    newEl = np.array([el[6], el[2], el[5], el[7], -1, -1, -1, -1])
                elif chstr == '11111110':
                    newEl = np.array([el[7], el[3], el[6], el[4], -1, -1, -1, -1])
            else:
                newEl = np.copy(el)
        else:
            newEl = np.copy(el)
            
        newEls.append(newEl)
    
    meshData['elements'] = np.array(newEls)
    
    return eliminateHangingNodes(meshData)
    
    
def projectMeshToSurface(meshData, surfaceData, searchRad, elSet='all', newSetName='projected'):
    meshFaces = getSurfaceFaces(meshData, elSet)
    
    meshNdSet = set() ## set of nodes on the surface of meshData
    for fk in meshFaces:
        for nd in meshFaces[fk]:
            meshNdSet.add(nd)
    
    surfGL = getMeshSpatialList(surfaceData['nodes'], surfaceData['elements'])  ## spatial grid for elements of surfaceData
    for i, el in enumerate(surfaceData['elements']):
        crd = getElCoord(el, surfaceData['nodes'])
        cent = getElCentroid(crd)
        surfGL.addEntry(i, cent)
    
    totNds = len(meshData['nodes'])
    ndProjFound = -1*np.ones(totNds, dtype=int)
    ndProjPt = np.zeros((totNds,3), dtype=float)
    for nd in meshNdSet:
        nearEls = surfGL.findInRadius(meshData['nodes'][nd], searchRad)
        nearPO = {'distance': 1.0e+100, 'nVec': np.zeros(8,dtype=float)}
        nearPCrd = np.zeros(3, dtype=float)
        for el in nearEls:
            ndLabs = surfaceData['elements'][el]
            elCrd = getElCoord(ndLabs, surfaceData['nodes'])
            if ndLabs[3] == -1:
                etp = 'shell3'
            else:
                etp = 'shell4'
            pO = getProjDist(elCrd, etp, meshData['nodes'][nd])
            if pO['distance'] < nearPO['distance']:
                nearPO = pO
                nearPCrd = np.matmul(elCrd, pO['nVec'])
        if nearPO['distance'] > 0.0 and nearPO['distance'] <= searchRad:
            ndProjFound[nd] = 1
            ndProjPt[nd] = nearPCrd
            
    newNds = list()
    newEls = list()
    for fk in meshFaces:
        flen = len(meshFaces[fk])
        allNdFound = True
        for nd in meshFaces[fk]:
            if ndProjFound[nd] == -1:
                allNdFound = False
        if allNdFound:
            newEl = -1*np.ones(8, dtype=int)
            newEl[0:flen] = meshFaces[fk]
            offSet = totNds + len(newNds)
            newEl[flen:2*flen] = np.array(range(0, flen)) + offSet
            newEls.append(newEl)
            for nd in meshFaces[fk]:
                newNds.append(ndProjPt[nd])
                
    allLen = totNds + len(newNds)
    allNds = np.zeros((allLen,3), dtype=float)
    allNds[0:totNds] = meshData['nodes']
    allNds[totNds:allLen] = np.array(newNds)
    
    totEls = len(meshData['elements'])
    allLen = totEls + len(newEls)
    allEls = np.zeros((allLen,8), dtype=int)
    allEls[0:totEls] = meshData['elements']
    allEls[totEls:allLen] = np.array(newEls)
    
    meshData['nodes'] = allNds
    meshData['elements'] = allEls
    
    badEls = checkAllJacobians(allNds, allEls, maxCond=50)
    if len(badEls) > 0:
        print("Warning: elements with negative or ill-conditioned jacobian in surface projected mesh:")
        print(badEls)
        
        elElim = -1*np.ones(len(allEls), dtype=int)
        for el in badEls:
            elElim[el] = 1
            
        ndElim = np.ones(len(allNds), dtype=int)
        for i, el in enumerate(allEls):
            if elElim[i] == -1:
                for nd in el:
                    if nd != -1:
                        ndElim[nd] = -1
                        
        ndNewLab = -1*np.ones(len(allNds), dtype=int)
        newNds = list()
        i = 0
        for j, nd in enumerate(allNds):
            if ndElim[j] == -1:
                ndNewLab[j] = i
                newNds.append(nd)
                i += 1
                
        newEls = list()
        for j, el in enumerate(allEls):
            if elElim[j] == -1:
                newEl = list()
                for nd in el:
                    if nd == -1:
                        newEl.append(-1)
                    else:
                        newEl.append(ndNewLab[nd])
                newEls.append(newEl)
        
        meshData['nodes'] = np.array(newNds)
        meshData['elements'] = np.array(newEls)
    
    setList = list(range(totNds, len(meshData['nodes'])))
    meshData = addNodeSet(meshData, {newSetName: setList})
    
    setList = list(range(totEls, len(meshData['elements'])))
    meshData = addElementSet(meshData, {newSetName: setList})
    
    return mergeDuplicateNodes(meshData)

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