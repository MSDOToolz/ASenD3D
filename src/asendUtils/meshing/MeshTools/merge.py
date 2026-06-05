# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 06:37:55 2026

@author: evans
"""

import numpy as np

from asendUtils.meshing.MeshTools.extract import *

## - Convert list of mesh objects into a single merged mesh, returning sets representing the elements/nodes from the original meshes     
def mergeDuplicateNodes(meshData,tolerance=None):
    allNds = meshData['nodes']
    allEls = meshData['elements']
    totNds = len(allNds)
    totEls = len(allEls)
    elDim = len(allEls[0])

    # avgSp = getAverageNodeSpacing(meshData['nodes'], meshData['elements'])
    # sp = 2*avgSp
    nodeGL = getMeshSpatialList(allNds,allEls)
    avgSp = 0.5*nodeGL.xGSz
    if(tolerance == None):
        tol = 1.0e-2*avgSp
    else:
        tol = tolerance
    
    i = 0
    for nd in allNds:
        nodeGL.addEntry(i,nd)
        i = i + 1
    
    ndElim = -np.ones(totNds,dtype=int)
    ndNewInd = -np.ones(totNds,dtype=int)
    for n1i in range(0,totNds):
        if(ndElim[n1i] == -1):
            nearNds = nodeGL.findInRadius(allNds[n1i],tol)
            for n2i in nearNds:
                if(n2i > n1i and ndElim[n2i] == -1):
                    proj = allNds[n2i] - allNds[n1i]
                    dist = np.linalg.norm(proj)
                    if(dist < tol):
                        ndElim[n2i] = n1i
    ndi = 0
    nodesFinal = list()
    for n1i in range(0,totNds):
        if(ndElim[n1i] == -1):
            nodesFinal.append(allNds[n1i])
            ndNewInd[n1i] = ndi
            ndi = ndi + 1
    nodesFinal = np.array(nodesFinal)
    for eli in range(0,totEls):
        for j in range(0,elDim):
            nd = allEls[eli,j]
            if(nd != -1):
                if(ndElim[nd] == -1):
                    allEls[eli,j] = ndNewInd[nd]
                else:
                    allEls[eli,j] = ndNewInd[ndElim[nd]]
    
    meshData['nodes'] = nodesFinal
    meshData['elements'] = allEls
    
    try:
        newSets = dict()
        ndSets = meshData['sets']['node']
        for ns in ndSets:
            newLabs = set()
            for nd in ndSets[ns]:
                if(ndElem[nd] == -1):
                    newLabs.add(ndNewInd[nd])
                else:
                    newLabs.add(ndNewInd[ndElim[nd]])
            newSets[ns] = list(newLabs)
        meshData['sets']['node'] = newSets
    except:
        pass
    
    return meshData

def mergeMeshes(mData1,mData2,tolerance=None):
    mergedData = dict()
    nds1 = mData1['nodes']
    nLen1 = len(nds1)
    nds2 = mData2['nodes']
    nLen2 = len(nds2)
    totNds = nLen1 + nLen2
    mrgNds = np.zeros((totNds,3),dtype=float)
    mrgNds[0:nLen1] = nds1
    mrgNds[nLen1:totNds] = nds2
    els1 = mData1['elements']
    eLen1 = len(els1)
    els2 = mData2['elements']
    eLen2 = len(els2)
    totEls = eLen1 + eLen2
    eCols = len(els1[0])
    mrgEls = -1*np.ones((totEls,eCols),dtype=int)
    mrgEls[0:eLen1] = els1
    for i, el in enumerate(els2,start=eLen1):
        addVec = np.zeros(eCols,dtype=int)
        for j, nd in enumerate(el):
            if(nd > -1):
                addVec[j] = nLen1
        mrgEls[i] = el + addVec
    mergedData['nodes'] = mrgNds
    mergedData['elements'] = mrgEls
    mergedData['sets'] = dict()
    mergedData['sets']['node'] = dict()
    mergedData['sets']['element'] = dict()
    try:
        for ns in mData1['sets']['node']:
            mergedData['sets']['node'][ns] = mData1['sets']['node'][ns]
    except:
        pass
    try:
        for es in mData1['sets']['element']:
            mergedData['sets']['element'][es] = mData1['sets']['element'][es]
    except:
        pass
    try:
        ndSets = mData2['sets']['node']
        for ns in ndSets:
            labs = list()
            for nd in ndSets[ns]:
                labs.append(nd+nLen1)
            try:
                mergedData['sets']['node'][ns].extend(labs)
            except:
                mergedData['sets']['node'][ns] = labs
    except:
        pass
    try:
        elSets = mData2['sets']['element']
        for es in elSets:
            labs = list()
            for el in elSets[es]:
                labs.append(el+eLen1)
            try:
                mergedData['sets']['element'][es].extend(labs)
            except:
                mergedData['sets']['element'][es] = labs
    except:
        pass
    return mergeDuplicateNodes(mergedData,tolerance)