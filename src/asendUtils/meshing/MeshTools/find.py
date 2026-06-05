# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 06:36:35 2026

@author: evans
"""

import numpy as np
import copy

from asendUtils.meshing.ElementUtils import *
from asendUtils.meshing.MeshTools.extract import *

def addNodeSet(meshData,newSet):
    try:
        for ns in newSet:
            meshData['sets']['node'][ns] = newSet[ns]
    except:
        try:
            meshData['sets']['node'] = newSet
        except:
            meshData['sets'] = {'node': newSet}
    return meshData

def addElementSet(meshData,newSet):
    try:
        for es in newSet:
            meshData['sets']['element'][es] = newSet[es]
    except:
        try:
            meshData['sets']['element'] = newSet
        except:
            meshData['sets'] = {'element': newSet}
    return meshData

def getNearestNodes(meshData,pt,numNds,setName):
    ptAr = np.array(pt)
    nearLst = list()
    for i, nd in enumerate(meshData['nodes']):
        vec = nd - ptAr
        dist = np.linalg.norm(vec)
        inserted = False
        if(len(nearLst) == 0):
            nearLst.append([i,dist])
            inserted = True
        else:
            for j, nrNd in enumerate(nearLst):
                if((not inserted) and dist < nrNd[1]):
                    nearLst.insert(j,[i,dist])
                    inserted = True
        if((not inserted) and (len(nearLst) < numNds)):
            nearLst.append([i,dist])
        if(len(nearLst) > numNds):
            nearLst.pop(numNds)
    minLab = list()
    for nd in nearLst:
        minLab.append(nd[0])
        
    newSet = {setName: minLab}
    return addNodeSet(meshData,newSet)

def getNodeSetInRadius(meshData,pt,rad,setName):
    labs = list()
    ptAr = np.array(pt)
    for i, nd in enumerate(meshData['nodes']):
        vec = nd - ptAr
        dist = np.linalg.norm(vec)
        if(dist < rad):
            labs.append(i)
    newSet = {setName: labs}
    return addNodeSet(meshData,newSet)

def getNodeSetNearLine(meshData,pt,dirVec,rad,setName):
    mag = np.linalg.norm(dirVec)
    unitDir = (1.0/mag)*np.array(dirVec)
    ptAr = np.array(pt)
    labs = list()
    for i, nd in enumerate(meshData['nodes']):
        ptond = nd - ptAr
        dp = np.dot(ptond,unitDir)
        normVec = ptond - dp*unitDir
        dist = np.linalg.norm(normVec)
        if(dist < rad):
            labs.append(i)
    newSet = {setName: labs}
    return addNodeSet(meshData,newSet)

def getNodeSetNearPlane(meshData,pt,normDir,dist,setName):
    ptAr = np.array(pt)
    nAr = np.array(normDir)
    mag = np.linalg.norm(nAr)
    unitNorm = (1.0/mag)*nAr
    labs = list()
    for i, nd in enumerate(meshData['nodes']):
        ptond = nd - ptAr
        dp = np.dot(ptond,unitNorm)
        if(abs(dp) < dist):
            labs.append(i)
    newSet = {setName: labs}
    return addNodeSet(meshData,newSet)

def getNodeSetInXYZRange(meshData,setName,xRange=None,yRange=None,zRange=None):
    if(xRange == None):
        xRng = [-1.0e+100,1.0e+100]
    else:
        xRng = xRange
    if(yRange == None):
        yRng = [-1.0e+100,1.0e+100]
    else:
        yRng = yRange
    if(zRange == None):
        zRng = [-1.0e+100,1.0e+100]
    else:
        zRng = zRange
    labs = list()
    for i, nd in enumerate(meshData['nodes']):
        if(nd[0] >= xRng[0] and nd[0] <= xRng[1]):
            if(nd[1] >= yRng[0] and nd[1] <= yRng[1]):
                if(nd[2] >= zRng[0] and nd[2] <= zRng[1]):
                    labs.append(i)
    newSet = {setName: labs}
    return addNodeSet(meshData,newSet)

def getConnectedNodeSet(meshData,setName,newSetName):
    
    ndConn = list()
    for i in range(0, len(meshData['nodes'])):
        ndConn.append(set())
    
    for el in meshData['elements']:
        for nd in el:
            if nd != -1:
                for nd2 in el:
                    if nd2 != -1:
                        ndConn[nd].add(nd2)
                        ndConn[nd2].add(nd)
                        
    prevSet = set()
    newSet = set(meshData['sets']['node'][setName])
    while len(prevSet) < len(newSet):
        prevSet = copy.deepcopy(newSet)
        for nd in prevSet:
            newSet = newSet.union(ndConn[nd])
    
    return addNodeSet(meshData, {newSetName: list(newSet)})

def getPeriodicSets(meshData,xDim,yDim,zDim,setNames=None):
    if(setNames is None):
        sN = ['periodicXMin','periodicXMax',
              'periodicYMin','periodicYMax',
              'periodicZMin','periodicZMax',
              'xMinRef','xMaxRef',
              'yMinRef','yMaxRef',
              'zMinRef','zMaxRef']
    else:
        sN = setNames
    nodes = meshData['nodes']
    
    xMin = np.min(nodes[:,0])
    xMax = np.max(nodes[:,0])
    xMid = 0.5*(xMin + xMax)
    yMin = np.min(nodes[:,1])
    yMax = np.max(nodes[:,1])
    yMid = 0.5*(yMin + yMax)
    zMin = np.min(nodes[:,2])
    zMax = np.max(nodes[:,2])
    zMid = 0.5*(zMin + zMax)
    meshData = getNearestNodes(meshData,[xMin,yMid,zMid],1,sN[6])
    meshData = getNearestNodes(meshData,[xMax,yMid,zMid],1,sN[7])
    meshData = getNearestNodes(meshData,[xMid,yMin,zMid],1,sN[8])
    meshData = getNearestNodes(meshData,[xMid,yMax,zMid],1,sN[9])
    meshData = getNearestNodes(meshData,[xMid,yMid,zMin],1,sN[10])
    meshData = getNearestNodes(meshData,[xMid,yMid,zMax],1,sN[11])
    
    # nSp = getAverageNodeSpacing(nodes,meshData['elements'])
    # gSp = 2.0*nSp
    gL = getMeshSpatialList(nodes,meshData['elements'])
    nSp = 0.5*gL.xGSz
    srcTol = 1.0e-4*nSp
    for i, nd in enumerate(nodes):
        gL.addEntry(i,nd)
    xMinSet = list()
    xMaxSet = list()
    yMinSet = list()
    yMaxSet = list()
    zMinSet = list()
    zMaxSet = list()
    xV = np.array([xDim,0.,0.])
    yV = np.array([0.,yDim,0.])
    zV = np.array([0.,0.,zDim])
    for i, nd in enumerate(nodes):
        srchPt = nd + xV
        nearNds = gL.findInRadius(srchPt,nSp)
        for nrNd in nearNds:
            dVec = srchPt - nodes[nrNd]
            dist = np.linalg.norm(dVec)
            if(dist < srcTol):
                xMinSet.append(i)
                xMaxSet.append(nrNd)
        srchPt = nd + yV
        nearNds = gL.findInRadius(srchPt,nSp)
        for nrNd in nearNds:
            dVec = srchPt - nodes[nrNd]
            dist = np.linalg.norm(dVec)
            if(dist < srcTol):
                yMinSet.append(i)
                yMaxSet.append(nrNd)
        srchPt = nd + zV
        nearNds = gL.findInRadius(srchPt,nSp)
        for nrNd in nearNds:
            dVec = srchPt - nodes[nrNd]
            dist = np.linalg.norm(dVec)
            if(dist < srcTol):
                zMinSet.append(i)
                zMaxSet.append(nrNd)
    meshData = addNodeSet(meshData,{sN[0]: xMinSet})
    meshData = addNodeSet(meshData,{sN[1]: xMaxSet})
    meshData = addNodeSet(meshData,{sN[2]: yMinSet})
    meshData = addNodeSet(meshData,{sN[3]: yMaxSet})
    meshData = addNodeSet(meshData,{sN[4]: zMinSet})
    meshData = addNodeSet(meshData,{sN[5]: zMaxSet})
    return meshData

def getNearestElements(meshData,pt,numEls,setName):
    nodes = meshData['nodes']
    ptAr = np.array(pt)
    nearLst = list()
    for i, el in enumerate(meshData['elements']):
        eCrd = getElCoord(el,nodes)
        eCent = getElCentroid(eCrd)
        vec = eCent - ptAr
        dist = np.linalg.norm(vec)
        inserted = False
        if(len(nearLst) == 0):
            nearLst.append([i,dist])
            inserted = True
        else:
            for j, nrEl in enumerate(nearLst):
                if((not inserted) and dist < nrEl[1]):
                    nearLst.insert(j,[i,dist])
                    inserted = True
        if((not inserted) and (len(nearLst) < numEls)):
            nearLst.append([i,dist])
        if(len(nearLst) > numEls):
            nearLst.pop(numEls)
    minLab = list()
    for el in nearLst:
        minLab.append(el[0])
        
    newSet = {setName: minLab}
    return addNodeSet(meshData,newSet)

def getElementSetInRadius(meshData,pt,rad,setName):
    allNds = meshData['nodes']
    ptAr = np.array(pt)
    labs = list()
    for i, eRow in enumerate(meshData['elements']):
        eCrd = getElCoord(eRow,allNds)
        cent = getElCentroid(eCrd)
        distVec = cent - ptAr
        dist = np.linalg.norm(distVec)
        if(dist < rad):
            labs.append(i)
    newSet = {setName: labs}
    return addElementSet(meshData,newSet)

def getElementSetNearLine(meshData,pt,dirVec,rad,setName):
    nodes = meshData['nodes']
    mag = np.linalg.norm(dirVec)
    unitDir = (1.0/mag)*np.array(dirVec)
    ptAr = np.array(pt)
    labs = list()
    for i, el in enumerate(meshData['elements']):
        eCrd = getElCoord(el,nodes)
        cent = getElCentroid(eCrd)
        ptoel = cent - ptAr
        dp = np.dot(ptoel,unitDir)
        normVec = ptoel - dp*unitDir
        dist = np.linalg.norm(normVec)
        if(dist < rad):
            labs.append(i)
    newSet = {setName: labs}
    return addElementSet(meshData,newSet)

def getElementSetNearPlane(meshData,pt,normDir,dist,setName):
    nodes = meshData['nodes']
    mag = np.linalg.norm(normDir)
    unitNorm = (1.0/mag)*normDir
    ptAr = np.array(pt)
    labs = list()
    for i, el in enumerate(meshData['elements']):
        eCrd = getElCoord(el,nodes)
        cent = getElCentroid(eCrd)
        ptoel = cent - pt
        dp = np.dot(ptoel,unitNorm)
        if(abs(dp) < dist):
            labs.append(i)
    newSet = {setName: labs}
    return addElementSet(meshData,newSet)

def getElementSetInXYZRange(meshData,setName,xRange=None,yRange=None,zRange=None):
    if(xRange == None):
        xRng = [-1.0e+100,1.0e+100]
    else:
        xRng = xRange
    if(yRange == None):
        yRng = [-1.0e+100,1.0e+100]
    else:
        yRng = yRange
    if(zRange == None):
        zRng = [-1.0e+100,1.0e+100]
    else:
        zRng = zRange
    nodes = meshData['nodes']
    labs = list()
    for i, el in enumerate(meshData['elements']):
        eCrd = getElCoord(el,nodes)
        cent = getElCentroid(eCrd)
        if(cent[0] >= xRng[0] and cent[0] <= xRng[1]):
            if(cent[1] >= yRng[0] and cent[1] <= yRng[1]):
                if(cent[2] >= zRng[0] and cent[2] <= zRng[1]):
                    labs.append(i)
    newSet = {setName: labs}
    return addElementSet(meshData,newSet)

def getSetInterfaceNodes(meshData,nodeSet1,nodeSet2,newSet1Name,newSet2Name,maxDist):
    nds = meshData['nodes']
    set1Labs = meshData['sets']['node'][nodeSet1]
    set2Labs = meshData['sets']['node'][nodeSet2]
    newLabs1 = set()
    newLabs2 = set()
    for s1 in set1Labs:
        crd1 = nds[s1]
        for s2 in set2Labs:
            crd2 = nds[s2]
            dVec = crd1 - crd2
            dist = np.linalg.norm(dVec)
            if(dist < maxDist):
                newLabs1.add(s1)
                newLabs2.add(s2)
    meshData['sets']['node'][newSet1Name] = list(newLabs1)
    meshData['sets']['node'][newSet2Name] = list(newLabs2)
    return meshData

def getMeshInterfaceNodes(mesh1Data,mesh2Data,nodeSet1,nodeSet2,newSet1Name,newSet2Name,maxDist):
    nds1 = mesh1Data['nodes']
    nds2 = mesh2Data['nodes']
    set1Labs = mesh1Data['sets']['node'][nodeSet1]
    set2Labs = mesh2Data['sets']['node'][nodeSet2]
    newLabs1 = set()
    newLabs2 = set()
    for s1 in set1Labs:
        crd1 = nds1[s1]
        for s2 in set2Labs:
            crd2 = nds2[s2]
            dVec = crd1 - crd2
            dist = np.linalg.norm(dVec)
            if(dist < maxDist):
                newLabs1.add(s1)
                newLabs2.add(s2)
    mesh1Data['sets']['node'][newSet1Name] = list(newLabs1)
    mesh2Data['sets']['node'][newSet2Name] = list(newLabs2)
    return [mesh1Data,mesh2Data]

def getSurfaceNodes(meshData,elSet,newSetName,normDir,normTol=5.0):
    nds = meshData['nodes']
    els = meshData['elements']
    mag = np.linalg.norm(normDir)
    unitNorm = (1.0/mag)*normDir
    cosTol = np.cos(normTol*np.pi/180.0)
    faceDic = getSurfaceFaces(meshData,elSet)
    surfSet = set()
    for fk in faceDic:
        glob = faceDic[fk]
        if(glob is not None):
            gLen = len(glob)
            if(gLen == 3):
                v1 = nds[glob[1]] - nds[glob[0]]
                v2 = nds[glob[2]] - nds[glob[1]]
            elif(gLen == 4):
                v1 = nds[glob[2]] - nds[glob[0]]
                v2 = nds[glob[3]] - nds[glob[1]]
            cp = crossProd(v1,v2)
            mag = np.linalg.norm(cp)
            fcNrm = (1.0/mag)*cp
            dp = np.dot(fcNrm,unitNorm)
            if(dp >= cosTol):
                for nd in glob:
                    surfSet.add(nd)
    newSet = {newSetName: list(surfSet)}
    return addNodeSet(meshData,newSet)

def getForceElementCloud(meshData,nodeSet1,nodeSet2,elSetName,maxDist):
    nds = meshData['nodes']
    gL = getMeshSpatialList(nds, meshData['elements'])
    for s2 in meshData['sets']['node'][nodeSet2]:
        crds = nds[s2]
        gL.addEntry(s2,crds)
    frcEls = list()
    for s1 in meshData['sets']['node'][nodeSet1]:
        crd1 = nds[s1]
        nearNds = gL.findInRadius(crd1,maxDist)
        for nrNd in nearNds:
            if(nrNd != s1):
                crd2 = nds[nrNd]
                dVec = crd2 - crd1
                dist = np.linalg.norm(dVec)
                if(dist < maxDist):
                    frcEls.append([s1,nrNd])
    try:
        stLen = len(meshData['elements'])
        newLen = stLen + len(frcEls)
        newEls = -1*np.ones((newLen,2),dtype=int)
        newEls[0:stLen] = meshData['elements']
        newEls[stLen:newLen] = np.array(frcEls)
        meshData['elements'] = newEls
        newSet = {elSetName: list(range(stLen,newLen))}
        return addElementSet(meshData,newSet)
    except:
        newLen = len(frcEls)
        meshData['elements'] = np.array(frcEls)
        newSet = {elSetName: list(range(0,newLen))}
        return addElementSet(meshData,newSet)

def getMatchingNodeSet(meshData,elSet,ndSetName):
    elements = meshData['elements']
    elSets = meshData['sets']['element']
    ns = set()
    for ei in elSets[elSet]:
        for elnd in elements[ei]:
            if(elnd > -1):
                ns.add(int(elnd))
    newSet = {ndSetName: list(ns)}
    return addNodeSet(meshData,newSet)
    # return meshData

def getMatchingElementSet(meshData,nodeSet,elSetName,optn='allNodes'):
    ## optn = 'allNodes' or 'anyNode'
    es = list()
    nsLabs = set(meshData['sets']['node'][nodeSet])
    for eli, el in enumerate(meshData['elements']):
        hit = 0
        ct = 0
        for elNd in el:
            if(elNd != -1):
                if(elNd in nsLabs):
                    hit += 1
                ct += 1
        if((optn == 'allNodes' and hit == ct) or (optn == 'anyNode' and hit > 0)):
            es.append(eli)
    newSet = {elSetName: es}
    return addElementSet(meshData,newSet)

def getAllMatchingNodeSets(meshData): ## Name changed from getMatchingNodeSets
    elements = meshData['elements']
    elSets = meshData['sets']['element']
    nodeSets = dict()
    for es in elSets:
        ns = set()
        for ei in elSets[es]:
            for elnd in elements[ei]:
                if(elnd > -1):
                    ns.add(int(elnd))
        nodeSets[es] = list(ns)
    try:
        for ns in nodeSets:
            meshData['sets']['node'][ns] = nodeSets[ns]
    except:
        meshData['sets']['node'] = nodeSets
        
    return meshData

def getExtrudedSets(meshData,numLayers):
    numEls = len(meshData['elements'])
    numNds = len(meshData['nodes'])
    extSets = dict()
    try:
        elSets = meshData['sets']['element']
        extES = dict()
        for es in elSets:
            labels = list()
            for lay in range(0,numLayers):
                for ei in elSets[es]:
                    newLab = ei + numEls*lay
                    labels.append(newLab)
            extES[es] = labels
        extSets['element'] = extES
    except:
        pass
    
    try:
        ndSets = meshData['sets']['node']
        extNS = dict()
        for ns in ndSets:
            labels = list()
            for lay in range(0,(numLayers + 1)):
                for ni in ndSets[ns]:
                    newLab = ni + numNds*lay
                    labels.append(newLab)
            extNS[ns] = labels
        extSets['node'] = extNS
    except:
        pass
        
    return extSets

def getNodeSetUnion(meshData,setList,newSetName):
    un = set()
    ndSets = meshData['sets']['node']
    for ns in ndSets:
        if(ns in setList):
            thisSet = set(ndSets[ns])
            un = un.union(thisSet)
    meshData['sets']['node'][newSetName] = list(un)
    return meshData

def getNodeSetIntersection(meshData,setList,newSetName):
    intsct = set(range(0,len(meshData['nodes'])))
    ndSets = meshData['sets']['node']
    for ns in ndSets:
        if(ns in setList):
            thisSet = set(ndSets[ns])
            intsct = intsct.intersection(thisSet)
    meshData['sets']['node'][newSetName] = list(intsct)
    return meshData

def subtractNodeSet(meshData,set1,set2,newSetName):
    labs = list()
    s1 = meshData['sets']['node'][set1]
    s2 = set(meshData['sets']['node'][set2])
    for nd in s1:
        if(nd not in s2):
            labs.append(nd)
    meshData['sets']['node'][newSetName] = labs
    return meshData

def getElementSetUnion(meshData,setList,newSetName):
    un = set()
    elSets = meshData['sets']['element']
    for es in elSets:
        if(es in setList):
            thisSet = set(elSets[es])
            un = un.union(thisSet)
    meshData['sets']['element'][newSetName] = list(un)
    return meshData

def getElementSetIntersection(meshData,setList,newSetName):
    intsct = set(range(0,len(meshData['elements'])))
    elSets = meshData['sets']['element']
    for es in elSets:
        if(es in setList):
            thisSet = set(elSets[es])
            intsct = intsct.intersection(thisSet)
    meshData['sets']['element'][newSetName] = list(intsct)
    return meshData

def subtractElementSet(meshData,set1,set2,newSetName):
    labs = list()
    s1 = meshData['sets']['element'][set1]
    s2 = set(meshData['sets']['element'][set2])
    for nd in s1:
        if(nd not in s2):
            labs.append(nd)
    meshData['sets']['element'][newSetName] = labs
    return meshData