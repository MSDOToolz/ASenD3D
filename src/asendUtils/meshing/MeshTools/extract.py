# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 06:35:58 2026

@author: evans
"""

import numpy as np
from asendUtils.meshing.ElementUtils import *
from asendUtils.meshing.SpatialGridList2D import *
from asendUtils.meshing.SpatialGridList3D import *

def getDirectionCosines(xDir,xyDir):
    mag = np.linalg.norm(xDir)
    a1 = (1.0/mag)*xDir
    zDir = cross_prod(xDir,xyDir)
    mag = np.linalg.norm(zDir)
    a3 = (1.0/mag)*zDir
    a2 = cross_prod(a3,a1)
    dirCos = np.array([a1,a2,a3])
    return dirCos

def getAverageNodeSpacing(nodes,elements):
    totDist = 0.0
    ct = 0
    for el in elements:
        for ndi in el:
            if(ndi > -1):
                nd1 = nodes[ndi]
                for ndi2 in el:
                    if(ndi2 > -1 and ndi2 != ndi):
                        nd2 = nodes[ndi2]
                        vec = nd1 - nd2
                        dist = np.linalg.norm(vec)
                        totDist = totDist + dist
                        ct = ct + 1
    return totDist/ct

def checkAllJacobians(nodes,elements):
    failedEls = set()
    for ei, el in enumerate(elements):
        elCrd = getElCoord(el, nodes)
        nn = len(elCrd[0])
        if(nn == 8):
            elType = 'brick8'
        elif(nn == 6):
            elType = 'wedge6'
        else:
            elType = ''
        passed = checkJacobian(elCrd,elType)
        if(not passed):
            failedEls.add(ei)
    return failedEls

def getMeshSpatialList(nodes,elements,xSpacing=0,ySpacing=0,zSpacing=0):
    totNds = len(nodes)
    spaceDim = len(nodes[0])

    maxX = np.amax(nodes[:,0])
    minX = np.amin(nodes[:,0])
    maxY = np.amax(nodes[:,1])
    minY = np.amin(nodes[:,1])
    # nto1_2 = np.power(totNds,0.5)
    # nto1_3 = np.power(totNds,0.3333333)
    avgSp = getAverageNodeSpacing(nodes, elements)
    if(spaceDim == 3):
        maxZ = np.amax(nodes[:,2])
        minZ = np.amin(nodes[:,2])
        dimVec = np.array([(maxX-minX),(maxY-minY),(maxZ-minZ)])
        meshDim = np.linalg.norm(dimVec)
        maxX = maxX + 0.01*meshDim
        minX = minX - 0.01*meshDim
        maxY = maxY + 0.01*meshDim
        minY = minY - 0.01*meshDim
        maxZ = maxZ + 0.01*meshDim
        minZ = minZ - 0.01*meshDim
        if(xSpacing == 0):
            xS = 2*avgSp
        else:
            xS = xSpacing
        if(ySpacing == 0):
            yS = 2*avgSp
        else:
            yS = ySpacing
        if(zSpacing == 0):
            zS = 2*avgSp
        else:
            zS = zSpacing
        meshGL = SpatialGridList3D(minX,maxX,minY,maxY,minZ,maxZ,xS,yS,zS)
        #tol = 1.0e-6*meshDim/nto1_3
    else:
        dimVec = np.array([(maxX-minX),(maxY-minY)])
        meshDim = np.linalg.norm(dimVec)
        maxX = maxX + 0.01*meshDim
        minX = minX - 0.01*meshDim
        maxY = maxY + 0.01*meshDim
        minY = minY - 0.01*meshDim
        if(xSpacing == 0):
            xS = 2*avgSp
        else:
            xS = xSpacing
        if(ySpacing == 0):
            yS = 2*avgSp
        else:
            yS = ySpacing
        meshGL = SpatialGridList2D(minX,maxX,minY,maxY,xS,yS)
        #tol = 1.0e-6*meshDim/nto1_2
    return meshGL

def getSurfaceFaces(meshData, elSet='all'):
    els = meshData['elements']
    if elSet == 'all':
        eset = list(range(0, len(els)))
    else:
        eset = meshData['sets']['element'][elSet]
    
    faceDic = dict()
    for ei in eset:
        fcStr, globFc = getSortedFaceStrings(els[ei])
        for fi, fk in enumerate(fcStr):
            if fk in faceDic:
                faceDic[fk] = None
            else:
                faceDic[fk] = globFc[fi]
    
    fcOut = dict()
    for fk in faceDic:
        fdat = faceDic[fk]
        if fdat != None:
            fcOut[fk] = fdat
            
    return fcOut

def getSurfaceMesh(meshData, elSet='all'):
    faces = getSurfaceFaces(meshData, elSet)
    surfNodes = set()
    for fk in faces:
        for nd in faces[fk]:
            surfNodes.add(nd)
    
    ndNewLab = -1*np.ones(len(meshData['nodes']), dtype=int)
    newNds = list()
    i = 0
    for j, nd in meshData['nodes']:
        if j in surfNodes:
            ndNewLab[j] = i
            newNds.append(nd)
            i += 1
    
    newEls = list()
    for fk in faces:
        newEl = np.array([-1,-1,-1,-1])
        for i, nd in enumerate(faces[fk]):
            newEl[i] = ndNewLab[nd]
        newEls.append(newEl)
        
    return {'nodes': np.array(newNds), 'elements': np.array(newEls)}

def tie2MeshesConstraints(tiedMesh,tgtMesh,maxDist):
    tiedNds = tiedMesh['nodes']
    tgtNds = tgtMesh['nodes']
    tgtEls = tgtMesh['elements']
    elGL = getMeshSpatialList(tgtNds,tgtEls)
    radius = elGL.xGSz
    if(radius < maxDist):
        radius = maxDist
    ei = 0
    for el in tgtEls:
        fstNd = tgtNds[el[0]]
        elGL.addEntry(ei,fstNd)
        ei = ei + 1
    ni = 0
    solidStr = 'tet4 wedge6 brick8'
    constraints = list()
    for nd in tiedNds:
        nearEls = elGL.findInRadius(nd,radius)
        minDist = 1.0e+100
        minPO = dict()
        minEi = -1
        for ei in nearEls:
            if(len(tgtEls[ei]) <= 4):
                if(tgtEls[ei,3] == -1):
                    elType = 'shell3'
                else:
                    elType = 'shell4'
            elif(len(tgtEls[ei]) <= 8):
                if(tgtEls[ei,4] == -1):
                    elType = 'tet4'
                elif(tgtEls[ei,6] == -1):
                    elType = 'wedge6'
                else:
                    elType = 'brick8'
            else:
                pstr = 'Warning: encountered unsupported element type in tie2MeshesConstraints'
            xC = []
            yC = []
            zC = []
            for en in tgtEls[ei]:
                if(en > -1):
                    xC.append(tgtNds[en,0])
                    yC.append(tgtNds[en,1])
                    zC.append(tgtNds[en,2])
            elCrd = np.array([xC,yC,zC])
            pO = getProjDist(elCrd,elType,nd)
            if(elType in solidStr):
                if(pO['distance'] > 0.0):
                    solidPO = getSolidSurfProj(elCrd,elType,nd)
                    if(solidPO['distance'] < minDist):
                        minDist = solidPO['distance']
                        minPO = solidPO
                        minEi = ei
                else:
                    minDist = 0.0
                    minPO = pO
                    minEi = ei
            else:
                if(pO['distance'] < minDist):
                    minDist = pO['distance']
                    minPO = pO
                    minEi = ei
        if(minDist < maxDist):
            newConst = dict()
            terms = list()
            newTerm = dict()
            newTerm['nodeSet'] = 'tiedMesh'
            newTerm['node'] = ni
            newTerm['coef'] = -1.0
            terms.append(newTerm)
            nVec = minPO['nVec']
            nVi = 0
            for en in tgtEls[minEi]:
                if(en > -1):
                    newTerm = dict()
                    newTerm['nodeSet'] = 'targetMesh'
                    newTerm['node'] = en
                    newTerm['coef'] = nVec[nVi]
                    terms.append(newTerm)
                    nVi = nVi + 1
            newConst['terms'] = terms
            newConst['rhs'] = 0.0
            constraints.append(newConst)
        ni = ni + 1
    return constraints

def tie2SetsConstraints(mesh,tiedSetName,tgtSetName,maxDist):
    try:
        elements = mesh['elements']
        nodes = mesh['nodes']
        elSets = mesh['sets']['element']
        ndSets = mesh['sets']['node']
        tgtSet = elSets[tgtSetName]
        tiedSet = ndSets[tiedSetName]
        tgtNdSet = ndSet[tgtSetName]
        fstEl = tgtSet[0]
        fstNd = tgtSet[0]
        fstTgtNd = tgtNdSet[0]
    except:
        raise Exception('There was a problem accessing the mesh data in tie2SetsConstraints().  Check the set names and make sure nodes, elements and sets exist in the input mesh')

    tgtNdCrd = []
    for ni in tgtNdSet:
        tgtNdCrd.append(nodes[ni])
    tgtNdCrd = np.array(tgtNdCrd)
    
    #radius = getAverageNodeSpacing(nodes, elements)
    elGL = getMeshSpatialList(tgtNdCrd, elements)
    radius = elGL.xGSz
    if(radius < maxDist):
        radius = maxDist
    for ei in tgtSet:
        fstNd = tgtNdCrd[elements[ei,0]]
        elGL.addEntry(ei,fstNd)
        ei = ei + 1    
    
    solidStr = 'tet4 wedge6 brick8'
    constraints = list()
    for ni in tiedSet:
        nd = nodes[ni]
        nearEls = elGL.findInRadius(nd,radius)
        minDist = 1.0e+100
        minPO = dict()
        minEi = -1
        for ei in nearEls:
            if(len(elements[ei]) <= 4):
                if(elements[ei,3] == -1):
                    elType = 'shell3'
                else:
                    elType = 'shell4'
            elif(len(elements[ei]) <= 8):
                if(elements[ei,4] == -1):
                    elType = 'tet4'
                elif(elements[ei,6] == -1):
                    elType = 'wedge6'
                else:
                    elType = 'brick8'
            else:
                pstr = 'Warning: encountered unsupported element type in tie2SetsConstraints'
                print(pstr)
            xC = []
            yC = []
            zC = []
            for en in elements[ei]:
                if(en > -1):
                    xC.append(nodes[en,0])
                    yC.append(nodes[en,1])
                    zC.append(nodes[en,2])
            elCrd = np.array([xC,yC,zC])
            pO = getProjDist(elCrd,elType,nd)
            if(elType in solidStr):
                if(pO['distance'] > 0.0):
                    solidPO = getSolidSurfProj(elCrd,elType,nd)
                    if(solidPO['distance'] < minDist):
                        minDist = solidPO['distance']
                        minPO = solidPO
                        minEi = ei
                else:
                    minDist = 0.0
                    minPO = pO
                    minEi = ei
            else:
                if(pO['distance'] < minDist):
                    minDist = pO['distance']
                    minPO = pO
                    minEi = ei
        if(minDist < maxDist and (ni not in elements[minEi])):
            newConst = dict()
            terms = list()
            newTerm = dict()
            newTerm['nodeSet'] = tiedSetName
            newTerm['node'] = ni
            newTerm['coef'] = -1.0
            terms.append(newTerm)
            nVec = minPO['nVec']
            nVi = 0
            for en in elements[minEi]:
                if(en > -1):
                    newTerm = dict()
                    newTerm['nodeSet'] = tgtSetName
                    newTerm['node'] = en
                    newTerm['coef'] = nVec[nVi]
                    terms.append(newTerm)
                    nVi = nVi + 1
            newConst['terms'] = terms
            newConst['rhs'] = 0.0
            constraints.append(newConst)
    
    return constraints

def getNodeFieldFunction(meshData,funType,params,nodeSet=None):
    nodes = meshData['nodes']
    if(nodeSet == None):
        ns = list(range(0,len(nodes)))
    else:
        ns = meshData['sets']['node'][nodeSet]
    if(funType == 'radialShift'):
        pt = np.array(params['pt'])
        maxR = params['radius']
    elif(funType == 'boxShift'):
        xRange = params['xRange']
        yRange = params['yRange']
        zRange = params['zRange']
        xMid = 0.5*(xRange[1] + xRange[0])
        xhL = 0.5*(xRange[1] - xRange[0])
        yMid = 0.5*(yRange[1] + yRange[0])
        yhL = 0.5*(yRange[1] - yRange[0])
        zMid = 0.5*(zRange[1] + zRange[0])
        zhL = 0.5*(zRange[1] - zRange[0])
    elif(funType == 'planeShift'):
        pt = np.array(params['pt'])
        vec = np.array(params['vec'])
        mag = np.linalg.norm(vec)
        vec = (1.0/mag)*vec
        maxR = params['radius']
    elif(funType == 'pointPolar'):
        pt = np.array(params['pt'])
        vec = np.array(params['vec'])
        mag = np.linalg.norm(vec)
        vec = (1.0/mag)*vec
        maxR = params['radius']
        coef = 3.4938562148434
    elif(funType == 'planePolar'):
        pt = np.array(params['pt'])
        vec = np.array(params['vec'])
        mag = np.linalg.norm(vec)
        vec = (1.0/mag)*vec
        maxR = params['radius']
        coef = 3.4938562148434
    fVals = list()
    for ndi in ns:
        nCrd = nodes[ndi]
        if(funType == 'radialShift'):
            dVec = nCrd - pt
            rad = np.linalg.norm(dVec)
            x = rad/maxR
            if(abs(x) <= 1.0):
                f = (x+1.0)*(x+1.0)*(1.0-x)*(1.0-x)
                fVals.append(f)
            else:
                fVals.append(0.0)
        elif(funType == 'boxShift'):
            x = (nCrd[0] - xMid)/xhL
            y = (nCrd[1] - yMid)/yhL
            z = (nCrd[2] - zMid)/zhL
            if(abs(x) <= 1.0 and abs(y) <= 1.0 and abs(z) <= 1.0):
                fx = (x+1.0)*(x+1.0)*(1.0-x)*(1.0-x)
                fy = (y+1.0)*(y+1.0)*(1.0-y)*(1.0-y)
                fz = (z+1.0)*(z+1.0)*(1.0-z)*(1.0-z)
                f = fx*fy*fz
                fVals.append(f)
            else:
                fVals.append(0.0)
        elif(funType == 'planeShift'):
            dVec = nCrd - pt
            x = np.dot(dVec,vec)/maxR
            if(abs(x) <= 1.0):
                f = (x+1.0)*(x+1.0)*(1.0-x)*(1.0-x)
                fVals.append(f)
            else:
                fVals.append(0.0)
        elif(funType == 'pointPolar'):
            dVec = nCrd - pt
            rad = np.linalg.norm(dVec)
            dp = np.dot(dVec,vec)
            x = rad/maxR
            if(abs(x) <= 1.0):
                f = (dp/rad)*coef*x*(x+1.0)*(x+1.0)*(1.0-x)*(1.0-x)
                fVals.append(f)
            else:
                fVals.append(0.0)
        elif(funType == 'planePolar'):
            dVec = nCrd - pt
            dp = np.dot(dVec,vec)
            x = dp/maxR
            if(abs(x) <= 1.0):
                f = coef*x*(x+1.0)*(x+1.0)*(1.0-x)*(1.0-x)
                fVals.append(f)
            else:
                fVals.append(0.0)
    return fVals
 
def getElementFieldFunction(meshData,funType,params,elementSet=None):
    nodes = meshData['nodes']
    elements = meshData['elements']
    if(elementSet == None):
        es = list(range(0,len(elements)))
    else:
        es = meshData['sets']['element'][elementSet]
    if(funType == 'radialShift'):
        pt = np.array(params['pt'])
        maxR = params['radius']
    elif(funType == 'boxShift'):
        xRange = params['xRange']
        yRange = params['yRange']
        zRange = params['zRange']
        xMid = 0.5*(xRange[1] + xRange[0])
        xhL = 0.5*(xRange[1] - xRange[0])
        yMid = 0.5*(yRange[1] + yRange[0])
        yhL = 0.5*(yRange[1] - yRange[0])
        zMid = 0.5*(zRange[1] + zRange[0])
        zhL = 0.5*(zRange[1] - zRange[0])
    elif(funType == 'planeShift'):
        pt = np.array(params['pt'])
        vec = np.array(params['vec'])
        mag = np.linalg.norm(vec)
        vec = (1.0/mag)*vec
        maxR = params['radius']
    elif(funType == 'pointPolar'):
        pt = np.array(params['pt'])
        vec = np.array(params['vec'])
        mag = np.linalg.norm(vec)
        vec = (1.0/mag)*vec
        maxR = params['radius']
        coef = 3.4938562148434
    elif(funType == 'planePolar'):
        pt = np.array(params['pt'])
        vec = np.array(params['vec'])
        mag = np.linalg.norm(vec)
        vec = (1.0/mag)*vec
        maxR = params['radius']
        coef = 3.4938562148434
    fVals = list()
    for eli in es:
        eCrd = getElCoord(elements[eli],nodes)
        eCent = getElCentroid(eCrd)
        if(funType == 'radialShift'):
            dVec = eCent - pt
            rad = np.linalg.norm(dVec)
            x = rad/maxR
            if(abs(x) <= 1.0):
                f = (x+1.0)*(x+1.0)*(1.0-x)*(1.0-x)
                fVals.append(f)
            else:
                fVals.append(0.0)
        elif(funType == 'boxShift'):
            x = (eCent[0] - xMid)/xhL
            y = (eCent[1] - yMid)/yhL
            z = (eCent[2] - zMid)/zhL
            if(abs(x) <= 1.0 and abs(y) <= 1.0 and abs(z) <= 1.0):
                fx = (x+1.0)*(x+1.0)*(1.0-x)*(1.0-x)
                fy = (y+1.0)*(y+1.0)*(1.0-y)*(1.0-y)
                fz = (z+1.0)*(z+1.0)*(1.0-z)*(1.0-z)
                f = fx*fy*fz
                fVals.append(f)
            else:
                fVals.append(0.0)
        elif(funType == 'planeShift'):
            dVec = eCent - pt
            x = np.dot(dVec,vec)/maxR
            if(abs(x) <= 1.0):
                f = (x+1.0)*(x+1.0)*(1.0-x)*(1.0-x)
                fVals.append(f)
            else:
                fVals.append(0.0)
        elif(funType == 'pointPolar'):
            dVec = eCent - pt
            rad = np.linalg.norm(dVec)
            dp = np.dot(dVec,vec)
            x = rad/maxR
            if(abs(x) <= 1.0):
                f = (dp/rad)*coef*x*(x+1.0)*(x+1.0)*(1.0-x)*(1.0-x)
                fVals.append(f)
            else:
                fVals.append(0.0)
        elif(funType == 'planePolar'):
            dVec = eCent - pt
            dp = np.dot(dVec,vec)
            x = dp/maxR
            if(abs(x) <= 1.0):
                f = coef*x*(x+1.0)*(x+1.0)*(1.0-x)*(1.0-x)
                fVals.append(f)
            else:
                fVals.append(0.0)
    return fVals