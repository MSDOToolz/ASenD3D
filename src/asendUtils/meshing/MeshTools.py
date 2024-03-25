import numpy as np
import plotly.graph_objects as go
from asendUtils.meshing.SpatialGridList2D import *
from asendUtils.meshing.SpatialGridList3D import *
from asendUtils.meshing.ElementUtils import *

def meshFromScratch(nodes,elements):
    meshData = dict()
    meshData['nodes'] = np.array(nodes)
    meshData['elements'] = np.array(elements)
    return meshData

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
    ei = 0
    for el in elements:
        xC = []
        yC = []
        zC = []
        for nd in el:
            if(nd > -1):
                xC.append(nodes[nd,0])
                yC.append(nodes[nd,1])
                zC.append(nodes[nd,2])
        elCrd = np.array([xC,yC,zC])
        nn = len(xC)
        if(nn == 8):
            elType = 'brick8'
        elif(nn == 6):
            elType = 'wedge6'
        else:
            elType = ''
        passed = checkJacobian(elCrd,elType)
        if(not passed):
            failedEls.add(ei)
        ei = ei + 1
    return failedEls

def getMeshSpatialList(nodes,xSpacing=0,ySpacing=0,zSpacing=0):
    totNds = len(nodes)
    spaceDim = len(nodes[0])

    maxX = np.amax(nodes[:,0])
    minX = np.amin(nodes[:,0])
    maxY = np.amax(nodes[:,1])
    minY = np.amin(nodes[:,1])
    nto1_2 = np.power(totNds,0.5)
    nto1_3 = np.power(totNds,0.3333333)
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
            xS = 0.5*(maxX - minX)/nto1_3
        else:
            xS = xSpacing
        if(ySpacing == 0):
            yS = 0.5*(maxY - minY)/nto1_3
        else:
            yS = ySpacing
        if(zSpacing == 0):
            zS = 0.5*(maxZ - minZ)/nto1_3
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
            xS = 0.5*(maxX - minX)/nto1_2
        else:
            xS = xSpacing
        if(ySpacing == 0):
            yS = 0.5*(maxY - minY)/nto1_2
        else:
            yS = ySpacing
        meshGL = SpatialGridList2D(minX,maxX,minY,maxY,xS,yS)
        #tol = 1.0e-6*meshDim/nto1_2
    return meshGL

## - Convert list of mesh objects into a single merged mesh, returning sets representing the elements/nodes from the original meshes     
def mergeDuplicateNodes(meshData,tolerance=None):
    allNds = meshData['nodes']
    allEls = meshData['elements']
    totNds = len(allNds)
    totEls = len(allEls)
    elDim = len(allEls[0])

    avgSp = getAverageNodeSpacing(meshData['nodes'], meshData['elements'])
    sp = 2*avgSp
    nodeGL = getMeshSpatialList(allNds,sp,sp,sp)
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
        newSets = list()
        for ns in meshData['sets']['node']:
            newSet = dict()
            newSet['name'] = ns['name']
            newLabs = set()
            for nd in ns['labels']:
                if(ndElem[nd] == -1):
                    newLabs.add(ndNewInd[nd])
                else:
                    newLabs.add(ndNewInd[ndElim[nd]])
            newSet['labels'] = list(newLabs)
            newSets.append(newSet)
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
    mergedData['sets']['node'] = list()
    mergedData['sets']['element'] = list()
    try:
        mergedData['sets']['node'].extend(mData1['sets']['node'])
    except:
        pass
    try:
        mergedData['sets']['element'].extend(mData1['sets']['element'])
    except:
        pass
    try:
        for ns in mData2['sets']['node']:
            newSet = dict()
            newSet['name'] = ns['name']
            labs = list()
            for nd in ns['labels']:
                labs.append(nd+nLen1)
            newSet['labels'] = labs
            mergedData['sets']['node'].append(newSet)
    except:
        pass
    try:
        for es in mData2['sets']['element']:
            newSet = dict()
            newSet['name'] = es['name']
            labs = list()
            for el in es['labels']:
                labs.append(el+eLen1)
            newSet['labels'] = labs
            mergedData['sets']['element'].append(newSet)
    except:
        pass
    return mergeDuplicateNodes(mergedData,tolerance)
    
def splitToTri(meshData):
    newEList = list()
    for el in meshData['elements']:
        if(el[3] == -1):
            newEList.append(el)
        else:
            newE = np.array([el[0],el[1],el[2],-1])
            newEList.append(newE)
            newE = np.array([el[0],el[2],el[3],-1])
            newEList.append(newE)
    meshData['elements'] = np.array(newEList)
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

def addNodeSet(meshData,newSet):
    try:
        meshData['sets']['node'].append(newSet)
    except:
        nSets = list()
        nSets.append(newSet)
        try:
            meshData['sets']['node'] = nSets
        except:
            sets = dict()
            sets['node'] = nSets
            meshData['sets'] = sets
    return meshData

def addElementSet(meshData,newSet):
    try:
        meshData['sets']['element'].append(newSet)
    except:
        elSets = list()
        elSets.append(newSet)
        try:
            meshData['sets']['element'] = elSets
        except:
            sets = dict()
            sets['element'] = elSets
            meshData['sets'] = sets
    return meshData

def addFreeNodes(meshData,ndList,setName):
    stLen = len(meshData['nodes'])
    newLen = len(ndList)
    totLen = stLen + newLen
    newNds = np.zeros((totLen,3),dtype=float)
    newNds[0:stLen] = meshData['nodes'].copy()
    newNds[stLen:totLen] = ndList.copy()
    newSet = dict()
    newSet['name'] = setName
    newSet['labels'] = list(range(stLen,totLen))
    meshData['nodes'] = newNds
    try:
        meshData['sets']['node'].append(newSet)
    except:
        try:
            meshData['sets']['node'] = [newSet]
        except:
            sets = dict()
            sets['node'] = [newSet]
            meshData['sets'] = sets
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

def getNearestNodes(meshData,pt,numNds,setName):
    newSet = dict()
    newSet['name'] = setName
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
        
    newSet['labels'] = minLab
    return addNodeSet(meshData,newSet)

def getNodeSetInRadius(meshData,pt,rad,setName):
    newSet = dict()
    newSet['name'] = setName
    labs = list()
    ptAr = np.array(pt)
    for i, nd in enumerate(meshData['nodes']):
        vec = nd - ptAr
        dist = np.linalg.norm(vec)
        if(dist < rad):
            labs.append(i)
    newSet['labels'] = labs
    return addNodeSet(meshData,newSet)

def getNodeSetNearLine(meshData,pt,dirVec,rad,setName):
    mag = np.linalg.norm(dirVec)
    unitDir = (1.0/mag)*dirVec
    ptAr = np.array(pt)
    newSet = dict()
    newSet['name'] = setName
    labs = list()
    for i, nd in enumerate(meshData['nodes']):
        ptond = nd - ptAr
        dp = np.dot(ptond,unitDir)
        normVec = ptond - dp*unitDir
        dist = np.linalg.norm(normVec)
        if(dist < rad):
            labs.append(i)
    newSet['labels'] = labs
    return addNodeSet(meshData,newSet)

def getNodeSetNearPlane(meshData,pt,normDir,dist,setName):
    ptAr = np.array(pt)
    nAr = np.array(normDir)
    mag = np.linalg.norm(nAr)
    unitNorm = (1.0/mag)*nAr
    newSet = dict()
    newSet['name'] = setName
    labs = list()
    for i, nd in enumerate(meshData['nodes']):
        ptond = nd - ptAr
        dp = np.dot(ptond,unitNorm)
        if(abs(dp) < dist):
            labs.append(i)
    newSet['labels'] = labs
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
    newSet = dict()
    newSet['name'] = setName
    labs = list()
    for i, nd in enumerate(meshData['nodes']):
        if(nd[0] >= xRng[0] and nd[0] <= xRng[1]):
            if(nd[1] >= yRng[0] and nd[1] <= yRng[1]):
                if(nd[2] >= zRng[0] and nd[2] <= zRng[1]):
                    labs.append(i)
    newSet['labels'] = labs
    return addNodeSet(meshData,newSet)

def getPeriodicSets(meshData,xDim,yDim,zDim,setNames=None):
    if(setNames is None):
        sN = ['periodicXMin','periodicXMax',
              'periodicYMin','periodicYMax',
              'periodicZMax','periodicZMax',
              'xMinRef','xMaxRef',
              'yMinRef','yMaxRef',
              'zMinRef','zMaxRef']
    else:
        sN = setNames
    nodes = meshData['nodes']
    nSp = getAverageNodeSpacing(nodes,meshData['elements'])
    gSp = 2.0*nSp
    gL = getMeshSpatialList(nodes,gSp,gSp,gSp)
    srcTol = 1.0e-4*nSp
    for i, nd in enumerate(nodes):
        gL.addEntry(i,nd)
    xMinSet = list()
    xMaxSet = list()
    yMinSet = list()
    yMaxSet = list()
    zMinSet = list()
    zMaxSet = list()
    xMinR = None
    xMaxR = None
    yMinR = None
    yMaxR = None
    zMinR = None
    zMaxR = None
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
                if(xMinR is None):
                    xMinR = i
                    xMaxR = nrNd
                else:
                    xMinSet.append(i)
                    xMaxSet.append(nrNd)
        srchPt = nd + yV
        nearNds = gL.findInRadius(srchPt,nSp)
        for nrNd in nearNds:
            dVec = srchPt - nodes[nrNd]
            dist = np.linalg.norm(dVec)
            if(dist < srcTol):
                if(yMinR is None):
                    yMinR = i
                    yMaxR = nrNd
                else:
                    yMinSet.append(i)
                    yMaxSet.append(nrNd)
        srchPt = nd + zV
        nearNds = gL.findInRadius(srchPt,nSp)
        for nrNd in nearNds:
            dVec = srchPt - nodes[nrNd]
            dist = np.linalg.norm(dVec)
            if(dist < srcTol):
                if(zMinR is None):
                    zMinR = i
                    zMaxR = nrNd
                else:
                    zMinSet.append(i)
                    zMaxSet.append(nrNd)
    meshData = addNodeSet(meshData,{'name': sN[0], 'labels': xMinSet})
    meshData = addNodeSet(meshData,{'name': sN[1], 'labels': xMaxSet})
    meshData = addNodeSet(meshData,{'name': sN[2], 'labels': yMinSet})
    meshData = addNodeSet(meshData,{'name': sN[3], 'labels': yMaxSet})
    meshData = addNodeSet(meshData,{'name': sN[4], 'labels': zMinSet})
    meshData = addNodeSet(meshData,{'name': sN[5], 'labels': zMaxSet})
    meshData = addNodeSet(meshData,{'name': sN[6], 'labels': [xMinR]})
    meshData = addNodeSet(meshData,{'name': sN[7], 'labels': [xMaxR]})
    meshData = addNodeSet(meshData,{'name': sN[8], 'labels': [yMinR]})
    meshData = addNodeSet(meshData,{'name': sN[9], 'labels': [yMaxR]})
    meshData = addNodeSet(meshData,{'name': sN[10], 'labels': [zMinR]})
    meshData = addNodeSet(meshData,{'name': sN[11], 'labels': [zMaxR]})
    return meshData

def getNearestElements(meshData,pt,numEls,setName):
    nodes = meshData['nodes']
    newSet = dict()
    newSet['name'] = setName
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
        
    newSet['labels'] = minLab
    return addNodeSet(meshData,newSet)

def getElementSetInRadius(meshData,pt,rad,setName):
    ## optn='allNodes': elements whose nodes are all in range
    ## optn='anyNode': elements with any node in range
    ## optn='centroid': elements with centroid in range
    allNds = meshData['nodes']
    ptAr = np.array(pt)
    newSet = dict()
    newSet['name'] = setName
    labs = list()
    # if(optn == 'allNodes'):
    #     for i, eRow in enumerate(meshData['elements']):
    #         inRng = True
    #         for nd in eRow:
    #             if(nd != -1):
    #                 distVec = allNds[nd] - ptAr
    #                 dist = np.linalg.norm(distVec)
    #                 if(dist > rad):
    #                     inRng = False
    #         if(inRng):
    #             labs.append(i)
    # elif(optn == 'anyNode'):
    #     for i, eRow in enumerate(meshData['elements']):
    #         inRng = False
    #         for nd in eRow:
    #             if(nd != -1):
    #                 distVec = allNds[nd] - ptAr
    #                 dist = np.linalg.norm(distVec)
    #                 if(dist < rad):
    #                     inRng = True
    #         if(inRng):
    #             labs.append(i)
    # else:
    for i, eRow in enumerate(meshData['elements']):
        # cent = np.zeros(3,dtype=float)
        # ct = 0
        # for nd in eRow:
        #     if(nd != -1):
        #         cent = cent + allNds[nd]
        #         ct = ct + 1
        # cent = (1.0/ct)*cent
        eCrd = getElCoord(eRow,allNds)
        cent = getElCentroid(eCrd)
        distVec = cent - ptAr
        dist = np.linalg.norm(distVec)
        if(dist < rad):
            labs.append(i)
    newSet['labels'] = labs
    return addElementSet(meshData,newSet)

def getElementSetNearLine(meshData,pt,dirVec,rad,setName):
    nodes = meshData['nodes']
    mag = np.linalg.norm(dirVec)
    unitDir = (1.0/mag)*dirVec
    ptAr = np.array(pt)
    newSet = dict()
    newSet['name'] = setName
    labs = list()
    for i, el in enumerate(meshData['elements']):
        eCrd = getElCoord(el,nodes)
        cent = getElCent(eCrd)
        ptoel = cent - ptAr
        dp = np.dot(ptoel,unitDir)
        normVec = ptoel - dp*unitDir
        dist = np.linalg.norm(normVec)
        if(dist < rad):
            labs.append(i)
    newSet['labels'] = labs
    return addElementSet(meshData,newSet)

def getElementSetNearPlane(meshData,pt,normDir,dist,setName):
    nodes = meshData['nodes']
    mag = np.linalg.norm(normDir)
    unitNorm = (1.0/mag)*normDir
    ptAr = np.array(pt)
    newSet = dict()
    newSet['name'] = setName
    labs = list()
    for i, el in enumerate(meshData['elements']):
        eCrd = getElCoord(el,nodes)
        cent = getElCentroid(eCrd)
        ptoel = cent - pt
        dp = np.dot(ptoel,unitNorm)
        if(abs(dp) < dist):
            labs.append(i)
    newSet['labels'] = labs
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
    newSet = dict()
    newSet['name'] = setName
    labs = list()
    for i, el in enumerate(meshData['elements']):
        eCrd = getElCoord(el,nodes)
        cent = getElCentroid(eCrd)
        if(cent[0] >= xRng[0] and cent[0] <= xRng[1]):
            if(cent[1] >= yRng[0] and cent[1] <= yRng[1]):
                if(cent[2] >= zRng[0] and cent[2] <= zRng[1]):
                    labs.append(i)
    newSet['labels'] = labs
    return addElementSet(meshData,newSet)

def getSetInterfaceNodes(meshData,nodeSet1,nodeSet2,newSet1Name,newSet2Name,maxDist):
    nds = meshData['nodes']
    set1Labs = []
    set2Labs = []
    newLabs1 = set()
    newLabs2 = set()
    for ns in meshData['sets']['node']:
        if(ns['name'] == nodeSet1):
            set1Labs = ns['labels']
        if(ns['name'] == nodeSet2):
            set2Labs = ns['labels']
    for s1 in set1Labs:
        crd1 = nds[s1]
        for s2 in set2Labs:
            crd2 = nds[s2]
            dVec = crd1 - crd2
            dist = np.linalg.norm(dVec)
            if(dist < maxDist):
                newLabs1.add(s1)
                newLabs2.add(s2)
    newSet = dict()
    newSet['name'] = newSet1Name
    newSet['labels'] = list(newLabs1)
    meshData['sets']['node'].append(newSet)
    newSet = dict()
    newSet['name'] = newSet2Name
    newSet['labels'] = list(newLabs2)
    meshData['sets']['node'].append(newSet)
    return meshData

def getMeshInterfaceNodes(mesh1Data,mesh2Data,nodeSet1,nodeSet2,newSet1Name,newSet2Name,maxDist):
    nds1 = mesh1Data['nodes']
    nds2 = mesh2Data['nodes']
    set1Labs = []
    set2Labs = []
    newLabs1 = set()
    newLabs2 = set()
    for ns in mesh1Data['sets']['node']:
        if(ns['name'] == nodeSet1):
            set1Labs = ns['labels']
    for ns in mesh2Data['sets']['node']:
        if(ns['name'] == nodeSet2):
            set2Labs = ns['labels']
    for s1 in set1Labs:
        crd1 = nds1[s1]
        for s2 in set2Labs:
            crd2 = nds2[s2]
            dVec = crd1 - crd2
            dist = np.linalg.norm(dVec)
            if(dist < maxDist):
                newLabs1.add(s1)
                newLabs2.add(s2)
    newSet = dict()
    newSet['name'] = newSet1Name
    newSet['labels'] = list(newLabs1)
    mesh1Data['sets']['node'].append(newSet)
    newSet = dict()
    newSet['name'] = newSet2Name
    newSet['labels'] = list(newLabs2)
    mesh2Data['sets']['node'].append(newSet)
    return [mesh1Data,mesh2Data]

def getMatchingNodeSet(meshData,elSet,ndSetName):
    elements = meshData['elements']
    elSets = meshData['sets']['element']
    for es in elSets:
        if(es['name'] == elSet):
            ns = set()
            for ei in es['labels']:
                for elnd in elements[ei]:
                    if(elnd > -1):
                        ns.add(elnd)
            newSet = dict()
            newSet['name'] = ndSetName
            newSet['labels'] = list(ns)
            return addNodeSet(meshData,newSet)
    return meshData

def getMatchingElementSet(meshData,nodeSet,elSetName,optn='allNodes'):
    ## optn = 'allNodes' or 'anyNode'
    for ns in meshData['sets']['node']:
        if(ns['name'] == nodeSet):
            es = list()
            nsLabs = set(ns['labels'])
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
            newSet = dict()
            newSet['name'] = elSetName
            newSet['labels'] = es
            return addElementSet(meshData,newSet)
    return meshData

def getAllMatchingNodeSets(meshData): ## Name changed from getMatchingNodeSets
    elements = meshData['elements']
    elSets = meshData['sets']['element']
    nodeSets = list()
    for es in elSets:
        ns = set()
        for ei in es['labels']:
            for elnd in elements[ei]:
                if(elnd > -1):
                    ns.add(elnd)
        newSet = dict()
        newSet['name'] = es['name']
        newSet['labels'] = list(ns)
        nodeSets.append(newSet)
    try:
        meshData['sets']['node'].extend(nodeSets)
    except:
        meshData['sets']['node'] = nodeSets
        
    return meshData

def getExtrudedSets(meshData,numLayers):
    numEls = len(meshData['elements'])
    numNds = len(meshData['nodes'])
    extSets = dict()
    try:
        elSets = meshData['sets']['element']
        extES = list()
        for es in elSets:
            labels = list()
            for lay in range(0,numLayers):
                for ei in es['labels']:
                    newLab = ei + numEls*lay
                    labels.append(newLab)
            newSet = dict()
            newSet['name'] = es['name']
            newSet['labels'] = labels
            extES.append(newSet)
        extSets['element'] = extES
    except:
        pass
    
    try:
        ndSets = meshData['sets']['node']
        extNS = list()
        for ns in ndSets:
            labels = list()
            for lay in range(0,(numLayers + 1)):
                for ni in ns['labels']:
                    newLab = ni + numNds*lay
                    labels.append(newLab)
            newSet = dict()
            newSet['name'] = ns['name']
            newSet['labels'] = labels
            extNS.append(newSet)
        extSets['node'] = extNS
    except:
        pass        
        
    return extSets

def getNodeSetUnion(meshData,setList,newSetName):
    un = set()
    for ns in meshData['sets']['node']:
        if(ns['name'] in setList):
            thisSet = set(ns['labels'])
            un = un.union(thisSet)
    newSet = dict()
    newSet['name'] = newSetName
    newSet['labels'] = list(un)
    meshData['sets']['node'].append(newSet)
    return meshData

def getNodeSetIntersection(meshData,setList,newSetName):
    intsct = set(range(0,len(meshData['nodes'])))
    for ns in meshData['sets']['node']:
        if(ns['name'] in setList):
            thisSet = set(ns['labels'])
            intsct = intsct.intersection(thisSet)
    newSet = dict()
    newSet['name'] = newSetName
    newSet['labels'] = list(intsct)
    meshData['sets']['node'].append(newSet)
    return meshData

def subtractNodeSet(meshData,set1,set2,newSetName):
    labs = list()
    for ns in meshData['sets']['node']:
        if(ns['name'] == set1):
            s1 = ns
        if(ns['name'] == set2):
            s2 = set(ns['labels'])
    for nd in s1['labels']:
        if(nd not in s2):
            labs.append(nd)
    newSet = dict()
    newSet['name'] = newSetName
    newSet['labels'] = labs
    meshData['sets']['node'].append(newSet)
    return meshData

def getElementSetUnion(meshData,setList,newSetName):
    un = set()
    for es in meshData['sets']['element']:
        if(es['name'] in setList):
            thisSet = set(es['labels'])
            un = un.union(thisSet)
    newSet = dict()
    newSet['name'] = newSetName
    newSet['labels'] = list(un)
    meshData['sets']['element'].append(newSet)
    return meshData

def getElementSetIntersection(meshData,setList,newSetName):
    intsct = set(range(0,len(meshData['elements'])))
    for es in meshData['sets']['element']:
        if(es['name'] in setList):
            thisSet = set(es['labels'])
            intsct = intsct.intersection(thisSet)
    newSet = dict()
    newSet['name'] = newSetName
    newSet['labels'] = list(intsct)
    meshData['sets']['element'].append(newSet)
    return meshData

def subtractElementSet(meshData,set1,set2,newSetName):
    labs = list()
    for ns in meshData['sets']['element']:
        if(ns['name'] == set1):
            s1 = ns
        if(ns['name'] == set2):
            s2 = set(ns['labels'])
    for nd in s1['labels']:
        if(nd not in s2):
            labs.append(nd)
    newSet = dict()
    newSet['name'] = newSetName
    newSet['labels'] = labs
    meshData['sets']['element'].append(newSet)
    return meshData
    
def make3D(meshData):
    numNodes = len(meshData['nodes'])
    nodes3D = np.zeros((numNodes,3))
    nodes3D[:,0:2] = meshData['nodes']
    dataOut = dict()
    dataOut['nodes'] = nodes3D
    dataOut['elements'] = meshData['elements']
    return dataOut

def tie2MeshesConstraints(tiedMesh,tgtMesh,maxDist):
    tiedNds = tiedMesh['nodes']
    tgtNds = tgtMesh['nodes']
    tgtEls = tgtMesh['elements']
    radius = getAverageNodeSpacing(tgtNds,tgtEls)
    elGL = getMeshSpatialList(tgtNds,radius,radius,radius)
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
        for es in elSets:
            if(es['name'] == tgtSetName):
                tgtSet = es['labels']
        for ns in ndSets:
            if(ns['name'] == tiedSetName):
                tiedSet = ns['labels']
            if(ns['name'] == tgtSetName):
                tgtNdSet = ns['labels']
        fstEl = tgtSet[0]
        fstNd = tgtSet[0]
        fstTgtNd = tgtNdSet[0]
    except:
        raise Exception('There was a problem accessing the mesh data in tie2SetsConstraints().  Check the set names and make sure nodes, elements and sets exist in the input mesh')

    tgtNdCrd = []
    for ni in tgtNdSet:
        tgtNdCrd.append(nodes[ni])
    tgtNdCrd = np.array(tgtNdCrd)
    
    radius = getAverageNodeSpacing(nodes, elements)
    elGL = getMeshSpatialList(tgtNdCrd,radius,radius,radius)
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
        ns = dict()
        ns['name'] = 'all'
        nLen = len(nodes)
        ns['labels'] = list(range(0,nLen))
    else:
        ns = dict()
        ns['labels'] = list()
        for ndset in meshData['sets']['node']:
            if(ndset['name'] == nodeSet):
                ns = ndset
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
    for ndi in ns['labels']:
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
        es = dict()
        es['name'] = 'all'
        eLen = len(elements)
        es['labels'] = list(range(0,eLen))
    else:
        es = dict()
        es['labels'] = list()
        for elset in meshData['sets']['element']:
            if(elset['name'] == elementSet):
                es = elset
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
    for eli in es['labels']:
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