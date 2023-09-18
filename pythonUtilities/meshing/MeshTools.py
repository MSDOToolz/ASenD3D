import numpy as np
import plotly.graph_objects as go
from SpatialGridList2DClass import *
from SpatialGridList3DClass import *
from ElementUtils import *

def getMeshSpatialList(nodes):
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
        xSpacing = 0.5*(maxX - minX)/nto1_3
        ySpacing = 0.5*(maxY - minY)/nto1_3
        zSpacing = 0.5*(maxZ - minZ)/nto1_3
        meshGL = SpatialGridList3D(minX,maxX,minY,maxY,minZ,maxZ,xSpacing,ySpacing,zSpacing)
        #tol = 1.0e-6*meshDim/nto1_3
    else:
        dimVec = np.array([(maxX-minX),(maxY-minY)])
        meshDim = np.linalg.norm(dimVec)
        maxX = maxX + 0.01*meshDim
        minX = minX - 0.01*meshDim
        maxY = maxY + 0.01*meshDim
        minY = minY - 0.01*meshDim
        xSpacing = 0.5*(maxX - minX)/nto1_3
        ySpacing = 0.5*(maxY - minY)/nto1_3
        meshGL = SpatialGridList2D(minX,maxX,minY,maxY,xSpacing,ySpacing)
        #tol = 1.0e-6*meshDim/nto1_2
    return meshGL

## - Convert list of mesh objects into a single merged mesh, returning sets representing the elements/nodes from the original meshes     
def mergeDuplicateNodes(meshData):
    allNds = meshData['nodes']
    allEls = meshData['elements']
    totNds = len(allNds)
    totEls = len(allEls)
    elDim = len(allEls[0])

    nodeGL = getMeshSpatialList(allNds)
    glDim = nodeGL.getDim()
    mag = np.linalg.norm(glDim)
    nto1_3 = np.power(len(allNds),0.3333333333)
    tol = 1.0e-6*mag/nto1_3
    
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
    
    return meshData

def getMatchingNodeSets(meshData):
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
    elGL = getMeshSpatialList(tgtNds)
    glDim = elGL.getDim()
    mag = np.linalg.norm(glDim)
    radius = mag/np.power(len(tgtNds),0.33333333)
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
    
    elGL = getMeshSpatialList(tgtNdCrd)
    glDim = elGL.getDim()
    mag = np.linalg.norm(glDim)
    radius = mag/np.power(len(tgtNdCrd),0.33333333)
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
        if(minDist < maxDist):
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
 
def plotShellMesh(meshData):
    xLst = meshData['nodes'][:,0]
    yLst = meshData['nodes'][:,1]
    try:
        zLst = meshData['nodes'][:,2]
    except:
        zLst = np.zeros(len(xLst))
    value = list()
    v1 = list()
    v2 = list()
    v3 = list()
    i = 0
    for el in meshData['elements']:
        v1.append(el[0])
        v2.append(el[1])
        v3.append(el[2])
        value.append(np.sin(i))
        if(el[3] != -1):
            v1.append(el[0])
            v2.append(el[2])
            v3.append(el[3])
            value.append(np.sin(i))
        i = i + 1
    fig = go.Figure(data=[
        go.Mesh3d(
            x=xLst,
            y=yLst,
            z=zLst,
            colorbar_title = '',
            colorscale=[[0.0, 'white'],
                        [0.5, 'gray'],
                        [1.0, 'black']],
            intensity=value,
            intensitymode='cell',
            i=v1,
            j=v2,
            k=v3,
            name='',
            showscale=True
        )
    ])

    fig.show()
    
def plotSolidMesh(meshData):
    xLst = meshData['nodes'][:,0]
    yLst = meshData['nodes'][:,1]
    zLst = meshData['nodes'][:,2]
    value = list()
    v1 = list()
    v2 = list()
    v3 = list()
    i = 0
    for el in meshData['elements']:
        si = np.sin(i)
        if(el[4] == -1):
            v1.append(el[0])
            v2.append(el[1])
            v3.append(el[2])
            value.append(si)
            v1.append(el[0])
            v2.append(el[1])
            v3.append(el[3])
            value.append(si)
            v1.append(el[0])
            v2.append(el[2])
            v3.append(el[3])
            value.append(si)
            v1.append(el[1])
            v2.append(el[2])
            v3.append(el[3])
            value.append(si)
        elif(el[6] == -1):
            v1.append(el[0])
            v2.append(el[1])
            v3.append(el[2])
            value.append(si)
            v1.append(el[3])
            v2.append(el[4])
            v3.append(el[5])
            value.append(si)
            v1.append(el[0])
            v2.append(el[1])
            v3.append(el[3])
            value.append(si)
            v1.append(el[1])
            v2.append(el[3])
            v3.append(el[4])
            value.append(si)
            
            v1.append(el[0])
            v2.append(el[2])
            v3.append(el[3])
            value.append(si)
            v1.append(el[2])
            v2.append(el[3])
            v3.append(el[5])
            value.append(si)
            v1.append(el[1])
            v2.append(el[2])
            v3.append(el[4])
            value.append(si)
            v1.append(el[2])
            v2.append(el[4])
            v3.append(el[5])
            value.append(si)
        else:
            v1.append(el[0])
            v2.append(el[3])
            v3.append(el[4])
            value.append(si)
            v1.append(el[3])
            v2.append(el[4])
            v3.append(el[7])
            value.append(si)
            v1.append(el[1])
            v2.append(el[2])
            v3.append(el[5])
            value.append(si)
            v1.append(el[2])
            v2.append(el[5])
            v3.append(el[6])
            value.append(si)
            
            v1.append(el[0])
            v2.append(el[1])
            v3.append(el[4])
            value.append(si)
            v1.append(el[1])
            v2.append(el[4])
            v3.append(el[5])
            value.append(si)
            v1.append(el[2])
            v2.append(el[3])
            v3.append(el[6])
            value.append(si)
            v1.append(el[3])
            v2.append(el[6])
            v3.append(el[7])
            value.append(si)
            
            v1.append(el[0])
            v2.append(el[1])
            v3.append(el[2])
            value.append(si)
            v1.append(el[0])
            v2.append(el[2])
            v3.append(el[3])
            value.append(si)
            v1.append(el[4])
            v2.append(el[5])
            v3.append(el[6])
            value.append(si)
            v1.append(el[4])
            v2.append(el[6])
            v3.append(el[7])
            value.append(si)
        i = i + 1
    fig = go.Figure(data=[
        go.Mesh3d(
            x=xLst,
            y=yLst,
            z=zLst,
            colorbar_title = '',
            colorscale=[[0.0, 'white'],
                        [0.5, 'gray'],
                        [1.0, 'black']],
            intensity=value,
            intensitymode='cell',
            i=v1,
            j=v2,
            k=v3,
            name='',
            showscale=True
        )
    ])

    fig.show()
## -Create node/element set within a spatial range or radius