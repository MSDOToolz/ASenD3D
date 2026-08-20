# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 09:35:47 2023

@author: evans
"""
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

def plotNodes(meshData):
    xLst = meshData['nodes'][:,0]
    yLst = meshData['nodes'][:,1]
    zLst = meshData['nodes'][:,2]
    
    fig = go.Figure(data=[go.Scatter3d(x=xLst,y=yLst,z=zLst,mode='markers')])
    
    fig.show()
    
def removeUnusedNodes(ndCoord,verts,ndVals=None):
    v1 = verts['v1']
    v2 = verts['v2']
    v3 = verts['v3']
    
    ndSet = set()
    for i, v in enumerate(v1):
        ndSet.add(v)
        ndSet.add(v2[i])
        ndSet.add(v3[i])
        
    return reduceToNodeSet(ndCoord, verts, ndSet, ndVals=ndVals)

def reduceToNodeSet(ndCoord, verts, ndSet, ndVals=None, fcVals=None):
    xLst = ndCoord['xLst']
    yLst = ndCoord['yLst']
    zLst = ndCoord['zLst']
    v1 = verts['v1']
    v2 = verts['v2']
    v3 = verts['v3']
    newX = list()
    newY = list()
    newZ = list()
    newVals = list()
    ndNewLab = -np.ones(len(xLst), dtype=int)
    j = 0
    for i, x in enumerate(xLst):
        if i in ndSet:
            newX.append(x)
            newY.append(yLst[i])
            newZ.append(zLst[i])
            ndNewLab[i] = j
            j += 1
            if ndVals != None:
                newVals.append(ndVals[i])
            
    newV1 = list()
    newV2 = list()
    newV3 = list()
    for i in range(0, len(v1)):
        if v1[i] in ndSet and v2[i] in ndSet and v3[i] in ndSet:
            newV1.append(ndNewLab[v1[i]])
            newV2.append(ndNewLab[v2[i]])
            newV3.append(ndNewLab[v3[i]])
            if fcVals != None:
                newVals.append(fcVals[i])
            
    return {'xLst': newX, 'yLst': newY, 'zLst': newZ}, {'v1': newV1, 'v2': newV2, 'v3': newV3}, newVals
            

def plotShellMesh(meshData):
    xLst = meshData['nodes'][:,0]
    yLst = meshData['nodes'][:,1]
    try:
        zLst = meshData['nodes'][:,2]
    except:
        zLst = np.zeros(len(xLst),dtype=float)
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
            colorscale='turbo',
            intensity=value,
            intensitymode='cell',
            i=v1,
            j=v2,
            k=v3,
            name='',
            showscale=True
        )
    ])
    
    # xMax = np.max(xLst)
    # xMin = np.min(xLst)
    # xLen = xMax - xMin
    # xMid = 0.5*(xMax+xMin)
    # yMax = np.max(yLst)
    # yMin = np.min(yLst)
    # yLen = yMax - yMin
    # yMid = 0.5*(yMax+yMin)
    # zMax = np.max(zLst)
    # zMin = np.min(zLst)
    # zLen = zMax - zMin
    # zMid = 0.5*(zMax+zMin)
    
    # maxLen = np.max([xLen,yLen,zLen])
    # hL = 0.5*maxLen
    # scn = {'xaxis': {'range': [(xMid-hL), (xMid+hL)]},
    #         'yaxis': {'range': [(yMid-hL), (yMid+hL)]},
    #         'zaxis': {'range': [(zMid-hL), (zMid+hL)]}}
    # fig.update_layout(scene=scn)

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
            colorscale='turbo',
            intensity=value,
            intensitymode='cell',
            i=v1,
            j=v2,
            k=v3,
            name='',
            showscale=True
        )
    ])
    
    xMax = np.max(xLst)
    xMin = np.min(xLst)
    xLen = xMax - xMin
    xMid = 0.5*(xMax+xMin)
    yMax = np.max(yLst)
    yMin = np.min(yLst)
    yLen = yMax - yMin
    yMid = 0.5*(yMax+yMin)
    zMax = np.max(zLst)
    zMin = np.min(zLst)
    zLen = zMax - zMin
    zMid = 0.5*(zMax+zMin)
    
    maxLen = np.max([xLen,yLen,zLen])
    hL = 0.5*maxLen
    scn = {'xaxis': {'range': [(xMid-hL), (xMid+hL)]},
            'yaxis': {'range': [(yMid-hL), (yMid+hL)]},
            'zaxis': {'range': [(zMid-hL), (zMid+hL)]}}
    fig.update_layout(scene=scn)

    fig.show()
    

def plotMeshSolution(nodeCrd,values,faceVerts,valMode='vertex',xRange=None,yRange=None,zRange=None,title=''):
    if xRange == None:
        xMax = np.max(nodeCrd['xLst'])
        xMin = np.min(nodeCrd['xLst'])
    else:
        xMax = xRange[1]
        xMin = xRange[0]
        
    if yRange == None:
        yMax = np.max(nodeCrd['yLst'])
        yMin = np.min(nodeCrd['yLst'])
    else:
        yMax = yRange[1]
        yMin = yRange[0]
    
    if zRange == None:
        zMax = np.max(nodeCrd['zLst'])
        zMin = np.min(nodeCrd['zLst'])
    else:
        zMax = zRange[1]
        zMin = zRange[0]
        
    xLen = xMax - xMin
    xMid = 0.5*(xMax+xMin)
    yLen = yMax - yMin
    yMid = 0.5*(yMax+yMin)
    zLen = zMax - zMin
    zMid = 0.5*(zMax+zMin)
    
    maxLen = np.max([xLen,yLen,zLen])
    hL = 0.75*maxLen
    
    xAug = [xMid - hL, xMid + hL, xMid, xMid, xMid, xMid]
    yAug = [yMid, yMid, yMid - hL, yMid + hL, yMid, yMid]
    zAug = [zMid, zMid, zMid, zMid, zMid - hL, zMid + hL]
    
    nodeCrd['xLst'].extend(xAug)
    nodeCrd['yLst'].extend(yAug)
    nodeCrd['zLst'].extend(zAug)
    
    fig = go.Figure(data=[
        go.Mesh3d(
            x=nodeCrd['xLst'],
            y=nodeCrd['yLst'],
            z=nodeCrd['zLst'],
            colorbar_title = title,
            colorscale='turbo',
            intensity=values,
            intensitymode=valMode,
            i=faceVerts['v1'],
            j=faceVerts['v2'],
            k=faceVerts['v3'],
            name='',
            showscale=True
        )
    ])
    
    scn = {'xaxis': {'showbackground': False},
           'yaxis': {'showbackground': False},
           'zaxis': {'showbackground': False}}
    fig.update_layout(scene=scn)

    fig.show()
    return

def animateMeshSolution(nodeCrd,values,faceVerts,valMode='vertex',xRange=None,yRange=None,zRange=None,title=''):
    if xRange == None:
        xMax = -1.0e-100
        xMin = 1.0e+100
        for nc in nodeCrd:
            xMaxi = np.max(nc['xLst'])
            if xMaxi > xMax:
                xMax = xMaxi
            xMini = np.min(nc['xLst'])
            if xMini < xMin:
                xMin = xMini
    else:
        xMax = xRange[1]
        xMin = xRange[0]
        
    if yRange == None:
        yMax = -1.0e-100
        yMin = 1.0e+100
        for nc in nodeCrd:
            yMaxi = np.max(nc['yLst'])
            if yMaxi > yMax:
                yMax = yMaxi
            yMini = np.min(nc['yLst'])
            if yMini < yMin:
                yMin = yMini
    else:
        yMax = yRange[1]
        yMin = yRange[0]
    
    if zRange == None:
        zMax = -1.0e-100
        zMin = 1.0e+100
        for nc in nodeCrd:
            zMaxi = np.max(nc['zLst'])
            if zMaxi > zMax:
                zMax = zMaxi
            zMini = np.min(nc['zLst'])
            if zMini < zMin:
                zMin = zMini
    else:
        zMax = zRange[1]
        zMin = zRange[0]
        
    xLen = xMax - xMin
    xMid = 0.5*(xMax+xMin)
    yLen = yMax - yMin
    yMid = 0.5*(yMax+yMin)
    zLen = zMax - zMin
    zMid = 0.5*(zMax+zMin)
    
    maxLen = np.max([xLen,yLen,zLen])
    hL = 0.75*maxLen
    
    xAug = [xMid - hL, xMid + hL, xMid, xMid, xMid, xMid]
    yAug = [yMid, yMid, yMid - hL, yMid + hL, yMid, yMid]
    zAug = [zMid, zMid, zMid, zMid, zMid - hL, zMid + hL]
        
    frameList = list()
    for ni, nC in enumerate(nodeCrd):
        if(ni > 0):
            nC['xLst'].extend(xAug)
            nC['yLst'].extend(yAug)
            nC['zLst'].extend(zAug)
            frm = go.Frame(data=[go.Mesh3d(
                x=nC['xLst'],
                y=nC['yLst'],
                z=nC['zLst'],
                colorbar_title = title,
                colorscale='turbo',
                intensity=values[ni],
                intensitymode=valMode,
                i=faceVerts['v1'],
                j=faceVerts['v2'],
                k=faceVerts['v3'],
                name='',
                showscale=True
                )],
                )
            frameList.append(frm)
    
    fig = go.Figure(
        data=[go.Mesh3d(
            x=nodeCrd[0]['xLst'],
            y=nodeCrd[0]['yLst'],
            z=nodeCrd[0]['zLst'],
            colorbar_title = title,
            colorscale='turbo',
            intensity=values[0],
            intensitymode=valMode,
            i=faceVerts['v1'],
            j=faceVerts['v2'],
            k=faceVerts['v3'],
            name='',
            showscale=True
            )],
        layout=go.Layout(
            updatemenus=[dict(
                type='buttons',
                buttons=[dict(label='Play',
                              method='animate',
                              args=[None])])]
            ),
        frames=frameList
        )
    
    scn = {'xaxis': {'showbackground': False},
            'yaxis': {'showbackground': False},
            'zaxis': {'showbackground': False}}
    
    fig.update_layout(scene=scn)
    
    fig.show()
    
    return

def plotTimeHistory(seriesData,timePts,title='',xTitle='Time',yTitle=''):
    fig = go.Figure()
    fig.update_layout(title=title)
    fig.update_xaxes(title=xTitle)
    fig.update_yaxes(title=yTitle)
    for s in seriesData:
        fig.add_trace(go.Scatter(x=timePts,y=seriesData[s],mode='lines',name=s))
    
    fig.show()
    
    return

def plotFrequencySpectrum(seriesData,frequencies,title='',xTitle='Frequency',yTitle='Amplitude'):
    fig = go.Figure()
    fig.update_layout(title=title)
    fig.update_xaxes(title=xTitle)
    fig.update_yaxes(title=yTitle)
    freqLab = list()
    for f in frequencies:
        freqLab.append(str(f))
    for s in seriesData:
        fig.add_trace(go.Bar(x=freqLab,y=seriesData[s],name=s))
    fig.show()