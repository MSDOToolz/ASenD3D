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
    

def plotMeshSolution(nodeCrd,values,faceVerts,valMode='vertex',title=''):
    ## valMode = 'vertex' or 'cell'
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
    
    xMax = np.max(nodeCrd['xLst'])
    xMin = np.min(nodeCrd['xLst'])
    xLen = xMax - xMin
    xMid = 0.5*(xMax+xMin)
    yMax = np.max(nodeCrd['yLst'])
    yMin = np.min(nodeCrd['yLst'])
    yLen = yMax - yMin
    yMid = 0.5*(yMax+yMin)
    zMax = np.max(nodeCrd['zLst'])
    zMin = np.min(nodeCrd['zLst'])
    zLen = zMax - zMin
    zMid = 0.5*(zMax+zMin)
    
    maxLen = np.max([xLen,yLen,zLen])
    hL = 0.75*maxLen
    scn = {'xaxis': {'range': [(xMid-hL), (xMid+hL)], 'showbackground': False},
           'yaxis': {'range': [(yMid-hL), (yMid+hL)], 'showbackground': False},
           'zaxis': {'range': [(zMid-hL), (zMid+hL)], 'showbackground': False}}
    fig.update_layout(scene=scn)

    fig.show()
    return

def animateMeshSolution(nodeCrd,values,faceVerts,valMode='vertex',frameDuration=1000,title=''):
    frameList = list()
    for ni, nC in enumerate(nodeCrd):
        if(ni > 0):
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
                # buttons=[dict(label='Play',
                #               method='animate',
                #               args=[None,{'frame': {'duration': frameDuration, 'redraw': False},
                #                           'fromcurrent': True, 'transition': {'duration': 0}}])])]
            ),
        frames=frameList
        )
    
    xMax = np.max(nodeCrd[0]['xLst'])
    xMin = np.min(nodeCrd[0]['xLst'])
    xLen = xMax - xMin
    xMid = 0.5*(xMax+xMin)
    yMax = np.max(nodeCrd[0]['yLst'])
    yMin = np.min(nodeCrd[0]['yLst'])
    yLen = yMax - yMin
    yMid = 0.5*(yMax+yMin)
    zMax = np.max(nodeCrd[0]['zLst'])
    zMin = np.min(nodeCrd[0]['zLst'])
    zLen = zMax - zMin
    zMid = 0.5*(zMax+zMin)
    
    maxLen = np.max([xLen,yLen,zLen])
    hL = 0.75*maxLen
    scn = {'xaxis': {'range': [(xMid-hL), (xMid+hL)], 'showbackground': False},
            'yaxis': {'range': [(yMid-hL), (yMid+hL)], 'showbackground': False},
            'zaxis': {'range': [(zMid-hL), (zMid+hL)], 'showbackground': False}}
    
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