# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 09:35:47 2023

@author: evans
"""
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

def plotMeshSolution(nodeCrd,values,faceVerts,title=''):
    fig = go.Figure(data=[
        go.Mesh3d(
            x=nodeCrd['xLst'],
            y=nodeCrd['yLst'],
            z=nodeCrd['zLst'],
            colorbar_title = title,
            colorscale=[[0.0, 'blue'],
                        [0.5, 'green'],
                        [1.0, 'red']],
            intensity=values,
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
    hL = 0.5*maxLen
    scn = {'xaxis': {'range': [(xMid-hL), (xMid+hL)]},
           'yaxis': {'range': [(yMid-hL), (yMid+hL)]},
           'zaxis': {'range': [(zMid-hL), (zMid+hL)]}}
    fig.update_layout(scene=scn)

    fig.show()
    return

def plotTimeHistory(seriesData,timePts,field='',title=''):
    fig = go.Figure()
    fig.update_layout(title=title)
    fig.update_xaxes(title='Time')
    fig.update_yaxes(title=field)
    for s in seriesData:
        fig.add_trace(go.Scatter(x=timePts,y=seriesData[s],mode='lines',name=s))
    
    fig.show()
    
    return