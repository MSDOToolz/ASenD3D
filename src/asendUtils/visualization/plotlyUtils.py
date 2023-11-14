# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 09:35:47 2023

@author: evans
"""
import plotly.graph_objects as go

def plotMeshSolution(xLst,yLst,zLst,values,v1,v2,v3,title=''):
    fig = go.Figure(data=[
        go.Mesh3d(
            x=xLst,
            y=yLst,
            z=zLst,
            colorbar_title = title,
            colorscale=[[0.0, 'blue'],
                        [0.5, 'green'],
                        [1.0, 'red']],
            intensity=values,
            intensitymode='cell',
            i=v1,
            j=v2,
            k=v3,
            name='',
            showscale=True
        )
    ])

    fig.show()
    return
