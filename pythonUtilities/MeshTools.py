import numpy as np
import plotly.graph_objects as go

## - Convert list of mesh objects into a single merged mesh, returning sets representing the elements/nodes from the original meshes     
def merge2DMeshes(meshList):
    return
    
def make3D(meshData):
    numNodes = len(meshData['nodes'])
    nodes3D = np.zeros((numNodes,3))
    nodes3D[:,0:2] = meshData['nodes']
    dataOut = dict()
    dataOut['nodes'] = nodes3D
    dataOut['faces'] = meshData['elements']
    return dataOut
 
def plot2DMesh(mesh):
    xLst = mesh.nodes[0:mesh.numNodes,0]
    yLst = mesh.nodes[0:mesh.numNodes,1]
    zLst = np.zeros(mesh.numNodes)
    value = list()
    v1 = list()
    v2 = list()
    v3 = list()
    for i in range(0,mesh.numTriEls):
        v1.append(mesh.triElements[i,0])
        v2.append(mesh.triElements[i,1])
        v3.append(mesh.triElements[i,2])
        value.append(np.sin(i))
    for i in range(0,mesh.numQuadEls):
        v1.append(mesh.quadElements[i,0])
        v2.append(mesh.quadElements[i,1])
        v3.append(mesh.quadElements[i,2])
        value.append(np.sin(i))
        v1.append(mesh.quadElements[i,0])
        v2.append(mesh.quadElements[i,2])
        v3.append(mesh.quadElements[i,3])
        value.append(np.sin(i))
    fig = go.Figure(data=[
        go.Mesh3d(
            x=xLst,
            y=yLst,
            z=zLst,
            colorbar_title = '',
            colorscale=[[0.0, 'blue'],
                        [0.5, 'yellow'],
                        [1.0, 'red']],
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