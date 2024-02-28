# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 13:27:48 2024

@author: evans
"""

import plotly.express as px
import yaml
from asendUtils.meshing.Mesh2D import *
from asendUtils.meshing.Mesh3D import *
from asendUtils.meshing.MeshTools import *
from asendUtils.syst.pathTools import *
from asendUtils.visualization.plotlyUtils import *

# airFoilBnd = Boundary2D()
# keyPts = [[1.0,-0.00147],
#           [0.9405,0.00785],
#           [0.7721,-0.01045],
#           [0.55321,-0.04074],
#           [0.23023,-0.04541],
#           [0.0263,-0.01901],
#           [0.0,0.0],
#           [0.01937,0.02044],
#           [0.17422,0.05118],
#           [0.56482,0.06324],
#           [0.80119,0.04186],
#           [0.9405,0.0194],
#           [1.0,0.00147]]
# airFoilBnd.addSegment('curve',keyPts,100)

# keyPts = [[1.0,0.00147],
#           [1.0,-0.00147]]
# airFoilBnd.addSegment('line', keyPts, 1)

# bndMesh = airFoilBnd.getBoundaryMesh()

# fig = px.scatter(x=bndMesh['nodes'][:,0],y=bndMesh['nodes'][:,1])
# fig.show()   

rtDir = getEnvPath('rootpath')

afPath = rtDir + '/examples/common/crmAfPtsMedium.yaml'
inFile = open(afPath,'r')
afData = yaml.safe_load(inFile)
inFile.close()


xSectnMesh = Mesh2D(afData['nodes'],afData['edges'])
masterMesh = xSectnMesh.createUnstructuredMesh('quad')

masterMesh = make3D(masterMesh)

plotShellMesh(masterMesh)

