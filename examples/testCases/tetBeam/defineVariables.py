# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 19:08:29 2024

@author: evans
"""

from asendUtils.variables.DesignVariables import *

dVars = DesignVariables()

dVars.addDesignVariable('modulus',component=1,layer=0,elementSet='all')
dVars.addDesignVariable('density',layer=0,elementSet='all')
dVars.addDesignVariable('nodeCoord',component=1,nodeSet='xMax')
dVars.addDesignVariable('nodeCoord',component=3,nodeSet='zMax')

dVars.writeInput('tetBeamDVars.yaml')