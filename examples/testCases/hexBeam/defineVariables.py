# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 11:22:57 2023

@author: evans
"""

from asendUtils.variables.DesignVariables import *

dVars = DesignVariables()

dVars.addDesignVariable('modulus',component=1,layer=0,elementSet='all')
dVars.addDesignVariable('density',layer=0,elementSet='all')
dVars.addDesignVariable('nodeCoord',component=1,nodeSet='xMax')

dVars.writeInput('hexBeamDVars.yaml')