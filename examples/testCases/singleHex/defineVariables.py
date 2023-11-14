# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 18:53:34 2023

@author: evans
"""

from asendUtils.variables.DesignVariables import *

dVars = DesignVariables()

dVars.addDesignVariable('modulus',component=1,elementSet='all')
dVars.addDesignVariable('density',elementSet='all')
dVars.addDesignVariable('nodeCoord',component=1,nodeSet='xMax')

dVars.writeInput('singleHexDVars.yaml')