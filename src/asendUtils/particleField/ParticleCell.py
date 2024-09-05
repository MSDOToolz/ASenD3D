# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 19:56:14 2024

@author: evans
"""

import numpy as np
from asendUtils.particleField.FreeParticle import *

class ParticleCell:
    
    def __init__(self,xMin,xMax,yMin,yMax,zMin,zMax):
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        self.zMin = zMin
        self.zMax = zMax
        self.center = 0.5*np.array([xMax+xMin,yMax+yMin,zMax+zMin])
        self.active = True
        self.presDen = None
        self.presVel = None
        self.presTemp = None
        self.initDen = None
        self.initVel = None
        self.initTemp = None
        self.refLev = 0
        self.particleList = list()
        
    def deactivate(self):
        self.active = False