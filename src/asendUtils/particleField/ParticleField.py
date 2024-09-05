# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 07:39:49 2024

@author: evans
"""

import numpy as np
from asendUtils.particleField.ParticleCell import *

class ParticleField:
    
    def __init__(self,xMin,xMax,yMin,yMax,zMin,zMax,xSpacing,ySpacing,zSpacing):
        self.xMin = xMin
        self.xMax = xMax
        self.xSpacing = xSpacing
        self.yMin = yMin
        self.yMax = yMax
        self.ySpacing = ySpacing
        self.zMin = zMin
        self.zMax = zMax
        self.zSpacing = zSpacing
        
        self.xRows = int(np.ceil((xMax-xMin)/xSpacing))
        self.yRows = int(np.ceil((yMax-yMin)/ySpacing))
        self.zRows = int(np.ceil((zMax-zMin)/zSpacing))
        
        self.cellList = list()
        
        for i in range(0,totalCells):
            xList = list()
            for j in range(0,self.yRows):
                yList = list()
                for k in range(0,self.zRows):
                    xL = xMin + i*self.xSpacing
                    xH = xL + self.xSpacing
                    yL = yMin + j*self.ySpacing
                    yH = yL + self.ySpacing
                    zL = zMin + k*self.zSpacing
                    zH = zL + self.zSpacing
                    yList.append(ParticleCell(xL,xH,yL,yH,zL,zH))
                xList.append(yList)
            self.cellList.append(xList)
            
        self.cellSets = dict()
        
    def getMeshOverlapSet(self,meshData,setName):
        setList = list()
        for nd in meshData['nodes']:
            xRow = int(np.floor((nd[0] - self.xMin)/self.xSpacing))
            if(xRow >= 0 and xRow < self.xRows):
               yRow = int(np.floor((nd[1] - self.yMin)/self.ySpacing))
               if(yRow >= 0 and yRow < self.yRows):
                   zRow = int(np.floor((nd[2] - self.zMin)/self.zSpacing))
                   if(zRow >= 0 and zRow < self.zRows):
                       setList.append([xRow,yRow,zRow])
        self.cellSets[setName] = setList
        
    def deactivateSet(self,setName):
        for ci in self.cellSets[setName]:
            xI, yI, zI = ci
            self.cellList[xI][yI][zI].deactivate()