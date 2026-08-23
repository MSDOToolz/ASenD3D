# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 10:50:20 2023

@author: evans
"""

import numpy as np

class Section:
    
    def __init__(self,secType):
        self.secData = dict()
        self.secData['type'] = secType
        
    def setMaterial(self,matName):
        self.secData['material'] = matName
        
    def setOrientation(self,xDir,xyVec,listAsStr=True):
        oriLst = list(xDir)
        oriLst.extend(list(xyVec))
        if listAsStr:
            self.secData['orientation'] = str(oriLst)
        else:
            self.secData['orientation'] = oriLst
        
    def setZOffset(self,offset):
        try:
            self.secData['layup']['zOffset'] = offset
        except:
            newLayup = dict()
            newLayup['zOffset'] = offset
            newLayup['layers'] = list()
            self.secData['layup'] = newLayup
            
    def addLayer(self,material,thickness,angle=0.0):
        newLay = dict()
        newLay['material'] = material
        newLay['thickness'] = thickness
        newLay['angle'] = angle
        try:
            self.secData['layup']['layers'].append(newLay)
        except:
            newLayup = dict()
            newLayup['zOffset'] = 0.0
            newLayup['layers'] = [newLay]
            self.secData['layup'] = newLayup
            
    def setXSArea(self,area):
        try:
            self.secData['beamProps']['area'] = area
        except:
            newProps = dict()
            newProps['area'] = area
            self.secData['beamProps'] = newProps
            
    def setAreaMoment(self,I2=0.0,I3=0.0,I22=0.0,I33=0.0,I23=0.0,listAsStr=True):
        if listAsStr:
            ILst = str([I2,I3,I22,I33,I23])
        else:
            ILst = [I2,I3,I22,I33,I23]
        try:
            self.secData['beamProps']['I'] = ILst
        except:
            newProps = dict()
            newProps['I'] = ILst
            self.secData['beamProps'] = newProps
            
    def setPolarMoment(self,J):
        try:
            self.secData['beamProps']['J'] = J
        except:
            newProps = dict()
            newProps['J'] = J
            self.secData['beamProps'] = newProps
            
    def setPotentialField(self,coef,exp):
        potField = dict()
        potField['coef'] = coef
        potField['exp'] = exp
        self.secData['potField'] = potField
        
    def setDampingField(self,coef,exp):
        dampField = dict()
        dampField['coef'] = coef
        dampField['exp'] = exp
        self.secData['dampField'] = dampField
        
    def setMassPerElement(self,elMass):
        self.secData['massPerEl'] = elMass
        
    def setSpecHeat(self, specHeat):
        self.secData['specHeat'] = specHeat
            
    def setElementSet(self,elsetName):
        self.secData['elementSet'] = elsetName