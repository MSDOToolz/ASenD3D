# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 11:43:42 2023

@author: evans
"""

class Material:
    
    def __init__(self,name):
        self.name = name
        self.matData = dict()
        
    def setDensity(self,density):
        self.matData['density'] = density
        
    def setOrthotropic(self,E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,listAsStr=True):
        elastic = dict()
        if listAsStr:
            elastic['E'] = str([E1,E2,E3])
            elastic['nu'] = str([nu12,nu13,nu23])
            elastic['G'] = str([G12,G13,G23])
        else:
            elastic['E'] = [E1,E2,E3]
            elastic['nu'] = [nu12,nu13,nu23]
            elastic['G'] = [G12,G13,G23]
        self.matData['elastic'] = elastic
        
    def setIsotropic(self,E,nu,listAsStr=True):
        G = E/(2.0*(1.0+nu))
        self.setOrthotropic(E,E,E,nu,nu,nu,G,G,G,listAsStr)
        
    def setThermalConductivity(self,k11,k22,k33,k12=0.0,k13=0.0,k23=0.0,listAsStr=True):
        if listAsStr:
            cnd = str([k11,k22,k33,k12,k13,k23])
        else:
            cnd = [k11,k22,k33,k12,k13,k23]
        try:
            self.matData['thermal']['conductivity'] = cnd
        except:
            thermal = dict()
            thermal['conductivity'] = cnd
            self.matData['thermal'] = thermal
            
    def setThermalExpansion(self,E11,E22,E33,E12=0.0,E13=0.0,E23=0.0,listAsStr=True):
        if listAsStr:
            te = str([E11,E22,E33,E12,E13,E23])
        else:
            te = [E11,E22,E33,E12,E13,E23]
        try:
            self.matData['thermal']['expansion'] = te
        except:
            thermal = dict()
            thermal['expansion'] = te
            self.matData['thermal'] = thermal
            
    def setSpecificHeat(self,specHeat):
        try:
            self.matData['thermal']['specHeat'] = specHeat
        except:
            thermal = dict()
            thermal['specHeat'] = specHeat
            self.matData['thermal'] = thermal
            
    def setMisesStrength(self, misesStrength):
        try:
            self.matData['custom']['misesStrength'] = misesStrength
        except:
            self.matData['custom'] = {'misesStrength': misesStrength}
            
    def setOrthoStrength(self,TS1,TS2,TS3,CS1,CS2,CS3,S12,S13,S23,listAsStr=True):
        if listAsStr:
            try:
                self.matData['custom']['tensileStrength'] = str([TS1,TS2,TS3])
                self.matData['custom']['compressiveStrength'] = str([CS1,CS2,CS3])
                self.matData['custom']['shearStrength'] = str([S12,S13,S23])
            except:
                self.matData['custom'] = {'tensileStrength': str([TS1,TS2,TS3]),
                                         'compressiveStrength': str([CS1,CS2,CS3]),
                                         'shearStrength': str([S12,S13,S23])}
        else:
            try:
                self.matData['custom']['tensileStrength'] = [TS1,TS2,TS3]
                self.matData['custom']['compressiveStrength'] = [CS1,CS2,CS3]
                self.matData['custom']['shearStrength'] = [S12,S13,S23]
            except:
                self.matData['custom'] = {'tensileStrength': [TS1,TS2,TS3],
                                         'compressiveStrength': [CS1,CS2,CS3],
                                         'shearStrength': [S12,S13,S23]}
            
    def addCustomProperty(self,propName,propVal):
        try:
            self.matData['custom'][propName] = propVal
        except:
            self.matData['custom'] = {propName: propVal}