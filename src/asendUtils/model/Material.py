# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 11:43:42 2023

@author: evans
"""

class Material:
    
    def __init__(self,name):
        self.matData = dict()
        self.matData['name'] = name
        
    def setDensity(self,density):
        self.matData['density'] = density
        
    def setOrthotropic(self,E1,E2,E3,nu12,nu13,nu23,G12,G13,G23):
        elastic = dict()
        elastic['E'] = str([E1,E2,E3])
        elastic['nu'] = str([nu12,nu13,nu23])
        elastic['G'] = str([G12,G13,G23])
        self.matData['elastic'] = elastic
        
    def setThermalConductivity(self,k11,k22,k33,k12=0.0,k13=0.0,k23=0.0):
        try:
            self.matData['thermal']['conductivity'] = str([k11,k22,k33,k12,k13,k23])
        except:
            thermal = dict()
            thermal['conductivity'] = str([k11,k22,k33,k12,k13,k23])
            self.matData['thermal'] = thermal
            
    def setThermalExpansion(self,E11,E22,E33,E12=0.0,E13=0.0,E23=0.0):
        try:
            self.matData['thermal']['expansion'] = str([E11,E22,E33,E12,E13,E23])
        except:
            thermal = dict()
            thermal['expansion'] = str([E11,E22,E33,E12,E13,E23])
            self.matData['thermal'] = thermal
            
    def setSpecificHeat(self,specHeat):
        try:
            self.matData['thermal']['specHeat'] = specHeat
        except:
            thermal = dict()
            thermal['specHeat'] = specHeat
            self.matData['thermal'] = thermal