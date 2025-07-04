# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 14:12:22 2025

@author: evaande
"""

class Fluid:
    
    def __init__(self,name):
        self.name = name
        self.flData = dict()
        
    def setViscosity(self,viscosity):
        self.flData['viscosity'] = viscosity
        
    def setIdealGasConstant(self,idealGas):
        self.flData['idealGas'] = idealGas
        
    def setThermalConductivity(self,thermCond):
        try:
            self.flData['thermal']['conductivity'] = thermCond
        except:
            self.flData['thermal'] = {'conductivity': thermCond}
    
    def setSpecificHeat(self,specHeat):
        try:
            self.flData['thermal']['specHeat'] = specHeat
        except:
            self.flData['thermal'] = {'specHeat': specHeat}