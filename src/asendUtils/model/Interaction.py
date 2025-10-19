# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 14:30:24 2025

@author: evaande
"""

class Interaction:
    def __init__(self, name=None, nodeSet1="", nodeSet2="", maxDistance=None, maxNeighbors=None, maxDistRatio=None, idealGasConstant=None, activeTime=None):
        self.name = name
        self.data = dict()
        self.data['nodeSet1'] = nodeSet1
        self.data['nodeSet2'] = nodeSet2
        if maxDistance != None:
            self.data['maxDistance'] = maxDistance
        if maxNeighbors != None:
            self.data['maxNeighbors'] = maxNeighbors
        if maxDistRatio != None:
            self.data['maxDistRatio'] = maxDistRatio
        if idealGasConstant != None:
            self.data['idealGasConstant'] = idealGasConstant
        if activeTime == None:
            self.data['activeTime'] = '[0., 1.0e+100]'
        else:
            self.data['activeTime'] = str(activeTime)
        self.data['potField'] = {'coef': ['[0., 0.]', '[1.0e+100, 0.]'], 'exp': 1.0}
        self.data['dampField'] = {'coef': ['[0., 0.]', '[1.0e+100, 0.]'], 'exp': 1.0}
        self.data['thermField'] = {'condCoef': 0., 'radCoef': 0., 'refTemp': 0.}
        
    def setPotentialField(self, coefficient, exponent, coefTime=None):
        try:
            clst = list()
            for i, c in enumerate(coefficient):
                clst.append(str([coefTime[i], c]))
            self.data['potField'] = {'coef': clst, 'exp': exponent}
        except:
            clst = list()
            clst.append(str([0., coefficient]))
            clst.append(str([1.0e+100, coefficient]))
            self.data['potField'] = {'coef': clst, 'exp': exponent}
            
    def setDampingField(self, coefficient, exponent, coefTime=None):
        try:
            clst = list()
            for i, c in enumerate(coefficient):
                clst.append(str([coefTime[i], c]))
            self.data['dampField'] = {'coef': clst, 'exp': exponent}
        except:
            clst = list()
            clst.append(str([0., coefficient]))
            clst.append(str([1.0e+100, coefficient]))
            self.data['dampField'] = {'coef': clst, 'exp': exponent}
            
    def setThermalField(self, conductionCoef=0., radiationCoef=0., referenceTemp=0.):
        self.data['thermField'] = {'condCoef': conductionCoef, 'radCoef': radiationCoef, 'refTemp': referenceTemp}
        