# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 09:32:52 2025

@author: evaande
"""

class ParticleSource:
    def __init__(self, elementSet="", refNodes=None, velInLocal=True, xRange=None, yRange=None, zRange=None, activeTime=None):
        self.data = dict()
        self.data['elementSet'] = elementSet
        self.data['coordinates'] = list()
        if refNodes != None:
            self.data['refNodes'] = str(refNodes)
        if not velInLocal:
            self.data['velInLocal'] = 'yes'
        self.data['meanVel'] = list()
        self.data['temperature'] = list()
        self.data['frequency'] = list()
        if xRange != None:
            self.data['boundXRange'] = str(xRange)
        else:
            self.data['boundXRange'] = '[-1.0e+100, 1.0e+100]'
        if yRange != None:
            self.data['boundYRange'] = str(yRange)
        else:
            self.data['boundYRange'] = '[-1.0e+100, 1.0e+100]'
        if zRange != None:
            self.data['boundZRange'] = str(zRange)
        else:
            self.data['boundZRange'] = '[-1.0e+100, 1.0e+100]'
        if activeTime != None:
            self.data['activeTime'] = str(activeTime)
        else:
            self.data['activeTime'] = '[0.0, 1.0e+100]'
            
    def setCoordinates(self, x, y, z, refNode=None, time=None):
        if refNode != None:
            self.data['refNode'] = refNode
        if time == None:
            s1 = str([0., x, y, z])
            s2 = str([1.0e+100, x, y, z])
            self.data['coordinates'] = [s1,s2]
        else:
            clst = list()
            for i, t in enumerate(time):
                s = str([t, x[i], y[i], z[i]])
                clst.append(s)
            self.data['coordinates'] = clst
            
    def setVelocity(self, vx, vy, vz, vRandom=None, time=None):
        if vRandom != None:
            self.data['randomVel'] = vRandom
        if time != None:
            s1 = str([0., vx, vy, vz])
            s2 = str([1.0e+100, vx, vy, vz])
            self.data['meanVel'] = [s1,s2]
        else:
            clst = list()
            for i, t in enumerate(time):
                s = str([t, vx[i], vy[i], vz[i]])
                clst.append(s)
            self.data['meanVel'] = clst
            
    def setTemperature(self, temperature, time=None):
        if time != None:
            s1 = str([0., temperature])
            s2 = str([1.0e+100, temperature])
            self.data['temperature'] = [s1,s2]
        else:
            clst = list()
            for i, t in enumerate(time):
                s = str([t, temperature[i]])
                clst.append(s)
            self.data['temperature'] = clst
            
    def setFrequency(self, frequency, time=None):
        if time != None:
            s1 = str([0., frequency])
            s2 = str([1.0e+100, frequency])
            self.data['frequency'] = [s1,s2]
        else:
            clst = list()
            for i, t in enumerate(time):
                s = str([t, frequency[i]])
                clst.append(s)
            self.data['frequency'] = clst