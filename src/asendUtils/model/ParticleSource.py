# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 09:32:52 2025

@author: evaande
"""

import numpy as np
from asendUtils.model.Section import *

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
            
    def setCoordinates(self, x, y, z, refNodes=None, time=None):
        if refNodes != None:
            self.data['refNodes'] = str(refNodes)
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

def sourceGroupFromMesh(meshData, elsPerSource, massPerEl, specHeat, resXRng, resYRng, resZRng, elementSet="", refNodes=None, velInLocal=True, xRange=None, yRange=None, zRange=None, activeTime=None):
    inNds = meshData['nodes']
    numSrc = len(inNds)
    totParts = elsPerSource*numSrc
    srcLst = list()
    elsets = dict()
    sectns = list()
    for s in range(0, numSrc):
        snm = elementSet + '_' + str(s)
        newSrc = ParticleSource(elementSet=snm, refNodes=refNodes, velInLocal=velInLocal, xRange=xRange, yRange=yRange, zRange=zRange, activeTime=activeTime)
        newSrc.setCoordinates(inNds[s,0], inNds[s,1], inNds[s,2])
        srcLst.append(newSrc)
        elsets[snm] = list(range(s*elsPerSource, (s+1)*elsPerSource))
        newSec = Section('mass')
        newSec.setElementSet(snm)
        newSec.setMassPerElement(massPerEl)
        newSec.setSpecHeat(specHeat)
        sectns.append(newSec)
        

    nodes = list()
    elements = list()
    xLen = resXRng[1] - resXRng[0]
    yLen = resYRng[1] - resYRng[0]
    zLen = resZRng[1] - resZRng[0]
    base = totParts*xLen*xLen/(zLen*yLen)
    xRowsFlt = np.pow(base, 0.3333333333)
    xRows = int(np.ceil(xRowsFlt))
    yRows = int(np.ceil(xRows*yLen/xLen))
    zRows = int(np.ceil(xRows*zLen/xLen))
    xInc = xLen/xRows
    yInc = yLen/yRows
    zInc = zLen/zRows
    ct = 0
    for i in range(0, xRows):
        x = resXRng[0] + i*xInc
        for j in range(0, yRows):
            y = resYRng[0] + j*yInc
            for k in range(0, zRows):
                z = resZRng[0] + k*zInc
                if ct < totParts:
                    nodes.append([x,y,z])
                    elements.append(ct)
                    ct += 1
    outMesh = dict()
    outMesh['nodes'] = np.array(nodes)
    outMesh['elements'] = np.array(elements)
    outMesh['sets'] = {'node': dict(), 'element': elsets}
    
    return srcLst, sectns, outMesh
    
    
    