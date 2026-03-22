# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 14:30:24 2025

@author: evaande
"""

import numpy as np

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
        if coefTime != None:
            clst = list()
            for i, t in enumerate(coefTime):
                clst.append(str([t, coefficient[i]]))
            self.data['potField'] = {'coef': clst, 'exp': exponent}
        else:
            clst = list()
            clst.append(str([0., coefficient]))
            clst.append(str([1.0e+100, coefficient]))
            self.data['potField'] = {'coef': clst, 'exp': exponent}
            
    def setDampingField(self, coefficient, exponent, coefTime=None):
        if coefTime != None:
            clst = list()
            for i, t in enumerate(coefTime):
                clst.append(str([t, coefficient[i]]))
            self.data['dampField'] = {'coef': clst, 'exp': exponent}
        else:
            clst = list()
            clst.append(str([0., coefficient]))
            clst.append(str([1.0e+100, coefficient]))
            self.data['dampField'] = {'coef': clst, 'exp': exponent}
            
    def setThermalField(self, conductionCoef=0., radiationCoef=0., referenceTemp=0.):
        self.data['thermField'] = {'condCoef': conductionCoef, 'radCoef': radiationCoef, 'refTemp': referenceTemp}

def calcPotForce(surfNodes, freeNode, coef, expnt):
    F = np.zeros(3, dtype=float)
    for n in surfNodes:
        dvec = freeNode - n
        dmag = np.linalg.norm(dvec)
        fmag = coef/np.power(dmag, expnt)
        F += (fmag/dmag)*dvec
    return F

def dPotFdX(surfNodes, freeNode, coef, expnt, elsz):
    dFdX = np.zeros((3,3), dtype=float)
    Fn = calcPotForce(surfNodes, freeNode, coef, expnt)
    dx = 0.01*elsz
    for i in range(0,3):
        xp = freeNode.copy()
        xp[i] += dx
        Fp = calcPotForce(surfNodes, xp, coef, expnt)
        dFdX[:,i] = (1.0/dx)*(Fp - Fn)
    return dFdX
        
def contactInteraction(maxNormalStress, frictionCoef, elementSize, exp=4.0, name=None, nodeSet1="", nodeSet2="", maxDistance=None, activeTime=None):
    hsz = 0.5*elementSize
    surfNodes = np.array([[-hsz,-hsz,0.], 
                          [hsz,-hsz,0.], 
                          [hsz,hsz,0.], 
                          [-hsz,hsz,0.]])
    maxF = maxNormalStress*elementSize*elementSize
    appF = np.array([frictionCoef*maxF, 0., -maxF])
    coef = 0.0
    #exp = 1.0
    #dexp = 1.0
    loopct = 0
    ht = elementSize
    hfact = 0.5
    dFdX = np.zeros((3,3), dtype=float)
    while hfact < 0.99 and loopct < 100:
        #freeNd = np.array([0.,0.,hsz])
        freeNd = np.array([0., 0., ht])
        PF = calcPotForce(surfNodes, freeNd, 1.0, exp)
        coef = maxF/np.linalg.norm(PF)
        
        nlit = 0
        dxmag = elementSize
        while dxmag > 1.0e-6*elementSize and nlit < 50:
            totF = appF + calcPotForce(surfNodes, freeNd, coef, exp)
            dFdX = dPotFdX(surfNodes, freeNd, coef, exp, elementSize)
            try:
                dx = np.linalg.solve(dFdX, -totF)
                freeNd += dx
                dxmag = np.linalg.norm(dx)
                nlit += 1
            except:
                nlit = 50
                
        if nlit == 50:
            ht *= hfact
        else:
            ht /= hfact
            hfact = np.sqrt(hfact)
            ht *= hfact
        
        # k = (exp/hsz)*maxF
        # omega = np.sqrt(k)
        # dt = 0.2*np.pi/omega
        # prevX = np.array([0.,0.,hsz])
        # ppX = np.array([0.,0.,hsz])
        # for i in range(0,100):
        #     totF = appF + calcPotForce(surfNodes, freeNd, coef, exp)
        #     xNext = (dt*dt)*totF + 2.0*prevX - ppX
        #     ppX = prevX
        #     prevX = freeNd
        #     freeNd = xNext
        # if np.linalg.norm(freeNd) > 2*ht:
        #     #exp += dexp
        #     ht *= hfact
        # else:
        #     #exp += (0.63 - 1.0)*dexp
        #     #dexp *= 0.63
        #     ht /= hfact
        #     hfact = np.sqrt(hfact)
        #     ht *= hfact
        
        loopct += 1
        
    if loopct == 100:
        print("Warning: did not converge to a set of contact interaction parameters")
        
    print('Interaction: ' + name + ' equilibrium gap distance: ' + str(ht/elementSize) + ' X (element size)')
        
    if maxDistance == None:
        mD = elementSize
    else:
        mD = maxDistance
    newInt = Interaction(name=name, nodeSet1=nodeSet1, nodeSet2=nodeSet2, maxDistance=mD, maxNeighbors=4, activeTime=activeTime)
    newInt.setPotentialField(-coef, exp)
    return newInt

def collisionInteraction(mass, velocity, nearDist, expnt=4, numNdPairs=1, name=None, nodeSet1="", nodeSet2="", maxDistance=None, maxNeighbors=None, maxDistRatio=None, idealGasConstant=None, activeTime=None):
    intExp = -expnt + 1.0
    ke = 0.5*mass*velocity*velocity
    coef = ke*intExp/(np.power(nearDist,intExp)*numNdPairs)
    newInt = Interaction(name=name, nodeSet1=nodeSet1, nodeSet2=nodeSet2, maxDistance=maxDistance, maxNeighbors=maxNeighbors, maxDistRatio=maxDistRatio, activeTime=activeTime)
    newInt.setPotentialField(coef, expnt)
    return newInt
        
def idealGasInteraction(idealGasConst, specificHeat, conductivity, viscosity, spacing, refTemp=0.0, name=None, nodeSet1="", nodeSet2="", maxDistance=None, maxNeighbors=12, maxDistRatio=1.8, activeTime=None):
    potExp = 3.0*idealGasConst/specificHeat + 1.0
    area = 0.08333333333333*np.pi*spacing*spacing
    dampCoef = area*viscosity
    condCoef = area*conductivity
    newInt = Interaction(name=name, nodeSet1=nodeSet1, nodeSet2=nodeSet2, maxDistance=maxDistance, maxNeighbors=maxNeighbors, maxDistRatio=maxDistRatio, activeTime=activeTime)
    newInt.setPotentialField(0.0, potExp)
    newInt.setDampingField(dampCoef, 1.0)
    newInt.setThermalField(conductionCoef=condCoef,referenceTemp=refTemp)
    return newInt

def evalIncompR(cnst, bulkModulus, libEnergy, spacing, peakDist):
    ## cnst = [c1, exp1, c2, exp2]
    rVec = np.zeros(4, dtype=float)
    rVec[0] = cnst[0]*np.power(spacing, -cnst[1]) + cnst[2]*np.power(spacing, -cnst[3])
    rVec[1] = -cnst[1]*cnst[0]*np.power(spacing, -(cnst[1] + 1.0)) - cnst[3]*cnst[2]*np.power(spacing, -(cnst[3] + 1.0)) - 0.25*bulkModulus*np.pi*spacing
    rVec[2] = -cnst[1]*cnst[0]*np.power(peakDist, -(cnst[1] + 1.0)) - cnst[3]*cnst[2]*np.power(peakDist, -(cnst[3] + 1.0))
    rVec[3] = cnst[0]/( (cnst[1] - 1.0)*np.power(spacing, -(cnst[1] - 1.0)) ) + cnst[2]/( (cnst[3] - 1.0)*np.power(spacing, -(cnst[3] - 1.0)) ) - libEnergy
    return rVec

def incompFluidInteraction(bulkModulus, libEnergy, conductivity, viscosity, spacing, peakDist, refTemp=0.0, name=None, nodeSet1="", nodeSet2="", maxDistance=None, maxNeighbors=12, maxDistRatio=1.8, activeTime=None):
    cnst = np.array([-1.0, 3.0, 1.0, 2.0])
    dcMag = 1.0
    it = 0
    dRdC = np.zeros((4,4), dtype=float)
    while dcMag > 1.0e-12 and it < 20:
        rVec = evalIncompR(cnst, bulkModulus, libEnergy, spacing, peakDist)
        for i in range(0, 4):
            dci = 0.01*np.abs(cnst[i])
            cnst[i] += dci
            dR = evalIncompR(cnst, bulkModulus, libEnergy, spacing, peakDist)
            dRdC[:,i] = (1.0/dci)*(dR - rVec)
            cnst[i] -= dci
        dc = np.linalg.solve(dRdC, -1.0*rVec)
        cnst += dc
        dcMag = np.linalg.norm(dc)
        it += 1
        
    int1 = Interaction(name=(name + '_1'), nodeSet1=nodeSet1, nodeSet2=nodeSet2, maxDistance=maxDistance, maxNeighbors=maxNeighbors, maxDistRatio=maxDistRatio, activeTime=activeTime)
    int1.setPotentialField(cnst[0], cnst[1])
    area = 0.08333333333333*np.pi*spacing*spacing
    dampCoef = area*viscosity
    condCoef = area*conductivity
    int1.setDampingField(dampCoef, 1.0)
    int1.setThermalField(conductionCoef=condCoef, referenceTemp=refTemp)
    
    int2 = Interaction(name=(name + '_2'), nodeSet1=nodeSet1, nodeSet2=nodeSet2, maxDistance=maxDistance, maxNeighbors=maxNeighbors, maxDistRatio=maxDistRatio, activeTime=activeTime)
    int2.setPotentialField(cnst[2], cnst[3])
        
    return [int1, int2]

def evalMatR(cnst, elasticModulus, spacing, peakDist, uts):
    rVec = np.zeros(4, dtype=float)
    rVec[0] = cnst[0]*np.power(spacing, -cnst[1]) + cnst[2]*np.power(spacing, -cnst[3])
    rVec[1] = -cnst[1]*cnst[0]*np.power(spacing, -(cnst[1] + 1.0)) - cnst[3]*cnst[2]*np.power(spacing, -(cnst[3] + 1.0)) - 0.0833333333333*elasticModulus*np.pi*spacing
    rVec[2] = -cnst[1]*cnst[0]*np.power(peakDist, -(cnst[1] + 1.0)) - cnst[3]*cnst[2]*np.power(peakDist, -(cnst[3] + 1.0))
    rVec[3] = cnst[0]*np.power(peakDist, -cnst[1]) + cnst[2]*np.power(peakDist, -cnst[3]) - uts*0.08333333333333*np.pi*peakDist*peakDist
    return rVec 

def materialMechInteraction(elasticModulus, ultimateStrength, ultimateStrain, conductivity, spacing, refTemp=0.0, name=None, nodeSet1="", nodeSet2="", maxDistance=None, maxNeighbors=12, maxDistRatio=1.8, activeTime=None):
    peakDist = spacing*(1.0 + ultimateStrain)
    
    cnst = np.array([-1.0, 3.0, 1.0, 2.0])
    dcMag = 1.0
    it = 0
    dRdC = np.zeros((4,4), dtype=float)
    while dcMag > 1.0e-12 and it < 20:
        rVec = evalMatR(cnst, elasticModulus, spacing, peakDist, ultimateStrength)
        for i in range(0, 4):
            dci = 0.01*np.abs(cnst[i])
            cnst[i] += dci
            dR = evalMatR(cnst, elasticModulus, spacing, peakDist, ultimateStrength)
            dRdC[:,i] = (1.0/dci)*(dR - rVec)
            cnst[i] -= dci
        dc = np.linalg.solve(dRdC, -1.0*rVec)
        cnst += dc
        dcMag = np.linalg.norm(dc)
        it += 1
        
    int1 = Interaction(name=(name + '_1'), nodeSet1=nodeSet1, nodeSet2=nodeSet2, maxDistance=maxDistance, maxNeighbors=maxNeighbors, maxDistRatio=maxDistRatio, activeTime=activeTime)
    int1.setPotentialField(cnst[0], cnst[1])
    area = 0.08333333333333*np.pi*spacing*spacing
    condCoef = area*conductivity
    int1.setThermalField(conductionCoef=condCoef, referenceTemp=refTemp)
    
    int2 = Interaction(name=(name + '_2'), nodeSet1=nodeSet1, nodeSet2=nodeSet2, maxDistance=maxDistance, maxNeighbors=maxNeighbors, maxDistRatio=maxDistRatio, activeTime=activeTime)
    int2.setPotentialField(cnst[2], cnst[3])
    
    return [int1, int2]