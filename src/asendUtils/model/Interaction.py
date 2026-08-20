# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 14:30:24 2025

@author: evaande
"""

import numpy as np

class Interaction:
    def __init__(self, name=None, nodeSet1="", nodeSet2="", maxDistance=None, maxNeighbors=None, maxDistRatio=None, idealGasConstant=None, activeTime=None, listAsStr=True):
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
        if listAsStr:
            if activeTime == None:
                self.data['activeTime'] = '[0., 1.0e+100]'
            else:
                self.data['activeTime'] = str(activeTime)
            self.data['potField'] = {'coef': ['[0., 0.]', '[1.0e+100, 0.]'], 'exp': 1.0}
            self.data['dampField'] = {'coef': ['[0., 0.]', '[1.0e+100, 0.]'], 'exp': 1.0}
            self.data['thermField'] = {'condCoef': 0., 'radCoef': 0., 'refTemp': 0.}
        else:
            if activeTime == None:
                self.data['activeTime'] = [0., 1.0e+100]
            else:
                self.data['activeTime'] = activeTime
            self.data['potField'] = {'coef': [[0., 0.], [1.0e+100, 0.]], 'exp': 1.0}
            self.data['dampField'] = {'coef': [[0., 0.], [1.0e+100, 0.]], 'exp': 1.0}
            self.data['thermField'] = {'condCoef': 0., 'radCoef': 0., 'refTemp': 0.}
        
    def setPotentialField(self, coefficient, exponent, coefTime=None, listAsStr=True):
        if coefTime != None:
            clst = list()
            if listAsStr:
                for i, t in enumerate(coefTime):
                    clst.append(str([t, coefficient[i]]))
            else:
                for i, t in enumerate(coefTime):
                    clst.append([t, coefficient[i]])
            self.data['potField'] = {'coef': clst, 'exp': exponent}
        else:
            clst = list()
            if listAsStr:
                clst.append(str([0., coefficient]))
                clst.append(str([1.0e+100, coefficient]))
            else:
                clst.append([0., coefficient])
                clst.append([1.0e+100, coefficient])
            self.data['potField'] = {'coef': clst, 'exp': exponent}
            
    def setDampingField(self, coefficient, distExponent, velExponent, coefTime=None, listAsStr=True):
        if coefTime != None:
            clst = list()
            if listAsStr:
                for i, t in enumerate(coefTime):
                    clst.append(str([t, coefficient[i]]))
            else:
                for i, t in enumerate(coefTime):
                    clst.append([t, coefficient[i]])
            self.data['dampField'] = {'coef': clst, 'distExp': distExponent, 'velExp': velExponent}
        else:
            clst = list()
            if listAsStr:
                clst.append(str([0., coefficient]))
                clst.append(str([1.0e+100, coefficient]))
            else:
                clst.append([0., coefficient])
                clst.append([1.0e+100, coefficient])
            self.data['dampField'] = {'coef': clst, 'distExp': distExponent, 'velExponent': velExponent}
            
    def setMagneticField(self, coefficient, distExponent, velExponent, coefTime=None, listAsStr=True):
        if coefTime != None:
            clst = list()
            if listAsStr:
                for i, t in enumerate(coefTime):
                    clst.append(str([t, coefficient[i]]))
            else:
                for i, t in enumerate(coefTime):
                    clst.append([t, coefficient[i]])
            self.data['magField'] = {'coef': clst, 'distExp': distExponent, 'velExp': velExponent}
        else:
            clst = list()
            if listAsStr:
                clst.append(str([0., coefficient]))
                clst.append(str([1.0e+100, coefficient]))
            else:
                clst.append([0., coefficient])
                clst.append([1.0e+100, coefficient])
            self.data['magField'] = {'coef': clst, 'distExp': distExponent, 'velExponent': velExponent}
            
    def setIdealGas(self, idealGasConst):
        self.data['idealGasConstant'] = idealGasConst
        
    def setIncompProps(self, bulkModulus, expansion, refDensity, refPressure):
        self.data['bulkModulus'] = bulkModulus
        self.data['expansion'] = expansion
        self.data['refDen'] = refDensity
        self.data['refPres'] = refPressure
            
    def setThermalField(self, conductionCoef=0., radiationCoef=0., referenceTemp=0.):
        self.data['thermField'] = {'condCoef': conductionCoef, 'radCoef': radiationCoef, 'refTemp': referenceTemp}

def calcPotForce(surfNodes, freeNode, coef, expnt):
    ## coef[0] = repulsive, coef[1] = attractive
    F = np.zeros(3, dtype=float)
    for n in surfNodes:
        dvec = freeNode - n
        dmag = np.linalg.norm(dvec)
        fmag = coef[0]/np.power(dmag, expnt[0]) - coef[1]/np.power(dmag, expnt[1])
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

def dPotFdCoef(surfNodes, freeNode, coef, expnt):
    dFdC = np.zeros((3,2), dtype=float)
    Fn = calcPotForce(surfNodes, freeNode, coef, expnt)
    coefCp = np.array(coef)
    for i in range(0, 2):
        dc = 0.01*coefCp[i]
        coefCp[i] += dc
        Fp = calcPotForce(surfNodes, freeNode, coefCp, expnt)
        dFdC[:,i] = (1.0/dc)*(Fp - Fn)
        coefCp[i] -= dc
    return dFdC

def contactInteraction(maxNormalStress, frictionCoef, elementSize, hiExp=2.0, distRat=1.5, name=None, nodeSet1="", nodeSet2="", maxDistance=None, activeTime=None):
    hsz = 0.5*elementSize
    surfNodes = np.array([[-hsz,-hsz,0.], 
                          [hsz,-hsz,0.], 
                          [hsz,hsz,0.], 
                          [-hsz,hsz,0.]])
    maxF = maxNormalStress*elementSize*elementSize
    appF = np.array([frictionCoef*maxF, 0., -maxF])
    
    exp = np.array([hiExp, 0.5*hiExp])
    coefMag = maxF*elementSize
    loopct = 0
    ht = elementSize
    hfact = 0.5
    dFdX = np.zeros((3,3), dtype=float)
    while hfact < 0.999 and loopct < 100:
        ## Find coefficients
        coef = np.array([coefMag, coefMag])
        dcMag = coefMag
        nearCrd = np.array([0., 0., ht])
        farCrd = distRat*nearCrd
        
        Fnear = calcPotForce(surfNodes, nearCrd, coef, exp)
        Ffar = calcPotForce(surfNodes, farCrd, coef, exp)
        R = np.array([Fnear[2] - maxF, 
                      Ffar[2]])
        dFnear = dPotFdCoef(surfNodes, nearCrd, coef, exp)
        dFfar = dPotFdCoef(surfNodes, farCrd, coef, exp)
        dRdC = np.array([dFnear[2,:], 
                          dFfar[2,:]])
        dc = np.linalg.solve(dRdC, -R)
        coef += dc
        
        ## Find equilibrium displacement
        freeNd = np.array([0., 0., ht])
        nlit = 0
        maxIt = 15
        dxmag = elementSize
        while dxmag > 1.0e-6*elementSize and nlit < maxIt:
            totF = appF + calcPotForce(surfNodes, freeNd, coef, exp)
            dFdX = dPotFdX(surfNodes, freeNd, coef, exp, elementSize)
            try:
                dx = np.linalg.solve(dFdX, -totF)
                freeNd += dx
                dxmag = np.linalg.norm(dx)
                nlit += 1
            except:
                nlit = maxIt
                
        if freeNd[0] < 0.0 or freeNd[0] > hsz or freeNd[2] < 0.0 or freeNd[2] > 2.0*ht:
            ht *= hfact
        else:
            ht /= hfact
            hfact = np.sqrt(hfact)
            ht *= hfact
        
        loopct += 1
        
    if loopct == 100:
        print("Warning: did not converge to a set of contact interaction parameters")
        
    if name == None:
        updtName = 'contact' + nodeSet1 + nodeSet2
    else:
        updtName = name
        
    print('Interaction: ' + updtName + ' equilibrium gap distance: ' + str(ht/elementSize) + ' X (element size)')
        
    if maxDistance == None:
        mD = distRat*ht
    else:
        mD = maxDistance

    repNm = updtName + '_rep'
    repInt = Interaction(name=repNm, nodeSet1=nodeSet1, nodeSet2=nodeSet2, maxDistance=mD, maxNeighbors=4, activeTime=activeTime)
    repInt.setPotentialField(-coef[0], exp[0])
    
    attNm = updtName + '_att'
    attInt = Interaction(name=attNm, nodeSet1=nodeSet1, nodeSet2=nodeSet2, maxDistance=mD, maxNeighbors=4, activeTime=activeTime)
    attInt.setPotentialField(coef[1], exp[1])
    
    return [repInt, attInt]

def collisionInteraction(mass, velocity, nearDist, expnt=4, numNdPairs=1, name=None, nodeSet1="", nodeSet2="", maxDistance=None, maxNeighbors=None, maxDistRatio=None, idealGasConstant=None, activeTime=None):
    intExp = -expnt + 1.0
    ke = 0.5*mass*velocity*velocity
    coef = ke*intExp/(np.power(nearDist,intExp)*numNdPairs)
    newInt = Interaction(name=name, nodeSet1=nodeSet1, nodeSet2=nodeSet2, maxDistance=maxDistance, maxNeighbors=maxNeighbors, maxDistRatio=maxDistRatio, activeTime=activeTime)
    newInt.setPotentialField(coef, expnt)
    return newInt
        
def idealGasInteraction(idealGasConst, specificHeat, conductivity, viscosity, spacing, 
                        refTemperature=0.0, name=None, nodeSet1="", nodeSet2="", maxDistance=None, maxNeighbors=12, maxDistRatio=1.8, activeTime=None):
    potExp = 3.0*idealGasConst/specificHeat + 1.0
    area = 0.08333333333333*np.pi*spacing*spacing
    dampCoef = area*viscosity
    condCoef = area*conductivity
    newInt = Interaction(name=name, nodeSet1=nodeSet1, nodeSet2=nodeSet2, maxDistance=maxDistance, maxNeighbors=maxNeighbors, maxDistRatio=maxDistRatio, activeTime=activeTime)
    newInt.setPotentialField(0.0, potExp)
    newInt.setDampingField(dampCoef, 1.0, 1.0)
    newInt.setIdealGas(idealGasConst)
    newInt.setThermalField(conductionCoef=condCoef,referenceTemp=refTemperature)
    return newInt

def incompFluidInteraction(bulkModulus, expansion, conductivity, viscosity, spacing, refDensity, refPressure, 
                           refTemperature=0, name=None, nodeSet1="", nodeSet2="", maxDistance=None, maxNeighbors=12, maxDistRatio=1.8, activeTime=None):
    area = 0.08333333333333*np.pi*spacing*spacing
    dampCoef = area*viscosity
    condCoef = area*conductivity
    newInt = Interaction(name=name, nodeSet1=nodeSet1, nodeSet2=nodeSet2, maxDistance=maxDistance, maxNeighbors=maxNeighbors, maxDistRatio=maxDistRatio, activeTime=activeTime)
    newInt.setPotentialField(0.0, 2.0)
    newInt.setDampingField(dampCoef, 1.0, 1.0)
    newInt.setIncompProps(bulkModulus, expansion, refDensity, refPressure)
    newInt.setThermalField(conductionCoef=condCoef,referenceTemp=refTemperature)
    return newInt

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