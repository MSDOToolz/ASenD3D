# -*- coding: utf-8 -*-
"""
Created on Sun Dec 31 14:01:19 2023

@author: evans
"""

import numpy as np

## define input
layers = list()

newLay = dict()
newLay['E11'] = 10000000.0
newLay['E22'] = 1000000.0
newLay['nu12'] = 0.0
newLay['G12'] = 5000000.0
newLay['thickness'] = 0.05
newLay['angle'] = 45

layers.append(newLay)

newLay = dict()
newLay['E11'] = 10000000.0
newLay['E22'] = 1000000.0
newLay['nu12'] = 0.0
newLay['G12'] = 5000000.0
newLay['thickness'] = 0.05
newLay['angle'] = -45

layers.append(newLay)

offset = 0.0

## calculate
totThk = 0.0
for l in layers:
    totThk = totThk + l['thickness']

ABD = np.zeros((6,6),dtype=float)
zCrd = -0.5*totThk*(1.0 + offset)
for l in layers:
    Smat = np.zeros((3,3),dtype=float)
    Smat[0,0] = 1.0/l['E11']
    Smat[0,1] = -l['nu12']/l['E11']
    Smat[1,0] = Smat[0,1]
    Smat[1,1] = 1.0/l['E22']
    Smat[2,2] = 1.0/l['G12']
    Qmat = np.linalg.inv(Smat)
    angRad = l['angle']*np.pi/180.0
    a11 = np.cos(angRad)
    a21 = np.sin(angRad)
    a12 = -a21
    a22 = a11
    Ts = np.array([[a11*a11,a12*a12,2.0*a11*a12],
                   [a21*a21,a22*a22,2.0*a22*a21],
                   [a11*a21,a12*a22,(a11*a22 + a12*a21)]])
    Te = np.array([[a11*a11,a12*a12,a11*a12],
                   [a21*a21,a22*a22,a22*a21],
                   [2.0*a11*a21,2.0*a12*a22,(a11*a22 + a12*a21)]])
    Teinv = np.linalg.inv(Te)
    TQ = np.matmul(Ts,Qmat)
    Qlam = np.matmul(TQ,Teinv)
    zNext = zCrd + l['thickness']
    ABD[0:3,0:3] = ABD[0:3,0:3] + (zNext - zCrd)*Qlam
    ABD[0:3,3:6] = ABD[0:3,3:6] + 0.5*(zNext*zNext - zCrd*zCrd)*Qlam
    ABD[3:6,0:3] = ABD[3:6,0:3] + 0.5*(zNext*zNext - zCrd*zCrd)*Qlam
    ABD[3:6,3:6] = ABD[3:6,3:6] + 0.333333333333*(zNext*zNext*zNext - zCrd*zCrd*zCrd)*Qlam
    zCrd = zNext

ABDInv = np.linalg.inv(ABD)

F = np.array([1.0,0.,0.,0.,0.,0.])

midStrn = np.matmul(ABDInv,F)

ux = 2.0*midStrn[0]
uy = 2.0*midStrn[1]
thetax = -midStrn[5]
uz = thetax