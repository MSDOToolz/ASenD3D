# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 14:56:54 2023

@author: evans
"""

dt = 0.01

t = 0.0
x = 1.0
v = 0.0
a = 0.0

xSoln = [x]
tSoln = [t]

while (x < 10):
    x = x + v*dt + 0.5*a*dt*dt
    t = t + dt
    xSoln.append(x)
    tSoln.append(t)
    v = v + a*dt
    a = (1.0-v)/(x*x)