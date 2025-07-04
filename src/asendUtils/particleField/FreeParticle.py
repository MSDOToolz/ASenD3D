# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 19:49:56 2024

@author: evans
"""

class FreeParticle:
    
    def __init__(self,position=None,velocity=None,mass=None,temp=None):
        self.position = position
        self.velocity = velocity 
        self.mass = mass
        self.temp = temp