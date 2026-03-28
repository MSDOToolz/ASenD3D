# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 06:37:06 2026

@author: evans
"""

import numpy as np

def meshFromScratch(nodes,elements):
    meshData = dict()
    meshData['nodes'] = np.array(nodes)
    meshData['elements'] = np.array(elements)
    return meshData