# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 09:01:50 2023

@author: evans
"""

import os
import sys

def getRootPath():
    aPath = os.path.abspath(sys.argv[0])
    pLst = aPath.split('src')
    rtPath = pLst[0]
    return rtPath