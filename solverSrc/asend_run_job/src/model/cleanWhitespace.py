# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 18:46:31 2023

@author: evans
"""

import os
import shutil

fileNames = ['element/element_equations.rs', 'element/element_fluid_eq.rs', 'element/element_fluid_fields.rs', 'element/element_meth.rs',
              'element/element_properties.rs','element/element_soln_fields.rs',
              'face/face_meth.rs','matrix_functions.rs',
              'node/node_meth.rs', 'interaction/interaction_meth.rs']

for fn in fileNames:
    inFile = open(fn,'r')
    outFile = open('temp.out','w')
    fileLine = inFile.readline()
    sinceText = 0
    while(fileLine != ''):
        if(fileLine.isspace()):
            sinceText = sinceText + 1
        else:
            sinceText = 0
        if(sinceText < 3):
            outFile.write(fileLine)
        fileLine = inFile.readline()
    inFile.close()
    outFile.close()

    if '/' in fn:
        lst = fn.split('/')
        os.rename('temp.out', lst[1])
        shutil.move(lst[1], lst[0])
    else:
        os.remove(fn)
        os.rename('temp.out', fn)

    print('copied ' + fn)