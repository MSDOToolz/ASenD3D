# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 18:46:31 2023

@author: evans
"""

import os
fileNames = ['element_equations.rs','element_meth.rs',
              'element_properties.rs','element_soln_fields.rs',
              'face_meth.rs','matrix_functions.rs',
              'node_meth.rs']

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
    os.remove(fn)
    os.rename('temp.out',fn)
    print('copied ' + fn)