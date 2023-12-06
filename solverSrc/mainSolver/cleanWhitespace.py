# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 18:46:31 2023

@author: evans
"""

import os

fileNames = ['ElementClass.cpp','ElementClass.h',
              'ElementEquations.cpp','ElementProperties.cpp',
              'ElementSolnFields.cpp','FaceClass.h',
              'FaceClass.cpp','matrixFunctions.cpp',
              'matrixFunctions.h','NodeClass.cpp','NodeClass.h']

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