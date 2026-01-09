# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 18:46:31 2023

@author: evans
"""

import os
import shutil

fileNames = [ 'model/element.rs',
              'model/element/element_equations.rs', 'model/element/element_meth.rs',
              'model/element/element_properties.rs','model/element/element_soln_fields.rs',
              'model/face/face_meth.rs','matrix_functions.rs',
              'model/node/node_meth.rs','model/interaction/interaction_meth.rs']

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
        lstLen = len(lst)
        nm = lst[lstLen-1]
        os.rename('temp.out', nm)
        os.remove(fn)
        pth = '/'.join(lst[0:lstLen-1])
        shutil.move(nm, pth)
    else:
        os.remove(fn)
        os.rename('temp.out', fn)

    print('copied ' + fn)