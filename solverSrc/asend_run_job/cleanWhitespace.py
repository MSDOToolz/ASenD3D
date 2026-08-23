# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 18:46:31 2023

@author: evans
"""

import os
import shutil
import time

fileNames = [ 'src/model/element.rs',
              'src/model/element/element_equations.rs', 'src/model/element/element_meth.rs',
              'src/model/element/element_properties.rs','src/model/element/element_soln_fields.rs',
              'src/model/face/face_meth.rs','src/matrix_functions.rs',
              'src/model/node/node_meth.rs','src/model/interaction/interaction_meth.rs']

for fn in fileNames:
    lst = fn.split('/')
    lstLen = len(lst)
    nm = lst[lstLen-1]
    
    inFile = open(fn,'r')
    outFile = open(nm,'w')
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

    success = False
    while not success:
        try:
            os.remove(fn)
            success = True
            print('removed file ' + fn)
        except:
            print('Problem encountered deleting file ' + fn)
            print('pausing to try again...')
            time.sleep(5)
             
    pth = '/'.join(lst[0:lstLen-1])
    success = False
    while not success:
        try:
            shutil.move(nm, pth)
            success = True
            print('copied new file ' + nm)
        except:
            print('Problem encountered moving file ' + nm)
            print('pausing to try again...')
            time.sleep(5)