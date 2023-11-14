# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 13:55:31 2023

@author: evans
"""
import os
import inspect
from asendUtils.syst.pathTools import *

thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
thisDir = thisDir.replace('\\','/')
solveDir = thisDir + '/bin/ASenD_runJob.exe'
setSolverPath(solveDir)

## Compile the core solver
# print('compiling core solver...')
# CC = 'g++'
# cmdLn = 'cd ' + rtPath + 'solverSrc/mainSolver'
# os.system('cd solverSrc/mainSolver')

# sourcFiles = ['ASenD_runJob.cpp',
#               'ConstraintClass.cpp','ConstraintClass.h']

# cmdLn = 'cd ' + wrkDir
# os.system(cmdLn)