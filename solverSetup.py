# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 13:55:31 2023

@author: evans
"""
import os
import inspect
from asendUtils.syst.pathTools import *
import subprocess

## Initialize solver path to the default
wkDir = os.getcwd()
thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
thisDir = thisDir.replace('\\','/')
solveDir = thisDir + '/bin/ASenDSolver/x64/Debug/ASenDSolver.exe' ##'/bin/ASenD_runJob.exe'
setSolverPath(solveDir)

## Compile the core solver
# print('compiling core solver...')
# CC = 'g++'

# binDir = thisDir + '/bin'
# if(not os.path.exists(binDir)):
#     os.mkdir(binDir)

# srcDir = thisDir + '/solverSrc/mainSolver'
# os.chdir(srcDir)

# sourceFiles = ['ASenD_runJob.cpp',
#               'ConstraintClass.cpp','ConstraintClass.h',
#               'DesignVariableClass.cpp','DesignVariableClass.h',
#               'DiffDoubClass.cpp','DiffDoubClass.h',
#               'ElementClass.cpp','ElementEquations.cpp','ElementProperties.cpp','ElementSolnFields.cpp','ElementClass.h',
#               'FaceClass.cpp','FaceClass.h',
#               'JobClass.cpp','JobClass.h',
#               'ListEntClass.cpp','ListEntClass.h',
#               'LoadClass.cpp','LoadClass.h',
#               'LowerTriMatClass.cpp','LowerTriMatClass.h',
#               'matrixFunctions.cpp','matrixFunctions.h',
#               'ModelClass.cpp','ModelAnalysis.cpp','ModelInput.cpp','ModelOutput.cpp','ModelClass.h',
#               'NodeClass.cpp','NodeClass.h',
#               'ObjectiveClass.cpp','ObjectiveClass.h',
#               'SectionClass.cpp','SectionClass.h',
#               'SetClass.cpp','SetClass.h']


# args = [CC,'-o',solveDir]

# args.extend(sourceFiles)

# subprocess.run(args,capture_output=True,text=True)


# os.chdir(wkDir)