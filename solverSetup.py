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
setEnvPath('rootpath',thisDir)

pf = os.sys.platform
if('win' in pf):
    solveDir = thisDir + '/VS/mainSolver/ASenD_runJob/x64/Release/ASenD_runJob.exe'
    setEnvPath('solverpath',solveDir)
    meshDir = thisDir + '/VS/mesher/unstrucTetMesh/x64/Release/unstrucTetMesh.exe'
    setEnvPath('mesherpath',meshDir)
else:
    solveDir = thisDir + '/bin/ASenD_runJob.exe'

    ## Compile the core solver
    print('compiling core solver...')
    CC = 'g++'
    
    binDir = thisDir + '/bin'
    if(not os.path.exists(binDir)):
        os.mkdir(binDir)
    
    srcDir = thisDir + '/solverSrc/mainSolver'
    os.chdir(srcDir)
    
    sourceFiles = ['ASenD_runJob.cpp',
                  'ConstraintClass.cpp','ConstraintClass.h',
                  'DesignVariableClass.cpp','DesignVariableClass.h',
                  'DiffDoubClass.cpp','DiffDoubClass.h',
                  'ElementClass.cpp','ElementEquations.cpp','ElementProperties.cpp','ElementSolnFields.cpp','ElementClass.h',
                  'FaceClass.cpp','FaceClass.h',
                  'JobClass.cpp','JobClass.h',
                  'ListEntClass.cpp','ListEntClass.h',
                  'LoadClass.cpp','LoadClass.h',
                  'LowerTriMatClass.cpp','LowerTriMatClass.h',
                  'matrixFunctions.cpp','matrixFunctions.h',
                  'ModelClass.cpp','ModelAnalysis.cpp','ModelInput.cpp','ModelOutput.cpp','ModelClass.h',
                  'NodeClass.cpp','NodeClass.h',
                  'ObjectiveClass.cpp','ObjectiveClass.h',
                  'SectionClass.cpp','SectionClass.h',
                  'SetClass.cpp','SetClass.h']
    
    
    args = [CC,'-o',solveDir]
    
    args.extend(sourceFiles)
    
    procRes = subprocess.run(args,capture_output=True,text=True)
    
    print(procRes.stdout)
    
    if(procRes.returncode != 0):
        errSt = 'Error: Problem compiling the core solver.  Check to make sure ' + CC + ' compiler is properly installed'
        raise Exception(errSt)
    
    setEnvPath('solverpath',solveDir)
    
    meshDir = thisDir + '/bin/unstrucTetMesh.exe'

    ## Compile the 3D unstructured mesher 
    print('compiling unstructured 3D mesher...')
    
    srcDir = thisDir + '/solverSrc/mesher'
    os.chdir(srcDir)
    
    sourceFiles = ['unstrucTetMesh.cpp',
                    'MeshElement.cpp', 'MeshElement.h',
                    'Mesher.cpp', 'Mesher.h',
                    'MeshFace.cpp', 'MeshFace.h',
                    'MeshNode.cpp', 'MeshNode.h',
                    'SpatialGrid.cpp', 'SpatialGrid.h',
                    'utilities.cpp', 'utilities.h']
    
    args = [CC,'-o',meshDir]
    
    args.extend(sourceFiles)
    
    procRes = subprocess.run(args,capture_output=True,text=True)
    
    print(procRes.stdout)
    
    if(procRes.returncode != 0):
        errSt = 'Error: Problem compiling the unstructured 3D mesh generator.  Check to make sure ' + CC + ' compiler is properly installed'
        raise Exception(errSt)
    
    setEnvPath('mesherpath',meshDir)
    
    os.chdir(wkDir)