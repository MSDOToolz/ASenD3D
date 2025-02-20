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

solveRt = thisDir + '/solverSrc/asend_run_job'
os.chdir(solveRt)

print('Setting up core solver...')
procRes = subprocess.run(['cargo', 'build', '--release'],capture_output=True,text=True)

if(procRes.returncode != 0):
    pf = os.sys.platform
    if('win' in pf):
        print('It looks like rustup, the tool for building the core solver locally is not installed.')
        print('You may attempt to use the included pre-build of the core solver, though it is recommended to build locally.')
        print('rustup can be installed from www.rust-lang.org/tools/install')
    else:
        print('It looks like rustup, the tool for building the core solver locally is not installed.')
        print('It can be installed from the command line, with')
        print("    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh")
        userIn = input('Would you like to install it now? (y/n)\n')
        userIn = userIn.strip()
        if(userIn == 'y' or userIn == 'Y'):
            rustRes = subprocess.run(["curl", "--proto", "'=https'", "--tlsv1.2", "-sSf", "https://sh.rustup.rs", "|", "sh"],capture_output=True,text=True)
            if(rustRes.returncode != 0):
                print("Problem installing rustup.  Please check www.rust-lang.org.tools.install for updates")
            else:
                print("runstup installed.  Retrying build of core solver...")
                procRes = subprocess.run(['cargo', 'build', '--release'],capture_output=True,text=True)
                if(procRes.returncode != 0):
                    print("Core solver build failed.  Please check www.rust-lang.org.tools.install for updates")
                else:
                    print("Core solver build successful")
                    
if(procRes.returncode == 0):
    solveDir = thisDir + '/solverSrc/asend_run_job/target/release/asend_run_job.exe'
    setEnvPath('solverpath', solveDir)
    # meshDir = thisDir + '/solverSrc/unstruc_tet_mesh/target/release/unstruc_tet_mesh.exe'
    # setEnvPath('mesherpath',meshDir)
                    
os.chdir(wkDir)

# pf = os.sys.platform
# if('win' in pf):
#     solveDir = thisDir + '/VS/mainSolver/ASenD_runJob/x64/Release/ASenD_runJob.exe'
#     setEnvPath('solverpath',solveDir)
#     meshDir = thisDir + '/VS/mesher/unstrucTetMesh/x64/Release/unstrucTetMesh.exe'
#     setEnvPath('mesherpath',meshDir)
# else:
#     solveDir = thisDir + '/bin/ASenD_runJob.exe'

#     ## Compile the core solver
#     print('compiling core solver...')
#     CC = 'g++'
    
#     binDir = thisDir + '/bin'
#     if(not os.path.exists(binDir)):
#         os.mkdir(binDir)
    
#     srcDir = thisDir + '/solverSrc/mainSolver'
#     os.chdir(srcDir)
    
#     sourceFiles = ['ASenD_runJob.cpp',
#                   'ConstraintClass.cpp','ConstraintClass.h',
#                   'DesignVariableClass.cpp','DesignVariableClass.h',
#                   'DiffDoubClass.cpp','DiffDoubClass.h',
#                   'ElementClass.cpp','ElementEquations.cpp','ElementProperties.cpp','ElementSolnFields.cpp','ElementClass.h',
#                   'FaceClass.cpp','FaceClass.h',
#                   'JobClass.cpp','JobClass.h',
#                   'ListEntClass.cpp','ListEntClass.h',
#                   'LoadClass.cpp','LoadClass.h',
#                   'LowerTriMatClass.cpp','LowerTriMatClass.h',
#                   'matrixFunctions.cpp','matrixFunctions.h',
#                   'ModelClass.cpp','ModelAnalysis.cpp','ModelInput.cpp','ModelOutput.cpp','ModelClass.h',
#                   'NodeClass.cpp','NodeClass.h',
#                   'ObjectiveClass.cpp','ObjectiveClass.h',
#                   'SectionClass.cpp','SectionClass.h',
#                   'SetClass.cpp','SetClass.h']
    
    
#     args = [CC,'-o',solveDir]
    
#     args.extend(sourceFiles)
    
#     procRes = subprocess.run(args,capture_output=True,text=True)
    
#     print(procRes.stdout)
    
#     if(procRes.returncode != 0):
#         errSt = 'Error: Problem compiling the core solver.  Check to make sure ' + CC + ' compiler is properly installed'
#         raise Exception(errSt)
    
#     setEnvPath('solverpath',solveDir)
    
#     meshDir = thisDir + '/bin/unstrucTetMesh.exe'

#     ## Compile the 3D unstructured mesher 
#     print('compiling unstructured 3D mesher...')
    
#     srcDir = thisDir + '/solverSrc/mesher'
#     os.chdir(srcDir)
    
#     sourceFiles = ['unstrucTetMesh.cpp',
#                     'MeshElement.cpp', 'MeshElement.h',
#                     'Mesher.cpp', 'Mesher.h',
#                     'MeshFace.cpp', 'MeshFace.h',
#                     'MeshNode.cpp', 'MeshNode.h',
#                     'SpatialGrid.cpp', 'SpatialGrid.h',
#                     'utilities.cpp', 'utilities.h']
    
#     args = [CC,'-o',meshDir]
    
#     args.extend(sourceFiles)
    
#     procRes = subprocess.run(args,capture_output=True,text=True)
    
#     print(procRes.stdout)
    
#     if(procRes.returncode != 0):
#         errSt = 'Error: Problem compiling the unstructured 3D mesh generator.  Check to make sure ' + CC + ' compiler is properly installed'
#         raise Exception(errSt)
    
#     setEnvPath('mesherpath',meshDir)
    
#     os.chdir(wkDir)