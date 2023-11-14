# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 23:08:45 2023

@author: evans
"""

import os
import sys
import platform

## Set up the root path
wrkDir = os.getcwd()
aPath = os.path.abspath(sys.argv[0])
pLst = aPath.split('env')
rtPath = pLst[0].replace('\\','/')
dest = rtPath + 'src/asendUtils/syst/getRootPath.py'
outFile = open(dest,'w')
outFile.write('\n')
outFile.write('def getRootPath():\n')
outLn = "    return '" + rtPath + "'\n"
outFile.write(outLn)
outFile.write('\n')
outFile.close()

## Initialize the environment file
envFN = rtPath + 'env/environment.yaml'
outFile = open(envFN,'w')
thisPlat = platform.uname()
thisOS = thisPlat.system
ln = 'platform: ' + thisOS + '\n'
outFile.write(ln)
ln = 'solverpath: ' + rtPath + 'bin/run_ASenDJob.exe\n'
outFile.write(ln)
outFile.close()

## Compile the core solver
print('compiling core solver...')
CC = 'g++'
cmdLn = 'cd ' + rtPath + 'solverSrc/mainSolver'
os.system('cd solverSrc/mainSolver')

sourcFiles = ['ASenD_runJob.cpp',
              'ConstraintClass.cpp','ConstraintClass.h']

cmdLn = 'cd ' + wrkDir
os.system(cmdLn)