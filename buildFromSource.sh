#!/bin/bash
_wkDir="$PWD"
cd sourceCode

gfortran -static createComplex.f90 -o createComplex.exe

## Create complex version of appropriate source files
./updateComplexSource.sh

## Compile all modules
gfortran -static -c ASenD_constantVals.f90
gfortran -static -c ASenD_globalData.f90
gfortran -static -c ASenD_input.f90
gfortran -static -c ASenD_r_overloadFunctions.f90
gfortran -static -c ASenD_c_overloadFunctions.f90
gfortran -static -c ASenD_solvers.f90
gfortran -static -c ASenD_r_designPropertyFunctions.f90
gfortran -static -c ASenD_c_designPropertyFunctions.f90
gfortran -static -c ASenD_r_elementEqns.f90
gfortran -static -c ASenD_c_elementEqns.f90
gfortran -static -c ASenD_output.f90
gfortran -static -c ASenD_bookKeeping.f90
gfortran -static -c ASenD_objective.f90
gfortran -static -c ASenD_commandFunctions.f90

## Compile main program for execution of ASenD3D job script
gfortran -static ASenD_runJob.f90 ASenD_constantVals.o ASenD_globalData.o ASenD_input.o ASenD_r_overloadFunctions.o ASenD_c_overloadFunctions.o ASenD_solvers.o ASenD_r_designPropertyFunctions.o ASenD_r_elementEqns.o ASenD_c_designPropertyFunctions.o ASenD_c_elementEqns.o ASenD_output.o ASenD_bookKeeping.o ASenD_objective.o ASenD_commandFunctions.o -o $_wkDir/bin/ASenD_runJob.exe

## Clean up unneeded residual files
./cleanUp.sh

cd ..
## End of script
