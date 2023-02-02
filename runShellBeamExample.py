## Run the example test number 1 for the shell model of a cantilever beam, which solves
## for the deflection under transverse loading in the Z-direction at the free tip,
## and visualize the transverse deflection distribution

import os
from pythonUtilities.ResultsProcessorClass import *

## Check to make sure output directory exist for shell canteliver beam model

if(not os.path.exists('examples/shellBeam/out')):
    os.system('mkdir examples/shellBeam/out')

## Execute the job, specified in the script runTest1.txt

os.system('./bin/ASenD_runJob.exe examples/shellBeam/runTest1.txt')

## Visualize the results

mFile = 'examples/shellBeam/shellBeamModelInput.yaml'
nFile = 'examples/shellBeam/out/test1Disp.yaml'

rp = ResultsProcessor(mFile,nodeResFile=nFile)
rp.plotNodeResults('displacement',3)

