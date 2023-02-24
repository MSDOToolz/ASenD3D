## Run the test cases in the examples directory, saving the results in the appropriate 'out' folder

import os

currDir = os.getcwd()

if('ASenD3D' in currDir):
    dirLst = currDir.split('ASenD3D')
    mainDir = dirLst[0] + 'ASenD3D'
    binDir = 'bin/ASenD_runJob.exe'
    exDir = 'examples/testModels'
    
    os.chdir(mainDir)
    
    outPath = exDir + '/compositePlate/out'
    if(not os.path.exists(outPath)):
        commandStr = 'mkdir ' + outPath
        os.system(commandStr)
    
    commandStr = './' + binDir + ' ' + exDir + '/compositePlate/runTest1.txt'
    os.system(commandStr)
    
    outPath = exDir + '/hexBeam/out'
    if(not os.path.exists(outPath)):
        commandStr = 'mkdir ' + outPath
        os.system(commandStr)
        
    for i in range(1,9):
        commandStr = './' + binDir + ' ' + exDir + '/hexBeam/runTest' + str(i) + '.txt'
        os.system(commandStr)
    
    outPath = exDir + '/shellBeam/out'
    if(not os.path.exists(outPath)):
        commandStr = 'mkdir ' + outPath
        os.system(commandStr)
        
    for i in range(1,9):
        commandStr = './' + binDir + ' ' + exDir + '/shellBeam/runTest' + str(i) + '.txt'
        os.system(commandStr)
        
    outPath = exDir + '/tetBeam/out'
    if(not os.path.exists(outPath)):
        commandStr = 'mkdir ' + outPath
        os.system(commandStr)
        
    for i in range(1,9):
        commandStr = './' + binDir + ' ' + exDir + '/tetBeam/runTest' + str(i) + '.txt'
        os.system(commandStr)
        
    os.chdir(currDir)
else:
    print('Error: unable to find examples directory.')
    print('Please navigate to your working directory for ASenD3D or set the directories for example tests')
    print('and binary executable in the variables binDir and exDir in runExampleTests.py')

##os.chdir('/tmp')