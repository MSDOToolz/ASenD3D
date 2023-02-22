## Run the test cases in the examples directory, saving the results in the appropriate 'out' folder

import os

currDir = os.getcwd()

if('ASenD3D' in currDir):
    dirLst = currDir.split('ASenD3D')
    mainDir = dirLst[0] + 'ASenD3D'
    binDir = 'bin/ASenD_runJob.exe'
    exDir = 'examples'
    
    os.chdir(mainDir)
    
    if(not os.path.exists('examples/compositePlate/out')):
        os.system('mkdir examples/compositePlate/out')
    
    commandStr = './' + binDir + ' ' + exDir + '/compositePlate/runStaticTest.txt'
    os.system(commandStr)
    
    if(not os.path.exists('examples/hexBeam/out')):
        os.system('mkdir examples/hexBeam/out')
        
    for i in range(1,9):
        commandStr = './' + binDir + ' ' + exDir + '/hexBeam/runTest' + str(i) + '.txt'
        os.system(commandStr)
    
    if(not os.path.exists('examples/shellBeam/out')):
        os.system('mkdir examples/shellBeam/out')
        
    for i in range(1,9):
        commandStr = './' + binDir + ' ' + exDir + '/shellBeam/runTest' + str(i) + '.txt'
        os.system(commandStr)
        
    if(not os.path.exists('examples/tetBeam/out')):
        os.system('mkdir examples/tetBeam/out')
        
    for i in range(1,9):
        commandStr = './' + binDir + ' ' + exDir + '/tetBeam/runTest' + str(i) + '.txt'
        os.system(commandStr)
        
    os.chdir(currDir)
    ##os.system("path/to/example.exe")
else:
    print('Error: unable to find examples directory.')
    print('Please navigate to your working directory for ASenD3D or set the directories for example tests')
    print('and binary executable in the variables binDir and exDir in runExampleTests.py')

##os.chdir('/tmp')