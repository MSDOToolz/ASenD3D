import os
import inspect
import yaml

def makeAbsolute(fileName):
    if(os.path.isabs(fileName)):
        absPath = fileName.replace('\\','/')
        return absPath
    else:
        currDir = os.getcwd()
        currDir = currDir.replace('\\','/')
        absPath = currDir + '/' + fileName
        return absPath 

def getSolverPath():
    thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    thisDir = thisDir.replace('\\','/')
    envFn = thisDir + '/environment.yaml'
    inFile = open(envFn,'r')
    envData = yaml.safe_load(inFile)
    inFile.close()
    return envData['solverpath']

def setSolverPath(newPath):
    pth = makeAbsolute(newPath)
    if(not os.path.exists(pth)):
        errStr = 'Error: the path ' + pth + ' does not exist.  Could not set solver path.'
        raise Exception(errStr)
    thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    thisDir = thisDir.replace('\\','/')
    envFn = thisDir + '/environment.yaml'
    outFile = open(envFn,'w')
    outLn = 'solverpath: ' + pth + '\n'
    outFile.write(outLn)
    outFile.close()