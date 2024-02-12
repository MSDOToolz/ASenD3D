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

def getEnvPath(tag):
    thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    thisDir = thisDir.replace('\\','/')
    envFn = thisDir + '/environment.yaml'
    try:
        inFile = open(envFn,'r')
        envData = yaml.safe_load(inFile)
        inFile.close()
        return envData[tag]
    except:
        errSt = 'Error: the ' + tag + ' path has not been set.  Set it with asendUtils.syst.pathTools.setEnvPath()'
        raise Exception(errSt)

def setEnvPath(tag,newPath):
    pth = makeAbsolute(newPath)
    if(not os.path.exists(pth)):
        errStr = 'Error: the path ' + pth + ' does not exist.  Could not set solver path.'
        raise Exception(errStr)
    thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    thisDir = thisDir.replace('\\','/')
    envFn = thisDir + '/environment.yaml'
    try:
        inFile = open(envFn,'r')
        envData = yaml.safe_load(inFile)
        inFile.close()
        envData[tag] = pth
        outFile = open(envFn,'w')
        yaml.dump(envData,outFile,sort_keys=False)
        outFile.close()
    except:
        envData = dict()
        envData[tag] = pth
        outFile = open(envFn,'w')
        yaml.dump(envData,outFile,sort_keys=False)
        outFile.close()