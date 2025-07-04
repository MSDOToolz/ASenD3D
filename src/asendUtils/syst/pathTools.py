import os
import inspect
import yaml
import shutil
import subprocess

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
        
def setUserUpdateFunction(path, rebuild_now=True):
    pth = makeAbsolute(path)
    if(not os.path.exists(pth)):
        errStr = 'Error: the path ' + pth + ' does not exist.  Could not set user-defined time-step update function.'
        raise Exception(errStr)
    if('.rs' not in path):
        print('Warning: The path provided for the user time step update function does not appear to be a Rust source file.')
        print('The function must be written in Rust following the example templates in ASenD3D/examples/userUpdateFunctions')
    thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    thisDir = thisDir.replace('\\','/')
    envFn = thisDir + '/environment.yaml'
    inFile = open(envFn,'r')
    envData = yaml.safe_load(inFile)
    inFile.close()
    dstPath = envData['rootpath'] + '/solverSrc/asend_run_job/src/user.rs'
    shutil.copy(path,dstPath)
    if rebuild_now:
        wkd = os.getcwd()
        slvrt = envData['rootpath'] + 'solverSrc/asend_run_job'
        os.chdir(slvrt)
        procRes = subprocess.run(['cargo', 'build', '--release'],capture_output=True,text=True)
        print('core solver re-build return status: ' + str(procRes.returncode))
        print('output:')
        print(procRes.stdout)
        print('err:')
        print(procRes.stderr)
        os.chdir(wkd)
    else:
        print('Warning: new user-defined time step update function will not take effect until core solver is rebuilt.')
    
def clearUserUpdateFunction(rebuild_now=True):
    thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    thisDir = thisDir.replace('\\','/')
    envFn = thisDir + '/environment.yaml'
    inFile = open(envFn,'r')
    envData = yaml.safe_load(inFile)
    inFile.close()
    dstPath = envData['rootpath'] + '/solverSrc/asend_run_job/src/user.rs'
    outFile = open(dstPath,'w')
    outFile.write('use crate::model::*;\n\n')
    outFile.write('impl Model {\n')
    outFile.write('    pub fn user_update_time_step(&mut self) {\n')
    outFile.write('        return;\n')
    outFile.write('    }\n')
    outFile.write('}\n')
    outFile.close()
    if rebuild_now:
        wkd = os.getcwd()
        slvrt = envData['rootpath'] + 'solverSrc/asend_run_job'
        os.chdir(slvrt)
        procRes = subprocess.run(['cargo', 'build', '--release'],capture_output=True,text=True)
        print('core solver re-build return status: ' + str(procRes.returncode))
        print('output:')
        print(procRes.stdout)
        print('err:')
        print(procRes.stderr)
        os.chdir(wkd)
    else:
        print('Warning: new user-defined time step update function will not take effect until core solver is rebuilt.')