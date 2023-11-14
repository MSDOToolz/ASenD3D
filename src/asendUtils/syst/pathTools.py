import os

def makeAbsolute(fileName):
    if(':' in fileName):
        return fileName
    else:
        currDir = os.getcwd()
        currDir = currDir.replace('\\','/')
        absPath = currDir + '/' + fileName
        return absPath 
