import os
import shutil
import time

fileNames = [ 'src/model/element.rs',
              'src/model/element/element_equations.rs', 'src/model/element/element_meth.rs',
              'src/model/element/element_properties.rs','src/model/element/element_soln_fields.rs',
              'src/model/face/face_meth.rs','src/matrix_functions.rs',
              'src/model/node/node_meth.rs','src/model/interaction/interaction_meth.rs']

##fileNames = ['matrixFunctions.cpp','matrixFunctions.h']

for fn in fileNames:
    lst = fn.split('/')
    lstLen = len(lst)
    nm = lst[lstLen-1]
    
    inFile = open(fn,'r')
    outFile = open(nm,'w')
    
    fileLine = inFile.readline()
    while(fileLine != ''):
        if('//dup1' in fileLine):
            lineList = list()
            while(('//end dup' not in fileLine) and fileLine != ''):
                lineList.append(fileLine)
                fileLine = inFile.readline()
            lineList.append(fileLine)
            
            for ln in lineList:
                outFile.write(ln)
                
            outFile.write(' \n')
            outFile.write('//skip \n')
            outFile.write(' \n')
            outFile.write('//DiffDoub1 versions: \n')
            
            preserve = False
            for ln in lineList:
                if('//preserve' in ln):
                    preserve = True
                elif('//end preserve' in ln):
                    preserve = False
                else:
                    if(preserve):
                        outFile.write(ln)
                    else:
                        newln = ln.replace('DiffDoub0','DiffDoub1')
                        newln = newln.replace('_dfd0','_dfd1')
                        outFile.write(newln)
                    
            outFile.write(' \n')
            outFile.write('//end skip \n')
            fileLine = inFile.readline()
        elif('//dup2' in fileLine):
            lineList = list()
            while(('//end dup' not in fileLine) and fileLine != ''):
                lineList.append(fileLine)
                fileLine = inFile.readline()
            lineList.append(fileLine)
            
            for ln in lineList:
                outFile.write(ln)
            
            outFile.write(' \n')            
            outFile.write('//skip \n')
            outFile.write(' \n')
            outFile.write('//DiffDoub1 versions: \n')
            
            preserve = False
            for ln in lineList:
                if('//preserve' in ln):
                    preserve = True
                elif('//end preserve' in ln):
                    preserve = False
                else:
                    if(preserve):
                        outFile.write(ln)
                    else:
                        newln = ln.replace('DiffDoub0','DiffDoub1')
                        newln = newln.replace('_dfd0','_dfd1')
                        outFile.write(newln)
                    
            outFile.write(' \n')
            outFile.write('//DiffDoub2 versions: \n')
                    
            preserve = False
            for ln in lineList:
                if('//preserve' in ln):
                    preserve = True
                elif('//end preserve' in ln):
                    preserve = False
                else:
                    if(preserve):
                        outFile.write(ln)
                    else:
                        newln = ln.replace('DiffDoub0','DiffDoub2')
                        newln = newln.replace('_dfd0','_dfd2')
                        outFile.write(newln)
            
            outFile.write(' \n')            
            outFile.write('//end skip \n')
            fileLine = inFile.readline()
        elif('//skip' in fileLine):
            while(('//end skip' not in fileLine) and fileLine != ''):
                fileLine = inFile.readline()
            fileLine = inFile.readline()
        else:
            outFile.write(fileLine)
            fileLine = inFile.readline()
    
    inFile.close()
    outFile.close()
    
    success = False
    while not success:
        try:
            os.remove(fn)
            success = True
            print('removed file ' + fn)
        except:
            print('Problem encountered deleting file ' + fn)
            print('pausing to try again...')
            time.sleep(5)
            
    pth = '/'.join(lst[0:lstLen-1])
    success = False
    while not success:
        try:
            shutil.move(nm, pth)
            success = True
            print('copied new file ' + nm)
        except:
            print('Problem encountered moving file ' + nm)
            print('pausing to try again...')
            time.sleep(5)
    