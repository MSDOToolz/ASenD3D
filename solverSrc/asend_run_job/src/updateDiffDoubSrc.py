import os
import shutil

fileNames = [ 'model/element.rs',
              'model/element/element_equations.rs', 'model/element/element_meth.rs',
              'model/element/element_properties.rs','model/element/element_soln_fields.rs',
              'model/face/face_meth.rs','matrix_functions.rs',
              'model/node/node_meth.rs','model/interaction/interaction_meth.rs']

##fileNames = ['matrixFunctions.cpp','matrixFunctions.h']

for fn in fileNames:
    inFile = open(fn,'r')
    outFile = open('temp.out','w')
    
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
    
    if '/' in fn:
        lst = fn.split('/')
        lstLen = len(lst)
        nm = lst[lstLen-1]
        os.rename('temp.out', nm)
        os.remove(fn)
        pth = '/'.join(lst[0:lstLen-1])
        shutil.move(nm, pth)
    else:
        os.remove(fn)
        os.rename('temp.out', fn)
    
    print('copied '+fn)