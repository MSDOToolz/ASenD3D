import os

fileNames = ['ElementClass.cpp','ElementClass.h',
              'ElementEquations.cpp','ElementProperties.cpp',
              'ElementSolnFields.cpp','FaceClass.h',
              'FaceClass.cpp','matrixFunctions.cpp',
              'matrixFunctions.h','NodeClass.cpp','NodeClass.h']

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
                        ln = ln.replace('DiffDoub0','DiffDoub1')
                        ln = ln.replace('_dfd0','_dfd1')
                        outFile.write(ln)
                    
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
                        ln = ln.replace('DiffDoub0','DiffDoub1')
                        ln = ln.replace('_dfd0','_dfd1')
                        outFile.write(ln)
                    
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
                        ln = ln.replace('DiffDoub0','DiffDoub2')
                        ln = ln.replace('_dfd0','_dfd2')
                        outFile.write(ln)
            
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
    
    os.remove(fn)
    os.rename('temp.out',fn)
    
    print('copied '+fn)