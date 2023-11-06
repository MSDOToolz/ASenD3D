import os

fileNames = ['ElementClass.cpp','ElementClass.h',
              'ElementEquations.cpp','ElementProperties.cpp',
              'ElementSolnFields.cpp','FaceClass.h',
              'FaceClass.cpp','matrixFunctions.cpp',
              'matrixFunctions.h','NodeClass.cpp','NodeClass.h']

##fileNames = ['matrixFunctions.cpp','matrixFunctions.h']

leadingStrings = [' ','(',',',':','\t']
trailingStrings = [' ','&','*','[','(','Stress']

strings = list()
repStrings1 = list()
repStrings2 = list()

for ls in leadingStrings:
    for ts in trailingStrings:
        lst = [ls,ts]
        strings.append('Doub'.join(lst))
        repStrings1.append('DiffDoub'.join(lst))
        repStrings2.append('Diff2Doub'.join(lst))

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
            outFile.write('//DiffDoub versions: \n')
            
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
                        for i in range(0,len(strings)):
                            ln = ln.replace(strings[i],repStrings1[i])
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
            outFile.write('//DiffDoub versions: \n')
            
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
                        for i in range(0,len(strings)):
                            ln = ln.replace(strings[i],repStrings1[i])
                        outFile.write(ln)
                    
            outFile.write(' \n')
            outFile.write('//Diff2Doub versions: \n')
                    
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
                        for i in range(0,len(strings)):
                            ln = ln.replace(strings[i],repStrings2[i])
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