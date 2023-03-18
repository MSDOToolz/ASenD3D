import os

fileNames = ['NodeClass.h']

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
            
            for ln in lineList:
                if('//' not in ln):
                    newLn = ln.replace('Doub','DiffDoub')
                    outFile.write(newLn)
                    
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
            
            for ln in lineList:
                if('//' not in ln):
                    newLn = ln.replace('Doub','DiffDoub')
                    outFile.write(newLn)
                    
            outFile.write(' \n')
            outFile.write('//Diff2Doub versions: \n')
                    
            for ln in lineList:
                if('//' not in ln):
                    newLn = ln.replace('Doub','Diff2Doub')
                    outFile.write(newLn)
            
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