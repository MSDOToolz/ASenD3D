class InputCoverter:

    def __init__(self):

    def covertInput(self,inputFile,inputFormat,outputFile):
        if(inputFormat == 'abaqus'):
            self.convertAbaqus(inputFile,outputFile):

    def convertAbaqus(self,inputFile,outputFile):
        inFile = open(inputFile,'r')

