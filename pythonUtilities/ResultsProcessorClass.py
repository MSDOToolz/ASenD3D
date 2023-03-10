from ruamel.yaml import YAML
import plotly.graph_objects as go

class ResultsProcessor:

    def __init__(self,modelFile,nodeResFile='none',elementResFile='none',modalResFile='none',objResFile='none'):
        self.yamlReader = YAML()

        inFile = open(modelFile,'r')
        self.modelData = self.yamlReader.load(inFile)
        inFile.close()

        if(not nodeResFile == 'none'):
            self.loadNodeResults(nodeResFile)
        else:
            self.nodeData = dict()

        if(not elementResFile == 'none'):
            self.loadElementResults(elementResFile)
        else:
            self.elementData = dict()
            
        if(not modalResFile == 'none'):
            self.loadModalResults(modalResFile)
        else:
            self.modalData = dict()
            
        if(not objResFile == 'none'):
            self.loadObjectiveResults(objResFile)
        else:
            self.objectiveData = dict()

    def loadNodeResults(self,nodeResFile):
        inFile = open(nodeResFile,'r')
        self.nodeData = self.yamlReader.load(inFile)
        inFile.close()

    def loadElementResults(self,elementResFile):
        inFile = open(elementResFile,'r')
        self.elementData = self.yamlReader.load(inFile)
        inFile.close()
        
    def loadModalResults(self,modalResFile):
        inFile = open(modalResFile,'r')
        self.modalData = self.yamlReader.load(inFile)
        inFile.close()
        
    def loadObjectiveResults(self,objResFile):
        inFile = open(objResFile,'r')
        self.objectiveData = self.yamlReader.load(inFile)
        inFile.close()

    def plotNodeResults(self,field,component=1):
        xLst = []
        yLst = []
        zLst = []
        for nd in self.modelData['nodes']:
            xLst.append(nd[1])
            yLst.append(nd[2])
            zLst.append(nd[3])
        value = []
        minVal = 10000000000000
        maxVal = -minVal
        for nd in self.nodeData[field]:
            val = nd[component]
            value.append(val)
            if(val < minVal):
                minVal = val
            elif(val > maxVal):
                maxVal = val
        midVal = 0.5*(minVal + maxVal)
        v1 = []
        v2 = []
        v3 = []
        for et in self.modelData['elements']:
            if('brick' in et['type']):
                for el in et['connectivity']:
                    n1 = el[1] - 1
                    n2 = el[2] - 1
                    n3 = el[3] - 1
                    n4 = el[4] - 1
                    n5 = el[5] - 1
                    n6 = el[6] - 1
                    n7 = el[7] - 1
                    n8 = el[8] - 1
                    
                    v1.append(n1)
                    v2.append(n4)
                    v3.append(n5)

                    v1.append(n4)
                    v2.append(n5)
                    v3.append(n8)

                    v1.append(n2)
                    v2.append(n3)
                    v3.append(n6)

                    v1.append(n3)
                    v2.append(n6)
                    v3.append(n7)

                    v1.append(n1)
                    v2.append(n2)
                    v3.append(n5)

                    v1.append(n2)
                    v2.append(n5)
                    v3.append(n6)

                    v1.append(n3)
                    v2.append(n4)
                    v3.append(n7)

                    v1.append(n4)
                    v2.append(n7)
                    v3.append(n8)

                    v1.append(n1)
                    v2.append(n2)
                    v3.append(n3)

                    v1.append(n1)
                    v2.append(n3)
                    v3.append(n4)

                    v1.append(n5)
                    v2.append(n6)
                    v3.append(n7)

                    v1.append(n5)
                    v2.append(n7)
                    v3.append(n8)
            elif('tet' in et['type']):
                for el in et['connectivity']:
                    n1 = el[1] - 1
                    n2 = el[2] - 1
                    n3 = el[3] - 1
                    n4 = el[4] - 1
                    
                    v1.append(n1)
                    v2.append(n2)
                    v3.append(n3)

                    v1.append(n1)
                    v2.append(n2)
                    v3.append(n4)

                    v1.append(n1)
                    v2.append(n3)
                    v3.append(n4)

                    v1.append(n2)
                    v2.append(n3)
                    v3.append(n4)
            elif('shell4' in et['type']):
                for el in et['connectivity']:
                    n1 = el[1] - 1
                    n2 = el[2] - 1
                    n3 = el[3] - 1
                    n4 = el[4] - 1
                    
                    v1.append(n1)
                    v2.append(n2)
                    v3.append(n3)

                    v1.append(n1)
                    v2.append(n3)
                    v3.append(n4)
            elif('shell3' in et['type']):
                for el in et['connectivity']:
                    n1 = el[1] - 1
                    n2 = el[2] - 1
                    n3 = el[3] - 1
                    
                    v1.append(n1)
                    v2.append(n2)
                    v3.append(n3)
        cbTitle = field + str(component)
        fig = go.Figure(data=[
            go.Mesh3d(
                x=xLst,
                y=yLst,
                z=zLst,
                colorbar_title = cbTitle,
                colorscale=[[0.0, 'blue'],
                            [0.5, 'yellow'],
                            [1.0, 'red']],
                intensity=value,
                i=v1,
                j=v2,
                k=v3,
                name='',
                showscale=True
            )
        ])

        fig.show()
