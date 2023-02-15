import numpy as np

class Segment2D():

    def __init__(self,segType,keyPts,numEls):
        self.segType = segType ## line, curve, arc
        self.keyPts = keyPts
        self.numEls = numEls
        
    def getNodesEdges(self):
        if(self.segType == 'line'):
            pt1 = np.array(self.keyPts[0])
            pt2 = np.array(self.keyPts[1])
            proj = pt2 - pt1
            steps = (1.0/self.numEls)*np.array(range(0,self.numEls+1))
            nds = list()
            for st in steps:
                nd = pt1 + st*proj
                nds.append(nd)
            nodes = np.array(nds)
            eN1 = np.array(range(0,self.numEls),dtype=int)
            eN2 = np.array(range(1,self.numEls+1),dtype=int)
            edges = np.transpose(np.array([eN1,eN2]))
            output = dict()
            output['nodes'] = nodes
            output['edges'] = edges
            return output
        elif(segType == 'curve'):
            pass
        elif(segType == 'arc'):
            pass