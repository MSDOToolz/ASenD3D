#from ruamel.yaml import YAML
import numpy as np
import pandas as pd
import yaml
from yaml import CLoader as Loader
from asendUtils.visualization.plotlyUtils import *

class ResultsProcessor:

    def __init__(self,modelFile,dVarFile=None,nodeResFile=None,elementResFile=None,modalResFile=None,objResFile=None):

        inFile = open(modelFile,'r')
        self.modelData = yaml.load(inFile,Loader=Loader)
        inFile.close()
        
        if(not dVarFile == None):
            self.loadDesignVariableInput(dVarFile)

        if(not nodeResFile == None):
            self.loadNodeResults(nodeResFile)
        else:
            self.nodeData = pd.DataFrame()

        if(not elementResFile == None):
            self.loadElementResults(elementResFile)
        else:
            self.elementData = pd.DataFrame()
            
        if(not modalResFile == None):
            self.loadModalVals(modalResFile)
        else:
            self.modalVals = pd.DataFrame()
        self.modalVec = pd.DataFrame()
            
        if(not objResFile == None):
            self.loadObjectiveResults(objResFile)
        else:
            self.objectiveVals = pd.DataFrame()
            self.objectiveGrad = pd.DataFrame()
            
    def loadDesignVariableInput(self,dVarFile):
        inFile = open(dVarFile)
        self.dVarData = yaml.load(inFile,Loader=Loader)
        inFile.close()

    def loadNodeResults(self,nodeResFile):
        self.nodeData = pd.read_csv(nodeResFile,index_col=0)
        
    def loadElementResults(self,elementResFile):
        self.elementData = pd.read_csv(elementResFile)
                
    def loadModalVals(self,modalResFile):
        self.modalVals = pd.read_csv(modalResFile)
                
    def loadModalVec(self,modalResFile,mode=0):
        rep = '_mode' + str(mode) + '.'
        f_name = modalResFile.replace('.',rep)
        self.modalVec = pd.read_csv(f_name,index_col=0)
        
    def loadObjectiveResults(self,objResFile):
        self.objectiveVals = pd.read_csv(objResFile,index_col=0)
        f_name = objResFile.replace('.','_grad.')
        self.objectiveGrad = pd.read_csv(f_name,index_col=0)
                
    def getPlotNdElSet(self,elementSet):
        if(elementSet == 'all'):
            elSet = set()
            for et in self.modelData['elements']:
                for el in et['connectivity']:
                    elSet.add(el[0])
        else:
            elSet = set(self.modelData['sets']['element'][elementSet])
        return elSet
        
    def buildNodalPlotCrd(self,elSet,deformed=False,defScaleFact=1.0,massElOptns=None):
        ## elSet a python set() of the desired element labels
        allNds = self.modelData['nodes']
        numNds = len(allNds)
        
        numMassNds = 0
        if massElOptns != None:
            if massElOptns['showAsDots']:
                for et in self.modelData['elements']:
                    if et['type'] == 'mass':
                        for el in et['connectivity']:
                            if el[0] in elSet:
                                numMassNds += 8
                                
        totNds = numNds + numMassNds
        crdAr = np.zeros((totNds,3),dtype=float)
        
        if(deformed):
            uar = defScaleFact*np.array(self.nodeData.loc[list(range(0,numNds)), ['U1','U2','U3']])
            crdAr[0:numNds] = uar
        
        for nd in allNds:
            lab = nd[0]
            crdAr[lab] += np.array(nd[1:4])
        
        if massElOptns != None:
            if massElOptns['showAsDots']:
                numEls = 0
                for et in self.modelData['elements']:
                    numEls += len(et['connectivity'])
                elMass = np.zeros(numEls, dtype=float)
                for sec in self.modelData['sections']:
                    if sec['type'] == 'mass':
                        for el in self.modelData['sets']['element'][sec['elementSet']]:
                            elMass[el] = sec['massPerEl']
                            
                szConst = massElOptns['refSize']/massElOptns['refMass']
                ind = numNds
                for et in self.modelData['elements']:
                    if et['type'] == 'mass':
                        for el in et['connectivity']:
                            if el[0] in elSet:
                                esz = szConst*np.power(elMass[el[0]], 0.333333333)
                                ndCrd = np.array(allNds[el[1]][1:4])
                                if deformed:
                                    ndCrd += uar[el[1]]
                                shft = 0.5*esz*np.array([[-1, -1, -1],
                                                         [1, -1, -1],
                                                         [1, 1, -1],
                                                         [-1, 1, -1],
                                                         [-1, -1, 1],
                                                         [1, -1, 1],
                                                         [1, 1, 1],
                                                         [-1, 1, 1]])
                                for s in shft:
                                    crdAr[ind] = ndCrd + s
                                    ind += 1
            
        xLst = list(crdAr[:,0])
        yLst = list(crdAr[:,1])
        zLst = list(crdAr[:,2])
            
        return {'xLst': xLst, 'yLst': yLst, 'zLst': zLst}

    def buildElementVertexList(self,elSet,massElOptns=None):
        # elSet = python set() with the desired element labels
        v1 = []
        v2 = []
        v3 = []
        mndInd = len(self.modelData['nodes'])
        for et in self.modelData['elements']:
            if('brick' in et['type'] or (et['type'] == 'mass' and massElOptns['showAsDots']) ):
                for el in et['connectivity']:
                    eli = el[0]
                    if(eli in elSet):
                        if 'brick' in et['type']:
                            n1 = el[1]
                            n2 = el[2]
                            n3 = el[3]
                            n4 = el[4]
                            n5 = el[5]
                            n6 = el[6]
                            n7 = el[7]
                            n8 = el[8]
                        else:
                            n1 = mndInd
                            n2 = mndInd + 1 
                            n3 = mndInd + 2 
                            n4 = mndInd + 3 
                            n5 = mndInd + 4 
                            n6 = mndInd + 5 
                            n7 = mndInd + 6 
                            n8 = mndInd + 7
                            mndInd += 8
                        
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
            elif('wedge' in et['type']):
                for el in et['connectivity']:
                    eli = el[0]
                    if(eli in elSet):
                        n1 = el[1]
                        n2 = el[2]
                        n3 = el[3]
                        n4 = el[4]
                        n5 = el[5]
                        n6 = el[6]
                        
                        v1.append(n1)
                        v2.append(n2)
                        v3.append(n3)
                        
                        v1.append(n4)
                        v2.append(n5)
                        v3.append(n6)
                        
                        v1.append(n1)
                        v2.append(n2)
                        v3.append(n4)
                        
                        v1.append(n2)
                        v2.append(n4)
                        v3.append(n5)
                        
                        v1.append(n1)
                        v2.append(n3)
                        v3.append(n4)
                        
                        v1.append(n3)
                        v2.append(n4)
                        v3.append(n6)
                        
                        v1.append(n2)
                        v2.append(n3)
                        v3.append(n5)
                        
                        v1.append(n3)
                        v2.append(n5)
                        v3.append(n6)
            elif('tet' in et['type']):
                for el in et['connectivity']:
                    eli = el[0]
                    if(eli in elSet):
                        n1 = el[1]
                        n2 = el[2]
                        n3 = el[3]
                        n4 = el[4]
                        
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
                    eli = el[0]
                    if(eli in elSet):
                        n1 = el[1]
                        n2 = el[2]
                        n3 = el[3]
                        n4 = el[4]
                        
                        v1.append(n1)
                        v2.append(n2)
                        v3.append(n3)
    
                        v1.append(n1)
                        v2.append(n3)
                        v3.append(n4)
            elif('shell3' in et['type']):
                for el in et['connectivity']:
                    eli = el[0]
                    if(eli in elSet):
                        n1 = el[1]
                        n2 = el[2]
                        n3 = el[3]
                        
                        v1.append(n1)
                        v2.append(n2)
                        v3.append(n3)
                    
        return {'v1': v1, 'v2': v2, 'v3': v3}
    
    def getFaceValues(self,elSet,elVal,massElOptns=None):
        fcVal = list()
        for et in self.modelData['elements']:
            eTp = et['type']
            for el in et['connectivity']:
                eli = el[0]
                if(eli in elSet):
                    si = elVal[eli]
                    if('brick' in eTp or ('mass' in eTp and massElOptns['showAsDots'])):
                        fcVal.extend([si,si,si,si,si,si,si,si,si,si,si,si])
                    elif('wedge' in eTp):
                        fcVal.extend([si,si,si,si,si,si,si,si])
                    elif('tet' in eTp):
                        fcVal.extend([si,si,si,si])
                    elif('shell4' in eTp):
                        fcVal.extend([si,si])
                    elif('shell3' in eTp):
                        fcVal.extend([si])
                        
        return fcVal
    
    def plotElementProperty(self,prop='section',elementSet='all',nodeSet='all',massElOptns=None):
        elSet = self.getPlotNdElSet(elementSet)
        ndCrd = self.buildNodalPlotCrd(elSet,massElOptns=massElOptns)
        verts = self.buildElementVertexList(elSet,massElOptns=massElOptns)
        
        if(prop == 'section'):
            numEls = 0
            for et in self.modelData['elements']:
                numEls = numEls + len(et['connectivity'])
            elVal = np.zeros(numEls,dtype=int)
            
            elSets = self.modelData['sets']['element']
            
            for si, sec in enumerate(self.modelData['sections']):
                setNm = sec['elementSet']
                for eli in elSets[setNm]:
                    elVal[eli] = si
        
        fcVals = self.getFaceValues(elSet,elVal)
        if nodeSet == 'all':
            ndCrd, verts, ndVals = removeUnusedNodes(ndCrd, verts)
        else:
            ndCrd, verts, fcVals = reduceToNodeSet(ndCrd, verts, set(self.modelData['sets']['node'][nodeSet]), fcVals=fcVals)
        cbTitle = prop
        plotMeshSolution(ndCrd,fcVals,verts,valMode='cell',title=cbTitle)

    def plotNodeResults(self,field,component=1,elementSet='all',nodeSet='all',deformed=False,defScaleFact=1.0,massElOptns=None):
        if massElOptns == None:
            massElOptns = {'showAsDots': False}
        if 'showAsDots' not in massElOptns:
            massElOptns['showAsDots'] = False
        
        elSet = self.getPlotNdElSet(elementSet)
        ndCrd = self.buildNodalPlotCrd(elSet,deformed,defScaleFact,massElOptns)
        numNds = len(self.modelData['nodes'])
        ndRng = list(range(0,numNds))
        
        abrv = {'displacement': ['U1','U2','U3','R1','R2','R3'],
                'velocity': ['V1','V2','V3','RV1','RV2','RV3'],
                'acceleration': ['A1','A2','A3','RA1','RA2','RA3'],
                'temperature': ['T'],
                'tdot': ['TDOT'],
                'reactionForce': ['RF1','RF2','RF3','RM1','RM2','RM3'],
                'reactionHeatGen': ['RHG']}
        
        if (component == 'mag'):
            cols = abrv[field][0:3]
            uAr = np.array(self.nodeData.loc[ndRng, cols])
            values = list()
            for u in uAr:
                values.append(np.linalg.norm(u))
                
            if massElOptns['showAsDots']:
                for et in self.modelData['elements']:
                    if et['type'] == 'mass':
                        for el in et['connectivity']:
                            if el[0] in elSet:
                                u = uAr[el[1]]
                                umag = np.linalg.norm(u)
                                for i in range(0,8):
                                    values.append(umag)
        else:
            fldLab = abrv[field][component-1]
            values = list(self.nodeData.loc[ndRng, fldLab])
            
            if massElOptns['showAsDots']:
                for et in self.modelData['elements']:
                    if et['type'] == 'mass':
                        for el in et['connectivity']:
                            if el[0] in elSet:
                                v = values[el[1]]
                                for i in range(0,8):
                                    values.append(v)
        
        verts = self.buildElementVertexList(elSet,massElOptns)
        if nodeSet == 'all':
            ndCrd, verts, values = removeUnusedNodes(ndCrd, verts, ndVals=values)
        else:
            ndCrd, verts, values = reduceToNodeSet(ndCrd, verts, set(self.modelData['sets']['node'][nodeSet]), ndVals=values)
        cbTitle = field + str(component)
        plotMeshSolution(ndCrd,values,verts,valMode='vertex',title=cbTitle)
        
    def plotElementResults(self,field,component=1,elementSet='all',nodeSet='all',layer=0,deformed=False,defScaleFact=1.0,massElOptns=None):
        if massElOptns == None:
            massElOptns = {'showAsDots': False}
        if 'showAsDots' not in massElOptns:
            massElOptns['showAsDots'] = False
        
        elSet = self.getPlotNdElSet(elementSet)
        
        abrv = {'stress': ['S11','S22','S33','S12','S13','S23','MISES','PS1','PS2','PS3'],
                'strain': ['E11','E22','E33','E12','E13','E23','PE1','PE2','PE3'],
                'strainEnergDen': ['SE'],
                'tsaiWu': ['TSAIWU'],
                'sectionFrcMom': ['SECT_F1','SECT_F2','SECT_F3','SECT_M1','SECT_M2','SECT_M3'],
                'sectionDef': ['SECT_E1','SECT_E2','SECT_E3','SECT_K1','SECT_K2','SECT_K3'],
                'heatFlux': ['HFLX1','HFLX2','HFLX3'],
                'tempGradient': ['TGRAD1','TGRAD2','TGRAD3']}
        
        allColLabs = set()
        for f in abrv:
            allColLabs = allColLabs.union(set(abrv[f]))
        
        try:
            fldLab = abrv[field][component-1]
        except:
            if(str(component) in allColLabs):
                fldLab = component
            else:
                print('Error: unrecognized element result component ' + str(component) + 'plotElementResults() failed')
                return
                                        
        ndCrd = self.buildNodalPlotCrd(elSet,deformed,defScaleFact,massElOptns)
        numNds = len(self.modelData['nodes'])
        
        numEls = 0
        for et in self.modelData['elements']:
            numEls = numEls + len(et['connectivity'])
        
        df1 = self.elementData[self.elementData['int_pt'] == 0]
        df2 = df1[df1['layer'] == layer]
        elValues = np.zeros(numEls,dtype=float)
        for r, ei in enumerate(df2['element']):
            elValues[ei] = df2.loc[r,fldLab]
        
        fcVals = self.getFaceValues(elSet,elValues,massElOptns)
        verts = self.buildElementVertexList(elSet,massElOptns)
        if nodeSet == 'all':
            ndCrd, verts, ndVals = removeUnusedNodes(ndCrd, verts)
        else:
            ndCrd, verts, fcVals = reduceToNodeSet(ndCrd, verts, set(self.modelData['sets']['node'][nodeSet]), fcVals=fcVals)
        cbTitle = field + str(component)
        plotMeshSolution(ndCrd,fcVals,verts,valMode='cell',title=cbTitle)
        
    def plotModalResults(self,elementSet='all',defScaleFact=1.0,massElOptns=None):
        nodeCopy = self.nodeData.copy()
        self.nodeData = self.modalVec
        self.plotNodeResults('displacement',component='mag',elementSet=elementSet,deformed=True,defScaleFact=defScaleFact,massElOptns=massElOptns)
        self.nodeData = nodeCopy
        
    def animateNodeResults(self,fileName,field,timeSteps,component=1,elementSet='all',nodeSet='all',deformed=False,defScaleFact=1.0,massElOptns=None):
        if massElOptns == None:
            massElOptns = {'showAsDots': False}
        if 'showAsDots' not in massElOptns:
            massElOptns['showAsDots'] = False
            
        elSet = self.getPlotNdElSet(elementSet)
        
        abrv = {'displacement': ['U1','U2','U3','R1','R2','R3'],
                'velocity': ['V1','V2','V3','RV1','RV2','RV3'],
                'acceleration': ['A1','A2','A3','RA1','RA2','RA3'],
                'temperature': ['T'],
                'tdot': ['TDOT'],
                'reactionForce': ['RF1','RF2','RF3','RM1','RM2','RM3'],
                'reactionHeatGen': ['RHG']}
        fldLab = abrv[field][component-1]
        
        allNdCrd = list()
        allNdValues = list()
        numNds = len(self.modelData['nodes'])
        ndRng = list(range(0,numNds))
        fnLst = fileName.split('.')
        verts = self.buildElementVertexList(elSet,massElOptns)
        valAr = np.zeros(numNds,dtype=float)
        firstStep = True
        for ts in timeSteps:
            print('animate time step: ' + str(ts))
            fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
            self.loadNodeResults(fn)
            ndCrd = self.buildNodalPlotCrd(elSet,deformed,defScaleFact,massElOptns)
            if (component == 'mag'):
                cols = abrv[field][0:3]
                uAr = np.array(self.nodeData.loc[ndRng, cols])
                ndValues = list()
                for u in uAr:
                    ndValues.append(np.linalg.norm(u))
                    
                if massElOptns['showAsDots']:
                    for et in self.modelData['elements']:
                        if et['type'] == 'mass':
                            for el in et['connectivity']:
                                if el[0] in elSet:
                                    u = uAr[el[1]]
                                    umag = np.linalg.norm(u)
                                    for i in range(0,8):
                                        ndValues.append(umag)
            else:
                ndValues = list(self.nodeData.loc[ndRng, fldLab])
                
                if massElOptns['showAsDots']:
                    for et in self.modelData['elements']:
                        if et['type'] == 'mass':
                            for el in et['connectivity']:
                                if el[0] in elSet:
                                    v = ndValues[el[1]]
                                    for i in range(0,8):
                                        ndValues.append(v)
            
            if nodeSet == 'all':
                ndCrd, redVerts, ndValues = removeUnusedNodes(ndCrd, verts, ndVals=ndValues)
            else:
                ndCrd, redVerts, ndValues = reduceToNodeSet(ndCrd, verts, set(self.modelData['sets']['node'][nodeSet]), ndVals=ndValues)
            allNdCrd.append(ndCrd)
            allNdValues.append(ndValues)
            if(firstStep):
                allNdCrd.append(ndCrd)
                allNdValues.append(ndValues)
                firstStep = False
        cbTitle = field + str(component)
        animateMeshSolution(allNdCrd,allNdValues,redVerts,'vertex',title=cbTitle)
        
    def animateElementResults(self,fileName,field,timeSteps,component=1,elementSet='all',nodeSet='all',layer=0,deformed=False,defScaleFact=1.0,nodeResFile=None,massElOptns=None):
        if massElOptns == None:
            massElOptns = {'showAsDots': False}
        if 'showAsDots' not in massElOptns:
            massElOptns['showAsDots'] = False
        
        elSet = self.getPlotNdElSet(elementSet)
        
        abrv = {'stress': ['S11','S22','S33','S12','S13','S23','MISES','PS1','PS2','PS3'],
                'strain': ['E11','E22','E33','E12','E13','E23','PE1','PE2','PE3'],
                'strainEnergDen': ['SE'],
                'tsaiWu': ['TSAIWU'],
                'sectionFrcMom': ['SECT_F1','SECT_F2','SECT_F3','SECT_M1','SECT_M2','SECT_M3'],
                'sectionDef': ['SECT_E1','SECT_E2','SECT_E3','SECT_K1','SECT_K2','SECT_K3'],
                'heatFlux': ['HFLX1','HFLX2','HFLX3'],
                'tempGradient': ['TGRAD1','TGRAD2','TGRAD3']}
        
        allColLabs = set()
        for f in abrv:
            allColLabs = allColLabs.union(set(abrv[f]))
        
        try:
            fldLab = abrv[field][component-1]
        except:
            if(str(component) in allColLabs):
                fldLab = component
            else:
                print('Error: unrecognized element result component ' + str(component) + 'animateElementResults() failed')
                return
        
        allNdCrd = list()
        allFcValues = list()
        numNds = len(self.modelData['nodes'])
        fnLst = fileName.split('.')
        if(nodeResFile != None):
            ndFnLst = nodeResFile.split('.')
        verts = self.buildElementVertexList(elSet,massElOptns)
        
        numEls = 0
        for et in self.modelData['elements']:
            numEls = numEls + len(et['connectivity'])
        elValues = np.zeros(numEls,dtype=float)
        
        for ts in timeSteps:
            if(deformed and nodeResFile != None):
                fn = ndFnLst[0] + '_timestep' + str(ts) + '.' + ndFnLst[1]
                self.loadNodeResults(fn)
            ndCrd = self.buildNodalPlotCrd(ndSet,deformed,defScaleFact,massElOptns)
            fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
            self.loadElementResults(fn)
            
            df1 = self.elementData[self.elementData['int_pt'] == 0]
            df2 = df1[df1['layer'] == layer]
            elValues = np.zeros(numEls,dtype=float)
            for r, ei in enumerate(df2['elements']):
                elValues[ei] = df2.loc[r,fldLab]
            
            fcVals = self.getFaceValues(elSet,elValues,massElOptns)
            if nodeSet == 'all':
                ndCrd, redVerts, ndVals = removeUnusedNodes(ndCrd, verts)
            else:
                ndCrd, redVerts, fcVals = reduceToNodeSet(ndCrd, verts, set(self.modelData['sets']['node'][nodeSet]), fcVals=fcVals)
            allNdCrd.append(ndCrd)
            allFcValues.append(fcVals)
        cbTitle = field + str(component)
        animateMeshSolution(allNdCrd,allFcValues,redVerts,'cell',title=cbTitle)
        
    def animateModalSolution(self,elementSet='all',nodeSet='all',defScaleFact=1.0,massElOptns=None):
        if massElOptns == None:
            massElOptns = {'showAsDots': False}
        if 'showAsDots' not in massElOptns:
            massElOptns['showAsDots'] = False
        
        elSet = self.getPlotNdElSet(elementSet)
        nodeCopy = self.nodeData.copy()
        self.nodeData = self.modalVec
        
        allNdCrd = list()
        allNdValues = list()
        numNds = len(self.modelData['nodes'])
        valAr = np.zeros(numNds,dtype=float)
        verts = self.buildElementVertexList(elSet, massElOptns)
        
        for theta in range(0,360,30):
            tRad = 0.0174533*theta
            sinTh = np.math.sin(tRad)
            sf = sinTh*defScaleFact
            ndCrd = self.buildNodalPlotCrd(ndSet,deformed=True,defScaleFact=sf,massElOptns=massElOptns)
            
            cols = ['U1','U2','U3']
            uAr = np.array(self.nodeData.loc[list(range(0, numNds)), cols])
            ndValues = list()
            for u in uAr:
                ndValues.append(np.linalg.norm(u))
                
            if massElOptns['showAsDots']:
                for et in self.modelData['elements']:
                    if et['type'] == 'mass':
                        for el in et['connectivity']:
                            if el[0] in elSet:
                                u = uAr[el[1]]
                                umag = np.linalg.norm(u)
                                for i in range(0,8):
                                    ndValues.append(umag)    
            
            if nodeSet == 'all':
                ndCrd, redVerts, ndValues = removeUnusedNodes(ndCrd, verts, ndValues)
            else:
                ndCrd, redVerts, ndValues = reduceToNodeSet(ndCrd, verts, set(self.modelData['sets']['node'][nodeSet]), ndVals=ndValues)
            allNdCrd.append(ndCrd)
            allNdValues.append(ndValues)
        cbTitle = 'displacement'
        animateMeshSolution(allNdCrd,allNdValues,redVerts,title=cbTitle)
        self.nodeData = nodeCopy
        
    def extractNodeHistory(self,fileName,field,timeSteps,nodeSet):
        fnLst = fileName.split('.')
        
        abrv = {'displacement': ['U1','U2','U3','R1','R2','R3'],
                'velocity': ['V1','V2','V3','RV1','RV2','RV3'],
                'acceleration': ['A1','A2','A3','RA1','RA2','RA3'],
                'temperature': ['T'],
                'tdot': ['TDOT'],
                'reactionForce': ['RF1','RF2','RF3','RM1','RM2','RM3'],
                'reactionHeatGen': ['RHG']}
        
        rescols = abrv[field]
        
        try:
            ndI = int(nodeSet)
            vals = list()
            timePts = list()
            for ts in timeSteps:
                fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
                self.loadNodeResults(fn)
                resrow = self.nodeData.loc[ndI]
                vals.append(list(resrow[rescols]))
                timepts.append(resrow['time'])
            series = dict()
            lab = 'node_' + str(ndI)
            series[lab] = vals
        except:
            series = dict()
            nsLabs = self.modelData['sets']['node'][nodeSet]
            for nd in nsLabs:
                lab = 'node_' + str(nd)
                series[lab] = list()
            timePts = list()
            for ts in timeSteps:
                fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
                self.loadNodeResults(fn)
                timePts.append(self.nodeData.loc[0,'time'])
                for nd in nsLabs:
                    lab = 'node_' + str(nd)
                    resrow = self.nodeData.loc[nd, rescols]
                    series[lab].append(list(resrow))
        return series, timePts
        
    def nodeHistorySeries(self,fileName,field,timeSteps,nodeSet,component=1):
        series, timePts = self.extractNodeHistory(fileName,field,timeSteps,nodeSet)
        newSeries = dict()
        for k in series:
            sList = list()
            for v in series[k]:
                if(component == 'mag'):
                    vec = np.array(v[0:3])
                    val = np.linalg.norm(vec)
                else:
                    val = v[component]
                sList.append(val)
            newSeries[k] = sList
        return newSeries, timePts
        
    def plotNodeHistory(self,fileName,field,timeSteps,nodeSet,component=1,xTitle='Time',yTitle=None):
        series, timePts = self.nodeHistorySeries(fileName,field,timeSteps,nodeSet,component)
        if(yTitle == None):
            ytitle = field + str(component)
        else:
            ytitle = yTitle
        plotTimeHistory(series,timePts,xTitle=xTitle,yTitle=ytitle)
        
    def plotElementHistory(self,fileName,field,timeSteps,elementSet,layer=0,component=1,xTitle='Time',yTitle=None):
        fnLst = fileName.split('.')
        
        abrv = {'stress': ['S11','S22','S33','S12','S13','S23','MISES','PS1','PS2','PS3'],
                'strain': ['E11','E22','E33','E12','E13','E23','PE1','PE2','PE3'],
                'strainEnergDen': ['SE'],
                'tsaiWu': ['TSAIWU'],
                'sectionFrcMom': ['SECT_F1','SECT_F2','SECT_F3','SECT_M1','SECT_M2','SECT_M3'],
                'sectionDef': ['SECT_E1','SECT_E2','SECT_E3','SECT_K1','SECT_K2','SECT_K3'],
                'heatFlux': ['HFLX1','HFLX2','HFLX3'],
                'tempGradient': ['TGRAD1','TGRAD2','TGRAD3']}
        
        allColLabs = set()
        for f in abrv:
            allColLabs = allColLabs.union(set(abrv[f]))
        
        try:
            rescol = abrv[field][component-1]
        except:
            if(str(component) in allColLabs):
                rescol = component
            else:
                print('Error: unrecognized element result component ' + str(component) + 'plotElementHistory() failed')
                return
        
        try:
            elI = int(elementSet)
            vals = list()
            timePts = list()
            for ts in timeSteps:
                fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
                self.loadElementResults(fn)
                
                df1 = self.elementData[self.elementData['element'] == elI]
                df2 = df1[df1['int_pt'] == 0]
                resrow = df2[df2['layer'] == layer]
                vals.append(resrow[rescol])
                timePts.append(resrow['time'])
            series = dict()
            lab = 'element_' + str(elI)
            series[lab] = vals
            if(yTitle == None):
                ytitle = field + str(component)
            else:
                ytitle = yTitle
            plotTimeHistory(series,timePts,xTitle=xTitle,yTitle=yTitle)
        except:
            series = dict()
            esLabs = self.modelData['sets']['element'][elementSet]
            for el in esLabs:
                lab = 'element_' + str(el)
                series[lab] = list()
            timePts = list()
            for ts in timeSteps:
                fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
                self.loadElementResults(fn)
                timePts.append(self.elementData.loc[0, 'time'])
                
                df1 = self.elementData[self.elementData['element'].isin(esLabs)]
                df2 = df1[df1['int_pt'] == 0]
                df3 = df2[df2['layer'] == layer]
                for r, ei in enumerate(df3['element']):
                    lab = 'element_' + str(ei)
                    series[lab].append(df3.loc[r,rescol])
            if(yTitle == None):
                ytitle = field + str(component)
            else:
                ytitle = yTitle
            plotTimeHistory(series,timePts,xTitle=xTitle,yTitle=ytitle)
            
    def plotElementSensitivity(self,elementSet='all',nodeSet='all',dVarSet='all',magnitude=False):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
        if(dVarSet == 'all'):
            lenD = len(self.dVarData['designVariables'])
            dSet = set(range(0,lenD))
        else:
            dSet = set(dVarSet)
        
        numEls = 0
        for et in self.modelData['elements']:
            numEls = numEls + len(et['connectivity'])
        elValues = np.zeros(numEls,dtype=float)
        
        modSets = self.modelData['sets']['element']
        #objGrad = self.objectiveData['objectiveGradient']
        objGrad = list(self.objectiveGrad['objGrad'])
        for di, d in enumerate(self.dVarData['designVariables']):
            if(di in dSet):
                try:
                    dES = d['elementSet']
                    try:
                        eli = int(dES)
                        if(magnitude):
                            elValues[eli] = abs(objGrad[di])
                        else:
                            elValues[eli] = objGrad[di]
                    except:
                        # for eS in modSets:
                        #     if(eS['name'] == dES):
                        for el in modSets[dES]:
                            if(magnitude):
                                elValues[el] = abs(objGrad[di])
                            else:
                                elValues[el] = objGrad[di]
                except:
                    pass
        
        ndCrd = self.buildNodalPlotCrd(ndSet)
        fcVals = self.getFaceValues(elSet,elValues)
        verts = self.buildElementVertexList(elSet)
        if nodeSet == 'all':
            ndCrd, verts, ndVals = removeUnusedNodes(ndCrd, verts)
        else:
            ndCrd, verts, fcVals = reduceToNodeSet(ndCrd, verts, set(self.modelData['sets']['node'][nodeSet]), fcVals=fcVals)
        plotMeshSolution(ndCrd,fcVals,verts,valMode='cell')
        
    def extractModalAmplitudes(self,nodeResFile,modalResFile,timeSteps,nodeSet,modeList):
        try:
            nSet = {int(nodeSet)}
        except:
            nSet = set(self.modelData['sets']['node'][nodeSet])
        
        normMdDisp = dict()
        for mi in modeList:
            self.loadModalVec(modalResFile,mi)
            uAr = np.array(self.modalVec[['U1','U2','U3']])
            mdStr = str(mi)
            maxMag = 0.0
            for u in uAr:
                mag = np.linalg.norm(u)
                if(mag > maxMag):
                    maxMag = mag
            mFact = 1.0/maxMag
            for ni in nSet:
                nu = np.array(self.modalVec.loc[ni, ['U1','U2','U3']])
                dk = str(ni) + ',' + mdStr
                normMdDisp[dk] = mFact*nu
        
        self.loadModalVals(modalResFile)
        freq = list(self.modalVals['frequency'])
        
                
        series, timePts = self.extractNodeHistory(fileName,'displacement',timeSteps,nodeSet)
        outDat = dict()
        nRows = 3*len(timeSteps)
        nCols = 2*len(modeList) + 1
        mat = np.zeros((nRows,nCols),dtype=float)
        uVec = np.zeros(nRows,dtype=float)
        for sk in series:
            kLst = sk.split('_')
            nStr = kLst[1]
            i = 0
            for j, nu in enumerate(series[sk]):
                uVec[i:i+3] = np.array(nu[0:3])
                t = timePts[j]
                i2 = 0
                for mi in modeList:
                    mStr = str(mi)
                    dk = nStr + ',' + mStr
                    omegaT = 2.0*np.pi*freq[mi]*t
                    mat[i:i+3,i2] = normMdDisp[dk]*np.sin(omegaT)
                    mat[i:i+3,i2+1] = normMdDisp[dk]*np.cos(omegaT)
                    i2 += 2
                i += 3
            mat[:,nCols-1] = 1.0
            Q, R = np.linalg.qr(mat)
            rhs = np.matmul(uVec,Q)
            soln = np.linalg.solve(R,rhs)
            ndDic = dict()
            i2 = 0
            for mi in modeList:
                vec = soln[i2:i2+2]
                amp = np.linalg.norm(vec)
                dk = 'mode_' + str(mi)
                ndDic[dk] = amp
                i2 += 2
            outDat[sk] = ndDic
            
        return outDat, freq
           
    def extractNodeFrequencies(self,fileName,field,timeSteps,nodeSet,freq,component=1):
        series, timePts = self.nodeHistorySeries(fileName,field,timeSteps,nodeSet,component)
        npts = len(timePts)
        tAr = np.array(timePts)
        pi2 = 2.0*np.pi
        ons = np.ones(npts,dtype=float)
        seriesAmp = dict()
        for sLab in series:
            vAr = np.array(series[sLab])
            seriesAmp[sLab] = list()
            mT = [ons]
            for f in freq:
                omega = pi2*f
                sVec = np.sin(omega*tAr)
                cVec = np.cos(omega*tAr)
                mT.append(sVec)
                mT.append(cVec)
            mat = np.transpose(mT)
            Q, R = np.linalg.qr(mat)
            rhs = np.matmul(vAr,Q)
            soln = np.linalg.solve(R,rhs)
            j = 1
            for f in freq:
                amp = np.linalg.norm(soln[j:j+2])
                seriesAmp[sLab].append(amp)
                j += 2
        return seriesAmp
           
    def plotNodeFrequencies(self,fileName,field,timeSteps,nodeSet,freq,component=1,xTitle='Frequency',yTitle='Amplitude'):
        seriesAmp = self.extractNodeFrequencies(fileName, field, timeSteps, nodeSet, freq, component)
        plotFrequencySpectrum(seriesAmp,freq,xTitle=xTitle,yTitle=yTitle)