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
        # inFile = open(nodeResFile,'r')
        # self.nodeData = yaml.load(inFile,Loader=Loader)
        # inFile.close()

    def loadElementResults(self,elementResFile):
        self.elementData = pd.read_csv(elementResFile)
        # inFile = open(elementResFile,'r')
        # self.elementData = yaml.load(inFile,Loader=Loader)
        # inFile.close()
        
    def loadModalVals(self,modalResFile):
        self.modalVals = pd.read_csv(modalResFile)
        # inFile = open(modalResFile,'r')
        # self.modalData = yaml.load(inFile,Loader=Loader)
        # inFile.close()
        
    def loadModalVec(self,modalResFile,mode=0):
        rep = '_mode' + str(mode) + '.'
        f_name = modalResFile.replace('.',rep)
        self.modalVec = pd.read_csv(f_name,index_col=0)
        
    def loadObjectiveResults(self,objResFile):
        self.objectiveVals = pd.read_csv(objResFile,index_col=0)
        f_name = objResFile.replace('.','_grad.')
        self.objectiveGrad = pd.read_csv(f_name,index_col=0)
        # inFile = open(objResFile,'r')
        # self.objectiveData = yaml.load(inFile,Loader=Loader)
        # inFile.close()
        
    def getPlotNdElSet(self,elementSet):
        numNds = len(self.modelData['nodes'])
        ndSet = set(range(0,numNds))
        if(elementSet == 'all'):
            #ndSet = set(range(0,numNds))
            elSet = set()
            for et in self.modelData['elements']:
                for el in et['connectivity']:
                    elSet.add(el[0])
        else:
            for es in self.modelData['sets']['element']:
                if(es['name'] == elementSet):
                    elSet = set(es['labels'])
                    # ndSet = set()
                    # for et in self.modelData['elements']:
                    #     for el in et['connectivity']:
                    #         lab = el[0]
                    #         if(lab in elSet):
                    #             for ni, nd in enumerate(el):
                    #                 if(ni > 0):
                    #                     ndSet.add(nd)
        return ndSet, elSet
        
    def buildNodalPlotCrd(self,ndSet,deformed=False,defScaleFact=1.0):
        ## ndSet a python set() of the desired node labels
        allNds = self.modelData['nodes']
        numNds = len(allNds)
        crdAr = np.zeros((numNds,3),dtype=float)
        for nd in allNds:
            lab = nd[0]
            if(lab in ndSet):
                crdAr[lab] = np.array(nd[1:4])
        if(deformed):
            crdAr +=  defScaleFact*np.array(self.nodeData.loc[list(range(0,numNds)), ['U1','U2','U3']])
        setCrd = crdAr[list(ndSet)]
        xLst = list(crdAr[:,0])
        yLst = list(crdAr[:,1])
        zLst = list(crdAr[:,2])
            
        return {'xLst': xLst, 'yLst': yLst, 'zLst': zLst}

    def buildElementVertexList(self,elSet):
        # elSet = python set() with the desired element labels
        v1 = []
        v2 = []
        v3 = []
        for et in self.modelData['elements']:
            if('brick' in et['type']):
                for el in et['connectivity']:
                    eli = el[0]
                    if(eli in elSet):
                        n1 = el[1]
                        n2 = el[2]
                        n3 = el[3]
                        n4 = el[4]
                        n5 = el[5]
                        n6 = el[6]
                        n7 = el[7]
                        n8 = el[8]
                        
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
    
    def getFaceValues(self,elSet,elVal):
        fcVal = list()
        for et in self.modelData['elements']:
            eTp = et['type']
            for el in et['connectivity']:
                eli = el[0]
                if(eli in elSet):
                    si = elVal[eli]
                    if('brick' in eTp):
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
    
    def plotElementProperty(self,prop='section',elementSet='all'):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
        ndCrd = self.buildNodalPlotCrd(ndSet)
        verts = self.buildElementVertexList(elSet)
        
        if(prop == 'section'):
            numEls = 0
            for et in self.modelData['elements']:
                numEls = numEls + len(et['connectivity'])
            elVal = np.zeros(numEls,dtype=int)
            
            elSets = self.modelData['sets']['element']
            setDic = dict()
            for si, es in enumerate(elSets):
                setDic[es['name']] = si
            
            for si, sec in enumerate(self.modelData['sections']):
                setNm = sec['elementSet']
                seti = setDic[setNm]
                for eli in elSets[seti]['labels']:
                    elVal[eli] = si
        
        fcVals = self.getFaceValues(elSet,elVal)
        cbTitle = prop
        plotMeshSolution(ndCrd,fcVals,verts,valMode='cell',title=cbTitle)

    def plotNodeResults(self,field,component=1,elementSet='all',deformed=False,defScaleFact=1.0):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
        ndCrd = self.buildNodalPlotCrd(ndSet,deformed,defScaleFact)
        numNds = len(self.modelData['nodes'])
        
        abrv = {'displacement': ['U1','U2','U3','R1','R2','R3'],
                'velocity': ['V1','V2','V3','RV1','RV2','RV3'],
                'acceleration': ['A1','A2','A3','RA1','RA2','RA3'],
                'temperature': ['T'],
                'tdot': ['TDOT'],
                'reactionForce': ['RF1','RF2','RF3','RM1','RM2','RM3'],
                'reactionHeatGen': ['RHG']}
        
        if (component == 'mag'):
            cols = abrv[field][0:3]
            uAr = np.array(self.nodeData.loc[list(ndSet), cols])
            values = list()
            for u in uAr:
                values.append(np.linalg.norm(u))
        else:
            fldLab = abrv[field][component-1]
            values = list(self.nodeData.loc[list(ndSet), fldLab])
        
        # valAr = np.zeros(numNds,dtype=float)
        # for nd in self.nodeData['nodeResults'][field]:
        #     lab = nd[0]
        #     if(component == 'mag'):
        #         vec = np.array(nd[1:4])
        #         val = np.linalg.norm(vec)
        #     else:
        #         val = nd[component]
        #     valAr[lab] = val
        # values = list()
        # for nd in self.modelData['nodes']:
        #     lab = nd[0]
        #     if(lab in ndSet):
        #         values.append(valAr[lab])
        
        verts = self.buildElementVertexList(elSet)
        cbTitle = field + str(component)
        plotMeshSolution(ndCrd,values,verts,valMode='vertex',title=cbTitle)
        
    def plotElementResults(self,field,component=1,elementSet='all',layer=0,deformed=False,defScaleFact=1.0):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
        
        abrv = {'stress': ['S11','S22','S33','S12','S13','S23'],
                'strain': ['E11','E22','E33','E12','E13','E23'],
                'strainEnergDen': ['SE'],
                'sectionFrcMom': ['SECT_F1','SECT_F2','SECT_F3','SECT_M1','SECT_M2','SECT_M3'],
                'sectionDef': ['SECT_E1','SECT_E2','SECT_E3','SECT_K1','SECT_K2','SECT_K3'],
                'heatFlux': ['HFLX1','HFLX2','HFLX3'],
                'tempGradient': ['TGRAD1','TGRAD2','TGRAD3']}
        
        fldLab = abrv[field][component-1]
                                        
        ndCrd = self.buildNodalPlotCrd(ndSet,deformed,defScaleFact)
        numNds = len(self.modelData['nodes'])
        
        numEls = 0
        for et in self.modelData['elements']:
            numEls = numEls + len(et['connectivity'])
        
        df1 = self.elementData[self.elementData['int_pt'] == 0]
        df2 = df1[df1['layer'] == layer]
        elValues = np.zeros(numEls,dtype=float)
        for r, ei in enumerate(df2['element']):
            elValues[ei] = df2.loc[r,fldLab]
        
        # elValues = np.zeros(numEls,dtype=float)
        # for el in self.elementData['elementResults'][field]:
        #     if(el[2] == layer):
        #         lab = el[0]
        #         i = component + 2
        #         val = el[i]
        #         elValues[lab] = val
        
        fcVals = self.getFaceValues(elSet,elValues)
        verts = self.buildElementVertexList(elSet)
        cbTitle = field + str(component)
        plotMeshSolution(ndCrd,fcVals,verts,valMode='cell',title=cbTitle)
        
    def plotModalResults(self,elementSet='all',defScaleFact=1.0):
        nodeCopy = self.nodeData.copy()
        # for md in self.modalData['modalResults']['modes']:
        #     if(md['mode'] == mode):
        #         self.nodeData['nodeResults'] = dict()
        #         self.nodeData['nodeResults']['displacement'] = md['displacement']
        self.nodeData = self.modalVec
        self.plotNodeResults('displacement',component='mag',elementSet=elementSet,deformed=True,defScaleFact=defScaleFact)
        self.nodeData = nodeCopy
        
    def animateNodeResults(self,fileName,field,timeSteps,component=1,elementSet='all',deformed=False,defScaleFact=1.0,frameDuration=1000):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
        
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
        fnLst = fileName.split('.')
        verts = self.buildElementVertexList(elSet)
        valAr = np.zeros(numNds,dtype=float)
        firstStep = True
        for ts in timeSteps:
            print('animate time step: ' + str(ts))
            fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
            self.loadNodeResults(fn)
            ndCrd = self.buildNodalPlotCrd(ndSet,deformed,defScaleFact)
            if (component == 'mag'):
                cols = abrv[field][0:3]
                uAr = np.array(self.nodeData.loc[list(ndSet), cols])
                ndValues = list()
                for u in uAr:
                    ndValues.append(np.linalg.norm(u))
            else:
                ndValues = list(self.nodeData.loc[list(ndSet), fldLab])
            # for nd in self.nodeData['nodeResults'][field]:
            #     lab = nd[0]
            #     if(component == 'mag'):
            #         vec = np.array(nd[1:4])
            #         val = np.linalg.norm(vec)
            #     else:
            #         val = nd[component]
            #     valAr[lab] = val
            # ndValues = list()
            # for nd in self.modelData['nodes']:
            #     lab = nd[0]
            #     if(lab in ndSet):
            #         ndValues.append(valAr[lab])
            allNdCrd.append(ndCrd)
            allNdValues.append(ndValues)
            if(firstStep):
                allNdCrd.append(ndCrd)
                allNdValues.append(ndValues)
                firstStep = False
        cbTitle = field + str(component)
        animateMeshSolution(allNdCrd,allNdValues,verts,'vertex',frameDuration,cbTitle)
        
    def animateElementResults(self,fileName,field,timeSteps,component=1,elementSet='all',layer=0,deformed=False,defScaleFact=1.0,nodeResFile=None,frameDuration=1000):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
        
        abrv = {'stress': ['S11','S22','S33','S12','S13','S23'],
                'strain': ['E11','E22','E33','E12','E13','E23'],
                'strainEnergDen': ['SE'],
                'sectionFrcMom': ['SECT_F1','SECT_F2','SECT_F3','SECT_M1','SECT_M2','SECT_M3'],
                'sectionDef': ['SECT_E1','SECT_E2','SECT_E3','SECT_K1','SECT_K2','SECT_K3'],
                'heatFlux': ['HFLX1','HFLX2','HFLX3'],
                'tempGradient': ['TGRAD1','TGRAD2','TGRAD3']}
        fldLab = abrv[field][component-1]
        
        allNdCrd = list()
        allFcValues = list()
        numNds = len(self.modelData['nodes'])
        fnLst = fileName.split('.')
        if(nodeResFile != None):
            ndFnLst = nodeResFile.split('.')
        verts = self.buildElementVertexList(elSet)
        
        numEls = 0
        for et in self.modelData['elements']:
            numEls = numEls + len(et['connectivity'])
        elValues = np.zeros(numEls,dtype=float)
        
        for ts in timeSteps:
            if(deformed and nodeResFile != None):
                fn = ndFnLst[0] + '_timestep' + str(ts) + '.' + ndFnLst[1]
                self.loadNodeResults(fn)
            ndCrd = self.buildNodalPlotCrd(ndSet,deformed,defScaleFact)
            fn = fnLst[0] + '_timestep' + str(ts) + '.' + fnLst[1]
            self.loadElementResults(fn)
            
            df1 = self.elementData[self.elementData['int_pt'] == 0]
            df2 = df1[df1['layer'] == layer]
            elValues = np.zeros(numEls,dtype=float)
            for r, ei in enumerate(df2['elements']):
                elValues[ei] = df2.loc[r,fldLab]
            
            # elValues = np.zeros(numEls,dtype=float)
            # for el in self.elementData['elementResults'][field]:
            #     if(el[2] == layer):
            #         lab = el[0]
            #         i = component + 2
            #         val = el[i]
            #         elValues[lab] = val
            
            fcVals = self.getFaceValues(elSet,elValues)
            allNdCrd.append(ndCrd)
            allFcValues.append(fcVals)
        cbTitle = field + str(component)
        animateMeshSolution(allNdCrd,allFcValues,verts,'cell',frameDuration,cbTitle)
        
    def animateModalSolution(self,elementSet='all',defScaleFact=1.0):
        ndSet, elSet = self.getPlotNdElSet(elementSet)
        nodeCopy = self.nodeData.copy()
        # for md in self.modalData['modalResults']['modes']:
        #     if(md['mode'] == mode):
        #         self.nodeData['nodeResults'] = dict()
        #         self.nodeData['nodeResults']['displacement'] = md['displacement']
        self.nodeData = self.modalVec
        
        allNdCrd = list()
        allNdValues = list()
        numNds = len(self.modelData['nodes'])
        valAr = np.zeros(numNds,dtype=float)
        verts = self.buildElementVertexList(elSet)
        
        for theta in range(0,360,30):
            tRad = 0.0174533*theta
            sinTh = np.math.sin(tRad)
            sf = sinTh*defScaleFact
            ndCrd = self.buildNodalPlotCrd(ndSet,deformed=True,defScaleFact=sf)
            
            cols = ['U1','U2','U3']
            uAr = np.array(self.nodeData.loc[list(ndSet), cols])
            ndValues = list()
            for u in uAr:
                ndValues.append(np.linalg.norm(u))
            # for nd in self.nodeData['nodeResults']['displacement']:
            #     lab = nd[0]
            #     if(lab in ndSet):
            #         vec = np.array(nd[1:4])
            #         valAr[lab] = np.linalg.norm(vec)
            # ndValues = list()
            # for nd in self.modelData['nodes']:
            #     lab = nd[0]
            #     if(lab in ndSet):
            #         ndValues.append(valAr[lab])
            allNdCrd.append(ndCrd)
            allNdValues.append(ndValues)
        cbTitle = 'displacement'
        animateMeshSolution(allNdCrd,allNdValues,verts,title=cbTitle)
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
                # for nd in self.nodeData['nodeResults'][field]:
                #     if(nd[0] == ndI):
                #         vals.append(nd[1:])
                #         timePts.append(self.nodeData['nodeResults']['time'])
            series = dict()
            lab = 'node_' + str(ndI)
            series[lab] = vals
        except:
            series = dict()
            nsLabs = []
            for ns in self.modelData['sets']['node']:
                if(ns['name'] == nodeSet):
                    nsLabs = ns['labels']
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
                # for nd in self.nodeData['nodeResults'][field]:
                #     try:
                #         lab = 'node_' + str(nd[0])
                #         series[lab].append(nd[1:])
                #     except:
                #         pass
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
        
        abrv = {'stress': ['S11','S22','S33','S12','S13','S23'],
                'strain': ['E11','E22','E33','E12','E13','E23'],
                'strainEnergDen': ['SE'],
                'sectionFrcMom': ['SECT_F1','SECT_F2','SECT_F3','SECT_M1','SECT_M2','SECT_M3'],
                'sectionDef': ['SECT_E1','SECT_E2','SECT_E3','SECT_K1','SECT_K2','SECT_K3'],
                'heatFlux': ['HFLX1','HFLX2','HFLX3'],
                'tempGradient': ['TGRAD1','TGRAD2','TGRAD3']}
        
        rescol = abrv[field][component-1]
        
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
                # for el in self.elementData['elementResults'][field]:
                #     if(el[0] == elI):
                #         vals.append(el[component+2])
                #         timePts.append(self.elementData['elementResults']['time'])
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
            esLabs = []
            for es in self.modelData['sets']['element']:
                if(es['name'] == elementSet):
                    esLabs = es['labels']
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
                # for el in self.elementData['elementResults'][field]:
                #     try:
                #         lab = 'element_' + str(el[0])
                #         series[lab].append(el[component+2])
                #     except:
                #         pass
            if(yTitle == None):
                ytitle = field + str(component)
            else:
                ytitle = yTitle
            plotTimeHistory(series,timePts,xTitle=xTitle,yTitle=ytitle)
            
    def plotElementSensitivity(self,elementSet='all',dVarSet='all',magnitude=False):
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
                        for eS in modSets:
                            if(eS['name'] == dES):
                                for el in eS['labels']:
                                    if(magnitude):
                                        elValues[el] = abs(objGrad[di])
                                    else:
                                        elValues[el] = objGrad[di]
                except:
                    pass
        
        ndCrd = self.buildNodalPlotCrd(ndSet)
        fcVals = self.getFaceValues(elSet,elValues)
        verts = self.buildElementVertexList(elSet)
        plotMeshSolution(ndCrd,fcVals,verts,valMode='cell')
        
    def extractModalAmplitudes(self,nodeResFile,modalResFile,timeSteps,nodeSet,modeList):
        try:
            nSet = {int(nodeSet)}
        except:
            for ns in self.modelData['sets']['node']:
                if(ns['name'] == nodeSet):
                    nSet = set(ns['labels'])
        
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
        
        # for md in self.modalData['modalResults']['modes']:
        #     mdStr = str(md['mode'])
        #     maxMag = 0.
        #     for nu in md['displacement']:
        #         vec = np.array(nu[1:4])
        #         mag = np.linalg.norm(vec)
        #         if(mag > maxMag):
        #             maxMag = mag
        #     mFact = 1.0/maxMag
        #     for nu in md['displacement']:
        #         lab = nu[0]
        #         if(lab in nSet):
        #             ndStr = str(lab)
        #             dk = ndStr + ',' + mdStr
        #             normMdDisp[dk] = mFact*np.array(nu[1:4])
                    
        # freq = self.modalData['modalResults']['frequencies']
        
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
        
    # def extractNodeFrequencies(self,fileName,field,timeSteps,nodeSet,component=1):
    #     series, timePts = self.nodeHistorySeries(fileName,field,timeSteps,nodeSet,component)
    #     npts = len(timePts)
    #     totalPeriod = timePts[npts-1] - timePts[0]
    #     lowFreq = 1.0/totalPeriod
    #     timeStep = totalPeriod/npts
    #     highFreq = 1.0/(2*timeStep)
    #     nFreq = int(np.ceil(highFreq/lowFreq))
    #     freq = list()
    #     for fi in range(0,nFreq):
    #         freq.append((fi+1)*lowFreq)
    #     pi2 = 2.0*np.pi
    #     tAr = np.array(timePts)
    #     ons = np.ones(npts,dtype=float)
    #     seriesAmp = dict()
    #     for sLab in series:
    #         vAr = np.array(series[sLab])
    #         seriesAmp[sLab] = list()
    #         for fi in range(0,nFreq):
    #             omega = pi2*(fi+1)*lowFreq
    #             sVec = np.sin(omega*tAr)
    #             cVec = np.cos(omega*tAr)
    #             mT = np.array([sVec,cVec,ons])
    #             mat = np.transpose(mT)
    #             Q, R = np.linalg.qr(mat)
    #             rhs = np.matmul(vAr,Q)
    #             soln = np.linalg.solve(R,rhs)
    #             amp = np.linalg.norm(soln[0:2])
    #             seriesAmp[sLab].append(amp)
    #     return seriesAmp, freq
    
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
        
        # npts = len(timePts)
        # totalPeriod = timePts[npts-1] - timePts[0]
        # lowFreq = 1.0/totalPeriod
        # timeStep = totalPeriod/npts
        # highFreq = 1.0/(2*timeStep)
        # nFreq = int(np.ceil(highFreq/lowFreq))
        # freq = list()
        # for fi in range(0,nFreq):
        #     freq.append((fi+1)*lowFreq)
        # pi2 = 2.0*np.pi
        # tAr = np.array(timePts)
        # ons = np.ones(npts,dtype=float)
        # seriesAmp = dict()
        # for sLab in series:
        #     vAr = np.array(series[sLab])
        #     seriesAmp[sLab] = list()
        #     for fi in range(0,nFreq):
        #         omega = pi2*(fi+1)*lowFreq
        #         sVec = np.sin(omega*tAr)
        #         cVec = np.cos(omega*tAr)
        #         mT = np.array([sVec,cVec,ons])
        #         mat = np.transpose(mT)
        #         Q, R = np.linalg.qr(mat)
        #         rhs = np.matmul(vAr,Q)
        #         soln = np.linalg.solve(R,rhs)
        #         amp = np.linalg.norm(soln[0:2])
        #         seriesAmp[sLab].append(amp)
        # return seriesAmp, freq
    
    def plotNodeFrequencies(self,fileName,field,timeSteps,nodeSet,freq,component=1,xTitle='Frequency',yTitle='Amplitude'):
        seriesAmp = self.extractNodeFrequencies(fileName, field, timeSteps, nodeSet, freq, component)
        plotFrequencySpectrum(seriesAmp,freq,xTitle=xTitle,yTitle=yTitle)