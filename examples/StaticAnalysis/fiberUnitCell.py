# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 09:49:49 2024

@author: evans
"""

import os
import sys
from asendUtils.syst.pathTools import *
from asendUtils.model.Model import *
from asendUtils.job.ASenDJob import *
from asendUtils.ResultsProcessor import *

constFile = 'fiberUnitCell/constraints.yaml'
loadFile = 'fiberUnitCell/loads.yaml'
nodeResFile = 'fiberUnitCell/results/nodeResults.yaml'
elResFile = 'fiberUnitCell/results/elementResults.yaml'
jobFile = 'fiberUnitCell/job.yaml'

if(not os.path.exists('fiberUnitCell')):
    os.mkdir('fiberUnitCell')
    
if(not os.path.exists('fiberUnitCell/results')):
    os.mkdir('fiberUnitCell/results')

rtDir = getEnvPath('rootpath')
modFile = rtDir + '/examples/common/fiberUnitCell.yaml'
modScr = rtDir + '/examples/modelGeneration'
if(not os.path.exists(modFile)):
    sys.path.append(modScr)
    import fiberCompositeUnitCell
    
constMod = Model()
constMod.fixDisplacement('xMinRef',ux=0.,uy=0.,uz=0.)
constMod.fixDisplacement('xMaxRef',uy=0.,uz=0.)
constMod.fixDisplacement('yMinRef',uz=0.)
constMod.periodicDisplacement()
constMod.writeModelInput(constFile)

loadMod = Model()
loadMod.addNodalForce('xMaxRef',F=[1000.,0.,0.],M=[0.,0.,0.] ,stTime=0.5, endTime=1.5)

loadMod.addNodalForce('yMinRef',F=[0.,-100.,0.],M=[0.,0.,0.] ,stTime=1.5,endTime=2.5)
loadMod.addNodalForce('yMaxRef',F=[0.,100.,0.],M=[0.,0.,0.] ,stTime=1.5,endTime=2.5)

loadMod.addNodalForce('zMinRef',F=[0.,0.,-100.],M=[0.,0.,0] ,stTime=2.5,endTime=3.5)
loadMod.addNodalForce('zMaxRef',F=[0.,0.,100.],M=[0.,0.,0.] ,stTime=2.5,endTime=3.5)

loadMod.addNodalForce('yMinRef',F=[-100.,0.,0.],M=[0.,0.,0.] ,stTime=3.5,endTime=4.5)
loadMod.addNodalForce('yMaxRef', F=[100.,0.,0.], M=[0.,0.,0.],stTime=3.5,endTime=4.5)

loadMod.addNodalForce('zMinRef',F=[-100.,0.,0.],M=[0.,0.,0.] , stTime=4.5,endTime=5.5)
loadMod.addNodalForce('zMaxRef',F=[100.,0.,0.],M=[0.,0.,0.] , stTime=4.5,endTime=5.5)

loadMod.addNodalForce('yMaxRef',F=[0.,0.,100.],M=[0.,0.,0] , stTime=5.5,endTime=6.5)
loadMod.addNodalForce('zMinRef',F=[0.,-57.735,0.],M=[0.,0.,0.] , stTime=5.5,endTime=6.5)
loadMod.addNodalForce('zMaxRef',F=[0.,57.735,0.],M=[0.,0.,0.] , stTime=5.5,endTime=6.5)

loadMod.writeModelInput(loadFile)

job = ASenDJob()
job.readModelInput(modFile)
job.readConstraints(constFile)
job.readLoads(loadFile)
job.solve(staticLoadTime=[1.,2.,3.,4.,5.,6.] ,saveSolnHist=True,solnHistDir='fiberUnitCell/results')
job.writeNodeResults(nodeResFile,['displacement'],timeSteps=[0,1,2,3,4,5])
job.writeElementResults(elResFile,['stress'],timeSteps=[0,1,2,3,4,5])
job.writeJobInput(jobFile)
job.executeJob()

rp = ResultsProcessor(modFile)
nrfLst = nodeResFile.split('.')
erfLst = elResFile.split('.')
for i in range(0,6):
    nrf = nrfLst[0] + '_timestep' + str(i) + '.yaml'
    erf = erfLst[0] + '_timestep' + str(i) + '.yaml'
    rp.loadNodeResults(nrf)
    rp.loadElementResults(erf)
    rp.plotElementResults('stress',component=(i+1),deformed=True,defScaleFact=100.)
#rp.plotElementResults('stress',component=1,elementSet='all', deformed=True,defScaleFact=100.)


