# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 08:00:15 2026

@author: evaande
"""

from asendUtils.model.Model import *
from asendUtils.model.Interaction import *
from asendUtils.model.ParticleSource import *
from asendUtils.meshing.MeshTools import *
from asendUtils.job.ASenDJob import *

nodes = [[1.0, 0.0, 0.0],
         [-1.0, 0.0, 0.0]]
srcMesh = meshFromScratch(nodes, [])

srcs, sectns, meshDat = sourceGroupFromMesh(srcMesh, 1, 1.0, 1.0, [-10,-9], [0, 1], [0, 1], elementSet='particles',
                                            xRange=[-5,5], yRange=[-5,5], zRange=[-5,5])

for s in srcs:
    s.setFrequency(1)
    if s.data['elementSet'] == 'particles_0':
        s.setVelocity(-1.0, 0.0, 0.0)
    else:
        s.setVelocity(1.0, 0.0, 0.0)
        
intctn = collisionInteraction(2.0, 1.0, 0.1, name='particle_int', nodeSet1='particles_0', nodeSet2='particles_1', maxDistance=0.5)

## Define model
myMod = Model()
myMod.addMeshData(meshDat, meshType='mass')
for s in sectns:
    myMod.addSection(s)
myMod.addInteraction(intctn)
for s in srcs:
    myMod.addParticleSource(s)
    
myMod.writeModelInput('collidingParticles.yaml')

