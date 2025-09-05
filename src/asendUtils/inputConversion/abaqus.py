# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 13:38:55 2025

@author: evaande
"""

import yaml
from yaml import CLoader as Cld
import numpy as np

def asendToAbaqus(inFile, outFile):
    inStr = open(inFile, 'r')
    modDat = yaml.load(inStr, Loader=Cld)
    inStr.close()
    
    outStr = open(outFile, 'w')
    outStr.write(f'** Abaqus conversion of ASenD3D model input {inFile}\n**\n')
    
    hasSolid = False
    solidNds = set()
    solidEls = set()
    solidTypes = {'tet4': 'C3D4', 'wedge6': 'C3D6', 'brick8': 'C3D8', 'brickIM': 'C3D8I', 'tet10': 'C3D10'}
    
    for et in modDat['elements']:
        if et['type'] in solidTypes:
            hasSolid = True
            
    if hasSolid:
        outStr.write('** Solid Part\n')
        outStr.write('*Part, name=solid\n')
        for et in modDat['elements']:
            tp = et['type']
            if tp in solidTypes:
                for el in et['connectivity']:
                    solidEls.add(el[0])
                    for i in range(1, len(el)):
                        solidNds.add(el[i])
        
        outStr.write('*Node\n')
        for ni in solidNds:
            nstr = str(modDat['nodes'][ni])
            indStr = f'[{ni}'
            repStr = f'{ni+1}'
            nstr = nstr.replace(indStr,repStr)
            nstr = nstr.replace(']', '\n')
            outStr.write(nstr)
        
        for et in modDat['elements']:
            tp = et['type']
            if tp in solidTypes:
                outStr.write(f'*Element, type={solidTypes[tp]}\n')
                for el in et['connectivity']:
                    strEl = []
                    for i in el:
                        strEl.append(str(i+1))
                    outStr.write(', '.join(strEl) + '\n')
        
        elsets = modDat['sets']['element']
        for es in elsets:
            insct = solidEls.intersection(set(elsets[es]))
            if len(insct) > 0:
                outStr.write(f'*Elset, elset={es}_solid\n')
                lnLen = 0
                for el in insct:
                    outStr.write(f'{el+1}, ')
                    lnLen += 1 
                    if lnLen == 16:
                        outStr.write('\n')
                        lnLen = 0
                outStr.write('\n')
        
        for sec in modDat['sections']:
            if sec['type'] == 'solid':
                secSet = sec['elementSet']
                secMat = sec['material']
                try:
                    ori = str(sec['orientation']).strip('[]')
                    outStr.write(f'*Orientation, name={secSet}\n')
                    outStr.write(ori + '\n')
                    outStr.write('1, 0.\n')
                    outStr.write(f'*Solid Section, elset={secSet}_solid, orientation={secSet}, material={secMat}\n, \n')
                except:
                    outStr.write(f'*Solid Section, elset={secSet}_solid, material={secMat}\n, \n')
                
        
        outStr.write('*End Part\n')
        
    ## (Same for shell part, etc.)
    
    outStr.write('*Assembly, name=Assembly\n')
    
    if hasSolid:
        outStr.write('*Instance, name=solid_inst, part=solid\n')
        outStr.write('*End Instance\n')
        nsets = modDat['sets']['node']
        for ns in nsets:
            insct = solidNds.intersection(set(nsets[ns]))
            if len(insct) > 0:
                outStr.write(f'*Nset, nset={ns}_solid, instance=solid_inst\n')
                lnLen = 0
                for nd in insct:
                    outStr.write(f'{nd+1}, ')
                    lnLen += 1 
                    if lnLen == 16:
                        outStr.write('\n')
                        lnLen = 0 
                outStr.write('\n')
    
    outStr.write('*End Assembly\n')
    
    allMats = modDat['materials']
    for mat in allMats:
        outStr.write(f'*Material, name={mat}\n')
        
        den = allMats[mat]['density']
        outStr.write(f'*Density\n{den},\n')
        
        allE = allMats[mat]['elastic']
        E = allE['E']
        nu = allE['nu']
        G = allE['G']
        outStr.write('*Elastic, type=ENGINEERING CONSTANTS\n')
        outStr.write(f'{E[0]}, {E[1]}, {E[2]}, {nu[0]}, {nu[1]}, {nu[2]}, {G[0]}, {G[1]},\n {G[2]},\n')
                
    outStr.close()