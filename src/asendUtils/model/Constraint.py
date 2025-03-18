# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 12:14:29 2023

@author: evans
"""

class Constraint:
    
    def __init__(self,constType):
        types = 'displacement, temperature'
        if(constType not in types):
            errstr = 'Error: ' + constType + ' is not a currently supported constraint type. Valid types: ' + types
            raise TypeError(errstr)
        self.constData = dict()
        self.constData['type'] = constType
        
    def setActiveTime(self,stTime=0.0,endTime=1e+100):
        self.constData['activeTime'] = str([stTime,endTime])
        
    def addTerm(self,nodeSet,dof,coef):
        newTrm = dict()
        newTrm['nodeSet'] = nodeSet
        newTrm['dof'] = dof
        newTrm['coef'] = coef
        try:
            self.constData['terms'].append(newTrm)
        except:
            terms = list()
            terms.append(newTrm)
            self.constData['terms'] = terms
        return
    
    def setRHS(self,rhs):
        self.constData['rhs'] = rhs