# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 09:18:03 2023

@author: evans
"""

class ASenDJob:
    
    def __init__(self):
        self.jobData = dict()
        self.jobData['jobCommands'] = list()
        
    def readModelInput(self,fileName):
        newCmd = dict()
        newCmd['command'] = 'readModelInput'
        newCmd['fileName'] = fileName
        self.jobData['jobCommands'].append(newCmd)