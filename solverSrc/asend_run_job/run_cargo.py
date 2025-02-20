# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 21:07:23 2025

@author: evans
"""

import subprocess

res = subprocess.run(['cargo','build'],capture_output=True,text=True)

outLst = res.stderr.split('\n')

sinceErr = 0
ldsp = "        "
for s in outLst:
    if("error[" in s or "error:" in s):
        sinceErr = 0
        ldsp = ""
    else:
        sinceErr += 1
        ldsp = "        "
    if(sinceErr < 15):
        print(ldsp + s)
        
lstLen = len(outLst)
for i in range(lstLen-5,lstLen):
    print(outLst[i])

# outFile = open('cargoOutput.txt','w')
# outFile.write('begin out\n\n')
# outFile.write(res.stdout)
# outFile.write('\nbegin err\n\n')
# outFile.write(res.stderr)
# outFile.close()
