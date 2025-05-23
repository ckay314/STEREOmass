import numpy as np
import os
import sys


global basePath
# Path to the corset catalog directory
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'

yrs = ['2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014']
sats = ['A', 'B']
mns = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']



def checkYr(yr, sat):
    fullPath = basePath + sat + '/' + yr
    for mn in mns:
        mnPath = fullPath + '/'+ mn + '/'
        if os.path.isdir(mnPath):
            allFiles = os.listdir(mnPath)
            allFolds = []
            for aFile in allFiles:
                if '.' not in aFile: 
                    allFolds.append(aFile)
            # Alphabetize for my own happiness
            allFolds = np.sort(np.array([aFold for aFold in allFolds]))
            for aFold in allFolds:
                moreFolds = allFiles = os.listdir(mnPath+'/'+aFold)
                for bFold in moreFolds:
                    if bFold[0] != '.':
                        #print (bFold, mnPath+aFold+'/')
                        theFiles = os.listdir(mnPath+'/'+aFold+'/'+bFold)
                        hasSav = False
                        hasFts = False
                        for aFile in theFiles:
                            if aFile[0] != '.':
                                if 'sav' in aFile:
                                    hasSav = True
                                elif 'fts' in aFile:
                                    hasFts = True
                        if hasFts:
                            ftsFiles = os.listdir(mnPath+'/'+aFold+'/'+bFold)
                            if len(ftsFiles) == 0:
                                hasFts = False
                        if hasSav and not hasFts:
                            print (aFold, bFold, hasSav, hasFts)
    

for aYr in yrs:   
    for aSat in sats: 
        checkYr(aYr, aSat)