import numpy as np
import os
import sys
from scipy.io import readsav


global basePath
# Path to the corset catalog directory
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'

yrs = ['2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014']
sats = ['A', 'B']
mns = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

u2l = {'A':'a', 'B':'b'}


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
                        savFile = mnPath+aFold+'/'+bFold+'/'+bFold+'.sav'
                        infoFile = mnPath+aFold+'/'+bFold+'/'+bFold+'.info'
                        if os.path.exists(savFile):   
                            print('Reading in CORSET files at:')
                            print(savFile)
                            sav_data = readsav(savFile)
                            masks = sav_data['cat'][0][10]
                            dates = sav_data['cat'][0][2]
                            times = sav_data['cat'][0][3]
                            
                            f1 = open(infoFile, 'r')
                            infos = []
                            baseIdx  = -1
                            fileIdx0 = -1
                            counter = -1
                            for x in f1:
                                counter += 1
                                if 'Base image' in x:
                                    baseIdx = counter +1
                                elif 'Files used' in x:
                                    fileIdx0 = counter + 1
                                infos.append(x)
                            f1.close
                            
                            
                            for i in range(len(masks)):
                                myPrefix = '/Volumes/SRP/vourla1/secchi/lz/L0/' + u2l[sat] + '/img/cor2/'
                                myFile = infos[i+fileIdx0].replace(myPrefix, '').replace('\n','')
                                f2.write(bFold + ' ' + myFile + ' '+str(np.sum(masks[i] == -3)).rjust(8) + '\n')

f2 = open('CORSETpixels.dat', 'w')
for aYr in yrs:   
    for aSat in sats:
        print (aYr, aSat) 
        checkYr(aYr, aSat)
f2.close()