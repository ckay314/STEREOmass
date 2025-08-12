import numpy as np
import os
import sys

global basePath, massPath, outPath
# Path to the corset catalog directory
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'
# Mass files in local directory (for now?)
massPath = '/Users/kaycd1/STEREO_Mass/MassFits/'
# Output path, storing with mass for now but keeping flexible
outPath = '/Users/kaycd1/STEREO_Mass/MassFits/'


def getBase(date, corID):
    # Break date string into parts
    yr = date[0:4]
    date2 = date[2:]
    mn = date[4:6]

    # Pull satellite from corID
    sat = 'A'
    satL = 'a'
    if 'B' in corID:
        sat = 'B'
        satL = 'b'

    # Build the full path
    fullPath = basePath + sat + '/' + yr + '/' + mn + '/' + date2 + '/' + corID + '/' + corID
    infoFile = fullPath +'.info'
    
       
    # Open up the info file
    if os.path.exists(infoFile):
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
        #toRem ='/Volumes/SRP/vourla1/secchi/lz/L0/' + satL +'/img/cor2/' + date + '/'
        lessName = infos[baseIdx].replace('_d4c2'+sat+'.fts','').replace('\n','')
        slIdx = lessName.rfind('/')
        print (corID, lessName[slIdx+1:])


def getAll(thisYr):
    yrFile = 'allGood'+ str(thisYr) + '.txt'
    corData = np.genfromtxt(yrFile, dtype=str)
    allIDs = corData[:,1]
    allDates = corData[:,2]
    nItems = len(allIDs)
    startIdx = 0
    for j in range(nItems-startIdx):
        i = j + startIdx
        #print ('On case ', i, 'of ', nItems)
        myID = allIDs[i].replace('c', '')
        if thisYr == 2014:
            tempDate = allDates[i].replace('/','')
            myDate = tempDate[4:] + tempDate[0:2] + tempDate[2:4]
        else:
            myDate = allDates[i].replace('-', '')
        #print (myDate, myID)
        getBase(myDate, myID)

for yr in [2007 + i for i in range(8)]:
    getAll(yr)