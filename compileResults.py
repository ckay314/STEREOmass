import numpy as np
import os
import sys


global basePath
# Path to the corset catalog directory
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'

corsetDict = {}
cdFile = np.genfromtxt('corsetDirectory.txt', dtype=None)
padKeys = []
for line in cdFile:
    corsetDict[line[0]] = line[1]
    padKeys.append(line[0].zfill(12))
sortKeys = np.sort(np.array(padKeys)) 
sortKeys = np.array([aKey.strip('0') for aKey in sortKeys])
print (sortKeys)

def makeMegaFile():
    # Open the text files
    f1 = open('CORSET_MaxMasses.dat', 'w')
    f2 = open('CORSET_EndMasses.dat', 'w')
    
    CORSETdata = np.genfromtxt('allGood.txt', dtype=str)
    CORnames = CORSETdata[:,1]
    CORnames = np.array([name.replace('c','') for name in CORnames])
    
    nItems = len(sortKeys)
    for i in range(nItems):
        key = sortKeys[i]
        myPath = corsetDict[key]
        myDate = '20'+myPath[-7:-1]
        massFile = myPath+key+'/'+key+'_'+myDate+'_CORSETmass.txt'
        #print (massFile)
        if os.path.exists(massFile):
            CORidx = np.where(CORnames == key)[0]
            myStuff = CORSETdata[CORidx[0],4:11]
            
            # get the mass info from the indiv file
            thisFile = np.genfromtxt(massFile, dtype=str)
            print (massFile)
            try:
                tags    = thisFile[:,0]
                masses  = thisFile[:,1].astype(float)
                heights = thisFile[:,2]
            except:
                print ('one line')
                tags    = np.array([thisFile[0]])
                masses  = np.array([thisFile[1].astype(float)])
                heights = np.array([thisFile[2]])
            maxID = np.where(masses == np.max(masses))[0]
            maxID = maxID[-1] # couple cases with more than 1 at max val, take latest
            maxTime = tags[maxID].replace(key+'_','').replace('_','T')
            endTime = tags[-1].replace(key+'_','').replace('_','T')
            outStrMax = key + ' ' + myDate+'T'+CORSETdata[i,3].replace(':','') + ' '+ maxTime + ' ' + heights[maxID].rjust(5) + ' ' + str(masses[maxID]).rjust(5)
            outStrEnd = key + ' ' + myDate+'T'+CORSETdata[i,3].replace(':','') + ' '+ endTime + ' ' + heights[-1].rjust(5) + ' ' + str(masses[-1]).rjust(5)
            for item in myStuff:
                    outStrMax = outStrMax + ' ' + item.rjust(7)
                    outStrEnd = outStrEnd + ' ' + item.rjust(7)
            f1.write(outStrMax+'\n')
            f2.write(outStrEnd+'\n')
    f1.close()
    f2.close()
    

def makeMegaErFile():
    # Open the text files
    f1 = open('CORSETprofilesNOPIX.dat', 'w')
    
    CORSETdata = np.genfromtxt('allGood.txt', dtype=str)
    CORnames = CORSETdata[:,1]
    CORnames = np.array([name.replace('c','') for name in CORnames])
    
    nItems = len(sortKeys)
    for i in range(nItems):
        print (i, ' out of ', nItems)
        key = sortKeys[i]
        myPath = corsetDict[key]
        myDate = '20'+myPath[-7:-1]
        massFile = myPath+key+'/'+key+'_'+myDate+'_CORSETmass.txt'
        #print (massFile)
        if os.path.exists(massFile):
            CORidx = np.where(CORnames == key)[0]
            myStuff = CORSETdata[CORidx[0],4:11]
            
            # get the mass info from the indiv file
            thisFile = np.genfromtxt(massFile, dtype=str)
            try:
                tags    = thisFile[:,0]
                masses  = thisFile[:,1].astype(float)
                heights = thisFile[:,2]
            except:
                print ('one line')
                tags    = np.array([thisFile[0]])
                masses  = np.array([thisFile[1].astype(float)])
                heights = np.array([thisFile[2]])
                
            for j in range(len(masses)):
                outline = tags[j].replace(key+'_', '') + ' ' + key + ' ' + str(masses[j]) + ' ' + str(heights[j]) +'\n'
                f1.write(outline)   

    f1.close()

def addPix():
    data  = np.genfromtxt('CORSETprofilesNOPIX.dat', dtype = str)
    data2 = np.genfromtxt('CORSETpixels.dat', dtype = str)
    
    f1 = open('CORSETprofiles.dat', 'w')
    
    
    otherIDs = data2[:,0]
    otherDates = data2[:,1] 
    for i in range(len(data[:,0])):
        myDate = data[i,0]
        myID   = data[i,1]
        myMass = data[i,2]
        myHeight = data[i,3]
        if 'A' in myID:
            suff = '_d4c2A.fts'
        else:
            suff = '_d4c2B.fts'
        longDate = myDate[:8]+'/'+myDate + suff
        matchIdx = np.where((otherIDs == myID) & (otherDates == longDate))[0]
        
        outline = myDate + ' ' + myID + ' ' + myMass + ' ' + myHeight + ' ' + data2[matchIdx[0],2]
        f1.write(outline+'\n')

makeMegaErFile()
addPix()