import numpy as np
import sunpy.map
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
from scipy.io import readsav
import sys, os
from scipy.spatial import ConvexHull
from scipy.interpolate import CubicSpline

sys.path.append(os.path.abspath('/Users/kaycd1/wombat/')) 
from pyGCS import getGCS

from GCSmass import makeGCSmask, makeOutputs

# Make sunpy/astropy shut up about info/warning for missing metadata
import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('sunpy')
slogger.setLevel(logging.ERROR)
alogger = logging.getLogger('astropy')
alogger.setLevel(logging.ERROR)


global basePath, massPath, outPath
# Path to the corset catalog directory
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'
# Path to the fits files
fitsPathA = '/Volumes/SRP/vourla1/secchi/lz/L0/a/img/cor2/'
fitsPathB = '/Volumes/SRP/vourla1/secchi/lz/L0/b/img/cor2/'
# Mass files in local directory (for now?)
massPath = '/Users/kaycd1/STEREO_Mass/MassFits/'
# Output path, storing with mass for now but keeping flexible
outPath = '/Users/kaycd1/STEREO_Mass/GCSimages/'

def makeCORSETdict():
    data = np.genfromtxt('corsetDirectory.txt', dtype = str)
    CORdict = {}
    for i in range(len(data[:,0])):
        CORdict[data[i,0]] = data[i,1]
    return CORdict


def getCORSETinfo(corfile):
    if 'A' in corfile:
        sat, satL = 'A', 'a'
        pref = fitsPathA
        
    elif 'B' in corfile:
        sat, satL = 'B', 'b'
        pref = fitsPathB
    
    corPath = CORdict[corfile]
    savFile = corPath + corfile + '/' + corfile +'.sav'
    if os.path.exists(savFile):
        sav_data = readsav(savFile)
    else:
        print ('CORSET sav does not exist')
        return [None], None, [None], [None]
    
    infoFile = savFile.replace('.sav','.info')
    hts = sav_data['cat'][0][8][0][0]
    masks = sav_data['cat'][0][10]
    
    # Open up the info file
    f1 = open(infoFile, 'r')
    infos = []
    baseIdx  = -1
    fileIdx0 = -1
    ccounter = -1
    for x in f1:
        ccounter += 1
        if 'Base image' in x:
            baseIdx = ccounter +1
        elif 'Files used' in x:
            fileIdx0 = ccounter + 1
        infos.append(x)
    f1.close
    fts0 = infos[baseIdx].replace('\n','')
    
    # This needs to point to the processed fits in their happy spots!!!!
    
    theFts = []
    i = fileIdx0
    while i <= ccounter:
        myL = infos[i]
        myL = myL.replace('\n','').replace(pref,'')
        theFts.append(myL)
        i+= 1
    # Make sure is sorted (not needed now?)        
    theFts = np.sort(np.array([a for a in theFts]))  
    return np.array(theFts), fts0, masks, hts


def processCase(aIm, aIm0, aUsr, myGCS, myCORmask, CORid):
    # Set up map for the background image 
    corPath = CORdict[CORid]
    mapPath = corPath+CORid+'/fts/'+CORid+'_'+aIm[9:].replace('d4c2A.fts', 'd4c2A_mass.fts').replace('d4c2B.fts', 'd4c2B_mass.fts')
    myMap = sunpy.map.Map(mapPath)

    if 'A' in CORid:
        sat, satL = 'A', 'a'
        pref = fitsPathA
        
    elif 'B' in CORid:
        sat, satL = 'B', 'b'
        pref = fitsPathB
    aIm = aIm.replace(pref, '')
    aIm0 = aIm0.replace(pref, '')
    
    
    # Get the appropriate mass file
    slashIdx = np.char.find(aIm, '/')
    shortName = aIm[slashIdx+1:].replace('.fts', '_mass.fts')
    shortestName = aIm[slashIdx+1:].replace('.fts', '')
    slashIdx0 = np.char.find(aIm0, '/')
    shortName0 = aIm0[slashIdx0+1:].replace('.fts', '_mass.fts')
    shortestName0 = aIm0[slashIdx0+1:].replace('.fts', '')
       
    # Use map and GCS to make masks
    # OO = outer only, ignores inner leg gap in GCS 
    lat, lon, tilt, AW, kap, ht = myGCS
    gMask, gMaskOO = makeGCSmask(lat, lon, tilt, AW, kap, ht, myMap)
    
    # Make all the outputs (numbers and fig)
    outstr = makeOutputs(myMap, myCORmask, gMask, gMaskOO, shortestName, shortestName0, aUsr, CORid=CORid)
    return outstr
    


def runAllCases(saveIt=False):
    # Open the files with info about GCSs
    dataCore = np.genfromtxt('coreGCScases.txt', dtype=str)
    dataIms = np.genfromtxt('haveMassImage.txt', dtype=str)
    
    if saveIt:
        f1 = open('GCSmassProfilesPOS2.txt', 'w')
        
    # Get list ocf unique COR time stamps
    timeIDs = np.unique(dataCore[:,0])
    counter =0
    n2do    = len(timeIDs)
    while counter < n2do:
    #for counter in [249]:
        counter += 1
        aTime = timeIDs[counter-1]
        print ('On event ', counter, ' out of ', n2do)
        
        # Get the matching idx and users
        myIdx  = np.where(dataCore[:,0] == aTime)[0] # same as np.where(dataIms[:,0] == aTime)[0])
        myUsrs = dataIms[myIdx, 1]
        
        # Get the CORSET match
        CORidA = dataIms[myIdx[0], 4]
        CORidB = dataIms[myIdx[0], 9]
        
        #print ('     ', aTime, CORidA, CORidB)
        #print ('     ', myUsrs)
        if CORidA != 'None':
            massFiles, fts0, masks, heights = getCORSETinfo(CORidA)
            print (len(myUsrs), 'fits and', len(heights), 'timesteps for A')            

            # Loop through the usrs for this event
            for i in range(len(myIdx)):
                idx = myIdx[i]
                usr = myUsrs[i]
                
                # Get the GCS parameters for this fit
                myLine = dataCore[idx, :]
                lat, lon, tilt, AW, kap = float(myLine[4]), float(myLine[5]), float(myLine[6]), float(myLine[7]), float(myLine[8])
                myGCS = [lat, lon, tilt, AW, kap, 0.0]
                
                # Loop through the time steps
                for j in range(len(heights)):
                    myGCS[-1] = heights[j]
                    myFile    = massFiles[j]
                    myMask    = masks[j]
                    
                    #print(j, usr, myGCS, myFile, fts0)
                    if myFile != None:
                        outstr = processCase(myFile, fts0, usr, myGCS, myMask, CORidA)
                        if outstr and saveIt:
                            f1.write(outstr+'{:8.2f}'.format(heights[j]) +'\n')
                    
                
        if CORidB != 'None':
            massFiles, fts0, masks, heights = getCORSETinfo(CORidB)
            print (len(myUsrs), 'fits and', len(heights), 'timesteps for B')
            
            # Loop through the usrs for this event
            for i in range(len(myIdx)):
                idx = myIdx[i]
                usr = myUsrs[i]
                
                # Get the GCS parameters for this fit
                myLine = dataCore[idx, :]
                lat, lon, tilt, AW, kap = float(myLine[4]), float(myLine[5]), float(myLine[6]), float(myLine[7]), float(myLine[8])
                myGCS = [lat, lon, tilt, AW, kap, 0.0]
                
                # Loop through the time steps
                for j in range(len(heights)):
                    myGCS[-1] = heights[j]
                    myFile    = massFiles[j]
                    myMask    = masks[j]
                    
                    #print(j, usr, myGCS, myFile, fts0)
                    if myFile != None:
                        outstr = processCase(myFile, fts0, usr, myGCS, myMask, CORidB)
                        if outstr and saveIt:
                            f1.write(outstr+'{:8.2f}'.format(heights[j]) +'\n')
    if saveIt:
        f1.close()

CORdict = makeCORSETdict()
runAllCases(saveIt=False)