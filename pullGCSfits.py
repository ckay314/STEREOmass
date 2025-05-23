import numpy as np
import sys, os
import datetime


# object that will hold all the results for an individual ensemble member
class ReconObj:
    "Container for all the reconstructions of a single CME"
    def __init__(self, name):
        self.nRecon = 0
        self.names = []
        self.strDate = name
        self.dateTime = datetime.datetime.strptime(name, "%Y-%m-%dT%H:%M:%S" )
        
        # Fit holders
        self.Lats  = {}
        self.Lons  = {}
        self.Tilts = {}
        self.AWs   = {}
        self.kaps  = {}
        self.vels  = {}
        self.Ms    = {}
        self.rs    = {}
        self.cLon  = None # for checking KINCAT to AFFECTS
        self.reconT = {}



coreDir = '/Users/kaycd1/LLAMAVERSE/CORE/'

files = {'AFFECTS':'AFFECTSgcs.dat', 'Gui':'Gui.dat', 'Kay':'Kay.dat', 'Majumdar':'Majumdar.dat', 'Martinic':'Martinic.dat', 'Sachdeva':'Sachdeva.dat', 'Temmer21':'Temmer21.dat'}

# OG file contents
# AFFECTS: 0 ID, 1 AFF_ID, 2 start, 3 carLon, 4 Lat, 5 tilt, 6 r, 7 kappa, 8 AW, 9 stoLon
# Gui:  0 ID, 1 time, 2 stonyLon, 3 Lat, 4 R, 5 k, 6 tilt, 7 AW, 8 v
# Kay: 0 ID, 1 Kay_ID, 2 start, 3 cor2 time, 4 lat, 5 carLon, 6 tilt, 7 R, 8 kappa, 9 AW, 10 vel, 11 M, 12 stoLon
# Majumdar: 0 ID, 1 time, 2 R, 3 CarLon, 4 Lat, 5 tilt, 6 k, 7 AW, 8 vMax, 9 vFinal, 10 stoLon
# Martinic: 0 ID, 1 MartinicID, 2 lat, 3 lon, 4 tilt, 5 AW, 6 kappa 
# Sachdeva: 0 ID, 1 time, 2 R, 3 v, 4 carLon, 5 Lat, 6 tilt, 7 k, 8 AW, 9 stoLon 
# Temmer21: 0 ID, 1 Tem Id, 2 time, 3 stoLon (?), 4 lat, 5 tilt, 6 AW, 7 kappa, 8 M, 9 vel, 10 source (who fit), 11 21.5 time - NS are dupes to Sachdeva?


# Build the object structure
allObjs = {}
for key in files.keys():
    afile = files[key]
    data = np.genfromtxt(coreDir+afile, dtype=str)
    nFits = len(data)
    myName = afile[:-4]
    
    for i in range(nFits):
        thisLine = data[i]
        myID = str(thisLine[0]) # name of the CME
        
        # Add an object for this event if needed
        if myID not in allObjs.keys():
            allObjs[myID] = ReconObj(thisLine[0])
        lat, lon, tilt, AW, kap, vel, M, rt, rfit = None, None, None, None, None, None, None, None, None
        
        allObjs[myID].names.append(myName)
        allObjs[myID].nRecon += 1
        
        # AFFECTS GCS 
        if afile == 'AFFECTSgcs.dat':
            lat, lon, tilt, AW, kap, rt, rfit = thisLine[4], thisLine[9], thisLine[5], thisLine[8], thisLine[7], thisLine[2], thisLine[6]
            cLon = float(thisLine[3])
            if cLon >180:
                cLon -= 360
            elif cLon < -180:
                cLon += 360.
            allObjs[myID].cLon = cLon
        
        # Gui    
        elif afile == 'Gui.dat':
            lat, lon, tilt, AW, kap, vel, rt, rfit = thisLine[3], thisLine[2], thisLine[6], thisLine[7], thisLine[5], thisLine[8], thisLine[1], thisLine[4]
            rt += ':00'
            
        # Kay
        elif afile == 'Kay.dat':
            lat, lon, tilt, AW, kap, vel, M, rt, rfit = thisLine[4], thisLine[12], thisLine[6], thisLine[9], thisLine[8], thisLine[10], thisLine[11], thisLine[3], thisLine[7]
            rt += ':00'
        
        # Majumdar
        elif afile == 'Majumdar.dat':     
            lat, lon, tilt, AW, kap, vel, rt, rfit = thisLine[4], thisLine[10], thisLine[5], thisLine[7], thisLine[6], thisLine[9], thisLine[1], thisLine[2]
        
        # Martinic
        elif afile == 'Martinic.dat':
            lat, lon, tilt, AW, kap, rt = thisLine[2], thisLine[3], thisLine[4], thisLine[5], thisLine[6], thisLine[1]
            rt += ':00'
        
        # Sachdeva
        elif afile == 'Sachdeva.dat':     
            lat, lon, tilt, AW, kap, vel, rt, rfit = thisLine[5], thisLine[9], thisLine[6], thisLine[8], thisLine[7], thisLine[3], thisLine[1], thisLine[2]
            rt += ':00'
            
        # Temmer 21  
        elif afile == 'Temmer21.dat':     
            lat, lon, tilt, AW, kap, vel, M, rt = thisLine[4], thisLine[3], thisLine[5], thisLine[6], thisLine[7], thisLine[9], float(thisLine[8])/1e15, thisLine[2]
            rt += ':00'  
                               
        # add everything into the objects    
        if lat:
            allObjs[myID].Lats[myName] = '{:.1f}'.format(float(lat))
        if lon:
            lon = float(lon)
            if lon > 180:
                lon -= 360.
            if lon < -180:
                lon += 360.
            allObjs[myID].Lons[myName] = '{:.1f}'.format(float(lon))
        if tilt:
            allObjs[myID].Tilts[myName] = '{:.1f}'.format(float(tilt))
        if AW:
            allObjs[myID].AWs[myName] = '{:.1f}'.format(float(AW))
        if kap:
            allObjs[myID].kaps[myName] = '{:.3f}'.format(float(kap))
        if vel:
            allObjs[myID].vels[myName] = '{:.1f}'.format(float(vel))
        if M:
            allObjs[myID].Ms[myName] = '{:.2f}'.format(float(M))
        if rfit:
            allObjs[myID].rs[myName] = '{:.2f}'.format(float(rfit))
        
        allObjs[myID].reconT[myName] = datetime.datetime.strptime(rt, "%Y-%m-%dT%H:%M:%S" )
            
# Pull in the file with the CoRe to CORSET associations
data = np.genfromtxt(coreDir + 'COR2times.txt', dtype=str)
corsetIDs = {}
for i in range(len(data)):
    thisLine = data[i]
    corsetIDs[str(thisLine[0])] = [thisLine[1], thisLine[4]]      
    
#matchKeys = corsetIDs.keys()        
keyList = np.sort(np.array([str(key) for key in allObjs.keys()]))
for key in keyList:
    myObj = allObjs[key]
    corsetA = corsetIDs[key][0]
    corsetB = corsetIDs[key][1]
    for name in myObj.names:
        r = ('None').rjust(6)
        v = ('None').rjust(8)
        M = ('None').rjust(6)
        reconT = ('None').rjust(16)
        
        if name in myObj.rs.keys():
            r = str(myObj.rs[name]).rjust(6)
        if name in myObj.vels.keys():
            v = str(myObj.vels[name]).rjust(8)
        if name in myObj.Ms.keys():
            M = str(myObj.Ms[name]).rjust(6)
        if name in myObj.reconT.keys():
            reconT = myObj.reconT[name]
            
            
        print (key, corsetA.ljust(14), corsetB.ljust(14), name.ljust(11), str(myObj.Lats[name]).rjust(6), str(myObj.Lons[name]).rjust(8), str(myObj.Tilts[name]).rjust(6), str(myObj.AWs[name]).rjust(6), str(myObj.kaps[name]).rjust(6), r, v, M, reconT)
    