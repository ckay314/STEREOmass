import numpy as np
import matplotlib.pyplot as plt
import datetime 
#from sunpy.coordinates import HeliocentricEarthEcliptic, get_body_heliographic_stonyhurst, get_horizons_coord
from sunpy.time import parse_time
import pickle

# make label text size bigger
plt.rcParams.update({'font.size':14})

dtor = np.pi / 180.

yrs = ['2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014']

t0 = datetime.datetime(2007,1,1)

saveIt = True
loadIt = True

class AllRes:
    "Container for all the results from mass party"
    def __init__(self, name):
            self.CMEtime       = name
            self.CMEtimeDT     = datetime.datetime.strptime(name, "%Y%m%d_%H%M%S" )
            # CORSET properties
            self.CORSETidA     = None
            self.CORSETidB     = None
            self.CORSETim0A    = None # prob not needed but dump in save structure anyway
            self.CORSETim0B    = None
            
            # CORSET max/end properites
            self.CORSETimMA    = None # im/time of max mass
            self.CORSETimMB    = None
            self.CORSETimEA    = None # im/time of end frame
            self.CORSETimEB    = None
            self.CORSETmassMA  = None
            self.CORSETmassMB  = None
            self.CORSETmassEA  = None
            self.CORSETmassEB  = None
            self.CORSEThtMA    = None
            self.CORSEThtMB    = None
            self.CORSEThtEA    = None
            self.CORSEThtEB    = None
            self.CORSETcpaMA   = None
            self.CORSETcpaMB   = None
            self.CORSETcpaEA   = None
            self.CORSETcpaEB   = None
            self.CORSETangMA   = None
            self.CORSETangMB   = None
            self.CORSETangEA   = None
            self.CORSETangEB   = None
            self.CORSETvelMA   = None
            self.CORSETvelMB   = None
            self.CORSETvelEA   = None
            self.CORSETvelEB   = None
            
            # CORSET profiles
            self.CORSETtimesA   = None
            self.CORSETheightsA = None
            self.CORSETtimesB   = None
            self.CORSETheightsB = None
            self.CORSETmassesA  = None
            self.CORSETmassesB  = None
            self.CORSETpixsA    = None
            self.CORSETpixsB    = None
            
            # CORSET CoM properties
            self.CCoMtimesA    = [None]
            self.CCoMtimesB    = [None]
            self.CCoMpixxA     = [None]
            self.CCoMpixyA     = [None]
            self.CCoMpixxB     = [None]
            self.CCoMpixyB     = [None]
            self.CCoMpixxA0    = [None]
            self.CCoMpixyA0    = [None]
            self.CCoMpixxB0    = [None]
            self.CCoMpixyB0    = [None]
            self.CCoMhtpixA    = [None]
            self.CCoMhtRsA     = [None]
            self.CCoMhtpixB    = [None]
            self.CCoMhtRsB     = [None]
            self.CCoM_CPAsA    = [None] 
            self.CCoM_CPAsB    = [None] 
            self.CCoM_hFsA     = [None] # front heights, or whatever
            self.CCoM_hFsB     = [None] # CORSET height actually is (but not CoM)
            self.CCoM_AWsA     = [None] 
            self.CCoM_AWsB     = [None]
            
            
            # CORSET Deproj Properties
            self.deprojTimes   = [None]
            self.deprojMasses  = [None]
            self.deprojHeightsA = [None]
            self.deprojHeightsB = [None]
            self.deprojLons    = [None]
            self.deprojSepA    = [None]
            self.deprojSepB    = [None]
            self.projMassesA   = [None]
            self.projMassesB   = [None]
            self.projHeightsA  = [None]
            self.projHeightsB  = [None]
            self.deprojLatA    = [None]
            self.deprojLatB    = [None]
            self.depSatLonsA   = [None]
            self.depSatLonsB   = [None]
            self.projCPAsA     = [None]
            self.projCPAsB     = [None]
            self.projCPAstdA   = None
            self.projCPAstdB   = None
            

            # CORSET Scaled Deproj Properties
            self.SCLdeprojTimes   = [None]
            self.SCLdeprojMasses  = [None]
            self.SCLdeprojHeightsA = [None]
            self.SCLdeprojHeightsB = [None]
            self.SCLdeprojLons    = [None]
            self.SCLdeprojSepA    = [None]
            self.SCLdeprojSepB    = [None]
            self.SCLprojMassesA   = [None]
            self.SCLprojMassesB   = [None]
            self.SCLprojHeightsA  = [None]
            self.SCLprojHeightsB  = [None]
            self.SCLdeprojLatA    = [None]
            self.SCLdeprojLatB    = [None]
            self.SCLdepSatLonsA   = [None]
            self.SCLdepSatLonsB   = [None]
            self.SCLprojCPAsA     = [None]
            self.SCLprojCPAsB     = [None]
            self.SCLprojCPAstdA   = None
            self.SCLprojCPAstdB   = None
            
            
            
            # GCS properties
            self.GCSnamesA     = None
            self.GCStimesA     = None
            self.GCSnamesB     = None
            self.GCStimesB     = None
            self.GCSvals       = {}
            self.GCSmassesA    = None # true masses using GCS height/time
            self.GCSmassesB    = None            
            self.GCSpixsA      = None
            self.GCSpixBs      = None
            self.GCSmassesOOA  = None # true masses using GCS height/time
            self.GCSmassesOOB  = None            
            self.GCSpixsOOA    = None
            self.GCSpixsOOB    = None
            self.GCS2CmassesA  = None # CORSET values at the same time as
            self.GCS2CmassesB  = None # the GCS recon           
            self.GCS2CpixsA    = None
            self.GCS2CpixBs    = None

            # GCS profiles
            self.GCSprofNamesA  = None
            self.GCSprofHtsA    = None
            self.GCSprofTimesA  = None
            self.GCSprofNamesB  = None
            self.GCSprofHtsB    = None
            self.GCSprofTimesB  = None
            self.GCSprofMassA   = None # approx masses using CORSET h/t
            self.GCSprofMassB   = None
            self.GCSprofPixA    = None
            self.GCSprofPixB    = None
            self.GCSprofMassOOA = None 
            self.GCSprofMassOOB = None
            self.GCSprofPixOOA  = None
            self.GCSprofPixOOB  = None
            
    def printMe(self):
        print ('')
        print ('|--------------------------------------------------------------|')
        print ('|--------------------------------------------------------------|')
        print ('                   Event at', self.CMEtimeDT)
        print ('|--------------------------------------------------------------|')
        print ('|--------------------------------------------------------------|')
        print ('')
        print ('')
        
        # Check if we have CORSET info
        doCORSET = False
        if self.CORSETidA or self.CORSETidB:
            doCORSET = True
        
        # Check if we have GCS info
        doGCS = False
        if len(self.GCSvals.keys()) > 0:
            doGCS = True
            
        if doCORSET:
            print ('|-------------------------------|')
            print ('       CORSET information')
            print ('|-------------------------------|')
            print ('   A id:', self.CORSETidA)
            print ('   B id:', self.CORSETidB)
            print ('')
            
            print ('Background Images (d4c2A/B):')
            print ('   A:', self.CORSETim0A)
            print ('   B:', self.CORSETim0B)
            print ('')
            
            print ('Profile information:')
            nA, nB = 0,0
            if type(self.CORSETtimesA) != type(None):
                nA = len(self.CORSETtimesA)
            if type(self.CORSETtimesB) != type(None):
                nB = len(self.CORSETtimesB)
            print ('   Times for STA:', nA)
            print ('   Times for STB:', nB)
            print ('')
            
            if nA > 0:
                print ('Time series in STA:')
                print ('   Time          Height (Rs) Mass(1e15 g) Pixels')
                for i in range(nA):
                    print('  ', self.CORSETtimesA[i], str(self.CORSETheightsA[i]).rjust(8), str(self.CORSETmassesA[i]).rjust(8), str(self.CORSETpixsA[i]).rjust(10))
                print ('')
            
            if nB > 0:
                print ('Time series in STB:')
                print ('   Time          Height (Rs) Mass(1e15 g) Pixels')
                for i in range(nB):
                    print('  ', self.CORSETtimesB[i], str(self.CORSETheightsB[i]).rjust(8), str(self.CORSETmassesB[i]).rjust(8), str(self.CORSETpixsB[i]).rjust(10))
                print ('')
            
            if nA > 0:
                print('Time of maximum mass (STA):')
                print('   Max reached at:'.ljust(19), self.CORSETimMA)
                print('   Height (Rs):'.ljust(19), self.CORSEThtMA)
                print('   Mass (1e15 g):'.ljust(19), self.CORSETmassMA)
                print('   CPA (deg):'.ljust(19), self.CORSETcpaMA)
                print('   Ang Wid (deg):'.ljust(19), self.CORSETangMA)
                print('   vel (km/s):'.ljust(19), self.CORSETvelMA)
                print('')
            if nB > 0:
                print('Time of maximum mass (STB):')
                print('   Max reached at:'.ljust(19), self.CORSETimMB)
                print('   Height (Rs):'.ljust(19), self.CORSEThtMB)
                print('   Mass (1e15 g):'.ljust(19), self.CORSETmassMB)
                print('   CPA (deg):'.ljust(19), self.CORSETcpaMB)
                print('   Ang Wid (deg):'.ljust(19), self.CORSETangMB)
                print('   vel (km/s):'.ljust(19), self.CORSETvelMB)
                print('')

            if nA > 0:
                print('End of time series (STA):')
                print('   End reached at:'.ljust(19), self.CORSETimEA)
                print('   Height (Rs):'.ljust(19), self.CORSEThtEA)
                print('   Mass (1e15 g):'.ljust(19), self.CORSETmassEA)
                print('   CPA (deg):'.ljust(19), self.CORSETcpaEA)
                print('   Ang Wid (deg):'.ljust(19), self.CORSETangEA)
                print('   vel (km/s):'.ljust(19), self.CORSETvelEA)
                print('')
            if nB > 0:
                print('End of time series (STB):')
                print('   End reached at:'.ljust(19), self.CORSETimEB)
                print('   Height (Rs):'.ljust(19), self.CORSEThtEB)
                print('   Mass (1e15 g):'.ljust(19), self.CORSETmassEB)
                print('   CPA (deg):'.ljust(19), self.CORSETcpaEB)
                print('   Ang Wid (deg):'.ljust(19), self.CORSETangEB)
                print('   vel (km/s):'.ljust(19), self.CORSETvelEB)
                print('')
                
            
        else:    
            print('No CORSET map for this event')
            
            
        if doGCS:
            print ('')
            print ('')
            print ('')
            print ('|-------------------------------|')
            print ('       GCS Reconstructions')
            print ('|-------------------------------|')
        
            print('')
            print('')
            print('GCS parameters from individual reconstructions')
            print('                     lat      lon     tilt      AW     kappa      ht    v    Fit Time ')
            print('                    (deg)    (deg)    (deg)    (deg)             (Rs)  (km/s)  ')
                        
            if type(self.GCSnamesA) != type(None):
                nGCS = len(self.GCSnamesA) # not necessarily same as profiles
                for i in range(nGCS):
                    myGCS = self.GCSvals[self.GCSnamesA[i]]
                    outstr = self.GCSnamesA[i].rjust(12) + ':     '
                    for item in myGCS:
                        outstr += str(item).rjust(6) + '   '
                    print (outstr + self.GCStimesA[i])
                print ('')
                print ('')    
            
                print('')
                print('GCS/CORSET Mass Comparison A')
                print ('                      CORSET vals           GCS values     GCS outer only')
                print ('    GCSFit        Mass(1e15 g) Pixels     Mass     Pixs     Mass     Pixs')
                for i in range(nGCS):
                    outstr = self.GCSnamesA[i].rjust(12) + '       '
                    print (outstr, str(self.GCS2CmassesA[i]).rjust(8), str(self.GCS2CpixsA[i]).rjust(8), str(self.GCSmassesA[i]).rjust(8), str(self.GCSpixsA[i]).rjust(8), str(self.GCSmassesOOA[i]).rjust(8), str(self.GCSpixsOOA[i]).rjust(8))
                print('')
            else:
                print('No GCS/CORSET match in A')
            
            
            if type(self.GCSnamesA) != type(None):
                nGCS = len(self.GCSnamesA)
                print('GCS/CORSET Mass Comparison B')
                print ('                      CORSET vals           GCS values     GCS outer only')
                print ('    GCSFit        Mass(1e15 g) Pixels     Mass     Pixs     Mass     Pixs')
                nGCS = len(self.GCSnamesB)
                for i in range(nGCS):
                    outstr = self.GCSnamesB[i].rjust(12) + '       '
                    print (outstr, str(self.GCS2CmassesB[i]).rjust(8), str(self.GCS2CpixsB[i]).rjust(8), str(self.GCSmassesB[i]).rjust(8), str(self.GCSpixsB[i]).rjust(8), str(self.GCSmassesOOB[i]).rjust(8), str(self.GCSpixsOOB[i]).rjust(8))
                print('')
            
            # profiles
            
            print ('')
            print ('')
            print ('')
            print ('Profiles using CORSET height/time and GCS params')
            if type(self.GCSprofNamesA) != type(None):
                nGCS = len(self.GCSprofNamesA) -1
                print ('Number of CORSET/GCS pairs in STA:', nGCS)
                allNames = '   '
                for i in range(nGCS):
                    allNames += self.GCSprofNamesA[i+1] + ' '
                print ('  Fit by:', allNames)
                print ('')
            
                for i in range(nGCS):
                    print ('Time Series in STA for', self.GCSprofNamesA[i+1])
                    print ('                             CORSET vals             GCS values     GCS outer only')
                    print ('    Time         Height (Rs) Mass(1e15 g) Pixels   Mass     Pixs     Mass     Pixs')
                    for j in range(len(self.GCSprofTimesA[i+1])):
                        print('   ', self.GCSprofTimesA[0][j], str(self.GCSprofHtsA[0][j]).rjust(8), str(self.GCSprofMassA[0][j]).rjust(8), str(self.GCSprofPixA[0][j]).rjust(8), str(self.GCSprofMassA[i+1][j]).rjust(8), str(self.GCSprofPixA[i+1][j]).rjust(8), str(self.GCSprofMassOOA[i+1][j]).rjust(8), str(self.GCSprofPixOOA[i+1][j]).rjust(8))
                    print ('')
                
                
            if type(self.GCSprofNamesB) != type(None):    
                nGCS = len(self.GCSprofNamesB) -1
                print ('Number of CORSET/GCS pairs in STA:', nGCS)
                allNames = '   '
                for i in range(nGCS):
                    allNames += self.GCSprofNamesB[i+1] + ' '
                print ('  Fit by:', allNames)
                print ('')
            
                print ('Profiles using CORSET height/time and GCS params:')
                for i in range(nGCS):    
                    print ('Time Series in STB for', self.GCSprofNamesB[i+1])
                    print ('                             CORSET vals             GCS values     GCS outer only')
                    print ('    Time         Height (Rs) Mass(1e15 g) Pixels   Mass     Pixs     Mass     Pixs')
                
                    for j in range(len(self.GCSprofHtsB[i+1])):
                        print('   ', self.GCSprofTimesB[0][j], str(self.GCSprofHtsB[0][j]).rjust(8), str(self.GCSprofMassB[0][j]).rjust(8), str(self.GCSprofPixB[0][j]).rjust(8), str(self.GCSprofMassB[i+1][j]).rjust(8), str(self.GCSprofPixB[i+1][j]).rjust(8), str(self.GCSprofMassOOB[i+1][j]).rjust(8), str(self.GCSprofPixOOB[i+1][j]).rjust(8))
                    print ('')
                       
        else:
            print('No GCS reconstructions for this event')
            
            
def makeThatRes():
    # Set up a dictionary that will contain the allres objects and
    # indexed by a time id
    res = {}
    
    # make a helper function to go from corset id to time id
    id2time = {}
    id2time[None] = None
    
    # |------------ allGood.txt ------------|
    # make a list of friends between CORSET a and b ids
    # sat, ID, date, time, CPA, wid, vrad, arad, vexp, aexp, friend!
    friendDict = {}
    data1 = np.genfromtxt('allGood.txt', dtype=str)
    for i in range(len(data1[:,1])):
        me = data1[i,1].replace('c','')
        if data1[i,10] != 'NULL':
            friend = data1[i,10].replace('c','')
        else:
            friend = None    
        friendDict[me] = friend
    
    # Have more corset cases than GCS ones, but not all GCS have a (working) CORSET match
    # -> need to sort out consistent indexing option
    # use earliest corset time if we have it, if not GCS fit time (or llamatime?)
    
    # |------------ CORSETpixels.dat ------------|
    # easiest way to get at all (good) event times
    # CORID, filename/time, pixels (does all the As then all the Bs)
    data2 = np.genfromtxt('CORSETpixels.dat', dtype=str)
    uniqIDs = set(data2[:,0])
    uniqIDs = np.sort(np.array([a for a in uniqIDs]))
    # collect the friends for pairs we've looked at already
    haveDone = []
    counter = 0
    for anId in uniqIDs:
        if anId not in haveDone:
            counter += 1
            haveDone.append(anId)
            myFriend = friendDict[anId]
            # find min corset time from pixel list 
            thisIdx = np.where(data2[:,0] == anId)[0]
            # probably already sorted but do anyway
            theTimes = np.sort(np.array([a for a in data2[thisIdx,1]])) 
            if theTimes[0][0] == '/':
                mintime1 = theTimes[0][10:25]
            else:
                mintime1 = theTimes[0][9:24]
            
            # Check the associated event to find earliest time
            if myFriend:
                haveDone.append(myFriend)
                thisIdx = np.where(data2[:,0] == myFriend)[0]
                # Check that it actually has pix/cormap
                if len(thisIdx) > 0: 
                    theTimes = np.sort(np.array([a for a in data2[thisIdx,1]])) 
                    if theTimes[0][0] == '/':
                        mintime2 = theTimes[0][10:25]
                    else:
                        mintime2 = theTimes[0][9:24]
                    
                    # compare with the other time
                    if mintime1 < mintime2:
                        mintime1 = mintime2
                        
                id2time[myFriend] = mintime1
            id2time[anId] = mintime1
                
            # add an entry to mega container
            # Fill in just the name, time and CORSET ids
            thisRes = AllRes(mintime1)
            if 'A' in anId:
                thisRes.CORSETidA = anId
            elif myFriend:
                if 'A' in myFriend:
                    thisRes.CORSETidA = myFriend
            if 'B' in anId:
                thisRes.CORSETidB = anId
            elif myFriend:
                if 'B' in myFriend:
                    thisRes.CORSETidB = myFriend
            
            res[mintime1] = thisRes
            
    
    
    # |------------ CORSET_MaxMasses.dat ------------|
    # CORID, baseIm time, max time, height, mass, CPA, wid, vrad, arad, vexp, aexp
    data3 = np.genfromtxt('CORSET_MaxMasses.dat', dtype=str)
    for i in range(len(data3[:,0])):
        myCORID = data3[i,0]
        myDateID = id2time[myCORID]
        myRes = res[myDateID]
        if 'A' in myCORID:
            myRes.CORSETim0A   = data3[i,1].replace('T', '_')
            myRes.CORSETimMA   = data3[i,2].replace('T', '_')
            myRes.CORSEThtMA   = data3[i,3]
            myRes.CORSETmassMA = float(data3[i,4])
            myRes.CORSETcpaMA  = float(data3[i,5])
            myRes.CORSETangMA  = float(data3[i,6])
            myRes.CORSETvelMA  = float(data3[i,7])
        elif 'B' in myCORID:
            myRes.CORSETim0B   = data3[i,1].replace('T', '_')
            myRes.CORSETimMB   = data3[i,2].replace('T', '_')
            myRes.CORSEThtMB   = float(data3[i,3])
            myRes.CORSETmassMB = float(data3[i,4])
            myRes.CORSETcpaMB  = float(data3[i,5])
            myRes.CORSETangMB  = float(data3[i,6])
            myRes.CORSETvelMB  = float(data3[i,7])                


    # |------------ CORSET_EndMasses.dat ------------|
    # CORID, baseIm time, max time, height, mass, CPA, wid, vrad, arad, vexp, aexp
    data4 = np.genfromtxt('CORSET_EndMasses.dat', dtype=str)
    for i in range(len(data4[:,0])):
        myCORID = data4[i,0]
        myDateID = id2time[myCORID]
        myRes = res[myDateID]
        if 'A' in myCORID:
            myRes.CORSETimEA   = data4[i,2].replace('T', '_')
            myRes.CORSEThtEA   = float(data4[i,3])
            myRes.CORSETmassEA = float(data4[i,4])
            myRes.CORSETcpaEA  = float(data4[i,5])
            myRes.CORSETangEA  = float(data4[i,6])
            myRes.CORSETvelEA  = float(data4[i,7])
        elif 'B' in myCORID:
            myRes.CORSETimEB   = data4[i,2].replace('T', '_')
            myRes.CORSEThtEB   = float(data4[i,3])
            myRes.CORSETmassEB = float(data4[i,4])
            myRes.CORSETcpaEB  = float(data4[i,5])
            myRes.CORSETangEB  = float(data4[i,6])
            myRes.CORSETvelEB  = float(data4[i,7])
    

    
    # |------------ GCStimeVals.txt ------------|
    # timeID, corIDA, corIDB, usr, lat, lon, tilt, AW, kap, ht, vel
    data4b = np.genfromtxt('GCStimeVals.txt',dtype=str)
    time2time = {}
    for i in range(len(data4b[:,0])):
        myida, myidb = data4b[i,1], data4b[i,2]
        mytime = data4b[i,0]
        myUsr  = data4b[i,3]
        myTidA, myTidB = None, None
        # pull time matchin a ID
        if myida in id2time.keys():
            myTidA = id2time[myida] 
        # pull time matching b ID
        if myidb in id2time.keys():
            myTidB = id2time[myidb] 
            
        # check if they are defined and match
        if (myTidA != None) & (myTidA == myTidB):
            resID = myTidA
            myRes = res[resID]
            myRes.CORSETidA = myida
            myRes.CORSETidB = myidb
        # check if one is defined
        elif (myTidA != None) or (myTidB != None):
            if (myTidA != None):
                resID = myTidA
                myRes = res[resID]
                myRes.CORSETidA = myida
                time2time[mytime] = myTidA
                if myida != 'None':
                    id2time[myida] = resID
                if myidb != 'None':
                    id2time[myidb] = resID
            else:
                resID = myTidB
                myRes = res[resID]
                myRes.CORSETidB = myidb
                time2time[mytime] = myTidB
                if myida != 'None':
                    id2time[myida] = resID
                if myidb != 'None':
                    id2time[myidb] = resID
                    
        # nobody defined, need to add
        else:
            myRes = AllRes(mytime)
            resID = mytime
            if myida != 'None':
                id2time[myida] = resID
            if myidb != 'None':
                id2time[myidb] = resID
        # add the GCS params to the res object
        myRes.GCSvals[myUsr] = data4b[i,4:]
        res[resID] = myRes  
              
    # |------------ GCSmassComparison.txt ------------|
    # Time ID string, CORID, COR mass, cor pix, user, GCS mass, GCSpix, outer only GCS mass, pix
    data5 = np.genfromtxt('GCSmassComparison.txt', dtype=str)
    for i in range(len(data5[:,0])):
        myCORid = data5[i,1]
        myDateID = None
        thisRes  = None # make sure not pulling in old one while testing this
        
        # Should already have res obj for each case, just have to find
        # option A from CORSET id
        if myCORid != 'None':
            if myCORid in id2time.keys():
                myDateID = id2time[myCORid]
                thisRes  = res[myDateID]
        else:
            # Option B check if the CORE time is in the obj keys
            fullTime = data5[i,0]
            myDateID = fullTime.replace('_d4c2A', '').replace('_d4c2B', '')
            # Option C convert from CORE time to CORSET time
            if myDateID not in res.keys():
                myDateID = time2time[myDateID]
        thisRes = res[myDateID]
        
        # Fill in the things! Need to check for appending to exist or at None
        # Time ID string, CORID, COR mass, cor pix, user, GCS mass, GCSpix, outer only GCS mass, pix
        myUsr = str(data5[i,4])
        if data5[i,2] != 'None':
            corMass, corPix = float(data5[i,2]), int(data5[i,3])
        else:
            corMass, corPix = None, None
        gcsMass, gcsPix = float(data5[i,5]), int(data5[i,6])
        gcsMassOO, gcsPixOO = float(data5[i,7]), int(data5[i,8])
        
        isA = False
        if 'A' in data5[i,0]:
            isA = True        
            myTime = data5[i,0].replace('_d4c2A', '')
        else:
            myTime = data5[i,0].replace('_d4c2B', '')
        
        # Check if we already have GCS arrays for this event
        if isA:
            if thisRes.GCSnamesA:
                if myUsr not in thisRes.GCSnamesA:
                    thisRes.GCSnamesA.append(myUsr)
                    thisRes.GCStimesA.append(myTime)
                    thisRes.GCSmassesA.append(gcsMass)
                    thisRes.GCSmassesOOA.append(gcsMassOO)
                    thisRes.GCSpixsA.append(gcsPix)
                    thisRes.GCSpixsOOA.append(gcsPixOO)
                    thisRes.GCS2CmassesA.append(corMass)
                    thisRes.GCS2CpixsA.append(corPix)
            else:
                thisRes.GCSnamesA    = [str(data5[i,4])]
                thisRes.GCStimesA    = [myTime]
                thisRes.GCSmassesA   = [gcsMass]
                thisRes.GCSmassesOOA = [gcsMassOO]
                thisRes.GCSmassesA   = [gcsMass]
                thisRes.GCSmassesOOA = [gcsMassOO]
                thisRes.GCSpixsA     = [gcsPix]
                thisRes.GCSpixsOOA   = [gcsPixOO]
                thisRes.GCS2CmassesA = [corMass]
                thisRes.GCS2CpixsA   = [corPix]
                
        else:
            if thisRes.GCSnamesB:
                if myUsr not in thisRes.GCSnamesB:
                    thisRes.GCSnamesB.append(myUsr)
                    thisRes.GCStimesB.append(myTime)
                    thisRes.GCSmassesB.append(gcsMass)
                    thisRes.GCSmassesOOB.append(gcsMassOO)
                    thisRes.GCSpixsB.append(gcsPix)
                    thisRes.GCSpixsOOB.append(gcsPixOO)
                    thisRes.GCS2CmassesB.append(corMass)
                    thisRes.GCS2CpixsB.append(corPix)
            else:
                thisRes.GCSnamesB    = [str(data5[i,4])]
                thisRes.GCStimesB    = [myTime]
                thisRes.GCSmassesB   = [gcsMass]
                thisRes.GCSmassesOOB = [gcsMassOO]
                thisRes.GCSmassesB   = [gcsMass]
                thisRes.GCSmassesOOB = [gcsMassOO]
                thisRes.GCSpixsB     = [gcsPix]
                thisRes.GCSpixsOOB   = [gcsPixOO]
                thisRes.GCS2CmassesB = [corMass]
                thisRes.GCS2CpixsB   = [corPix]
        
        # Add to array (needed for new cases)  
        res[myDateID] = thisRes
                      
            
    # |------------ GCSmassProfiles.txt ------------|
    # Time string, COR ID, COR mass, COR mix, usr, GCS mass, GCSpix, outer only GCS mass, pix, heights 
    # Shouldn't need to add anyone since already ran GCS given-height cases
    data6 = np.genfromtxt('GCSmassProfiles.txt', dtype=str)
    allCORids = data6[:,1]
    unCORids = set(allCORids)
    unCORids = np.sort(np.array([a for a in unCORids]))
        
    for aID in unCORids:
        tId = id2time[aID]
        myRes = res[tId]
        
        isA = False
        if 'A' in aID:
            isA = True        
            toRep = '_d4c2A'
        else:
           toRep = '_d4c2B'
        
        
        thisEv = np.where(allCORids == aID)[0]
        myUsrs = data6[thisEv, 4]
        unUsrs = set(myUsrs)
        for aUsr in unUsrs:
            subIdx = thisEv[np.where(myUsrs == aUsr)[0]]
            corMass, corPix = data6[subIdx, 2].astype(float), data6[subIdx, 3].astype(int)
            gcsMass, gcsPix = data6[subIdx, 5].astype(float), data6[subIdx, 6].astype(int)
            gcsMassOO, gcsPixOO = data6[subIdx, 7].astype(float), data6[subIdx, 8].astype(int)
            hts = data6[subIdx, 9]
            
            if isA:
                times = [a.replace(toRep, '') for a in data6[subIdx,0]]
                if myRes.GCSprofTimesA:
                    myRes.GCSprofNamesA.append(str(aUsr))
                    myRes.GCSprofTimesA.append(times)
                    myRes.GCSprofHtsA.append(hts)
                    myRes.GCSprofMassA.append(gcsMass)
                    myRes.GCSprofPixA.append(gcsPix)
                    myRes.GCSprofMassOOA.append(gcsMassOO)
                    myRes.GCSprofPixOOA.append(gcsPixOO)
                else:
                    myRes.GCSprofNamesA  = ['CORSET', str(aUsr)]
                    myRes.GCSprofTimesA  = [times, times]
                    myRes.GCSprofHtsA    = [hts, hts]
                    myRes.GCSprofMassA   = [corMass, gcsMass]
                    myRes.GCSprofPixA    = [corPix, gcsPix]
                    myRes.GCSprofMassOOA = [corMass, gcsMassOO]
                    myRes.GCSprofPixOOA  = [corPix, gcsPixOO]
                
            else:
                times = [a.replace(toRep, '') for a in data6[subIdx,0]]
                if myRes.GCSprofTimesB:
                    myRes.GCSprofNamesB.append(str(aUsr))
                    myRes.GCSprofTimesB.append(times)
                    myRes.GCSprofHtsB.append(hts)
                    myRes.GCSprofMassB.append(gcsMass)
                    myRes.GCSprofPixB.append(gcsPix)
                    myRes.GCSprofMassOOB.append(gcsMassOO)
                    myRes.GCSprofPixOOB.append(gcsPixOO)
                else:
                    myRes.GCSprofNamesB  = ['CORSET', str(aUsr)]
                    myRes.GCSprofTimesB  = [times, times]
                    myRes.GCSprofHtsB    = [hts, hts]
                    myRes.GCSprofMassB   = [corMass, gcsMass]
                    myRes.GCSprofPixB    = [corPix, gcsPix]
                    myRes.GCSprofMassOOB = [corMass, gcsMassOO]
                    myRes.GCSprofPixOOB  = [corPix, gcsPixOO]
  
  
                  
    # |------------ CORSETprofiles.dat ------------|
    # All the CORSET profiles (not just those with GCS friends)
    data7 = np.genfromtxt('CORSETprofiles.dat', dtype = str)
    allCORids = data7[:,1]
    unCORids = set(allCORids)
    unCORids = np.sort(np.array([a for a in unCORids]))
    
    for aID in unCORids:
        tId = id2time[aID]
        myRes = res[tId]
        
        isA = False
        if 'A' in aID:
            isA = True   
        
        thisEv = np.where(allCORids == aID)[0]
        myMass = data7[thisEv, 2]
        myPix  = data7[thisEv, 4]
        myHeight = data7[thisEv, 3]
        myTimes = data7[thisEv, 0]
        
        if isA:
            myRes.CORSETheightsA = myHeight.astype(float)
            myRes.CORSETtimesA   = myTimes
            myRes.CORSETmassesA  = myMass.astype(float)
            myRes.CORSETpixsA    = myPix.astype(int)
            
        else:
            myRes.CORSETheightsB = myHeight.astype(float)
            myRes.CORSETtimesB   = myTimes
            myRes.CORSETmassesB  = myMass.astype(float)
            myRes.CORSETpixsB    = myPix.astype(int)
            

    
    # |------------ CORSET_CoMs.dat ------------|
    # ID, mass fits file, sun x, suny, CoM x, CoM y, CoM height (pix), CoM height Rs, CORSET h, CPA, AW

    data8 = np.genfromtxt('CORSET_CoMs.dat', dtype=str)
    allCORids = data8[:,0]    
    uniqIDs = set(data8[:,0])
    uniqIDs = np.sort(np.array([a for a in uniqIDs]))
    for aID in uniqIDs:
        tId = id2time[aID]
        myRes = res[tId]
        
        isA = False
        if 'A' in aID:
            isA = True   
        
        thisEv = np.where(allCORids == aID)[0]
        sunx, suny = [], []
        comx, comy = [], []
        times, hpix, hrs = [], [], []
        hFs, CPAs, AWs = [], [], []
        if len(thisEv) > 0:
            for idx in thisEv:
                sunx.append(float(data8[idx,2]))
                suny.append(float(data8[idx,3]))
                comx.append(float(data8[idx,4]))
                comy.append(float(data8[idx,5]))
                hpix.append(float(data8[idx,6]))
                hrs.append(float(data8[idx,7]))
                hFs.append(float(data8[idx,8]))
                CPAs.append(float(data8[idx,9]))
                AWs.append(float(data8[idx,10]))
                times.append(data8[idx,1].replace('_d4c2A_mass.fts', '').replace('_d4c2B_mass.fts', ''))
        if isA:
            myRes.CCoMpixxA  = np.array(comx)
            myRes.CCoMpixyA  = np.array(comy)
            myRes.CCoMpixxA0 = np.array(sunx)
            myRes.CCoMpixyA0 = np.array(suny)
            myRes.CCoMhtpixA = np.array(hpix)
            myRes.CCoMhtRsA  = np.array(hrs)
            myRes.CCoMtimesA = np.array(times)
            myRes.CCoM_CPAsA = np.array(CPAs)
            myRes.CCoM_hFsA  = np.array(hFs)
            myRes.CCoM_AWsA  = np.array(AWs)

        else:
            myRes.CCoMpixxB  = np.array(comx)
            myRes.CCoMpixyB  = np.array(comy)
            myRes.CCoMpixxB0 = np.array(sunx)
            myRes.CCoMpixyB0 = np.array(suny)
            myRes.CCoMhtpixB = np.array(hpix)
            myRes.CCoMhtRsB  = np.array(hrs)
            myRes.CCoMtimesB = np.array(times)
            myRes.CCoM_CPAsB = np.array(CPAs)
            myRes.CCoM_hFsB  = np.array(hFs)
            myRes.CCoM_AWsB  = np.array(AWs)
            
            
    # |------------ CORSET_deProj.dat ------------|
    # output is key, IDA, IDB, time, lonCME, mass [5]
    # lonA, sepA, cpaA, latA, proj mA, projCOM htA, deprojCOM htA, p hFA, dp hFA [14]
    # lonB, sepB, cpaB, latB, proj mB, projCOM htB, deprojCOM htB, p hFB, dp hFB,
    data9 = np.genfromtxt('CORSET_deProj.dat', dtype=str)
    allKeys = data9[:,0]   
    uniqIDs = set(allKeys)
    uniqIDs = np.sort(np.array([a for a in uniqIDs]))
    counter = 0
    for key in uniqIDs:
        counter += 1
        myRes = res[key]
        thisEv = np.where(allKeys == key)[0]
        Ms, ts, HAs, HBs, dirs = [], [], [], [], []
        pMAs, pMBs, pHAs, pHBs = [], [], [], []
        sepAs, sepBs = [], []
        latAs, latBs = [], []
        cpaAs, cpaBs = [], []
        satLonAs, satLonBs = [], []
        for idx in thisEv:
            Ms.append(float(data9[idx,5]))
            ts.append(data9[idx,3])
            # Deprojected front height
            HAs.append(data9[idx,14])
            HBs.append(data9[idx,23])
            # Deprojected CME lon
            depLon = float(data9[idx,4])%360
            dirs.append(depLon)
            # Projected masses
            pMAs.append(float(data9[idx,10]))
            pMBs.append(float(data9[idx,19]))
            # Projected front heights
            pHAs.append(float(data9[idx,13]))
            pHBs.append(float(data9[idx,22]))
            # Separations (from plane of sky)
            sepAs.append(float(data9[idx,7]))
            sepBs.append(float(data9[idx,16]))
            # Calculated latitudes
            latAs.append(float(data9[idx,9]))
            latBs.append(float(data9[idx,18]))
            # Satellite longitude
            aLon = (float(data9[idx,6])+360)%360
            satLonAs.append(aLon)
            bLon = (float(data9[idx,15])+360)%360
            satLonBs.append(bLon)
            # track CPA to get spreaed bc sometimes CORSET is wonky
            cpaAs.append(float(data9[idx,8]))
            cpaBs.append(float(data9[idx,17]))
             #print (key, aLon, bLon, myRes.CORSETcpaMA, myRes.CORSETcpaMB)
        
        cpaAalt = np.copy(cpaAs)
        cpaAalt[np.where(cpaAalt > 180)] -= 360
        cpaBalt = np.copy(cpaBs)
        cpaBalt[np.where(cpaBalt > 180)] -= 360
        #errA1, errA2 = np.std(cpaAs), np.std(cpaAalt)
        myRes.projCPAstdA = np.min([np.std(cpaAs), np.std(cpaAalt)])
        myRes.projCPAstdB = np.min([np.std(cpaBs), np.std(cpaBalt)])
                    
        myRes.deprojTimes   = np.array(ts)
        myRes.deprojMasses  = np.array(Ms)
        myRes.deprojHeightsA = np.array(HAs)
        myRes.deprojHeightsB = np.array(HBs)
        myRes.deprojLons    = np.array(dirs)
        myRes.deprojSepA    = np.array(sepAs)
        myRes.deprojSepB    = np.array(sepBs)
        myRes.projMassesA   = np.array(pMAs)
        myRes.projMassesB   = np.array(pMBs)
        myRes.projHeightsA  = np.array(pHAs)
        myRes.projHeightsB  = np.array(pHBs)
        myRes.deprojLatA    = np.array(latAs)
        myRes.deprojLatB    = np.array(latBs)
        myRes.depSatLonsA   = np.array(satLonAs)
        myRes.depSatLonsB   = np.array(satLonBs)
        myRes.projCPAsA     = np.array(cpaAs)
        myRes.projCPAsB     = np.array(cpaBs)
              
    

    #for key in res.keys():
    #    print (key)
    #    print ('  ', res[key].CORSETidA, res[key].CORSETidB)       
    #print (len(res.keys()))
    return res, id2time, time2time


if __name__ == '__main__':
    if saveIt:
        res, id2time, time2time = makeThatRes()
        f =  open('allRes.pkl', 'wb')
        pickle.dump([res, id2time, time2time], f)
        f.close()
    
    if loadIt:
        f =  open('allRes.pkl', 'rb')
        res, id2time, time2time = pickle.load(f)
        f.close()    
    

    keys = np.sort(np.array([key for key in res.keys()]))
    #res[id2time['2087.1938_0A']].printMe()
    #for key in keys:
    #    res[key].printMe()