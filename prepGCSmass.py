import numpy as np
import matplotlib.pyplot as plt
from scipy.io import readsav
import sys, os
import datetime

global basePath, massPath, outPath
# Path to the corset catalog directory
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'
# Mass files in local directory (for now?)
massPath = '/Users/kaycd1/STEREO_Mass/MassFits/'
# Output path, storing with mass for now but keeping flexible
outPath = '/Users/kaycd1/STEREO_Mass/MassFits/'
# Path to raw data
fitsPath = '/Volumes/SRP/vourla1/secchi/lz/L0/'




dataGCS = np.genfromtxt('coreGCScases.txt', dtype=str)

def getFits():
    # Want to make a file with the base and gcs time fits files for every CoRe GCS case
    # and check if we already made those mass images for corset
    # Output: CoRe_Date, Catalog, GCS_Time, CORSET_ID, Sat, Base_fits, GCS_fits, doneInCORSET
    
    # Loop through the CoRe GCS file
    for aLine in dataGCS:
        # Add CoRe_Date
        myDate = aLine[0]
        outLine = myDate
        
        # Add Catalog
        outLine += aLine[3].rjust(14)
        
        # Add GCS time 
        myGCSdate = aLine[12] + 'T' + aLine[13]
        outLine += '   ' + myGCSdate
        
        # Reformat the CoRe Date
        yr = myDate[:4]
        mn = myDate[5:7]
        dy = myDate[8:10]
        cleanDate = yr+mn+dy
        lilDate = yr[2:]+mn+dy
        
        # Reformat the GCS Date
        yrG = myGCSdate[:4]
        mnG = myGCSdate[5:7]
        dyG = myGCSdate[8:10]
        cleanDateG = yrG+mnG+dyG
        lilDateG = yrG[2:]+mnG+dyG
        
        # Pull in the CORSET ID 
        myCORidA = aLine[1].replace('c','').replace('?','')
        myCORidB = aLine[2].replace('c','').replace('?','')
        corFold = '  None'
        baseIm = '                              None'
        GCSIm  = '                              None'
        
        # Run for STA 
        if myCORidA != 'Null':
            outLine += ' A ' + myCORidA
            corFold, baseIm, GCSIm = checkCOR(myCORidA, cleanDate, myGCSdate)
            baseIm = baseIm.replace('\n','')
            GCSIm = GCSIm.replace('\n','')
            if '/' not in GCSIm:
                friend = findFriend(myGCSdate, 'A')
                if friend:
                    GCSIm = friend
        else:
            outLine += ' A        None' 
            
        if corFold == '  None':
            friend = findFriend(myGCSdate, 'A')
        
        outLine += ' ' + ' ' + corFold + ' ' + baseIm + ' ' + GCSIm
        #print (outLine) 
        
        # Run for STB 
        if myCORidB != 'Null':
            outLine += ' B ' + myCORidB
            corFold, baseIm, GCSIm = checkCOR(myCORidB, cleanDate, myGCSdate)
            baseIm = baseIm.replace('\n','')
            GCSIm = GCSIm.replace('\n','')
            if '/' not in GCSIm:
                friend = findFriend(myGCSdate, 'B')
                if friend:
                    GCSIm = friend
        else:
            outLine += ' B        None' 
            
        if corFold == '  None':
            friend = findFriend(myGCSdate, 'B')
        
        outLine += ' ' + ' ' + corFold + ' ' + baseIm + ' ' + GCSIm
        if aLine[3] == 'Majumdar':
            print (outLine) 
        

def checkCOR(CORID, coret, GCSt):
    baseIm = '                              None'
    GCSIm  = '                              None'
    
    mySat = CORID[-1]
    yr = coret[:4]
    mn = coret[4:6]
    lilDate = coret[2:]
    
    # Figure out if STA or STB
    if mySat not in ['A', 'B']:
        print ('Error in determining satellite')
        return baseIm, GCSIm 
        
    # Find the corset info file
    infoFile = mySat+'/' +yr + '/' + mn + '/' + lilDate + '/' + CORID + '/' + CORID + '.info'
    theInfos = None
    # Does the CORSET file exist where expected
    if os.path.isfile(basePath+infoFile):
        fullInfoPath = basePath+infoFile
        corFold = lilDate 
    else:
        corFold = '  None'
        # If not check if that folder exists. If it does then just no fit
        upDir = mySat+'/' +yr + '/' + mn + '/' + lilDate + '/' + CORID + '/'
        # Otherwise check the day before to see if hiding there
        if not (os.path.isdir(basePath+upDir)):
            # check if maybe in the date before
            myDT   = datetime.datetime(int(yr), int(mn), int(coret[6:]))
            prevDT = myDT-datetime.timedelta(days=1)
            pmn    = str(prevDT.month).zfill(2)
            pdy    = str(prevDT.day).zfill(2)
            pLilDate = yr[2:] + pmn + pdy
            infoFile = mySat+'/' +yr + '/' + pmn + '/' + pLilDate + '/' + CORID + '/' + CORID + '.info'
            if os.path.isfile(basePath+infoFile):
                fullInfoPath = basePath+infoFile 
                corFold = pLilDate
    if corFold == '  None':
        return corFold, baseIm, GCSIm 
    
    # Use the info file to get base and prev fitted mass images
    # Open up the info file
    f1 = open(fullInfoPath, 'r')
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
    
    # Get the base image
    if mySat == 'A': 
        toRem = '/Volumes/SRP/vourla1/secchi/lz/L0/a/img/cor2/'
    elif mySat == 'B':
        toRem = '/Volumes/SRP/vourla1/secchi/lz/L0/b/img/cor2/'
    baseIm = infos[baseIdx].replace(toRem, '')
    
    # Get the (closest within reason) GCS match
    allIms = []
    i = fileIdx0
    strGCSt = GCSt.replace('-','').replace('_','').replace('T','').replace(':','')
    
    while i <= counter:
        # Some of the .infos have // typos that are messing this up so have to check and clean
        myInfo = infos[i].replace('//','/') 
        subName = myInfo.replace(toRem, '').replace('_d4c2'+mySat+'.fts','')[:-1].replace('_','')
        allIms.append(subName)
        if strGCSt in subName:
            GCSIm = myInfo.replace(toRem, '')
            return corFold, baseIm, GCSIm 
        i+= 1
    
    # Didn't have an exact match but check if anyone was close
    if len(allIms) > 1:
        diffs = []
        for aIm in allIms:
            ImInt = int(aIm[9:])
            diffs.append(np.abs(int(strGCSt) - ImInt))
        minDiff = np.min(diffs)
        if minDiff < 3000: # 30 mins in the int/str format
            gIdx = np.where(diffs == minDiff)[0]
            GCSIm = infos[fileIdx0+gIdx[0]].replace(toRem,'')
        
    return corFold, baseIm, GCSIm    


def findFriend(dateIn, sat):
    friend = None
    tag = {'a':'d4c2A', 'b':'d4c2B'} 
    if sat == 'A':
        sat = 'a'
    elif sat == 'B':
        sat = 'b'

    stripDate = dateIn[:10].replace('-','')
    stripTime = dateIn[11:].replace(':', '')
    
    fitsFold = sat+'/img/cor2/'+stripDate
    myTag = tag[sat]
    
    if os.path.isdir(fitsPath + fitsFold):
        # Get all the files within this folder
        theseFiles = os.listdir(fitsPath + fitsFold)
        
        # Double check that it found files
        if len(theseFiles) > 0:
            fitsTimes = []
            for aFile in theseFiles:
                if myTag in aFile:
                   fitsTimes.append(aFile[9:15])
            fitsTimes = np.array(fitsTimes)
            timeDiffs = np.abs(fitsTimes.astype(int) - int(stripTime))
            closest = np.where(timeDiffs == np.min(timeDiffs))[0]
            cFits = stripDate+'_'+fitsTimes[closest[0]]+'_'+myTag+'.fts'
            friend = stripDate +'/'+ cFits

    return friend








# These are old and can delete
def checkForCORSET():
    for aLine in dataGCS:
        myDate = aLine[0]
        myGCSdate = aLine[12] + 'T' + aLine[13]
        yr = myDate[:4]
        mn = myDate[5:7]
        dy = myDate[8:10]
        cleanDate = yr+mn+dy
        lilDate = yr[2:]+mn+dy
        
        myCORidA = aLine[1].replace('c','').replace('?','')
        myCORidB = aLine[2].replace('c','').replace('?','')
    
        #print (myDate, myCORidA, myCORidB)
        # look for file of form 2779.1417_0A.corset.ht in basepath + A/2014/08/140815/2851.2458_0A/ 
        Afile = 'A/' +yr + '/' + mn + '/' + lilDate + '/' + myCORidA + '/' + myCORidA + '.corset.ht'
        Bfile = 'B/' +yr + '/' + mn + '/' + lilDate + '/' + myCORidB + '/' + myCORidB + '.corset.ht'
        #print (myCORidA == 'Null', myCORidB == 'Null')
        haveA = os.path.isfile(basePath+Afile)
        haveB = os.path.isfile(basePath+Bfile)
        
        # If we have a CORSET file match need to double check that we actually have 
        # the specific GCS time already
        # Build the full path
        if haveA:
            haveA = False
            Afile = 'A/' +yr + '/' + mn + '/' + lilDate + '/' + myCORidA + '/' + myCORidA + '.info'
            # Open up the info file
            f1 = open(basePath+Afile, 'r')
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
            
            
            # Get the date/time of the individual files
            toRem ='/Volumes/SRP/vourla1/secchi/lz/L0/' + 'a' +'/img/cor2/' + cleanDate + '/'
            baseIm = infos[baseIdx].replace(toRem, '').replace('_d4c2'+'A'+'.fts','')
            print (infos[baseIdx])
            print (toRem)
            print (baseIm)
            print (sd)
            allIms = []
            i = fileIdx0
            GCSt = aLine[13].replace(':','')
            while i <= counter:
                # Some of the .infos have // typos that are messing this up so have to check and clean
                myInfo = infos[i].replace('//','/') 
                subName = myInfo.replace(toRem, '').replace('_d4c2'+'A'+'.fts','')[:-1]
                if GCSt in subName:
                    haveA = True
                    print ('found one')
                i+= 1
        
                    
        print (myDate, aLine[3].rjust(14), myCORidA.rjust(16), str(haveA).rjust(10), myCORidB.rjust(16), str(haveB).rjust(10), myGCSdate)
        # save this in haveMassImage.txt



def checkForFits():
    firstPass = np.genfromtxt('haveMassImage.txt', dtype=str)
    toCheck = {}
    # Find the lines with False for having mass images already
    # Make sure we a
    for aLine in firstPass:
        date = str(aLine[0])
        dateG = str(aLine[6])
        if aLine[3] == 'False':
            if date not in toCheck.keys():
                toCheck[date] = []
            toCheck[date].append(['a', dateG])
        if aLine[5] == 'False':
            if date not in toCheck.keys():
                toCheck[date] = []
            if 'b' not in toCheck[date]:
                toCheck[date].append(['b', dateG])
    
    tag = {'a':'d4c2A', 'b':'d4c2B'}        
    for key in toCheck.keys():
        stripDate = key[:10].replace('-','')
        stripTime = key[11:].replace(':', '')
        
        for item in toCheck[key]:
            sat = item[0]
            dateG = item[1]
            stripDateG = dateG[:10].replace('-','')
            stripTimeG = dateG[11:].replace(':', '')
            
            # Base img and fits img might be diff days
            fitsFold = sat+'/img/cor2/'+stripDate
            fitsFoldG = sat+'/img/cor2/'+stripDateG
            myTag = tag[sat]
            
            # Get the file for the base difference
            confirmBase = False
            if os.path.isdir(fitsPath + fitsFold):
                # Get all the files within this folder
                theseFiles = os.listdir(fitsPath + fitsFold)

                # Double check that it found files
                if len(theseFiles) > 0:
                    fitsTimes = []
                    for aFile in theseFiles:
                        if myTag in aFile:
                           fitsTimes.append(aFile[9:15])
                    fitsTimes = np.array(fitsTimes)
                    timeDiffs = np.abs(fitsTimes.astype(int) - int(stripTime))
                    closest = np.where(timeDiffs == np.min(timeDiffs))[0]
                    cFits = stripDate+'_'+fitsTimes[closest[0]]+'_'+myTag+'.fts'
                    confirmBase = os.path.isfile(fitsPath+fitsFold+'/'+cFits)
                    
            
            # Get the file at the GCS time
            confirmGCS = False
            if os.path.isdir(fitsPath + fitsFoldG):
                # Get all the files within this folder
                theseFiles = os.listdir(fitsPath + fitsFoldG)

                # Double check that it found files
                if len(theseFiles) > 0:
                    fitsTimes = []
                    for aFile in theseFiles:
                        if myTag in aFile:
                           fitsTimes.append(aFile[9:15])
                    fitsTimes = np.array(fitsTimes)
                    timeDiffs = np.abs(fitsTimes.astype(int) - int(stripTimeG))
                    closest = np.where(timeDiffs == np.min(timeDiffs))[0]
                    cFitsG = stripDateG+'_'+fitsTimes[closest[0]]+'_'+myTag+'.fts'
                    confirmGCS = os.path.isfile(fitsPath+fitsFoldG+'/'+cFitsG)
                                
            # List the files we need to get for each key
            if not confirmBase:
                cFits = '                     None'
            if not confirmGCS:
                cFitsG = '                     None'
                
            print (key, cFits, cFitsG)
            
            '''if os.path.isdir(fitsPath + fitsFold):
                theseFiles = os.listdir(fitsPath + fitsFold)
                myTag = tag[sat]
                # pick the subfiles nearest the time
                
                #print (key, len(theseFiles),  stripTime)
                if len(theseFiles) > 0:
                    fitsTimes = []
                    for aFile in theseFiles:
                        if myTag in aFile:
                           fitsTimes.append(aFile[9:15])
                    fitsTimes = np.array(fitsTimes)
                    timeDiffs = np.abs(fitsTimes.astype(int) - int(stripTime))
                    closest = np.where(timeDiffs == np.min(timeDiffs))[0]
                    cFits = stripDate+'_'+fitsTimes[closest[0]]+'_'+myTag+'.fts'
                    
                    
                    if os.path.isfile(fitsPath+fitsFold+'/'+cFits):
                        print (key, fitsPath+fitsFold+'/'+cFits)
                        print ('    ', stripDateG, stripTimeG)
                        print (sd)
                else:
                    print (key, 'NoFits')
            else:
                print (key, 'NoFits')'''
    # b/img/cor2/20080409/20080409_233754_d4c2B.fts

        
getFits()
#checkForCORSET()
#checkForFits()