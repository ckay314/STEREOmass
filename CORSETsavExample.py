from scipy.io import readsav
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from astropy.io import fits
import datetime 


global basePath
# Path to the corset catalog directory
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'


# |----- Set up a dictionary to go from ID to date -----|
def makeCORSETdict():
    data = np.genfromtxt('corsetDirectory.txt', dtype = str)
    CORdict = {}
    for i in range(len(data[:,0])):
        CORdict[data[i,0]] = data[i,1]
    return CORdict



# |----- Helper script to see what the ID tags are for a given day -----|
def getCORIDs(date):
    # date in format YYYY-MM-DD
    dateDT = datetime.datetime.strptime(date, "%Y-%m-%d" )
    yr, mn, dy = str(dateDT.year), str(dateDT.month).zfill(2),  str(dateDT.day).zfill(2)
    
    # Setup the correct path (separate ones for STEREO A & B)
    Apath = basePath + 'A/' + yr + '/' + mn + '/' + yr[2:]+mn+dy + '/'
    # Look for the CORSET ID folders and check to see if each one has a sav file (not all do)
    try:
        aFiles = os.listdir(Apath)
        for aFold in aFiles:
            if aFold[0] != '.':
                subFiles = os.listdir(Apath+aFold)
                # files save as ID.sav
                if (aFold+'.sav' in subFiles):
                    print (aFold)
    except:
        print ('No date folder for STA')
    # Repeat for B
    Bpath = basePath + 'B/' + yr + '/' + mn + '/' + yr[2:]+mn+dy + '/'
    try:
        bFiles = os.listdir(Bpath)
        for bFold in bFiles:
            if bFold[0] != '.':
                subFiles = os.listdir(Bpath+bFold)
                # files save as ID.sav
                if (bFold+'.sav' in subFiles):
                    print (bFold)
    except:
        print ('No date folder for STB')
        
def getCORSETmask(CORid):
    # Pull satellite from corID
    sat = 'A'
    satL = 'a'
    if 'B' in CORid:
        sat = 'B'
        satL = 'b'
    
    # Make the path to the save files
    path = CORdict[CORid] + CORid + '/'
    
    # Open the IDL save file
    sav_data = readsav(path+CORid+'.sav')
    
    # Pull out the masks -> array of [nTimes, 1024, 1024]
    masks = sav_data['cat'][0][10]
    
    
    
    # How to get the times for each mask
    # (Saved in the info file on the drive)
    infoFile =  path+CORid+'.info'
    
    # Open the file and read into array
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

    # Get the date/time of the individual files for each mask
    i = fileIdx0
    imgs = []
    while i <= counter:
        # Some of the .infos have // typos that are messing this up so have to check and clean
        myInfo = infos[i].replace('//','/') 
        print (myInfo.replace('\n',''))
        imgs.append(myInfo.replace('\n',''))
        i +=1
    baseIm = infos[baseIdx].replace('//','/').replace('\n','') 
        
    fig, ax = plt.subplots(1,1,  figsize=(8,8))
    # show example for mask 2
    maskIdx = 2 
    # Open up base time and desired step
    # Base
    hdul = fits.open(baseIm) 
    imData0 = hdul[0].data
    hdul.close()
    # Desired Time
    hdul = fits.open(imgs[maskIdx]) 
    imData = hdul[0].data
    hdul.close()
    
    diffIm = imData - imData0
    
    # Very lazy rebinning of obs data from 2048 to 1024 to match mask
    smallIm = diffIm[::2,::2]
    
    plt.imshow(smallIm, cmap='binary', origin='lower')
    plt.contour(masks[maskIdx,:,:,], levels=[-3], c='r')
    plt.show()

    
getCORIDs('2012-07-12')
print ('')
global CORdict
CORdict = makeCORSETdict()
getCORSETmask('2087.1938_0A')