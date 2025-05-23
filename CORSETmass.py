from scipy.io import readsav
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from astropy.io import fits


global basePath, massPath, outPath
# Path to the corset catalog directory
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'
# Mass files in local directory (for now?)
massPath = '/Users/kaycd1/STEREO_Mass/MassFits/'
# Output path, storing with mass for now but keeping flexible
outPath = '/Users/kaycd1/STEREO_Mass/MassFits/'


def getCORSETmass(date, corID, savepng=True, silent=False):
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
    savFile = fullPath +'.sav'
    
    if os.path.exists(savFile):   
        if not silent:
            print('Reading in CORSET files at:')
            print(savFile)
            sav_data = readsav(savFile)
    else:
        print('Cannot find ' + savFile)
        return
        
    # Find the masks
    masks = sav_data['cat'][0][10]
    hts = sav_data['cat'][0][8][0][0] 
    PAs = sav_data['cat'][0][4][0][0]
    wids = sav_data['cat'][0][5][0][0]

    # Open up the info file
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

    # Get the date/time of the individual files
    toRem ='/Volumes/SRP/vourla1/secchi/lz/L0/' + satL +'/img/cor2/' + date + '/'
    baseIm = infos[baseIdx].replace(toRem, '').replace('_d4c2'+sat+'.fts','')
    allIms = []
    i = fileIdx0
    while i <= counter:
        # Some of the .infos have // typos that are messing this up so have to check and clean
        myInfo = infos[i].replace('//','/') 
        subName = myInfo.replace(toRem, '').replace('_d4c2'+sat+'.fts','')[:-1]
        if len(subName) == 15:
            allIms.append(corID+'_'+myInfo.replace(toRem, '').replace('_d4c2'+sat+'.fts','')[:-1])
        else:
            cut1 = myInfo.replace('/Volumes/SRP/vourla1/secchi/lz/L0/' + satL +'/img/cor2/','')
            nextDay = cut1[:8]
            allIms.append(corID+'_'+cut1.replace(nextDay+'/','').replace('_d4c2'+sat+'.fts','')[:-1])
        i+= 1
    allIms = np.array(allIms)


    outFile = massPath+corID+'_'+date+'_CORSETmass.txt'
    f1 = open(outFile, 'w')

    for nowIdx in range(len(allIms)):
        # Open a mass file
        imFile = massPath + allIms[nowIdx]+'_d4c2'+sat+'_mass.fts'
        if not silent:
            print ('Reading in mass image file ', imFile)
        
        hdul = fits.open(imFile) 
        imData = hdul[0].data
        hdul.close()
        
        # Get the appropriate mask from corset results
        nowMask =   np.array(masks[nowIdx])
        isCME = np.where(nowMask == -3)  
        # simple sum within mask
        totMass = np.sum(imData*(nowMask == -3))/1e15
        
        # Print output to screen and file
        if not silent:
            print (allIms[nowIdx], '{:4.2f}'.format(totMass), '{:4.2f}'.format(hts[nowIdx]), PAs[nowIdx], wids[nowIdx])
            print ('')
        f1.write(allIms[nowIdx] + '{:8.2f}'.format(totMass) + '{:8.2f}'.format(hts[nowIdx]) + '{:10.2f}'.format(PAs[nowIdx]) + '{:10.2f}'.format(wids[nowIdx]) + '\n')
        
        # Make a figure if set to save it
        if savepng:
            # Scale it for plotting
            maxMass = np.percentile(np.abs(imData),95)
            rngIt = imData / maxMass
            rngIt[np.where(rngIt < -0.99)] = -0.99
            rngIt[np.where(rngIt > 0.99)] = 0.99
            
            # And make the plot
            fig, ax = plt.subplots(1,1,  figsize=(8,8))
            plt.imshow(rngIt, cmap='binary_r', origin='lower')
            plt.contour(nowMask, levels=[-3], colors=['red'])
            # Add im/base time and total mass as text
            ax.text(0.97, 0.04, '{:4.2f}'.format(totMass)+'x10$^{15}$ g', transform=ax.transAxes, color='w', horizontalalignment='right', verticalalignment='center', fontsize=12)
            ax.text(0.97, 0.07, allIms[nowIdx], transform=ax.transAxes, color='w', horizontalalignment='right', verticalalignment='center', fontsize=12)
            ax.text(0.97, 0.0, 'base: '+baseIm, transform=ax.transAxes, color='w', horizontalalignment='right', verticalalignment='center', fontsize=12)
            # Make the plot pretty
            ax.axis('off')
            plt.subplots_adjust(left=0.,right=1.,top=1,bottom=0)
            plt.savefig(massPath+allIms[nowIdx]+'_'+sat+'_CORSETmass.png')  
            plt.close()
        
    f1.close()  

def getAll(thisYr):
    yrFile = 'allGood'+ str(thisYr) + '.txt'
    corData = np.genfromtxt(yrFile, dtype=str)
    allIDs = corData[:,1]
    allDates = corData[:,2]
    nItems = len(allIDs)
    startIdx = 0
    for j in range(nItems-startIdx):
        i = j + startIdx
        print ('On case ', i, 'of ', nItems)
        myID = allIDs[i].replace('c', '')
        if thisYr == 2014:
            tempDate = allDates[i].replace('/','')
            myDate = tempDate[4:] + tempDate[0:2] + tempDate[2:4]
        else:
            myDate = allDates[i].replace('-', '')
        print (myDate, myID)
        getCORSETmass(myDate, myID)
    
    

#getAll(2014)    
#getCORSETmass('20120712', '2087.1938_0A')
#getCORSETmass('20120712', '2087.2042_0B')
#getCORSETmass('20120713', '2088.4542_0A')
#getCORSETmass('20070530', '216.5049_0B')