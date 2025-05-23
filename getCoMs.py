import numpy as np
import os
import sys
from scipy.io import readsav

import sunpy.map
from astropy.coordinates import SkyCoord

import matplotlib.pyplot as plt


# Make sunpy/astropy shut up about info/warning for missing metadata
import logging
logging.basicConfig(level='INFO')
slogger = logging.getLogger('sunpy')
slogger.setLevel(logging.ERROR)
alogger = logging.getLogger('astropy')
alogger.setLevel(logging.ERROR)


global basePath
# Path to the corset catalog directory
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'

#yrs = ['2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014']
yrs = ['2014']
sats = ['A', 'B']
mns = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

u2l = {'A':'a', 'B':'b'}


xpix = np.array(range(1024))
ypix = np.array(range(1024))
global xx, yy
xx, yy = np.meshgrid(xpix, ypix)

def calcCoM(ftsFile, mask):
    myMap = sunpy.map.Map(ftsFile)
    skyPt = SkyCoord(x=0, y=0, z=0, unit='R_sun', representation_type='cartesian', frame='heliographic_stonyhurst')
    sx, sy = myMap.wcs.world_to_pixel(skyPt)
    rSunInPix = (myMap.rsun_obs/myMap.scale[0]).to_value()
        
    imData = myMap.data
    mask0 = (mask == -3)
    
    # Get pix of center of mass
    totMass = np.sum(imData*mask0)
    avgx    = np.sum(imData*xx*mask0) / totMass
    avgy    = np.sum(imData*yy*mask0) / totMass
    
    # Get height, convert to pixels
    heightPix = np.sqrt((avgx-sx)**2 + (avgy-sy)**2)
    heightRs  = heightPix / rSunInPix
    
    '''maxMass = np.percentile(np.abs(imData),95)
    rngIt = imData / maxMass
    rngIt[np.where(rngIt < -0.99)] = -0.99
    rngIt[np.where(rngIt > 0.99)] = 0.99

    
    # And make the plot
    fig, ax = plt.subplots(1,1,  figsize=(8,8))
    plt.imshow(rngIt, cmap='binary_r', origin='lower')
    plt.contour(mask, levels=[-3], colors=['red'])
    ax.axis('off')
    plt.subplots_adjust(left=0.,right=1.,top=1,bottom=0)
    plt.plot(sx,sy,'ro')
    plt.plot(avgx,avgy,'ro')
    plt.show()'''
    
    return sx, sy, avgx, avgy, heightPix, heightRs
    

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
                            CPAs  = sav_data['cat'][0][4][0][0] 
                            AWs  = sav_data['cat'][0][5][0][0] 
                            hts = sav_data['cat'][0][8][0][0]
                            
                            
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
                                
                                # Get the center of mass
                                
                                lastBit = savFile.split('/')[-1]
                                endIdx  = savFile.find(lastBit)
                                ftsFold = savFile[:endIdx]+'fts/'
                                ftsFile = lastBit.replace('.sav', '')+'_'+myFile.split('/')[-1].replace('.fts','_mass.fts')
                                if os.path.exists(ftsFold+ftsFile):
                                    sx, sy, avgx, avgy, heightPix, heightRs = calcCoM(ftsFold+ftsFile, masks[i])
                                    stuff = [lastBit.replace('.sav',''), myFile.split('/')[-1].replace('.fts','_mass.fts'), sx, sy, avgx, avgy, heightPix, heightRs, hts[i], CPAs[i], AWs[i]]
                                    outstr = ''
                                    for item in stuff:
                                        outstr += str(item) + ' '
                                    #print (outstr)
                                    f2.write(outstr+'\n')
                                else:
                                     print ('Cannot find ', ftsFold+ftsFile)
                                #print(sd)
                                #f2.write(bFold + ' ' + myFile + ' '+str(np.sum(masks[i] == -3)).rjust(8) + '\n')
                                
def checkOne(savPath, corID, sat):
    savFile = savPath+corID+'.sav'
    infoFile = savPath+corID+'.info'
    print (infoFile)
    if os.path.exists(savFile):   
        print('Reading in CORSET files at:')
        print(savFile)
        sav_data = readsav(savFile)
        masks = sav_data['cat'][0][10]
        dates = sav_data['cat'][0][2]
        times = sav_data['cat'][0][3]
        CPAs  = sav_data['cat'][0][4][0][0] 
        AWs  = sav_data['cat'][0][5][0][0] 
        hts = sav_data['cat'][0][8][0][0]
        
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
            myFile = infos[i+fileIdx0].replace(myPrefix, '').replace('\n','').replace('//','/')
            
            # Get the center of mass
            
            lastBit = savFile.split('/')[-1]
            endIdx  = savFile.find(lastBit)
            ftsFold = savFile[:endIdx]+'fts/'
            ftsFile = lastBit.replace('.sav', '')+'_'+myFile.split('/')[1].replace('.fts','_mass.fts')
            
            if os.path.exists(ftsFold+ftsFile):
                sx, sy, avgx, avgy, heightPix, heightRs = calcCoM(ftsFold+ftsFile, masks[i])
                stuff = [lastBit, myFile.split('/')[1].replace('.fts','_mass.fts'), sx, sy, avgx, avgy, heightPix, heightRs, hts[i], CPAs[i], AWs[i]]
                outstr = ''
                for item in stuff:
                    outstr += str(item) + ' '
                #print (outstr)    
                f2.write(outstr+'\n')
            else:
                print ('Cannot find ', ftsFold+ftsFile)
                                    


#checkOne('/Volumes/SRP/vourla1/laura/web_database/catalog/A/2014/02/140212/2666.9021_0A/', '2666.9021_0A', 'A')

#checkOne('/Volumes/SRP/vourla1/laura/web_database/catalog/B/2014/05/140505/2749.4333_0B/', '2749.4333_0B', 'B')


''' 
 missing 20140505_225400 2749.3917_0A 2749.4333_0B'''

f2 = open('tempCORSET_CoMs.dat', 'w')
for aYr in yrs:   
    for aSat in sats:
        print (aYr, aSat) 
        checkYr(aYr, aSat)
f2.close()