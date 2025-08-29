from astropy.io import fits
import os
import numpy as np


fitsPath = '/Users/kaycd1/STEREO_Mass/MassFits/'

haveDone = []

if False:
    # Print total mass and n good pix per image
    allFiles = os.listdir(fitsPath)

    for aFile in allFiles:
        if ('POS' not in aFile) & ('d4c2B' not in aFile):
            bFile = aFile.replace('_d4c2A_mass.fts', '_d4c2B_mass.fts')
            if bFile in allFiles:
                try:
                    hdul = fits.open(fitsPath+aFile) 
                    imDataA = hdul[0].data
                    hdul.close()
            
                    hdul = fits.open(fitsPath+bFile) 
                    imDataB = hdul[0].data
                    hdul.close()
            
                    goodIdxA =  np.where(imDataA != 0)
                    goodIdxB = np.where(imDataB != 0)
            
                    print ('   ', np.sum(imDataA[goodIdxA]), len(goodIdxA[0]), np.sum(imDataB[goodIdxB]), len(goodIdxB[0]))
                except:
                    pass
                    
if True:
    data = np.genfromtxt('fullImComp.txt', dtype=str)
    totA = data[:,0].astype(float)
    totB = data[:,2].astype(float)
    pixA = data[:,1].astype(float)
    pixB = data[:,3].astype(float)
    denA = totA / pixA
    denB = totB / pixB
    print (np.mean(totA), np.median(totA), np.mean(denA), np.median(denA))
    print (np.mean(totB), np.median(totB), np.mean(denB), np.median(denB))
    print(np.mean(totA/totB), np.median(totA/totB), np.mean(denA/denB), np.median(denA/denB))